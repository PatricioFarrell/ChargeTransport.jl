
#=
# GaAs diode (1D) (but with Discontinuous Quantities)
([source code](SOURCE_URL))

We simulate charge transport in a GaAs pin diode, where we we assume that the electric
charge carriers may not be necessarily continuous.

We may infer the existence of interface species at the interior boundaries.

(DA: note to myself, check again the facts here!!!)

was noch offen ist:
[ ]: für visualisierung von Interface species aus voronoifvm sachen exportieren.
     mit den sachen, die ich benutze, funktioniert visualisierung nicht.
[ ]: exportiere alle sachen die ich brauche aus VoronoiFVM wie norm, subgrids und views

=#

module Ex401_PIN_DiscontqF

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot
using DelimitedFiles

## This function is used to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_intrinsic, h_pdoping)
    coord_ndoping    = collect(range(0.0, stop = h_ndoping, length = 3 * refinementfactor))
    coord_intrinsic  = collect(range(h_ndoping, stop = (h_ndoping + h_intrinsic), length = 3 * refinementfactor))
    coord_pdoping    = collect(range((h_ndoping + h_intrinsic), stop = (h_ndoping + h_intrinsic + h_pdoping), length = 3 * refinementfactor))
    coord            = glue(coord_ndoping, coord_intrinsic)
    coord            = glue(coord, coord_pdoping)

    return coord
end

function main(;n = 6, plotting = false, verbose = false, test = false, interfaceSpecies = true)

    PyPlot.close("all")
    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionAcceptor          = 1          # p doped region
    regionIntrinsic         = 2          # intrinsic region
    regionDonor             = 3          # n doped region
    regions                 = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions         = length(regions)

    ## boundary region numbers
    bregionAcceptor         = 1
    bregionDonor            = 2
    bregionJunction1        = 3
    bregionJunction2        = 4

    ## grid
    h_pdoping               = 2.0 * μm
    h_intrinsic             = 2.0 * μm
    h_ndoping               = 2.0 * μm
    refinementfactor        = 2^(n-1)
    h_pdoping               = 2.0    * μm
    coord                   = initialize_pin_grid(refinementfactor, h_pdoping,
                                                  h_intrinsic, h_ndoping)

    grid                    = simplexgrid(coord)

    ## cellmask! for defining the subregions and assigning region number
    cellmask!(grid, [0.0 * μm],                 [h_pdoping],                           regionAcceptor)   # p-doped   region = 1
    cellmask!(grid, [h_pdoping],                [h_pdoping + h_intrinsic],             regionIntrinsic)  # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic],  [h_pdoping + h_intrinsic + h_ndoping], regionDonor)      # n-doped   region = 3

    bfacemask!(grid, [h_pdoping],               [h_pdoping],                           bregionJunction1) # first  inner interface
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic],             bregionJunction2) # second inner interface

    if plotting
        gridplot(grid, Plotter = PyPlot, legend=:lt, fignumber = 1)
        PyPlot.title("Grid")
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    ## set indices of the quasi Fermi potentials
    iphin                = 1 # electron quasi Fermi potential
    iphip                = 2 # hole quasi Fermi potential

    if interfaceSpecies
        iphinb           = 3
        iphipb           = 4
        numberOfCarriers = 4
    else
        numberOfCarriers = 2
    end

    # We define the physical data.
    Ec                   = 1.424                *  eV
    Ev                   = 0.0                  *  eV
    Nc                   = 4.351959895879690e17 / (cm^3)
    Nv                   = 9.139615903601645e18 / (cm^3)
    mun                  = 8500.0               * (cm^2) / (V * s)
    mup                  = 400.0                * (cm^2) / (V * s)
    εr                   = 12.9                 *  1.0    # relative dielectric permittivity of GAs
    T                    = 300.0                *  K

    ## recombination parameters
    Auger                = 1.0e-29              * cm^6 / s
    SRH_TrapDensity      = 1.0e10               / cm^3
    SRH_LifeTime         = 1.0                  * ns
    Radiative            = 1.0e-10              * cm^3 / s

    ## doping
    dopingFactorNd       = 1.0
    dopingFactorNa       = 0.46
    Nd                   = dopingFactorNd * Nc
    Na                   = dopingFactorNa * Nv

    ni                   = sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * kB * T)) # intrinsic concentration
    voltageAcceptor      = 1.5                  * V

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    # We initialize the Data instance and fill in predefined data.
    data                                = Data(grid, numberOfCarriers)
    data.modelType                      = Stationary
    data.F                             .= Boltzmann
    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                  bulk_recomb_Auger = true,
                                                                  bulk_recomb_radiative = true,
                                                                  bulk_recomb_SRH = true)
    # These are the crucial lines where we impose the discontinuity!!!
    data.isContinuous[iphin]            = false
    data.isContinuous[iphip]            = false
    data.boundaryType[bregionAcceptor]  = OhmicContact
    # DA: hä warte mal, wieso verschiebe ich das nicht alles einfach nach intern? Mindmap für mich selbst erstellen.
    # Ich such mir einfach mein InterfaceModel aus
    # Weil das logischste ist doch, dass User hier InterfaceModelNone schreibt. Wir intern wissen aber, dass wegen isContinuous = false, wir andere breaction nehmen ....
    # other possibility: InterfaceModelDiscontqFInterfaceSpecies
    data.boundaryType[bregionJunction2] = InterfaceModelDiscontqF
    data.boundaryType[bregionDonor]     = OhmicContact
    data.fluxApproximation             .= ScharfetterGummel

    if interfaceSpecies == false
        data.boundaryType[bregionJunction1] = InterfaceModelDiscontqF
    else
        data.boundaryType[bregionJunction1] = InterfaceModelDiscontqFInterfaceSpecies
        # wäre schöner, wenn pro iphinb nur iphin, das wäre toll.
        # DA: wir müssen auf jeden fall iphin und iphinb in einer methode gemeinsam haben, damit
        # wir intern wissen welche interface species zu welcher bulk species gehört ...
        enable_interface_carriers!(data, bulkSpecies = [iphin, iphip], interfaceSpecies = [iphinb, iphipb], boundaryRegion = bregionJunction1)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                                = Params(grid, numberOfCarriers)

    params.temperature                                    = T
    params.UT                                             = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                           = -1
    params.chargeNumbers[iphip]                           =  1
    # DA: sollen wir interface carriers eigene chargeNumbers geben? Ja eigtl schon oder ...
    # weil dann können wir auch dieses enable_Interface_carriers umschreiben (glaube ich)

    for ibreg in 1:2   # boundary region data
        params.bDensityOfStates[iphin, ibreg]             = Nc
        params.bDensityOfStates[iphip, ibreg]             = Nv
        params.bBandEdgeEnergy[iphin, ibreg]              = Ec
        params.bBandEdgeEnergy[iphip, ibreg]              = Ev
    end

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg]                   = εr * ε0

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]               = Nc
        params.densityOfStates[iphip, ireg]               = Nv
        params.bandEdgeEnergy[iphin, ireg]                = Ec
        params.bandEdgeEnergy[iphip, ireg]                = Ev
        params.mobility[iphin, ireg]                      = mun
        params.mobility[iphip, ireg]                      = mup

        ## recombination parameters
        params.recombinationRadiative[ireg]               = Radiative
        params.recombinationSRHLifetime[iphin, ireg]      = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg]      = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg]   = SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg]   = SRH_TrapDensity
        params.recombinationAuger[iphin, ireg]            = Auger
        params.recombinationAuger[iphip, ireg]            = Auger

    end

    ## interior doping
    params.doping[iphin, regionDonor]                     = Nd
    params.doping[iphin, regionIntrinsic]                 = ni
    params.doping[iphip, regionIntrinsic]                 = 0.0
    params.doping[iphip, regionAcceptor]                  = Na

    ## boundary doping
    params.bDoping[iphin, bregionDonor]                   = Nd
    params.bDoping[iphip, bregionAcceptor]                = Na


    if interfaceSpecies
        δn                                                =  0.2 * eV
        δp                                                = -0.2 * eV
        data.d                                            = 6.28 * 10e-7 * cm
        params.bDensityOfStates[iphinb, bregionJunction1] = data.d * params.densityOfStates[iphin, regionIntrinsic]
        params.bDensityOfStates[iphipb, bregionJunction1] = data.d * params.densityOfStates[iphip, regionIntrinsic]

        params.bBandEdgeEnergy[iphinb, bregionJunction1]  = params.bandEdgeEnergy[iphin, regionIntrinsic] + δn
        params.bBandEdgeEnergy[iphipb, bregionJunction1]  = params.bandEdgeEnergy[iphip, regionIntrinsic] + δp

    else
        params.bReactDiscont[iphin, bregionJunction1]     = 1.0e15
        params.bReactDiscont[iphip, bregionJunction1]     = 1.0e15
    end

    # If you decrease these values you will observe discontinuity in the qFs.
    params.bReactDiscont[iphin, bregionJunction2]         = 1.0e15
    params.bReactDiscont[iphip, bregionJunction2]         = 1.0e15

    data.params                                           = params
    ctsys                                                 = System(grid, data, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Define outerior boundary conditions")
    end
    ################################################################################

    ## We set zero voltage ohmic contacts for each charge carrier at all outerior boundaries
    ## for the equilibrium calculations.
    set_contact!(ctsys, bregionAcceptor, Δu = 0.0)
    set_contact!(ctsys, bregionDonor,    Δu = 0.0)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control         = NewtonControl()
    control.verbose = verbose

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    ## initialize solution and starting vectors
    initialGuess  = unknowns(ctsys)
    solution      = unknowns(ctsys)
    solution      = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    initialGuess .= solution

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    data.calculationType = OutOfEquilibrium

    biasValues           = range(0, stop = voltageAcceptor, length = 32)
    IV                   = zeros(0)

    for Δu in biasValues

        if test == false
            println("Δu  = ", Δu )
        end

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        solve!(solution, initialGuess, ctsys, control = control, tstep = Inf)
        #solution = VoronoiFVM.solve(initialGuess, ctsys.fvmsys, [0.0, 1e1],control = control)

        initialGuess .= solution

        ## get I-V data
        val = get_current_val(ctsys, solution)

        push!(IV,  val)

    end # bias loop

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Some plotting")
    end
    ################################################################################

    function compute_densities(icc, ireg, solicc, psi)

        eta = data.params.chargeNumbers[icc] ./ data.params.UT .* ( (solicc .- psi) .+ data.params.bandEdgeEnergy[icc, ireg] ./ q )

        return data.params.densityOfStates[icc, ireg] .* data.F[icc].(eta)
    end

    if plotting

        sol_ref  = readdlm("data/reference-sol-PIN.dat") # [coord sol_iphin sol_iphip sol_ipsi]
        IV_ref   = readdlm("data/reference-IV-PIN.dat")
        vis1     = GridVisualizer(Plotter = PyPlot, layout=(2,1), size = (600,670), xlabel = "space [m]", ylabel = "potential [V]", fignumber=2)
        vis2     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "space [m]", ylabel = "density [m\$^{-3}\$]", fignumber=3)
        vis3     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "voltage [V]", ylabel = "current density [Am\$^{-2}\$]",fignumber=4)

        subgridB = subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
        phin_sol = views(solution, data.chargeCarrierList[iphin], subgridB, ctsys.fvmsys)
        phip_sol = views(solution, data.chargeCarrierList[iphip], subgridB, ctsys.fvmsys)
        psi_sol  = views(solution, data.index_psi, subgridB, ctsys.fvmsys)

        if interfaceSpecies

            # this is unfortunately not working soooo good ...
            bgrid    = subgrids(data.chargeCarrierList[iphinb], ctsys.fvmsys)

            phinb_sol = views(solution, data.chargeCarrierList[iphinb], bgrid, ctsys.fvmsys)
            phipb_sol = views(solution, data.chargeCarrierList[iphipb], bgrid, ctsys.fvmsys)
        end

        ###############################################################################
        ##########                         Potentials                        ##########
        ###############################################################################
        for i in eachindex(phin_sol)
            scalarplot!(    vis1[1, 1], subgridB[i], psi_sol[i], clear = false, color=:blue, linewidth = 5)
            if i == 3
                scalarplot!(vis1[1, 1], subgridB[i], psi_sol[i], clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")
            end
        end
        scalarplot!(vis1[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                    legend =:best)

        ##########################
        for i in eachindex(phin_sol)
            scalarplot!(vis1[2, 1], subgridB[i], phin_sol[i], clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis1[2, 1], subgridB[i], phip_sol[i], clear = false, color=:red)

            if i == 3
                scalarplot!(vis1[2, 1], subgridB[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                scalarplot!(vis1[2, 1], subgridB[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$", color=:red)
            end
        end

        # DA: current way out, when waiting for changes within ExtendableGrids and GridVisualize
        if interfaceSpecies
            PyPlot.figure(2)
            #scalarplot!(vis1[2, 1], bgrid, phinb_sol, clear = false, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{\\varphi}_n \$")
            #scalarplot!(vis1[2, 1], bgrid, phipb_sol, clear = false, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{\\varphi}_p \$")
            PyPlot.plot(coord[3*refinementfactor], phinb_sol, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{\\varphi}_n \$")
            PyPlot.plot(coord[3*refinementfactor], phipb_sol, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{\\varphi}_p \$")
            println("value phin_b = ", phinb_sol)
            println("value phip_b = ", phipb_sol)
        end


        scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
        scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                    legend =:best, show = true)

        ###############################################################################
        ##########                         Densities                         ##########
        ###############################################################################

        for i in eachindex(phin_sol)
            scalarplot!(vis2, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis2, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, color=:red)
            if i == 3
            scalarplot!(vis2, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
            scalarplot!(vis2, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, label ="\$ n_p \$", color=:red)
            end

        end

        if interfaceSpecies
            PyPlot.figure(3)
            eta_nb = -1/ data.params.UT * ( (phinb_sol[1] - psi_sol[1][end]) + data.params.bBandEdgeEnergy[iphinb, bregionJunction1]/q )
            eta_pb =  1/ data.params.UT * ( (phipb_sol[1] - psi_sol[1][end]) + data.params.bBandEdgeEnergy[iphipb, bregionJunction1]/q )

            # DA: divide by d such that it is three dimensional again?
            nb     = data.params.bDensityOfStates[iphinb, bregionJunction1] * data.F[iphin](eta_nb)#./data.d
            pb     = data.params.bDensityOfStates[iphipb, bregionJunction1] * data.F[iphip](eta_pb)#./data.d

            println("value n_b = ", nb)
            println("value p_b = ", pb)

            PyPlot.semilogy(coord[3*refinementfactor], nb, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{n}_n \$")
            PyPlot.semilogy(coord[3*refinementfactor], pb, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{n}_p \$")

            #scalarplot!(vis2, coord[3*refinementfactor], nb, clear = false, marker = "x", markersize = 12,  color =:darkgreen, label = "\$ n_\\bar{n} \$")
            #scalarplot!(vis2, coord[3*refinementfactor], pb, clear = false, marker = "x", markersize = 12,  color =:darkred, label = "\$ n_\\bar{p} \$")

        end

        # since we have a homogeneous set of parameters, region does not matter
        n = compute_densities(iphin, 1, sol_ref[:, 2], sol_ref[:, 4])
        p = compute_densities(iphip, 1, sol_ref[:, 3], sol_ref[:, 4])

        scalarplot!(vis2, sol_ref[:, 1], n, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
        scalarplot!(vis2, sol_ref[:, 1], p, clear = false, label = "ref sol", legend =:best, show = true)

        ###############################################################################
        ##########                            IV                             ##########
        ###############################################################################
        scalarplot!(vis3, biasValues,   abs.(IV),           clear = false, color=:green, linewidth = 5)
        scalarplot!(vis3, IV_ref[:, 1], abs.(IV_ref[:, 2]), clear = false, color=:black, linestyle=:dot)

    end # plotting

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solution))/length(solution) # when using sparse storage, we get NaN values in solution
    return testval

    if test == false
        println("*** done\n")
    end
end #  main

function test()
    testvalinterfaceSpecies = 0.32708816536853264; testvalDiscont = 0.4256408342166736
    main(test = true, interfaceSpecies = true) ≈ testvalinterfaceSpecies && main(test = true, interfaceSpecies = false) ≈ testvalDiscont
end

if test == false
    println("This message should show when the PIN module has successfully recompiled.")
end

end # module
