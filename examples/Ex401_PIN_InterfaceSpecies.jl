
#=
# GaAs diode (1D) (but with Discontinuous Quantities)
([source code](SOURCE_URL))

We simulate charge transport in a GaAs pin diode, where we we assume that the electric
charge carriers may not be necessarily continuous.

We infer the existence of interface species at one of the interior boundaries.


=#

module Ex401_PIN_InterfaceSpecies

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

function main(;n = 6, plotting = false, verbose = false, test = false, interfaceSpecies = true, leftInterface = true)

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
        ipsi             = 5
    else
        numberOfCarriers = 2
        ipsi             = 3
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
    data.boundaryType[bregionAcceptor]  = OhmicContact
    data.boundaryType[bregionDonor]     = OhmicContact
    data.fluxApproximation             .= ScharfetterGummel

    if leftInterface == true
        bregActive = bregionJunction1
        bregDeact  = bregionJunction2
        icoordJ    = 3*refinementfactor

        regl       = regionAcceptor
        regr       = regionIntrinsic
    else
        bregActive = bregionJunction2
        bregDeact  = bregionJunction1
        icoordJ    = 2 * 3*refinementfactor

        regl       = regionIntrinsic
        regr       = regionDonor
    end

    if interfaceSpecies
        enable_interface_carrier!(data, bulkCarrier = iphin, interfaceCarrier = iphinb, bregions = [bregActive])
        enable_interface_carrier!(data, bulkCarrier = iphip, interfaceCarrier = iphipb, bregions = [bregActive])
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
    if interfaceSpecies
        params.chargeNumbers[iphinb]                      = -1
        params.chargeNumbers[iphipb]                      =  1
    end


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

        data.d                                         = 6.28 * 10e-7 * cm

        if leftInterface
            dopingN                                    = 0.0
            dopingP                                    = data.d * Na
        else
            dopingN                                    = data.d * Nd
            dopingP                                    = 0.0
        end

        δn                                             =  -0.1  * eV
        δp                                             =   0.1  * eV
        params.bDoping[iphinb, bregActive]             = dopingN
        params.bDoping[iphipb, bregActive]             = dopingP

        params.bDensityOfStates[iphinb, bregActive]    = data.d * Nc
        params.bDensityOfStates[iphipb, bregActive]    = data.d * Nv

        params.bBandEdgeEnergy[iphinb, bregActive]     = Ec + δn
        params.bBandEdgeEnergy[iphipb, bregActive]     = Ev + δp

        params.bReactionCoefficient[iphin, bregActive] = params.UT * mun * Nc/data.d
        params.bReactionCoefficient[iphip, bregActive] = params.UT * mup * Nv/data.d

        # For the other interface, where we do not have interface species, we infer a high
        # reaction coefficient such that we observe continuity. If you still observe discontinuity
        # at the other interface without interface species, increase this value.
        params.bReactionCoefficient[iphin, bregDeact] = 1.0e15
        params.bReactionCoefficient[iphip, bregDeact] = 1.0e15
    end

    data.params                                           = params
    ctsys                                                 = System(grid, data, unknown_storage=:sparse)

    if test == false
        println("*** done\n")
    end

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

    #writedlm("data/reference-sol-PIN-EQ.dat", [coord solution'])

    if plotting
        function compute_densities(icc, ireg, solicc, psi)

            eta = data.params.chargeNumbers[icc] ./ data.params.UT .* ( (solicc .- psi) .+ data.params.bandEdgeEnergy[icc, ireg] ./ q )

            return data.params.densityOfStates[icc, ireg] .* data.F[icc].(eta)
        end

        sol_ref  = readdlm("data/reference-sol-PIN-EQ.dat") # [coord sol_iphin sol_iphip sol_ipsi]
        vis1     = GridVisualizer(Plotter = PyPlot, layout=(2,1), size = (600,670), xlabel = "space [m]", ylabel = "potential [V]", fignumber=2)
        vis2     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "space [m]", ylabel = "density [m\$^{-3}\$]", fignumber=3)

        if interfaceSpecies
            subgridB = subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
            phin_sol = views(solution, data.chargeCarrierList[iphin], subgridB, ctsys.fvmsys)
            phip_sol = views(solution, data.chargeCarrierList[iphip], subgridB, ctsys.fvmsys)
            psi_sol  = views(solution, data.index_psi, subgridB, ctsys.fvmsys)

            # this is unfortunately not working soooo good ...
            bgrid    = subgrids(data.chargeCarrierList[iphinb], ctsys.fvmsys)

            phinb_sol = views(solution, data.chargeCarrierList[iphinb], bgrid, ctsys.fvmsys)
            phipb_sol = views(solution, data.chargeCarrierList[iphipb], bgrid, ctsys.fvmsys)
        else
            subgridp  = subgrid(grid, [1])
            subgridi  = subgrid(grid, [2])
            subgridn  = subgrid(grid, [3])
            # actually, we do not need this cases, since for this set-up we have continuity
            # in the densities with no effects on interfaces, but still for demonstrational
            # purpose.
            psi_solp  = view(solution[ipsi, :],  subgridp)
            psi_soli  = view(solution[ipsi, :],  subgridi)
            psi_soln  = view(solution[ipsi, :],  subgridn)

            phin_solp = view(solution[iphin, :], subgridp)
            phin_soli = view(solution[iphin, :], subgridi)
            phin_soln = view(solution[iphin, :], subgridn)

            phip_solp = view(solution[iphip, :], subgridp)
            phip_soli = view(solution[iphip, :], subgridi)
            phip_soln = view(solution[iphip, :], subgridn)
        end

        ###############################################################################
        ##########                         Potentials                        ##########
        ###############################################################################
        if interfaceSpecies
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
            PyPlot.figure(2)
            PyPlot.plot(coord[icoordJ], phinb_sol, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{\\varphi}_n \$")
            PyPlot.plot(coord[icoordJ], phipb_sol, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{\\varphi}_p \$")
            println("value phin_b = ", phinb_sol)
            println("value phip_b = ", phipb_sol)

            scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
            scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best, show = true)
        else
            scalarplot!(vis1[1, 1], subgridp, psi_solp,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis1[1, 1], subgridi, psi_soli,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis1[1, 1], subgridn, psi_soln,  clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")
            scalarplot!(vis1[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best)
            scalarplot!(vis1[2, 1], subgridp, phin_solp, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis1[2, 1], subgridp, phip_solp, clear = false, color=:red)
            scalarplot!(vis1[2, 1], subgridi, phin_soli, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis1[2, 1], subgridi, phip_soli, clear = false, color=:red)
            scalarplot!(vis1[2, 1], subgridn, phin_soln, clear = false, label = "\$ \\varphi_n \$", color=:green)
            scalarplot!(vis1[2, 1], subgridn, phip_soln, clear = false, label = "\$ \\varphi_p \$", color=:red)
            scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
            scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best, show = true)

        end

         ###############################################################################
        ##########                         Densities                         ##########
        ###############################################################################

        if interfaceSpecies
            for i in eachindex(phin_sol)
                scalarplot!(vis2, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, linewidth = 5, yscale=:log)
                scalarplot!(vis2, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, color=:red)
                if i == 3
                    scalarplot!(vis2, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
                    scalarplot!(vis2, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, label ="\$ n_p \$", color=:red)
                end

            end

            PyPlot.figure(3)

            if leftInterface
                psib         = psi_sol[1][end]
            else
                psib         = psi_sol[2][end]
            end

            eta_nb = -1/ data.params.UT * ( (phinb_sol[1] - psib) + data.params.bBandEdgeEnergy[iphinb, bregActive]/q )
            eta_pb =  1/ data.params.UT * ( (phipb_sol[1] - psib) + data.params.bBandEdgeEnergy[iphipb, bregActive]/q )
            # DA: divide by d such that it is three dimensional again?
            nb     = data.params.bDensityOfStates[iphinb, bregActive] * data.F[iphin](eta_nb)#./data.d
            pb     = data.params.bDensityOfStates[iphipb, bregActive] * data.F[iphip](eta_pb)#./data.d

            println("value n_b = ", nb)
            println("value p_b = ", pb)

            PyPlot.semilogy(coord[icoordJ], nb, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{n}_n \$")
            PyPlot.semilogy(coord[icoordJ], pb, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{n}_p \$")

            # since we have a homogeneous set of parameters, region does not matter
            n = compute_densities(iphin, 1, sol_ref[:, 2], sol_ref[:, 4])
            p = compute_densities(iphip, 1, sol_ref[:, 3], sol_ref[:, 4])

            scalarplot!(vis2, sol_ref[:, 1], n, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
            scalarplot!(vis2, sol_ref[:, 1], p, clear = false, label = "ref sol", legend =:best, show = true)

        else
            scalarplot!(vis2, subgridp, compute_densities(iphin, regionAcceptor,  phin_solp, psi_solp), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis2, subgridp, compute_densities(iphip, regionAcceptor,  phip_solp, psi_solp), clear = false, color=:red)
            scalarplot!(vis2, subgridi, compute_densities(iphin, regionIntrinsic, phin_soli, psi_soli), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis2, subgridi, compute_densities(iphip, regionIntrinsic, phip_soli, psi_soli), clear = false, color=:red)
            scalarplot!(vis2, subgridn, compute_densities(iphin, regionDonor,     phin_soln, psi_soln), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
            scalarplot!(vis2, subgridn, compute_densities(iphip, regionDonor,     phip_soln, psi_soln), clear = false, label ="\$ n_p \$", color=:red)

            # since we have a homogeneous set of parameters, region does not matter
            n = compute_densities(iphin, 1, sol_ref[:, 2], sol_ref[:, 4])
            p = compute_densities(iphip, 1, sol_ref[:, 3], sol_ref[:, 4])

            scalarplot!(vis2, sol_ref[:, 1], n, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
            scalarplot!(vis2, sol_ref[:, 1], p, clear = false, label = "ref sol", legend =:best, show = true)

        end

    end

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

    if plotting

        sol_ref  = readdlm("data/reference-sol-PIN.dat") # [coord sol_iphin sol_iphip sol_ipsi]
        IV_ref   = readdlm("data/reference-IV-PIN.dat")
        vis3     = GridVisualizer(Plotter = PyPlot, layout=(2,1), size = (600,670), xlabel = "space [m]", ylabel = "potential [V]", fignumber=4)
        vis4     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "space [m]", ylabel = "density [m\$^{-3}\$]", fignumber=5)
        vis5     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "voltage [V]", ylabel = "current density [Am\$^{-2}\$]",fignumber=6)

        if interfaceSpecies
            subgridB = subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
            phin_sol = views(solution, data.chargeCarrierList[iphin], subgridB, ctsys.fvmsys)
            phip_sol = views(solution, data.chargeCarrierList[iphip], subgridB, ctsys.fvmsys)
            psi_sol  = views(solution, data.index_psi, subgridB, ctsys.fvmsys)

            # this is unfortunately not working soooo good ...
            bgrid    = subgrids(data.chargeCarrierList[iphinb], ctsys.fvmsys)

            phinb_sol = views(solution, data.chargeCarrierList[iphinb], bgrid, ctsys.fvmsys)
            phipb_sol = views(solution, data.chargeCarrierList[iphipb], bgrid, ctsys.fvmsys)
        else
            subgridp  = subgrid(grid, [1])
            subgridi  = subgrid(grid, [2])
            subgridn  = subgrid(grid, [3])
            # actually, we do not need this cases, since for this set-up we have continuity
            # in the densities with no effects on interfaces, but still for demonstrational
            # purpose.
            psi_solp  = view(solution[ipsi, :],  subgridp)
            psi_soli  = view(solution[ipsi, :],  subgridi)
            psi_soln  = view(solution[ipsi, :],  subgridn)

            phin_solp = view(solution[iphin, :], subgridp)
            phin_soli = view(solution[iphin, :], subgridi)
            phin_soln = view(solution[iphin, :], subgridn)

            phip_solp = view(solution[iphip, :], subgridp)
            phip_soli = view(solution[iphip, :], subgridi)
            phip_soln = view(solution[iphip, :], subgridn)


        end

        ###############################################################################
        ##########                         Potentials                        ##########
        ###############################################################################
        if interfaceSpecies
            for i in eachindex(phin_sol)
                scalarplot!(    vis3[1, 1], subgridB[i], psi_sol[i], clear = false, color=:blue, linewidth = 5)
                if i == 3
                    scalarplot!(vis3[1, 1], subgridB[i], psi_sol[i], clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")
                end
            end
            scalarplot!(vis3[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best)

            ##########################
            for i in eachindex(phin_sol)
                scalarplot!(vis3[2, 1], subgridB[i], phin_sol[i], clear = false, color=:green, linestyle=:solid, linewidth = 5)
                scalarplot!(vis3[2, 1], subgridB[i], phip_sol[i], clear = false, color=:red)

                if i == 3
                    scalarplot!(vis3[2, 1], subgridB[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                    scalarplot!(vis3[2, 1], subgridB[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$", color=:red)
                end
            end

            # DA: current way out, when waiting for changes within ExtendableGrids and GridVisualize
            PyPlot.figure(4)
            PyPlot.plot(coord[icoordJ], phinb_sol, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{\\varphi}_n \$")
            PyPlot.plot(coord[icoordJ], phipb_sol, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{\\varphi}_p \$")
            println("value phin_b = ", phinb_sol)
            println("value phip_b = ", phipb_sol)

            scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
            scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best, show = true)
        else
            scalarplot!(vis3[1, 1], subgridp, psi_solp,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis3[1, 1], subgridi, psi_soli,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis3[1, 1], subgridn, psi_soln,  clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")
            scalarplot!(vis3[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best)
            scalarplot!(vis3[2, 1], subgridp, phin_solp, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis3[2, 1], subgridp, phip_solp, clear = false, color=:red)
            scalarplot!(vis3[2, 1], subgridi, phin_soli, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis3[2, 1], subgridi, phip_soli, clear = false, color=:red)
            scalarplot!(vis3[2, 1], subgridn, phin_soln, clear = false, label = "\$ \\varphi_n \$", color=:green)
            scalarplot!(vis3[2, 1], subgridn, phip_soln, clear = false, label = "\$ \\varphi_p \$", color=:red)
            scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
            scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best, show = true)

        end
        ###############################################################################
        ##########                         Densities                         ##########
        ###############################################################################

        if interfaceSpecies
            for i in eachindex(phin_sol)
                scalarplot!(vis4, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, linewidth = 5, yscale=:log)
                scalarplot!(vis4, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, color=:red)
                if i == 3
                    scalarplot!(vis4, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
                    scalarplot!(vis4, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, label ="\$ n_p \$", color=:red)
                end

            end

            PyPlot.figure(5)


            if leftInterface
                psib         = psi_sol[1][end]
            else
                psib         = psi_sol[2][end]
            end

            eta_nb = -1/ data.params.UT * ( (phinb_sol[1] - psib) + data.params.bBandEdgeEnergy[iphinb, bregActive]/q )
            eta_pb =  1/ data.params.UT * ( (phipb_sol[1] - psib) + data.params.bBandEdgeEnergy[iphipb, bregActive]/q )

            # DA: divide by d such that it is three dimensional again?
            nb     = data.params.bDensityOfStates[iphinb, bregActive] * data.F[iphin](eta_nb)#./data.d
            pb     = data.params.bDensityOfStates[iphipb, bregActive] * data.F[iphip](eta_pb)#./data.d

            println("value n_b = ", nb)
            println("value p_b = ", pb)

            PyPlot.semilogy(coord[icoordJ], nb, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{n}_n \$")
            PyPlot.semilogy(coord[icoordJ], pb, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{n}_p \$")

            # since we have a homogeneous set of parameters, region does not matter
            n = compute_densities(iphin, 1, sol_ref[:, 2], sol_ref[:, 4])
            p = compute_densities(iphip, 1, sol_ref[:, 3], sol_ref[:, 4])

            scalarplot!(vis4, sol_ref[:, 1], n, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
            scalarplot!(vis4, sol_ref[:, 1], p, clear = false, label = "ref sol", legend =:best, show = true)

        else
            scalarplot!(vis4, subgridp, compute_densities(iphin, regionAcceptor,  phin_solp, psi_solp), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis4, subgridp, compute_densities(iphip, regionAcceptor,  phip_solp, psi_solp), clear = false, color=:red)
            scalarplot!(vis4, subgridi, compute_densities(iphin, regionIntrinsic, phin_soli, psi_soli), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis4, subgridi, compute_densities(iphip, regionIntrinsic, phip_soli, psi_soli), clear = false, color=:red)
            scalarplot!(vis4, subgridn, compute_densities(iphin, regionDonor,     phin_soln, psi_soln), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
            scalarplot!(vis4, subgridn, compute_densities(iphip, regionDonor,     phip_soln, psi_soln), clear = false, label ="\$ n_p \$", color=:red)

            # since we have a homogeneous set of parameters, region does not matter
            n = compute_densities(iphin, 1, sol_ref[:, 2], sol_ref[:, 4])
            p = compute_densities(iphip, 1, sol_ref[:, 3], sol_ref[:, 4])

            scalarplot!(vis4, sol_ref[:, 1], n, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
            scalarplot!(vis4, sol_ref[:, 1], p, clear = false, label = "ref sol", legend =:best, show = true)

        end

        ###############################################################################
        ##########                            IV                             ##########
        ###############################################################################
        scalarplot!(vis5, biasValues,   abs.(IV),           clear = false, color=:green, linewidth = 5)
        scalarplot!(vis5, IV_ref[:, 1], abs.(IV_ref[:, 2]), clear = false, color=:black, linestyle=:dot)

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
    testvalinterfaceSpecies = 0.32708816536853264; testval = 0.9896905415004197
    main(test = true, interfaceSpecies = true) ≈ testvalinterfaceSpecies && main(test = true, interfaceSpecies = false) ≈ testval
end

if test == false
    println("This message should show when the PIN module has successfully recompiled.")
end

end # module
