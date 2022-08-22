
module PIN_InterfaceSpecies

using VoronoiFVM       # PDE solver with a FVM spatial discretization
using ChargeTransport  # drift-diffusion solver
using ExtendableGrids  # grid initializer
using GridVisualize    # grid visualizer
using PyPlot           # solution visualizer
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

function main(;n = 5, Plotter = PyPlot, plotting = false, verbose = false, test = true, unknown_storage=:dense)

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
    bregions                = [bregionAcceptor, bregionDonor, bregionJunction1, bregionJunction2]
    numberOfBoundaryRegions = length(bregions)

    ## grid
    h_pdoping               = 2.0 * μm
    h_intrinsic             = 2.0 * μm
    h_ndoping               = 2.0 * μm
    refinementfactor        = 2^(n-1)
    h_pdoping               = 2.0    * μm
    coord                   = initialize_pin_grid(refinementfactor, h_ndoping,
                                                  h_intrinsic, h_pdoping)

    grid                    = simplexgrid(coord)

    ## cellmask! for defining the subregions and assigning region number
    cellmask!(grid, [0.0 * μm],                [h_pdoping],                           regionAcceptor)  # p-doped region = 1
    cellmask!(grid, [h_pdoping],               [h_pdoping + h_intrinsic],             regionIntrinsic) # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor)     # n-doped region = 3

    bfacemask!(grid, [h_pdoping],               [h_pdoping],                           bregionJunction1)  # first  inner interface
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic],             bregionJunction2)  # second inner interface

    if plotting
        gridplot(grid, Plotter = Plotter, legend=:lt)
        Plotter.title("Grid")
        Plotter.figure()
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
    iphin              = 1 # electron quasi Fermi potential
    iphip              = 2 # hole quasi Fermi potential
    iphin_b1           = 3
    iphip_b1           = 4
    numberOfCarriers   = 4

    # We define the physical data.
    Ec                 = 1.424                *  eV
    Ev                 = 0.0                  *  eV
    Nc                 = 4.351959895879690e17 / (cm^3)
    Nv                 = 9.139615903601645e18 / (cm^3)
    mun                = 8500.0               * (cm^2) / (V * s)
    mup                = 400.0                * (cm^2) / (V * s)
    εr                 = 12.9                 *  1.0    # relative dielectric permittivity of GAs
    T                  = 300.0                *  K

    ## recombination parameters
    Auger             = 1.0e-29              * cm^6 / s
    SRH_TrapDensity   = 1.0e10               / cm^3
    SRH_LifeTime      = 1.0                  * ns
    Radiative         = 1.0e-10              * cm^3 / s

    ## doping
    dopingFactorNd    = 1.0
    dopingFactorNa    = 0.46
    Nd                = dopingFactorNd * Nc
    Na                = dopingFactorNa * Nv

    ## intrinsic concentration
    ni                = sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * kB * T))

    ## contact voltages: we impose an applied voltage only on one boundary.
    ## At the other boundary the applied voltage is zero.
    voltageAcceptor   = 1.5                  * V

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
                                                                  bulk_recomb_Auger = false,
                                                                  bulk_recomb_radiative = false,
                                                                  bulk_recomb_SRH = false)

    data.isContinuous[iphin]            = false
    data.isContinuous[iphip]            = false

    data.boundaryType[bregionAcceptor]  = OhmicContact
    data.boundaryType[bregionJunction1] = InterfaceModelDiscontqF
    data.boundaryType[bregionJunction2] = InterfaceModelDiscontqFNoReaction
    data.boundaryType[bregionDonor]     = OhmicContact

    # wäre schöner, wenn pro iphin_b1 nur iphin, das wäre toll.
    enable_interface_carriers!(data, bulkSpecies = [iphin, iphip], interfaceSpecies = [iphin_b1, iphip_b1], boundaryRegion = bregionJunction1)

    data.fluxApproximation             .= ScharfetterGummel

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                              = Params(grid, numberOfCarriers)

    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1

    for ibreg in 1:2   # boundary region data
        params.bDensityOfStates[iphin, ibreg]           = Nc
        params.bDensityOfStates[iphip, ibreg]           = Nv
        params.bBandEdgeEnergy[iphin, ibreg]            = Ec
        params.bBandEdgeEnergy[iphip, ibreg]            = Ev
    end

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg]                 = εr * ε0

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nc
        params.densityOfStates[iphip, ireg]             = Nv
        params.bandEdgeEnergy[iphin, ireg]              = Ec
        params.bandEdgeEnergy[iphip, ireg]              = Ev
        params.mobility[iphin, ireg]                    = mun
        params.mobility[iphip, ireg]                    = mup

        ## recombination parameters
        params.recombinationRadiative[ireg]             = Radiative
        params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

    end

    ## inner boundary region data
    data.d                                              = 6.28 * 10e-7 * cm  # lattice size of perovskite (from Eames et al.):6.28 * 10e-8 * cm
    params.bDensityOfStates[iphin_b1, bregionJunction1] = data.d * params.densityOfStates[iphin, regionIntrinsic]
    params.bDensityOfStates[iphip_b1, bregionJunction1] = data.d * params.densityOfStates[iphip, regionIntrinsic]

    params.bBandEdgeEnergy[iphin_b1, bregionJunction1]  = params.bandEdgeEnergy[iphin, regionIntrinsic]
    params.bBandEdgeEnergy[iphip_b1, bregionJunction1]  = params.bandEdgeEnergy[iphip, regionIntrinsic]

    params.bReactDiscont[iphin, bregionJunction2]       = 1.0e15
    params.bReactDiscont[iphip, bregionJunction2]       = 1.0e15

    ## interior doping
    params.doping[iphin, regionDonor]                   = Nd
    params.doping[iphin, regionIntrinsic]               = ni
    params.doping[iphip, regionIntrinsic]               = 0.0
    params.doping[iphip, regionAcceptor]                = Na

    ## boundary doping
    params.bDoping[iphin, bregionDonor]                 = Nd
    params.bDoping[iphip, bregionAcceptor]              = Na

    # Region dependent params is now a substruct of data which is again a substruct of the system and will be parsed
    # in next step.
    data.params                                         = params

    # In the last step, we initialize our system with previous data which is likewise dependent on the parameters.
    # It is important that this is in the end, otherwise our VoronoiFVMSys is not dependent on the data we initialized
    # but rather on default data.
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    ################################################################################
    if test == false
        println("Define outerior boundary conditions")
    end
    ################################################################################

    ## We set zero voltage ohmic contacts for each charge carrier at all outerior boundaries
    ## for the equilibrium calculations.
    set_contact!(ctsys, bregionAcceptor, Δu = 0.0)
    set_contact!(ctsys, bregionDonor, Δu = 0.0)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control                   = NewtonControl()
    control.verbose           = verbose
    control.damp_initial      = 0.5
    control.damp_growth       = 1.21
    control.max_iterations    = 250
    control.tol_absolute      = 1.0e-12
    control.tol_relative      = 1.0e-12
    control.handle_exceptions = true
    control.tol_round         = 1.0e-8
    control.max_round         = 6

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    control.damp_initial  = 0.5
    control.damp_growth   = 1.2 # >= 1
    control.max_round     = 3

    ## initialize solution and starting vectors
    initialGuess          = unknowns(ctsys)
    solution              = unknowns(ctsys)

    solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    initialGuess         .= solution

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    data.calculationType      = OutOfEquilibrium


    maxBias    = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
    biasValues = range(0, stop = voltageAcceptor, length = 21)

    IV         = zeros(0)

    i = 0
    for Δu in biasValues

        i = i+1
        println("Δu  = ", Δu )

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        solve!(solution, initialGuess, ctsys, control = control, tstep = Inf)
        #solution = VoronoiFVM.solve(initialGuess, ctsys.fvmsys, [0.0, 1e1],control = control)

        initialGuess .= solution

        ## get I-V data

        factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
        tf      = testfunction(factory, [1], [2])
        I       = integrate(ctsys.fvmsys, tf, initialGuess)

        val = 0.0
        for ii = 1:length(I)-1
            val = val + abs(I[ii])
        end

        push!(IV,  val)



    end # bias loop

    function compute_densities(icc, ireg, phin, psi)
        eta = data.params.chargeNumbers[icc] ./ data.params.UT .* ( (phin .- psi) .+ data.params.bandEdgeEnergy[icc, ireg] ./ q )

        return data.params.densityOfStates[icc, ireg] .* data.F[icc].(eta)
    end
    # writedlm("reference-sol-PIN.dat", [coord solution'])
    # res = [biasValues IV]
    # writedlm("reference-IV-PIN.dat", res)

    vis = GridVisualizer(Plotter = PyPlot, layout=(3,1))

    subgrids = VoronoiFVM.subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
    phin_sol = VoronoiFVM.views(initialGuess, data.chargeCarrierList[iphin], subgrids, ctsys.fvmsys)
    phip_sol = VoronoiFVM.views(initialGuess, data.chargeCarrierList[iphip], subgrids, ctsys.fvmsys)
    psi_sol  = VoronoiFVM.views(initialGuess, data.index_psi, subgrids, ctsys.fvmsys)

    for i = 1:length(phin_sol)
        scalarplot!(vis[1, 1], subgrids[i], phin_sol[i], clear = false, color=:green)
        scalarplot!(vis[1, 1], subgrids[i], phip_sol[i], clear = false, color=:red)
        scalarplot!(vis[1, 1], subgrids[i], psi_sol[i],  clear = false, color=:blue)
        if i == 3
            scalarplot!(vis[1, 1], subgrids[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
            scalarplot!(vis[1, 1], subgrids[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$",  color=:red)
            scalarplot!(vis[1, 1], subgrids[i], psi_sol[i],  clear = false, label = "\$ \\psi \$",color=:blue)
        end
    end

    sol_ref = readdlm("data/reference-sol-PIN.dat")
    PyPlot.plot(sol_ref[:, 1], sol_ref[:, 2], linestyle="--", color = "black")
    PyPlot.plot(sol_ref[:, 1], sol_ref[:, 3], linestyle="--", color = "black")
    PyPlot.plot(sol_ref[:, 1], sol_ref[:, 4], linestyle="--", color = "black")
    Plotter.legend(fancybox = true, loc = "best", fontsize=11)
    Plotter.title("Solution with Bias")

    for i = 1:length(phin_sol)
        scalarplot!(vis[2, 1], subgrids[i], log.(compute_densities(iphin, subgrids[i][CellRegions][1], phin_sol[i], psi_sol[i])), clear = false, color=:green)
        scalarplot!(vis[2, 1], subgrids[i], log.(compute_densities(iphip, subgrids[i][CellRegions][1], phip_sol[i], psi_sol[i])), clear = false, color=:red)
    end
    ##########################################################
    scalarplot!(vis[3, 1], biasValues, log.(IV), clear = false, color=:green)
    IV_ref         = readdlm("data/reference-IV-PIN.dat")
    PyPlot.plot(IV_ref[:, 1], log.(IV_ref[:, 2]), linestyle="--", color = "black")

    testval = VoronoiFVM.norm(ctsys.fvmsys, solution, 2)
    return testval

    if test == false
        println("*** done\n")
    end
end #  main

function test()
    testval = 24.984716235891845
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when the PIN module has successfully recompiled.")
end

end # module
