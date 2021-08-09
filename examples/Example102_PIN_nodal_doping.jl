#=
# 102: 1D GaAs p-i-n diode with spacially varying doping.
([source code](SOURCE_URL))

Simulating charge transport in a GaAs pin diode. This means
the corresponding PDE problem corresponds to the van Roosbroeck
system of equations.
The simulations are performed out of equilibrium and for the
stationary problem.
=#

module Example102_PIN_nodal_doping

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize

function main(;Plotter = nothing, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # region numbers
    regionAcceptor          = 1                           # p doped region
    regionIntrinsic         = 2                           # intrinsic region
    regionDonor             = 3                           # n doped region
    regions                 = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions         = length(regions)

    # boundary region numbers
    bregionAcceptor         = 1
    bregionDonor            = 2
    bregions                = [bregionAcceptor, bregionDonor]
    numberOfBoundaryRegions = length(bregions)

    h_pdoping              = 0.1 * μm
    h_intrinsic            = 0.1 * μm
    h_ndoping              = 0.1 * μm

    coord                  = range(0.0, stop = h_ndoping + h_intrinsic + h_pdoping, length = 25)
    coord                  = collect(coord)
    grid                   = simplexgrid(coord)
    numberOfNodes          = length(coord)

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],                [h_pdoping],                           regionAcceptor,  tol = 1.0e-15)    # p-doped region = 1
    cellmask!(grid, [h_pdoping],                [h_pdoping + h_intrinsic],            regionIntrinsic, tol = 1.0e-15)    # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor,     tol = 1.0e-15)    # n-doped region = 3

    if plotting
        gridplot(grid, Plotter = Plotter)
        Plotter.title("Grid")
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # set indices of the quasi Fermi potentials
    iphin              = 1 # electron quasi Fermi potential
    iphip              = 2 # hole quasi Fermi potential
    numberOfCarriers   = 2 

    # physical data
    Ec                 = 1.424                *  eV
    Ev                 = 0.0                  *  eV
    Nc                 = 4.351959895879690e17 / (cm^3)
    Nv                 = 9.139615903601645e18 / (cm^3)
    mun                = 8500.0               * (cm^2) / (V * s)
    mup                = 400.0                * (cm^2) / (V * s)
    εr                 = 12.9                 *  1.0              # relative dielectric permittivity of GAs
    T                  = 300.0                *  K


    # recombination model
    bulk_recombination = bulk_recomb_model_trap_assisted

    # recombination parameters
    SRH_TrapDensity_n  = 4.760185435081902e5    / cm^3       
    SRH_TrapDensity_p  = 9.996936448738406e6    / cm^3
    SRH_LifeTime       = 1.0                    * ps   

    # contact voltages
    voltageAcceptor    = 1.4 * V

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define ChargeTransportSystem and fill in information about model")
    end
    ################################################################################

    # initialize ChargeTransportData instance and fill in data
    data                                = ChargeTransportData(grid, numberOfCarriers)

    #### declare here all necessary information concerning the model ###

    # Following variable declares, if we want to solve stationary or transient problem
    data.model_type                     = model_stationary

    # Following choices are possible for F: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                             .= Boltzmann

    #Here the user can specify, if they assume continuous or discontinuous charge carriers. We note that for a surface recombination model,
    # we encourage to use discontinuous electron and hole quasi Fermi potentials.
    data.isContinuous[iphin]            = true
    data.isContinuous[iphip]            = true

    # Following choices are possible for recombination model: bulk_recomb_model_none, bulk_recomb_model_trap_assisted, bulk_recomb_radiative, bulk_recomb_full <: bulk_recombination_model 
    data.bulk_recombination             = set_bulk_recombination(iphin = iphin, iphip = iphip, bulk_recombination_model = bulk_recombination)

    # Following choices are possible for boundary model: For contacts currently only ohmic_contact and schottky_contact are possible.
    # For inner boundaries we have interface_model_none, interface_model_surface_recombination, interface_model_ion_charge
    # (distinguish between left and right).
    data.boundary_type[bregionAcceptor] = ohmic_contact                       
    data.boundary_type[bregionDonor]    = ohmic_contact   
    
    # Following choices are possible for the flux_discretization scheme: ScharfetterGummel, ScharfetterGummel_Graded,
    # excessChemicalPotential, excessChemicalPotential_Graded, diffusionEnhanced, generalized_SG
    data.flux_approximation             = ScharfetterGummel

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define ChargeTransportParams and fill in physical parameters")
    end
    ################################################################################

    # Params is a struct which contains all necessary physical parameters. If one wants to simulate
    # space-dependent variable, one additionally needs to generate a ParamsNodal struct as done here for the doping.
    params                                              = ChargeTransportParams(grid, numberOfCarriers)
    paramsnodal                                         = ChargeTransportParamsNodal(grid, numberOfCarriers)

    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1

    for ibreg in 1:numberOfBoundaryRegions   # boundary region data

        params.bDensityOfStates[iphin, ibreg]           = Nc
        params.bDensityOfStates[iphip, ibreg]           = Nv
        params.bBandEdgeEnergy[iphin, ibreg]            = Ec
        params.bBandEdgeEnergy[iphip, ibreg]            = Ev
    end

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg]                 = εr

        # effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nc
        params.densityOfStates[iphip, ireg]             = Nv
        params.bandEdgeEnergy[iphin, ireg]              = Ec
        params.bandEdgeEnergy[iphip, ireg]              = Ev
        params.mobility[iphin, ireg]                    = mun
        params.mobility[iphip, ireg]                    = mup

        # recombination parameters
        params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity_n
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity_p

    end

    # initialize the space dependent doping
    # doping
    NDoping           =   1.0e17  / cm^3
    κ = 500.0
    for icoord = 1:numberOfNodes
        paramsnodal.doping[icoord] = NDoping * 0.5 * ( 1.0  +  tanh( (0.1 - coord[icoord]/μm) *κ )  - ( 1.0 + tanh( (coord[icoord]/μm - 0.2) * κ )) )
    end

    # Region dependent params is now a substruct of data which is again a substruct of the system and will be parsed 
    # in next step.
    data.params                                         = params
    # same holds for space dependent params, i.e. the doping
    data.paramsnodal                                    = paramsnodal

    # in the last step, we initialize our system with previous data which is likewise dependent on the parameters. 
    # important that this is in the end, otherwise our VoronoiFVMSys is not dependent on the data we initialized
    # but rather on default data.
    ctsys                                               = ChargeTransportSystem(grid, data, unknown_storage=unknown_storage)

    if test == false
        println("*** done\n")
    end

    # show region dependent physical parameters. show_params() only supports region dependent parameters, but, if one wishes to
    # print nodal dependent parameters, currently this is possible with println(ctsys.data.paramsnodal). We neglected here, since
    # in most applications where the numberOfNodes is >> 10 this would results in a large output in the terminal.
    if test == false
        show_params(ctsys)
    end

    if plotting == true
        ################################################################################
        println("Plot doping")
        ################################################################################
        Plotter.figure()
        plot_doping(Plotter, grid, paramsnodal)
        println("*** done\n")
    end
    
    ################################################################################
    if test == false
        println("Define outerior boundary conditions")
    end
    ################################################################################

    # set ohmic contacts for each charge carrier at all outerior boundaries. First, 
    # we compute equilibrium solutions. Hence the boundary values at the ohmic contacts
    # are zero.
    set_ohmic_contact!(ctsys, iphin, bregionAcceptor, 0.0)
    set_ohmic_contact!(ctsys, iphip, bregionAcceptor, 0.0)
    set_ohmic_contact!(ctsys, iphin, bregionDonor, 0.0)
    set_ohmic_contact!(ctsys, iphip, bregionDonor, 0.0)

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
    control.damp_growth       = 1.21
    control.max_iterations    = 250
    control.tol_absolute      = 1.0e-14
    control.tol_relative      = 1.0e-14
    control.handle_exceptions = true
    control.tol_round         = 1.0e-8
    control.max_round         = 5

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium for Boltzmann")
    end
    ################################################################################

    # control.damp_initial      = 0.001
    # control.damp_growth       = 1.21 # >= 1
    # control.max_round         = 4

    # initialize solution and starting vectors
    initialGuess          = unknowns(ctsys)
    solution              = unknowns(ctsys)

    solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

    initialGuess         .= solution 

    if plotting
        Plotter.figure()
        plot_energies(Plotter, grid, data, solution, "Equilibrium")
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution, "Equilibrium")
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "Equilibrium")
        Plotter.figure()
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    ctsys.data.calculation_type      = outOfEquilibrium
     
    maxBias                          = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
    biasValues                       = range(0, stop = maxBias, length = 41)
    IV                               = zeros(0)

    w_device                         = 0.1    * μm  # width of device
    z_device                         = 1.0e-5 * cm  # depth of device

    for Δu in biasValues

        # set non equilibrium boundary conditions
        set_ohmic_contact!(ctsys, iphin, bregionAcceptor, Δu)
        set_ohmic_contact!(ctsys, iphip, bregionAcceptor, Δu)

        solve!(solution, initialGuess, ctsys, control = control, tstep = Inf)

        initialGuess .= solution

        # get IV curve
        factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)

        # testfunction zero in bregionAcceptor and one in bregionDonor
        tf     = testfunction(factory, [bregionAcceptor], [bregionDonor])
        I      = integrate(ctsys.fvmsys, tf, solution)

        push!(IV,  abs.(w_device * z_device * (I[iphin] + I[iphip])))

    end # bias loop

    
    if plotting # plot solution and IV curve
        plot_energies(Plotter, grid, data, solution, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = false)
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = true)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = true)
        Plotter.figure()
        plot_IV(Plotter, biasValues,IV, biasValues[end], plotGridpoints = true)
    end

    testval = solution[15]
    return testval

end #  main

function test()
    testval = 1.4676876302354516
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when the PIN module is successfully recompiled.")
end

end # module
