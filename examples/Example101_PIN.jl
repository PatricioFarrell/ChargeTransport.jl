#=
# 101: 1D GaAs p-i-n diode.
([source code](SOURCE_URL))

Simulating charge transport in a GaAs pin diode. This means
the corresponding PDE problem corresponds to the van Roosbroeck
system of equations.
The simulations are performed out of equilibrium and for the
stationary problem.
=#

module Example101_PIN

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize

# function for initializing the grid for a possble extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_intrinsic, h_pdoping)
    coord_ndoping    = collect(range(0.0, stop = h_ndoping, length = 3 * refinementfactor))
    coord_intrinsic  = collect(range(h_ndoping, stop = (h_ndoping + h_intrinsic), length = 3 * refinementfactor))
    coord_pdoping    = collect(range((h_ndoping + h_intrinsic), stop = (h_ndoping + h_intrinsic + h_pdoping), length = 3 * refinementfactor))
    coord            = glue(coord_ndoping, coord_intrinsic)
    coord            = glue(coord, coord_pdoping)

    return coord
end


function main(;n = 3, Plotter = nothing, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

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

    # grid
    refinementfactor        = 2^(n-1)
    h_pdoping               = 2 * μm
    h_intrinsic             = 2 * μm
    h_ndoping               = 2 * μm
    coord                   = initialize_pin_grid(refinementfactor,
                                                  h_pdoping,
                                                  h_intrinsic,
                                                  h_ndoping)

    grid                    = ExtendableGrids.simplexgrid(coord)

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm], [h_pdoping], regionAcceptor)        # p-doped region = 1
    cellmask!(grid, [h_pdoping], [h_pdoping + h_intrinsic], regionIntrinsic)    # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor)     # n-doped region = 3

    if plotting
        GridVisualize.gridplot(grid, Plotter = Plotter, legend=:lt)
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

    numberOfCarriers  = 2 # electrons and holes

    # physical data
    Ec                = 1.424                *  eV
    Ev                = 0.0                  *  eV
    Nc                = 4.351959895879690e17 / (cm^3)
    Nv                = 9.139615903601645e18 / (cm^3)
    mun               = 8500.0               * (cm^2) / (V * s)
    mup               = 400.0                * (cm^2) / (V * s)
    εr                = 12.9                 *  1.0              # relative dielectric permittivity of GAs
    T                 = 300.0                *  K


    # recombination model
    bulk_recombination = bulk_recombination_full

    # recombination parameters
    Auger             = 1.0e-29   * cm^6 / s          # 1.0e-41
    SRH_TrapDensity   = 1.0e10    / cm^3              # 1.0e16
    SRH_LifeTime      = 1.0       * ns                # 1.0e10
    Radiative         = 1.0e-10   * cm^3 / s          # 1.0e-16

    # doping
    dopingFactorNd    =   1.0
    dopingFactorNa    =   0.46
    Nd                =   dopingFactorNd * Nc
    Na                =   dopingFactorNa * Nv

    # intrinsic concentration (not doping!)
    ni                =   sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * kB * T)) 

    # contact voltages: we impose an applied voltage only on one boundary.
    # At the other boundary the applied voltage is zero.
    voltageAcceptor   = 1.5 * V

    # interface model
    interface_reaction = interface_model_none

    # set the correct indices for each species (this is needed for giving the user the correct index set)
    # but likewise it is possible to define one owns index set, i.e. iphin, iphin, ipsi = 1:3, but
    # one needs to be aware of the remarks within the documentary.
    indexSet         = set_indices!(grid, numberOfCarriers, interface_reaction)

    iphin           = indexSet["iphin"]
    iphip           = indexSet["iphip"]
    ipsi            = indexSet["ipsi" ]


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

    # Following choices are possible for recombination model: bulk_recombination_model_none, bulk_recombination_model_trap_assisted, bulk_recombination_radiative, bulk_recombination_full <: bulk_recombination_model 
    data.bulk_recombination_model       = bulk_recombination

    # Following choices are possible for boundary model: For contacts currently only ohmic_contact and schottky_contact are possible.
    # For inner boundaries we have interface_model_none, interface_model_surface_recombination, interface_model_ion_charge
    # (distinguish between left and right).
    data.boundary_type[bregionAcceptor] = ohmic_contact                       
    data.boundary_type[bregionDonor]    = ohmic_contact   
    
    # Following choices are possible for the flux_discretization scheme: ScharfetterGummel, ScharfetterGummel_Graded,
    # excessChemicalPotential, excessChemicalPotential_Graded, diffusionEnhanced, generalized_SG
    data.flux_approximation             = excessChemicalPotential
    
    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define ChargeTransportParams and fill in physical parameters")
    end
    ################################################################################

    # Params is a struct which contains all necessary physical parameters. If one wants to simulate
    # space-dependent variable, one additionally needs to generate a ParamsNodal struct, see Example_103.
    params                                              = ChargeTransportParams(grid, numberOfCarriers)

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
        params.recombinationRadiative[ireg]             = Radiative
        params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
        params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

    end

    # interior doping
    params.doping[iphin, regionDonor]                   = Nd        # data.doping   = [0.0  Na;
    params.doping[iphin, regionIntrinsic]               = ni        #                  ni   0.0;
    params.doping[iphip, regionIntrinsic]               = 0.0        #                  Nd  0.0]
    params.doping[iphip, regionAcceptor]                = Na

    # boundary doping
    params.bDoping[iphin, bregionDonor]                 = Nd        # data.bDoping  = [0.0  Na;
    params.bDoping[iphip, bregionAcceptor]              = Na        #                  Nd  0.0]

    # Region dependent params is now a substruct of data which is again a substruct of the system and will be parsed 
    # in next step.
    data.params                                         = params

    # in the last step, we initialize our system with previous data which is likewise dependent on the parameters. 
    # important that this is in the end, otherwise our VoronoiFVMSys is not dependent on the data we initialized
    # but rather on default data.
    ctsys                                               = ChargeTransportSystem(grid, data, unknown_storage=unknown_storage)

    if test == false
        # show region dependent physical parameters. show_params() only supports region dependent parameters, but, if one wishes to
        # print nodal dependent parameters, currently this is possible with println(ctsys.data.paramsnodal). We neglected here, since
        # in most applications where the numberOfNodes is >> 10 this would results in a large output in the terminal.
        show_params(ctsys)
        println("*** done\n")
    end

    if plotting == true
        ################################################################################
        println("Plot electroneutral potential, band-edge energies and doping")
        ################################################################################
        psi0 = electroNeutralSolution!(grid, data)
        plot_energies(Plotter, grid, data)
        Plotter.figure()
        plot_doping(Plotter, grid, data)
        Plotter.figure()
        plot_electroNeutralSolutionBoltzmann(Plotter, grid, psi0, ;plotGridpoints=true)
        Plotter.figure()
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define outerior boundary conditions and enabled layers")
    end
    ################################################################################

    # set ohmic contacts for each charge carrier at all outerior boundaries. First, 
    # we compute equilibrium solutions. Hence the boundary values at the ohmic contacts
    # are zero.
    set_ohmic_contact!(ctsys, iphin, bregionAcceptor, 0.0)
    set_ohmic_contact!(ctsys, iphip, bregionAcceptor, 0.0)
    set_ohmic_contact!(ctsys, iphin, bregionDonor, 0.0)
    set_ohmic_contact!(ctsys, iphip, bregionDonor, 0.0)

    # enable all three species in all regions
    # entweder ct_enable_species! oder ChargeTransportInSolids.enable_species!
    ct_enable_species!(ctsys, ipsi,  regions)
    ct_enable_species!(ctsys, iphin, regions)
    ct_enable_species!(ctsys, iphip, regions)

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control                   = VoronoiFVM.NewtonControl()
    control.verbose           = verbose
    control.damp_initial      = 0.5
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
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    control.damp_initial  = 0.5
    control.damp_growth   = 1.2 # >= 1
    control.max_round     = 3

    initialGuess          = ct_unknowns(ctsys)
    solution              = ct_unknowns(ctsys)

    solution              = ct_equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

    initialGuess         .= solution 


    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    ctsys.data.calculation_type      = outOfEquilibrium

    if !(data.F == Boltzmann) # adjust control, when not using Boltzmann
        control.damp_initial      = 0.5
        control.damp_growth       = 1.2
        control.max_iterations    = 30
    end

    maxBias    = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 16)
    IV         = zeros(0)

    w_device = 0.5 * μm     # width of device
    z_device = 1.0e-4 * cm  # depth of device

    for Δu in biasValues

        # set non equilibrium boundary conditions
        set_ohmic_contact!(ctsys, iphin, bregionAcceptor, Δu)
        set_ohmic_contact!(ctsys, iphip, bregionAcceptor, Δu)

        ct_solve!(solution, initialGuess, ctsys, control = control, tstep = Inf)

        initialGuess .= solution

        # get IV curve
        factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)

        # testfunction zero in bregionAcceptor and one in bregionDonor
        tf     = testfunction(factory, [bregionAcceptor], [bregionDonor])
        I      = integrate(ctsys.fvmsys, tf, solution)

        push!(IV,  abs.(w_device * z_device * (I[iphin] + I[iphip])))

    end # bias loop

    if test == false
        println("*** done\n")
    end

    # plot solution and IV curve
    if verbose
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

    if test == false
        println("*** done\n")
    end

end #  main

function test()
    testval = 1.5068426773059806
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when the PIN module has successfully recompiled.")
end

end # module
