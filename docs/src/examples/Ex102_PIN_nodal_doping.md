# GaAs diode with spatially varying doping (1D).
([source code](https://github.com/PatricioFarrell/ChargeTransport.jl/tree/master/examplesEx102_PIN_nodal_doping.jl))

Simulating charge transport in a GaAs pin diode. This means the PDE problem corresponds to the
van Roosbroeck system of equations. The simulations are performed out of equilibrium and for
the stationary problem. A special feature here is that the doping is node-dependent.

````julia
module Ex102_PIN_nodal_doping

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot

function main(;Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

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

    h_pdoping               = 0.1    * μm
    h_intrinsic             = 0.1    * μm
    h_ndoping               = 0.1    * μm
    w_device                = 0.1    * μm  # width of device
    z_device                = 1.0e-5 * cm  # depth of device

    coord                   = range(0.0, stop = h_ndoping + h_intrinsic + h_pdoping, length = 25)
    coord                   = collect(coord)
    grid                    = simplexgrid(coord)
    numberOfNodes           = length(coord)

    # set different regions in grid
    cellmask!(grid, [0.0 * μm],                [h_pdoping],                           regionAcceptor,  tol = 1.0e-15)    # p-doped region = 1
    cellmask!(grid, [h_pdoping],               [h_pdoping + h_intrinsic],             regionIntrinsic, tol = 1.0e-15)    # intrinsic region = 2
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

    # Define the physical data.
    Ec                 = 1.424                *  eV
    Ev                 = 0.0                  *  eV
    Nc                 = 4.351959895879690e17 / (cm^3)
    Nv                 = 9.139615903601645e18 / (cm^3)
    mun                = 8500.0               * (cm^2) / (V * s)
    mup                = 400.0                * (cm^2) / (V * s)
    εr                 = 12.9                 *  1.0              # relative dielectric permittivity of GAs
    T                  = 300.0                *  K

    # recombination parameters
    SRH_TrapDensity_n  = 4.760185435081902e5    / cm^3
    SRH_TrapDensity_p  = 9.996936448738406e6    / cm^3
    SRH_LifeTime       = 1.0                    * ps

    # contact voltage
    voltageAcceptor    = 1.4 * V

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################
````

We initialize the Data instance and fill in predefined data.

````julia
    data                               = Data(grid, numberOfCarriers)

    # Possible choices: Stationary, Transient
    data.modelType                     = Stationary

    # Possible choices for F: Boltzmann, FermiDiracOneHalfBednarczyk,
    # FermiDiracOneHalfTeSCA, FermiDiracMinusOne, Blakemore
    data.F                            .= Boltzmann

    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                bulk_recomb_Auger = false,
                                                                bulk_recomb_radiative = false,
                                                                bulk_recomb_SRH = true)

    # Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceModelNone,
    # InterfaceModelSurfaceReco (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor]    = OhmicContact

    # Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    # ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation            .= ScharfetterGummel

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################
````

Define the Params and ParamsNodal struct.

````julia
    params                                              = Params(grid, numberOfCarriers)
    paramsnodal                                         = ParamsNodal(grid, numberOfCarriers)

    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1

    for ibreg in 1:numberOfBoundaryRegions  # boundary region data

        params.bDensityOfStates[iphin, ibreg]           = Nc
        params.bDensityOfStates[iphip, ibreg]           = Nv
        params.bBandEdgeEnergy[iphin, ibreg]            = Ec
        params.bBandEdgeEnergy[iphip, ibreg]            = Ev
    end

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg]                 = εr * ε0

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

    # initialize the space dependent doping (see FarrellPeschka2019, Computers & Mathematics with Applications, 2019).
    NDoping  = 1.0e17  / cm^3
    κ        = 500.0
    for icoord = 1:numberOfNodes
        paramsnodal.doping[icoord] = NDoping * 0.5 * ( 1.0  +  tanh( (0.1 - coord[icoord]/μm) *κ )  - ( 1.0 + tanh( (coord[icoord]/μm - 0.2) * κ )) )
    end

    data.params      = params
    data.paramsnodal = paramsnodal

    ctsys            = System(grid, data, unknown_storage=unknown_storage)

    if test == false
        println("*** done\n")
    end

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
        println("Define outer boundary conditions")
    end
    ################################################################################

    # set zero voltage ohmic contacts for each charge carrier at all outer boundaries.
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

    # initialize solution and starting vectors
    initialGuess  = unknowns(ctsys)
    solution      = unknowns(ctsys)

    solution      = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

    initialGuess .= solution

    if plotting
        # set legend for plotting routines. Either you can use the predefined labels or write your own.
        label_solution, label_density, label_energy = set_plotting_labels(data)

        Plotter.figure()
        plot_energies(Plotter,  grid, data, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution, "Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter,  grid, data, solution, "Equilibrium", label_solution)
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
````

Set calculationType to OutOfEquilibrium for starting with respective simulation.

````julia
    data.calculationType = OutOfEquilibrium

    maxBias              = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    biasValues           = range(0, stop = maxBias, length = 41)
    IV                   = zeros(0)

    for Δu in biasValues

        # set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

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
        plot_energies(Plotter, grid, data, solution, "Applied voltage Δu = $(biasValues[end])",  label_energy)
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "Applied voltage Δu = $(biasValues[end])",  label_solution, plotGridpoints = true)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution, "Applied voltage Δu = $(biasValues[end])", label_density,  plotGridpoints = true)
        Plotter.figure()
        plot_IV(Plotter, biasValues,IV, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = true)
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
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

