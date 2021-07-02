# 101: 1D GaAs p-i-n diode.
([source code](https://github.com/PatricioFarrell/ChargeTransportInSolids.jl/tree/master/examplesExample101_PIN.jl))

Simulating charge transport in a GaAs pin diode. This means
the corresponding PDE problem corresponds to the van Roosbroeck
system of equations.
The simulations are performed out of equilibrium and for the
stationary problem.

```julia
module Example101_PIN

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize
```

function for initializing the grid for a possble extension to other p-i-n devices.

```julia
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
```

region numbers

```julia
    regionAcceptor   = 1                           # p doped region
    regionIntrinsic  = 2                           # intrinsic region
    regionDonor      = 3                           # n doped region
    regions          = [regionAcceptor, regionIntrinsic, regionDonor]
```

boundary region numbers

```julia
    bregionAcceptor  = 1
    bregionDonor     = 2
    bregions         = [bregionAcceptor, bregionDonor]
```

grid

```julia
    refinementfactor = 2^(n-1)
    h_pdoping        = 2 * μm
    h_intrinsic      = 2 * μm
    h_ndoping        = 2 * μm
    coord            = initialize_pin_grid(refinementfactor,
                       h_pdoping,
                       h_intrinsic,
                       h_ndoping)

    grid             = ExtendableGrids.simplexgrid(coord)
    numberOfNodes    = length(coord)
```

set different regions in grid, doping profiles do not intersect

```julia
    cellmask!(grid, [0.0 * μm], [h_pdoping], regionAcceptor)        # p-doped region = 1
    cellmask!(grid, [h_pdoping], [h_pdoping + h_intrinsic], regionIntrinsic)    # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor)     # n-doped region = 3

    if plotting
        GridVisualize.gridplot(grid, Plotter = Plotter)
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
```

indices

```julia
    iphin                   = 1
    iphip                   = 2
    ipsi                    = 3
    species                 = [iphin, iphip, ipsi]
```

number of (boundary) regions and carriers

```julia
    numberOfRegions         = length(regions)
    numberOfBoundaryRegions = length(bregions)
    numberOfCarriers        = length(species) - 1
```

physical data

```julia
    Ec                = 1.424                *  eV
    Ev                = 0.0                  *  eV
    Nc                = 4.351959895879690e17 / (cm^3)
    Nv                = 9.139615903601645e18 / (cm^3)
    mun               = 8500.0               * (cm^2) / (V * s)
    mup               = 400.0                * (cm^2) / (V * s)
    εr                = 12.9                 *  1.0              # relative dielectric permittivity of GAs
    T                 = 300.0                *  K
```

recombination model

```julia
    recombinationOn = true
```

recombination parameters

```julia
    Auger             = 1.0e-29   * cm^6 / s          # 1.0e-41
    SRH_TrapDensity   = 1.0e10    / cm^3              # 1.0e16
    SRH_LifeTime      = 1.0       * ns                # 1.0e10
    Radiative         = 1.0e-10   * cm^3 / s          # 1.0e-16
```

doping

```julia
    dopingFactorNd    =   1.0
    dopingFactorNa    =   0.46
    Nd                =   dopingFactorNd * Nc
    Na                =   dopingFactorNa * Nv
```

intrinsic concentration (not doping!)

```julia
    ni                =   sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * kB * T))
```

contact voltages: we impose an applied voltage only on one boundary.
At the other boundary the applied voltage is zero.

```julia
    voltageAcceptor   = 1.5 * V

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define ChargeTransport data and fill in previously defined data")
    end
    ################################################################################
```

initialize ChargeTransport instance

```julia
    data      = ChargeTransportInSolids.ChargeTransportData(numberOfNodes,
                                                            numberOfRegions,
                                                            numberOfBoundaryRegions,
                                                            numberOfSpecies = numberOfCarriers + 1)
```

region independent data

```julia
    data.F                                           .= Boltzmann # Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.temperature                                  = T
    data.UT                                           = (kB * data.temperature) / q
    data.chargeNumbers[iphin]                         = -1
    data.chargeNumbers[iphip]                         =  1

    data.recombinationOn                              = recombinationOn
```

boundary region data

```julia
    for ibreg in 1:numberOfBoundaryRegions

        data.bDensityOfStates[iphin, ibreg]           = Nc
        data.bDensityOfStates[iphip, ibreg]           = Nv
        data.bBandEdgeEnergy[iphin, ibreg]            = Ec
        data.bBandEdgeEnergy[iphip, ibreg]            = Ev
    end
```

interior region data

```julia
    for ireg in 1:numberOfRegions

        data.dielectricConstant[ireg]                 = εr
```

dos, band edge energy and mobilities

```julia
        data.densityOfStates[iphin, ireg]             = Nc
        data.densityOfStates[iphip, ireg]             = Nv
        data.bandEdgeEnergy[iphin, ireg]              = Ec
        data.bandEdgeEnergy[iphip, ireg]              = Ev
        data.mobility[iphin, ireg]                    = mun
        data.mobility[iphip, ireg]                    = mup
```

recombination parameters

```julia
        data.recombinationRadiative[ireg]             = Radiative
        data.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
        data.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
        data.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
        data.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity
        data.recombinationAuger[iphin, ireg]          = Auger
        data.recombinationAuger[iphip, ireg]          = Auger

    end
```

interior doping

```julia
    data.doping[iphin, regionDonor]                   = Nd        # data.doping   = [0.0  Na;
    data.doping[iphin, regionIntrinsic]               = ni        #                  ni   ni;
    data.doping[iphip, regionIntrinsic]               = 0.0        #                  Nd  0.0]
    data.doping[iphip, regionAcceptor]                = Na
```

boundary doping

```julia
    data.bDoping[iphin, bregionDonor]                 = Nd        # data.bDoping  = [0.0  Na;
    data.bDoping[iphip, bregionAcceptor]              = Na        #                  Nd  0.0]

    if test == false
```

print data

```julia
        println(data)
        println("*** done\n")
    end

    if plotting == true
        ################################################################################
        println("Plot electroneutral potential, band-edge energies and doping")
        ################################################################################
        psi0 = ChargeTransportInSolids.electroNeutralSolution!(data, grid)
        ChargeTransportInSolids.plotEnergies(Plotter, grid, data)
        Plotter.figure()
        ChargeTransportInSolids.plotDoping(Plotter, grid, data)
        Plotter.figure()
        ChargeTransportInSolids.plotElectroNeutralSolutionBoltzmann(Plotter, grid, psi0, ;plotGridpoints=true)
        Plotter.figure()
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define physics and system")
    end
    ################################################################################

    # initializing physics environment ##
    physics = VoronoiFVM.Physics(
    data        = data,
    num_species = numberOfCarriers + 1,
    flux        = ChargeTransportInSolids.Sedan!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
    reaction    = ChargeTransportInSolids.reaction!,
    breaction   = ChargeTransportInSolids.breactionOhmic!
    )

    sys         = VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
```

enable all three species in all regions

```julia
    enable_species!(sys, ipsi,  regions)
    enable_species!(sys, iphin, regions)
    enable_species!(sys, iphip, regions)

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control = VoronoiFVM.NewtonControl()
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

    data.inEquilibrium             = true
```

initialize solution and starting vectors

```julia
    initialGuess                   = unknowns(sys)
    solution                       = unknowns(sys)
    @views initialGuess[ipsi,  :] .= 0.0
    @views initialGuess[iphin, :] .= 0.0
    @views initialGuess[iphip, :] .= 0.0

    function pre(u,lambda)
        sys.physics.data.λ1                         = lambda
        sys.boundary_values[iphin, bregionAcceptor] = 0.0
        sys.boundary_values[iphip, bregionAcceptor] = 0.0
    end

    control.damp_initial      = 0.01
    control.damp_growth       = 1.2 # >= 1
    control.max_round         = 3
```

set Dirichlet boundary conditions (Ohmic contacts), in Equilibrium we impose homogeneous Dirichlet conditions,
i.e. the boundary values at outer boundaries are zero.

```julia
    sys.boundary_factors[iphin, bregionDonor]    = VoronoiFVM.Dirichlet
    sys.boundary_values[iphin,  bregionDonor]    = 0.0 * V
    sys.boundary_factors[iphin, bregionAcceptor] = VoronoiFVM.Dirichlet
    sys.boundary_values[iphin,  bregionAcceptor] = 0.0 * V

    sys.boundary_factors[iphip, bregionDonor]    = VoronoiFVM.Dirichlet
    sys.boundary_values[iphip,  bregionDonor]    = 0.0 * V
    sys.boundary_factors[iphip, bregionAcceptor] = VoronoiFVM.Dirichlet
    sys.boundary_values[iphip,  bregionAcceptor] = 0.0 * V

    I = collect(20.0:-1:0.0)
    LAMBDA = 10 .^ (-I)
    prepend!(LAMBDA,0.0)

    for i in 1:length(LAMBDA)
        if test == false
            println("λ1 = $(LAMBDA[i])")
        end
        sys.physics.data.λ1 = LAMBDA[i]
        solve!(solution, initialGuess, sys, control = control, tstep=Inf)
        initialGuess = solution
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    data.inEquilibrium = false

    if !(data.F == ChargeTransportInSolids.Boltzmann) # adjust control, when not using Boltzmann
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
```

set non equilibrium boundary conditions

```julia
        sys.boundary_values[iphin, bregionAcceptor] = Δu
        sys.boundary_values[iphip, bregionAcceptor] = Δu

        solve!(solution, initialGuess, sys, control = control, tstep = Inf)

        initialGuess .= solution
```

get IV curve

```julia
        factory = VoronoiFVM.TestFunctionFactory(sys)
```

testfunction zero in bregionAcceptor and one in bregionDonor

```julia
        tf     = testfunction(factory, [bregionAcceptor], [bregionDonor])
        I      = integrate(sys, tf, solution)

        push!(IV,  abs.(w_device * z_device * (I[iphin] + I[iphip])))

    end # bias loop
```

plot solution and IV curve

```julia
    if plotting
        ChargeTransportInSolids.plotEnergies(Plotter, grid, data, solution, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = false)
        Plotter.figure()
        ChargeTransportInSolids.plotSolution(Plotter, coord, solution, 0.0, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = true)
        Plotter.figure()
        ChargeTransportInSolids.plotDensities(Plotter, grid, data, solution, "Applied voltage Δu = $(biasValues[end])", plotGridpoints = true)
        Plotter.figure()
        ChargeTransportInSolids.plotIV(Plotter, biasValues,IV, biasValues[end], plotGridpoints = true)
    end
    testval = solution[15]
    return testval

    if test == false
        println("*** done\n")
    end

end #  main

function test()
    testval=1.5068426773059806
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when the PIN module has successfully recompiled.")
end

end # module
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

