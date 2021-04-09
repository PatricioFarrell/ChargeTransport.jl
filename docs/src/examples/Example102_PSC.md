# 102: 1D PSC p-i-n device without mobile ions.
([source code](https://github.com/PatricioFarrell/ChargeTransportInSolids.jl/tree/master/examplesExample102_PSC.jl))

Simulating a three layer PSC device without mobile ions.
The simulations are performed out of equilibrium and with
abrupt interfaces.

This simulation coincides with the one made in Section 4.3
of Calado et al. (https://arxiv.org/abs/2009.04384).
The paramters can be found here:
https://github.com/barnesgroupICL/Driftfusion/blob/Methods-IonMonger-comparison/Input_files/IonMonger_default_noIR.csv.

```julia
module Example102_PSC

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using Printf
using GridVisualize

function main(;n = 8, Plotter = nothing, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################
```

region numbers

```julia
    regionDonor     = 1                           # n doped region
    regionIntrinsic = 2                           # intrinsic region
    regionAcceptor  = 3                           # p doped region
    regions         = [regionDonor, regionIntrinsic, regionAcceptor]
```

boundary region numbers

```julia
    bregionDonor    = 1
    bregionAcceptor = 2
    bregions        = [bregionAcceptor, bregionDonor]
```

grid
NB: Using geomspace to create uniform mesh is not a good idea. It may create virtual duplicates at boundaries.

```julia
    h_ndoping       = 9.90e-6 * cm
    h_intrinsic     = 4.00e-5 * cm + 2.0e-7 * cm # add 2.e-7 cm to this layer for agreement with grid of Driftfusion
    h_pdoping       = 1.99e-5 * cm

    x0              = 0.0 * cm
    δ               = 2*n        # the larger, the finer the mesh
    t               = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k               = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_n_u       = collect(range(x0, h_ndoping/2, step=h_ndoping/(0.8*δ)))
    coord_n_g       = geomspace(h_ndoping/2,
                                h_ndoping,
                                h_ndoping/(1.1*δ),
                                h_ndoping/(1.1*δ),
                                tol=t)
    coord_i_g1      = geomspace(h_ndoping,
                                h_ndoping+h_intrinsic/k,
                                h_intrinsic/(2.8*δ),
                                h_intrinsic/(2.8*δ),
                                tol=t)
    coord_i_g2      = geomspace(h_ndoping+h_intrinsic/k,
                                h_ndoping+h_intrinsic,
                                h_intrinsic/(2.8*δ),
                                h_intrinsic/(2.8*δ),
                                tol=t)
    coord_p_g       = geomspace(h_ndoping+h_intrinsic,
                                h_ndoping+h_intrinsic+h_pdoping/2,
                                h_pdoping/(1.6*δ),
                                h_pdoping/(1.6*δ),
                                tol=t)
    coord_p_u       = collect(range(h_ndoping+h_intrinsic+h_pdoping/2, h_ndoping+h_intrinsic+h_pdoping, step=h_pdoping/(1.3*δ)))

    coord           = glue(coord_n_u,coord_n_g,  tol=10*t)
    coord           = glue(coord,    coord_i_g1, tol=10*t)
    coord           = glue(coord,    coord_i_g2, tol=10*t)
    coord           = glue(coord,    coord_p_g,  tol=10*t)
    coord           = glue(coord,    coord_p_u,  tol=10*t)
    grid            = ExtendableGrids.simplexgrid(coord)
    numberOfNodes   = length(coord)
```

set different regions in grid, doping profiles do not intersect

```julia
    cellmask!(grid, [0.0 * μm],                [h_ndoping],                           regionDonor)     # n-doped region   = 1
    cellmask!(grid, [h_ndoping],               [h_ndoping + h_intrinsic],             regionIntrinsic) # intrinsic region = 2
    cellmask!(grid, [h_ndoping + h_intrinsic], [h_ndoping + h_intrinsic + h_pdoping], regionAcceptor)  # p-doped region   = 3

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
    iphin, iphip, ipsi      = 1:3
    species                 = [iphin, iphip, ipsi]
```

number of (boundary) regions and carriers

```julia
    numberOfRegions         = length(regions)
    numberOfBoundaryRegions = length(bregions) #+ length(iregions)
    numberOfCarriers        = length(species) - 1
```

temperature

```julia
    T               = 300.0                 *  K
```

band edge energies

```julia
    Eref            = 0.0        # reference energy

    Ec_d            = -4.0                  *  eV
    Ev_d            = -6.0                  *  eV

    Ec_i            = -3.7                  *  eV
    Ev_i            = -5.4                  *  eV

    Ec_a            = -3.1                  *  eV
    Ev_a            = -5.1                  *  eV

    EC              = [Ec_d, Ec_i, Ec_a]
    EV              = [Ev_d, Ev_i, Ev_a]
```

effective densities of state

```julia
    Nc_d            = 5.0e19                / (cm^3)
    Nv_d            = 5.0e19                / (cm^3)

    Nc_i            = 8.1e18                / (cm^3)
    Nv_i            = 5.8e18                / (cm^3)

    Nc_a            = 5.0e19                / (cm^3)
    Nv_a            = 5.0e19                / (cm^3)

    NC              = [Nc_d, Nc_i, Nc_a]
    NV              = [Nv_d, Nv_i, Nv_a]
```

mobilities

```julia
    μn_d            = 3.89                  * (cm^2) / (V * s)
    μp_d            = 3.89                  * (cm^2) / (V * s)

    μn_i            = 6.62e1                * (cm^2) / (V * s)
    μp_i            = 6.62e1                * (cm^2) / (V * s)

    μn_a            = 3.89e-1               * (cm^2) / (V * s)
    μp_a            = 3.89e-1               * (cm^2) / (V * s)

    μn              = [μn_d, μn_i, μn_a]
    μp              = [μp_d, μp_i, μp_a]
```

relative dielectric permittivity

```julia
    ε_d             = 10.0                  *  1.0
    ε_i             = 24.1                  *  1.0
    ε_a             = 3.0                   *  1.0

    ε               = [ε_d, ε_i, ε_a]
```

recombination model

```julia
    recombinationOn = true
```

radiative recombination

```julia
    r0_d            = 0.0e+0               * cm^3 / s
    r0_i            = 1.0e-12              * cm^3 / s
    r0_a            = 0.0e+0               * cm^3 / s

    r0              = [r0_d, r0_i, r0_a]
```

life times and trap densities

```julia
    τn_d            = 1.0e100              * s
    τp_d            = 1.0e100              * s

    τn_i            = 3.0e-10              * s
    τp_i            = 3.0e-8               * s
    τn_a            = τn_d
    τp_a            = τp_d

    τn              = [τn_d, τn_i, τn_a]
    τp              = [τp_d, τp_i, τp_a]
```

SRH trap energies (needed for calculation of recombinationSRHTrapDensity)

```julia
    Ei_d            = -5.0                 * eV
    Ei_i            = -4.55                * eV
    Ei_a            = -4.1                 * eV

    EI              = [Ei_d, Ei_i, Ei_a]
```

Auger recombination

```julia
    Auger           = 0.0
```

doping (doping values are from Phils paper, not stated in the parameter list online)

```julia
    Nd              =   1.03e18             / (cm^3)
    Na              =   1.03e18             / (cm^3)
    Ni_acceptor     =   8.32e7              / (cm^3)
```

contact voltages

```julia
    voltageAcceptor =  1.2                  * V
    voltageDonor    =  0.0                  * V

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
    data                                 = ChargeTransportInSolids.ChargeTransportData(numberOfNodes,
                                                                                       numberOfRegions,
                                                                                       numberOfBoundaryRegions,
                                                                                       numberOfSpecies = numberOfCarriers + 1)
```

region independent data

```julia
    data.F                                       .= Boltzmann # Boltzmann, FermiDiracMinusOne,FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA, Blakemore
    data.temperature                              = T
    data.UT                                       = (kB * data.temperature) / q
    data.contactVoltage[bregionAcceptor]          = voltageAcceptor
    data.contactVoltage[bregionDonor]             = voltageDonor
    data.chargeNumbers[iphin]                     = -1
    data.chargeNumbers[iphip]                     =  1
    data.Eref                                     =  Eref

    data.recombinationOn                          = recombinationOn
```

boundary region data

```julia
    data.bDensityOfStates[iphin, bregionDonor]    = Nc_d
    data.bDensityOfStates[iphip, bregionDonor]    = Nv_d

    data.bDensityOfStates[iphin, bregionAcceptor] = Nc_a
    data.bDensityOfStates[iphip, bregionAcceptor] = Nv_a

    data.bBandEdgeEnergy[iphin, bregionDonor]     = Ec_d + data.Eref
    data.bBandEdgeEnergy[iphip, bregionDonor]     = Ev_d + data.Eref

    data.bBandEdgeEnergy[iphin, bregionAcceptor]  = Ec_a + data.Eref
    data.bBandEdgeEnergy[iphip, bregionAcceptor]  = Ev_a + data.Eref
```

interior region data

```julia
    for ireg in 1:numberOfRegions

        data.dielectricConstant[ireg]                 = ε[ireg]
```

dos, band edge energy and mobilities

```julia
        data.densityOfStates[iphin, ireg]             = NC[ireg]
        data.densityOfStates[iphip, ireg]             = NV[ireg]

        data.bandEdgeEnergy[iphin, ireg]              = EC[ireg] + data.Eref
        data.bandEdgeEnergy[iphip, ireg]              = EV[ireg] + data.Eref

        data.mobility[iphin, ireg]                    = μn[ireg]
        data.mobility[iphip, ireg]                    = μp[ireg]
```

recombination parameters

```julia
        data.recombinationRadiative[ireg]             = r0[ireg]
        data.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        data.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        data.recombinationSRHTrapDensity[iphin, ireg] = ChargeTransportInSolids.trapDensity(iphin, ireg, data, EI[ireg])
        data.recombinationSRHTrapDensity[iphip, ireg] = ChargeTransportInSolids.trapDensity(iphip, ireg, data, EI[ireg])
        data.recombinationAuger[iphin, ireg]          = Auger
        data.recombinationAuger[iphip, ireg]          = Auger
    end
```

interior doping

```julia
    data.doping[iphin, regionDonor]               = Nd
    data.doping[iphip, regionIntrinsic]           = Ni_acceptor
    data.doping[iphip, regionAcceptor]            = Na
```

boundary doping

```julia
    data.bDoping[iphip, bregionAcceptor]          = Na        # data.bDoping  = [Na  0.0;
    data.bDoping[iphin, bregionDonor]             = Nd        #                  0.0  Nd]
```

print data

```julia
    if test == false
        println(data)
    end

    if test == false
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
    flux        = ChargeTransportInSolids.ScharfetterGummel!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
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

    sys.boundary_values[iphin,  bregionAcceptor] = data.contactVoltage[bregionAcceptor]
    sys.boundary_factors[iphin, bregionAcceptor] = VoronoiFVM.Dirichlet

    sys.boundary_values[iphin,  bregionDonor]    = data.contactVoltage[bregionDonor]
    sys.boundary_factors[iphin, bregionDonor]    = VoronoiFVM.Dirichlet

    sys.boundary_values[iphip,  bregionAcceptor] = data.contactVoltage[bregionAcceptor]
    sys.boundary_factors[iphip, bregionAcceptor] = VoronoiFVM.Dirichlet

    sys.boundary_values[iphip,  bregionDonor]    = data.contactVoltage[bregionDonor]
    sys.boundary_factors[iphip, bregionDonor]    = VoronoiFVM.Dirichlet

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
    control.max_iterations    = 300
    control.tol_absolute      = 1.0e-13
    control.tol_relative      = 1.0e-13
    control.handle_exceptions = true
    control.tol_round         = 1.0e-13
    control.max_round         = 5

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium for Boltzmann")
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

    control.damp_initial      = 0.1
    control.damp_growth       = 1.61 # >= 1
    control.max_round         = 5

    sys.boundary_values[iphin, bregionAcceptor] = 0.0 * V
    sys.boundary_values[iphip, bregionAcceptor] = 0.0 * V
    sys.physics.data.contactVoltage             = 0.0 * sys.physics.data.contactVoltage

    I = collect(20.0:-1:0.0)
    LAMBDA = 10 .^ (-I)
    prepend!(LAMBDA,0.0)
    for i in 1:length(LAMBDA)
        if test == false
            println("λ1 = $(LAMBDA[i])")
        end
        sys.physics.data.λ1 = LAMBDA[i]
        solve!(solution, initialGuess, sys, control = control, tstep=Inf)
        initialGuess .= solution
    end

    if plotting
```

ChargeTransportInSolids.plotEnergies(Plotter, grid, data, solution, "EQULIBRIUM (NO illumination)")
Plotter.figure()
ChargeTransportInSolids.plotDensities(Plotter, grid, data, solution, "EQULIBRIUM (NO illumination)")
Plotter.figure()

```julia
        ChargeTransportInSolids.plotSolution(Plotter, coord, solution, data.Eref, "EQULIBRIUM (NO illumination)")
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
    data.inEquilibrium = false

    control.damp_initial                             = 0.005
    control.damp_growth                              = 1.21 # >= 1
    control.max_round                                = 7
```

set non equilibrium boundary conditions

```julia
    sys.physics.data.contactVoltage[bregionDonor]    = voltageDonor
    sys.physics.data.contactVoltage[bregionAcceptor] = voltageAcceptor
    sys.boundary_values[iphin, bregionAcceptor]      = data.contactVoltage[bregionAcceptor]
    sys.boundary_values[iphip, bregionAcceptor]      = data.contactVoltage[bregionAcceptor]

    maxBias    = data.contactVoltage[bregionAcceptor]
    biasValues = range(0, stop = maxBias, length = 13)

    for Δu in biasValues
        if test == false
            println("Bias value: Δu = $(Δu) (no illumination)")
        end

        data.contactVoltage[bregionAcceptor]         = Δu
        sys.boundary_values[iphin, bregionAcceptor]  = Δu
        sys.boundary_values[iphip, bregionAcceptor]  = Δu

        solve!(solution, initialGuess, sys, control  = control, tstep = Inf)

        initialGuess .= solution

    end # bias loop

    #plotting
    if plotting
```

ChargeTransportInSolids.plotEnergies(Plotter, grid, data, solution, "Applied voltage Δu = $maxBias")
Plotter.figure()
ChargeTransportInSolids.plotDensities(Plotter, grid, data, solution, "Applied voltage Δu = $maxBias")
Plotter.figure()

```julia
        ChargeTransportInSolids.plotSolution(Plotter, coord, solution, data.Eref, "Applied voltage Δu = $maxBias")
    end
    if test == false
        println("*** done\n")
    end

    testval = solution[ipsi, 20]
    return testval

end #  main

function test()
    testval=-4.052650626421281
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

