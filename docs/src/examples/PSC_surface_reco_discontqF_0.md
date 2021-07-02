# xxx
([source code](https://github.com/PatricioFarrell/ChargeTransportInSolids.jl/tree/master/examplesPSC_surface_reco_discontqF_0.jl))

Simulating a three layer PSC device Pedot| MAPI | PCBM.

The paramters can be found here and are from
Calado et al.:
https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/pedotpss_mapi_pcbm.csv.
(with adjustments on layer lengths)

here: discont qF, i.e. three times phin and phip (on each subdomain one)

```julia
module PSC_surface_reco_discontqF_0

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize
using PyPlot
using DelimitedFiles

function main(;n = 6, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:dense)

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
```

inner boundary regions

```julia
    bregionJunction1 = 3
    bregionJunction2 = 4
    bregions         = [bregionAcceptor, bregionDonor, bregionJunction1, bregionJunction2]
    innerInterfaces  = true # put here true, when wanting additional effects on inner inerfaces.
```

grid
NB: Using geomspace to create uniform mesh is not a good idea. It may create virtual duplicates at boundaries.

```julia
    h_pdoping       = 3.00e-6 * cm + 1.0e-7 *cm# add 1.e-7 cm to this layer for agreement with grid of Driftfusion
    h_intrinsic     = 3.00e-5 * cm
    h_ndoping       = 8.50e-6 * cm + 1.0e-7 *cm# add 1.e-7 cm to this layer for agreement with grid of Driftfusion

    x0              = 0.0 * cm
    δ               = 4*n        # the larger, the finer the mesh
    t               = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k               = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_p_u       = collect(range(x0, h_pdoping/2, step=h_pdoping/(1.0*δ)))
    coord_p_g       = geomspace(h_pdoping/2,
                                h_pdoping,
                                h_pdoping/(1.2*δ),
                                h_pdoping/(0.6*δ),
                                tol=t)
    coord_i_g1      = geomspace(h_pdoping,
                                h_pdoping+h_intrinsic/k,
                                h_intrinsic/(7.1*δ),
                                h_intrinsic/(6.5*δ),
                                tol=t)
    coord_i_g2      = geomspace(h_pdoping+h_intrinsic/k,
                                h_pdoping+h_intrinsic,
                                h_intrinsic/(6.5*δ),
                                h_intrinsic/(7.1*δ),
                                tol=t)
    coord_n_g       = geomspace(h_pdoping+h_intrinsic,
                                h_pdoping+h_intrinsic+h_ndoping/2,
                                h_ndoping/(1.6*δ),
                                h_ndoping/(1.5*δ),
                                tol=t)
    coord_n_u       = collect(range(h_pdoping+h_intrinsic+h_ndoping/2, h_pdoping+h_intrinsic+h_ndoping, step=h_pdoping/(0.6*δ)))

    coord           = glue(coord_p_u,coord_p_g,  tol=10*t)
    coord           = glue(coord,    coord_i_g1, tol=10*t)
    coord           = glue(coord,    coord_i_g2, tol=10*t)
    coord           = glue(coord,    coord_n_g,  tol=10*t)
    coord           = glue(coord,    coord_n_u,  tol=10*t)
    grid            = ExtendableGrids.simplexgrid(coord)
    numberOfNodes   = length(coord)
```

set different regions in grid, doping profiles do not intersect

```julia
    cellmask!(grid, [0.0 * μm],                [h_pdoping],                           regionAcceptor, tol = 1.0e-15)     # n-doped region   = 1
    cellmask!(grid, [h_pdoping],               [h_pdoping + h_intrinsic],             regionIntrinsic, tol = 1.0e-15) # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor, tol = 1.0e-15)  # p-doped region   = 3
```

inner boundary regions

```julia
    bfacemask!(grid, [h_pdoping],               [h_pdoping],               bregionJunction1, tol = 1.0e-15)
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic], bregionJunction2, tol = 1.0e-15)

    if plotting
        GridVisualize.gridplot(grid, Plotter = Plotter, legend=:lt)
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
```

indices

```julia
    iphin1, iphin2, iphin3   = 1:3
    iphip1, iphip2, iphip3   = 4:6
    iphia                    = 7
    ipsi                     = 8
    speciesBulk              = [iphin1, iphin2, iphin3,
                                iphip1, iphip2, iphip3,
                                iphia, ipsi]
```

number of (boundary) regions and carriers

```julia
    numberOfRegions          = length(regions)
    numberOfBoundaryRegions  = length(bregions)
    numberOfCarriers         = length(speciesBulk) - 1
```

temperature

```julia
    T               = 300.0                 *  K
```

band edge energies

```julia
    Ec_a            = -3.0                  *  eV
    Ev_a            = -5.1                  *  eV

    Ec_i            = -3.8                  *  eV
    Ev_i            = -5.4                  *  eV
    ###################### adjust Na, Ea here #####################
    Nanion          = 1.21e22                / (cm^3)
    Ea_i            = -5.175                 *  eV
```

for the labels in the figures

```julia
    textEa          = Ea_i./eV
    textNa          = Nanion.*cm^3
    ###################### adjust Na, Ea here #####################
    Ec_d            = -3.8                  *  eV
    Ev_d            = -6.2                  *  eV

    EC              = [Ec_a, Ec_i, Ec_d]
    EV              = [Ev_a, Ev_i, Ev_d]
    EA              = [0.0,  Ea_i,  0.0]
```

effective densities of state

```julia
    Nc_a            = 1.0e20                / (cm^3)
    Nv_a            = 1.0e20                / (cm^3)

    Nc_i            = 1.0e19                / (cm^3)
    Nv_i            = 1.0e19                / (cm^3)

    Nc_d            = 1.0e19                / (cm^3)
    Nv_d            = 1.0e19                / (cm^3)

    NC              = [Nc_a, Nc_i, Nc_d]
    NV              = [Nv_a, Nv_i, Nv_d]
    NAnion          = [0.0,  Nanion, 0.0]
```

mobilities

```julia
    μn_a            = 0.1                   * (cm^2) / (V * s)
    μp_a            = 0.1                   * (cm^2) / (V * s)

    μn_i            = 2.00e1                * (cm^2) / (V * s)
    μp_i            = 2.00e1                * (cm^2) / (V * s)
    μa_i            = 1.00e-10              * (cm^2) / (V * s)

    μn_d            = 1.0e-3                * (cm^2) / (V * s)
    μp_d            = 1.0e-3                * (cm^2) / (V * s)

    μn              = [μn_a, μn_i, μn_d]
    μp              = [μp_a, μp_i, μp_d]
    μa              = [0.0,  μa_i, 0.0 ]
```

relative dielectric permittivity

```julia
    ε_a             = 4.0                   *  1.0
    ε_i             = 23.0                  *  1.0
    ε_d             = 3.0                   *  1.0

    ε               = [ε_a, ε_i, ε_d]
```

recombination model

```julia
    recombinationOn = true
```

radiative recombination

```julia
    r0_a            = 6.3e-11               * cm^3 / s
    r0_i            = 3.6e-12               * cm^3 / s
    r0_d            = 6.8e-11               * cm^3 / s

    r0              = [r0_a, r0_i, r0_d]
```

life times and trap densities

```julia
    τn_a            = 1.0e-6              * s
    τp_a            = 1.0e-6              * s

    τn_i            = 1.0e-7              * s
    τp_i            = 1.0e-7              * s
    τn_d            = τn_a
    τp_d            = τp_a

    τn              = [τn_a, τn_i, τn_d]
    τp              = [τp_a, τp_i, τp_d]
```

SRH trap energies (needed for calculation of recombinationSRHTrapDensity)

```julia
    Ei_a            = -4.05              * eV
    Ei_i            = -4.60              * eV
    Ei_d            = -5.00              * eV

    EI              = [Ei_a, Ei_i, Ei_d]
```

Auger recombination

```julia
    Auger           = 0.0
```

doping (doping values are from Phils paper, not stated in the parameter list online)

```julia
    Nd              =   2.089649130192123e17 / (cm^3)
    Na              =   4.529587947185444e18 / (cm^3)
    C0              =   1.0e18               / (cm^3)
```

contact voltages: we impose an applied voltage only on one boundary.
At the other boundary the applied voltage is zero.

```julia
    voltageAcceptor =  1.2                  * V

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
                                                                                       discontqF = true,
                                                                                       numberOfSpecies = numberOfCarriers + 1
                                                                                       ) # we need to add ipsi
```

region independent data

```julia
    #speciesBulk              = [iphin1, iphin2, iphin3, iphip1, iphip2, iphip3, iphia, ipsi]
    data.F                                        = [Boltzmann, Boltzmann, Boltzmann,
                                                     Boltzmann, Boltzmann, Boltzmann,
                                                     FermiDiracMinusOne]
    data.temperature                              = T
    data.UT                                       = (kB * data.temperature) / q
    data.chargeNumbers[iphin1]                    = -1
    data.chargeNumbers[iphin2]                    = -1
    data.chargeNumbers[iphin3]                    = -1
    data.chargeNumbers[iphip1]                    =  1
    data.chargeNumbers[iphip2]                    =  1
    data.chargeNumbers[iphip3]                    =  1
    data.chargeNumbers[iphia]                     =  1
    data.recombinationOn                          = recombinationOn
    data.innerInterfaces                          = false # innerInterfaces
```

boundary region data

```julia
    data.bDensityOfStates[iphin1, bregionAcceptor]   = Nc_a
    data.bDensityOfStates[iphip1, bregionAcceptor]   = Nv_a

    data.bBandEdgeEnergy[iphin1, bregionAcceptor]    = Ec_a
    data.bBandEdgeEnergy[iphip1, bregionAcceptor]    = Ev_a
    #####
    data.bDensityOfStates[iphin3, bregionDonor]      = Nc_d
    data.bDensityOfStates[iphip3, bregionDonor]      = Nv_d

    data.bBandEdgeEnergy[iphin3, bregionDonor]       = Ec_d
    data.bBandEdgeEnergy[iphip3, bregionDonor]       = Ev_d
```

HTL

```julia
    data.densityOfStates[iphin1, regionAcceptor]              = NC[regionAcceptor]
    data.densityOfStates[iphip1, regionAcceptor]              = NV[regionAcceptor]
    data.bandEdgeEnergy[iphin1, regionAcceptor]               = EC[regionAcceptor]
    data.bandEdgeEnergy[iphip1, regionAcceptor]               = EV[regionAcceptor]
    data.mobility[iphin1, regionAcceptor]                     = μn[regionAcceptor]
    data.mobility[iphip1, regionAcceptor]                     = μp[regionAcceptor]
    data.recombinationSRHLifetime[iphin1, regionAcceptor]     = τn[regionAcceptor]
    data.recombinationSRHLifetime[iphip1, regionAcceptor]     = τp[regionAcceptor]
    data.recombinationSRHTrapDensity[iphin1, regionAcceptor]  = ChargeTransportInSolids.trapDensity(iphin1, regionAcceptor, data, EI[regionAcceptor])
    data.recombinationSRHTrapDensity[iphip1, regionAcceptor]  = ChargeTransportInSolids.trapDensity(iphip1, regionAcceptor, data, EI[regionAcceptor])
    data.recombinationAuger[iphin1, regionAcceptor]           = Auger
    data.recombinationAuger[iphip1, regionAcceptor]           = Auger
```

active perovskite layer

```julia
    data.densityOfStates[iphin2, regionIntrinsic]             = NC[regionIntrinsic]
    data.densityOfStates[iphip2, regionIntrinsic]             = NV[regionIntrinsic]
    data.densityOfStates[iphia, regionIntrinsic]              = NAnion[regionIntrinsic]
    data.bandEdgeEnergy[iphin2, regionIntrinsic]              = EC[regionIntrinsic]
    data.bandEdgeEnergy[iphip2, regionIntrinsic]              = EV[regionIntrinsic]
    data.bandEdgeEnergy[iphia, regionIntrinsic]               = EA[regionIntrinsic]
    data.mobility[iphin2, regionIntrinsic]                    = μn[regionIntrinsic]
    data.mobility[iphip2, regionIntrinsic]                    = μp[regionIntrinsic]
    data.mobility[iphia, regionIntrinsic]                     = μa[regionIntrinsic]
    data.recombinationSRHLifetime[iphin2, regionIntrinsic]    = τn[regionIntrinsic]
    data.recombinationSRHLifetime[iphip2, regionIntrinsic]    = τp[regionIntrinsic]
    data.recombinationSRHTrapDensity[iphin2, regionIntrinsic] = ChargeTransportInSolids.trapDensity(iphin2, regionIntrinsic, data, EI[regionIntrinsic])
    data.recombinationSRHTrapDensity[iphip2, regionIntrinsic] = ChargeTransportInSolids.trapDensity(iphip2, regionIntrinsic, data, EI[regionIntrinsic])
    data.recombinationAuger[iphin2, regionIntrinsic]          = Auger
    data.recombinationAuger[iphip2, regionIntrinsic]          = Auger
```

ETL

```julia
    data.densityOfStates[iphin3, regionDonor]                 = NC[regionDonor]
    data.densityOfStates[iphip3, regionDonor]                 = NV[regionDonor]
    data.bandEdgeEnergy[iphin3, regionDonor]                  = EC[regionDonor]
    data.bandEdgeEnergy[iphip3, regionDonor]                  = EV[regionDonor]
    data.mobility[iphin3, regionDonor]                        = μn[regionDonor]
    data.mobility[iphip3, regionDonor]                        = μp[regionDonor]
    data.recombinationSRHLifetime[iphin3, regionDonor]        = τn[regionDonor]
    data.recombinationSRHLifetime[iphip3, regionDonor]        = τp[regionDonor]
    data.recombinationSRHTrapDensity[iphin3, regionDonor]     = ChargeTransportInSolids.trapDensity(iphin3, regionDonor, data, EI[regionDonor])
    data.recombinationSRHTrapDensity[iphip3, regionDonor]     = ChargeTransportInSolids.trapDensity(iphip3, regionDonor, data, EI[regionDonor])
    data.recombinationAuger[iphin3, regionDonor]              = Auger
    data.recombinationAuger[iphip3, regionDonor]              = Auger
```

interior region data

```julia
    for ireg in 1:numberOfRegions
        data.dielectricConstant[ireg]                 = ε[ireg]
        data.recombinationRadiative[ireg]             = r0[ireg]
    end
```

interior doping

```julia
    data.doping[iphip1, regionAcceptor]            = Na
    data.doping[iphia, regionIntrinsic]            = C0
    data.doping[iphin3, regionDonor]               = Nd
```

boundary doping

```julia
    data.bDoping[iphip1, bregionAcceptor]          = Na
    data.bDoping[iphin3, bregionDonor]             = Nd
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
    num_species = (numberOfCarriers + 1),
    flux        = ChargeTransportInSolids.Sedan!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
    reaction    = ChargeTransportInSolids.reaction!,
    breaction   = ChargeTransportInSolids.breactionOhmic!,
    storage     = ChargeTransportInSolids.storage!
    )

    sys         = VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
```

enable bulk species

```julia
    enable_species!(sys, iphin1, [regionAcceptor])
    enable_species!(sys, iphip1, [regionAcceptor])
    enable_species!(sys, iphin2, [regionIntrinsic])
    enable_species!(sys, iphip2, [regionIntrinsic])
    enable_species!(sys, iphia,  [regionIntrinsic])
    enable_species!(sys, iphin3, [regionDonor])
    enable_species!(sys, iphip3, [regionDonor])
    enable_species!(sys, ipsi, regions)

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
    control.tol_absolute      = 1.0e-12
    control.tol_relative      = 1.0e-12
    control.handle_exceptions = true
    control.tol_round         = 1.0e-12
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
    @views initialGuess           .= 0.0

    control.damp_initial      = 0.5
    control.damp_growth       = 1.21 # >= 1
    control.max_round         = 5
```

set Dirichlet boundary conditions (Ohmic contacts), in Equilibrium we impose homogeneous Dirichlet conditions,
i.e. the boundary values at outer boundaries are zero.

```julia
    sys.boundary_factors[iphin1, bregionAcceptor] = VoronoiFVM.Dirichlet
    sys.boundary_values[iphin1,  bregionAcceptor] = 0.0 * V

    sys.boundary_factors[iphip1, bregionAcceptor] = VoronoiFVM.Dirichlet
    sys.boundary_values[iphip1,  bregionAcceptor] = 0.0 * V
    ######
    sys.boundary_factors[iphin3, bregionDonor]    = VoronoiFVM.Dirichlet
    sys.boundary_values[iphin3,  bregionDonor]    = 0.0 * V

    sys.boundary_factors[iphip3, bregionDonor]    = VoronoiFVM.Dirichlet
    sys.boundary_values[iphip3,  bregionDonor]    = 0.0 * V


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

    dfusion_grid          = 1.0e-2.*readdlm("data-Na-1p21e22/Driftfusion-pedotpss-grid.dat")

    dfusion_n_interface   = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-n-EQ.dat")
    dfusion_p_interface   = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-p-EQ.dat")
    dfusion_a_interface   = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-a-EQ.dat")
    dfusion_psi_interface = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-psi-EQ.dat")

    dfusion_p             = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-p-EQ.dat")
    dfusion_n             = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-n-EQ.dat")
    dfusion_a             = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-a-EQ.dat")
    dfusion_psi           = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-psi-EQ.dat")


    if plotting
        Plotter.figure()
        ChargeTransportInSolids.plotDensities(Plotter, grid, data, solution,"Equilibrium; limited ions; \$ E_a =\$$(textEa)eV; \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$")
        PyPlot.semilogy(dfusion_grid', (dfusion_n[1,:]), linewidth = 3, linestyle= ":", color="black" )
        PyPlot.semilogy(dfusion_grid', (dfusion_p[1,:]), linewidth = 3, linestyle= ":", color="black" )
        PyPlot.semilogy(dfusion_grid', (dfusion_a[1,:]), linewidth = 3, linestyle= ":", color="black" )
```

```julia
        PyPlot.semilogy(dfusion_grid', (dfusion_n_interface[1,:]), linewidth = 3, linestyle= ":", color="grey" )
        PyPlot.semilogy(dfusion_grid', (dfusion_p_interface[1,:]), linewidth = 3, linestyle= ":", color="grey" )
        PyPlot.semilogy(dfusion_grid', (dfusion_a_interface[1,:]), linewidth = 3, linestyle= ":", color="grey" )
        PyPlot.axvline(h_pdoping, color="black", linestyle="solid")
        PyPlot.axvline(h_pdoping + h_intrinsic, color="black", linestyle="solid")
        ################
        Plotter.figure()
        ChargeTransportInSolids.plotSolution(Plotter, coord, solution, 0.0, "Equilibrium; limited ions; \$ E_a =\$$(textEa)eV; \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$")
        PyPlot.plot(dfusion_grid', (dfusion_psi[1,:]- 5.02*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="black")
```

```julia
        PyPlot.plot(dfusion_grid', (dfusion_psi_interface[1,:]- 5.02*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="grey")
        PyPlot.axvline(h_pdoping, color="black", linestyle="solid")
        PyPlot.axvline(h_pdoping + h_intrinsic, color="black", linestyle="solid")
    end

    if test == false
        println("*** done\n")
    end

   ################################################################################
    if test == false
        println("I-V Measurement Loop for turning electrochemical reaction on")
    end
    ################################################################################
    data.inEquilibrium = false

    control.damp_initial    = 0.5
    control.damp_growth     = 1.21 # >= 1
    control.max_round       = 5
    control.tol_absolute    = 1.0e-9
    control.tol_relative    = 1.0e-9
    control.tol_round       = 1.0e-9
```

there are different way to control timestepping
Here we assume these primary data

```julia
    scanrate                = 0.04 * V/s
    ntsteps                 = 40
    vend                    = voltageAcceptor
    v0                      = 0.0
```

The end time then is calculated here:

```julia
    tend                    = vend/scanrate
```

with fixed timestep sizes we can calculate the times
a priori

```julia
    tvalues                 = range(0, stop = tend, length = ntsteps)

    IVForward          = zeros(0) # for IV values
    biasValuesForward  = zeros(0) # for bias values

    for istep = 2:ntsteps
        t                   = tvalues[istep] # Actual time
        Δu                  = v0 + t*scanrate # Applied voltage
        Δt                  = t - tvalues[istep-1] # Time step size
```

Apply new voltage

```julia
        sys.boundary_values[iphin, bregionAcceptor]      = Δu
        sys.boundary_values[iphip, bregionAcceptor]      = Δu

        if test == false
            println("time value: t = $(t)")
        end
```

Solve time step problems with timestep Δt. initialGuess plays the role of the solution
from last timestep

```julia
        solve!(solution, initialGuess, sys, control = control, tstep = Δt)
```

get IV curve

```julia
        factory = VoronoiFVM.TestFunctionFactory(sys)
```

testfunction zero in bregionAcceptor and one in bregionDonor

```julia
        tf1     = testfunction(factory, [bregionDonor], [bregionAcceptor])
        I1      = integrate(sys, tf1, solution, initialGuess, Δt)

        push!(IVForward, (I1[ipsi] + I1[iphin] + I1[iphip] + I1[iphia]) )
        push!(biasValuesForward, Δu)

        initialGuess .= solution
    end # time loop


    #resForward = [biasValuesForward IVForward]
    #writedlm("jl-IV-forward-Na-$(textNa)-Ea-$(textEa).dat",resForward)
    #writedlm("jl-sol-Na-$(textNa)-Ea-$(textEa).dat", [coord solution'])

    dfusion_n             = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-n-t-30.dat")
    dfusion_p             = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-p-t-30.dat")
    dfusion_a             = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-a-t-30.dat")
    dfusion_psi           = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-psi-t-30.dat")

    dfusion_phin   = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-Efn-t-30.dat")
    dfusion_phip   = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-Efp-t-30.dat")
```

```julia
    dfusion_n_interface   = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-n-t-30.dat")
    dfusion_p_interface   = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-p-t-30.dat")
    dfusion_a_interface   = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-a-t-30.dat")
    dfusion_psi_interface = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-psi-t-30.dat")

    dfusion_phin_interface   = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-Efn-t-30.dat")
    dfusion_phip_interface   = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-Efp-t-30.dat")

    IV_Driftfusion           = readdlm("data-Na-1p21e22/Driftfusion-pedotpss-Na-1p21e22-IV-forward.dat")
    IV_Driftfusion_interface = readdlm("data-Na-1p21e22/with-interface/Driftfusion-pedotpss-Na-1p21e22-interface-reco-IV-forward.dat")
    IV_measured              = readdlm("data-Na-1p21e22/Driftfusion-IV-measurement-pcb-forward.dat")

    if plotting
        Plotter.figure()
        ChargeTransportInSolids.plotDensities(Plotter, grid, data, solution, "\$ \\Delta u = $(biasValuesForward[end])\$; \$ E_a =\$$(textEa)eV;  \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$")
        PyPlot.semilogy(dfusion_grid', (dfusion_n[1,:]), linewidth = 3, linestyle= ":", color="black" )
        PyPlot.semilogy(dfusion_grid', (dfusion_p[1,:]), linewidth = 3, linestyle= ":", color="black" )
        PyPlot.semilogy(dfusion_grid', (dfusion_a[1,:]), linewidth = 3, linestyle= ":", color="black" )
```

```julia
        PyPlot.semilogy(dfusion_grid', (dfusion_n_interface[1,:]), linewidth = 3, linestyle= ":", color="grey" )
        PyPlot.semilogy(dfusion_grid', (dfusion_p_interface[1,:]), linewidth = 3, linestyle= ":", color="grey" )
        PyPlot.semilogy(dfusion_grid', (dfusion_a_interface[1,:]), linewidth = 3, linestyle= ":", color="grey" )
        PyPlot.axvline(h_pdoping, color="black", linestyle="solid")
        PyPlot.axvline(h_pdoping + h_intrinsic, color="black", linestyle="solid")
        ################
        Plotter.figure()
        ChargeTransportInSolids.plotSolution(Plotter, coord, solution, data.Eref, "\$ \\Delta u = $(biasValuesForward[end])\$; \$ E_a =\$$(textEa)eV;  \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$")
        PyPlot.plot(dfusion_grid', (dfusion_psi[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="black")
        PyPlot.plot(dfusion_grid', (-dfusion_phin[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="black")
        PyPlot.plot(dfusion_grid', (-dfusion_phip[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="black")
```

```julia
        PyPlot.plot(dfusion_grid', (dfusion_psi_interface[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="grey")
        PyPlot.plot(dfusion_grid', (-dfusion_phin_interface[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="grey")
        PyPlot.plot(dfusion_grid', (-dfusion_phip_interface[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="grey")
        PyPlot.axvline(h_pdoping, color="black", linestyle="solid")
        PyPlot.axvline(h_pdoping + h_intrinsic, color="black", linestyle="solid")
        ################
        Plotter.figure()
        PyPlot.plot(IV_measured[:, 1], IV_measured[:, 2]*1.0e-3, label = "measurement",  linestyle="--", color = "black")
        Plotter.plot(biasValuesForward, IVForward.*(cm^2), label = "(without internal BC)",  linewidth= 3, linestyle="--", color="red")
        PyPlot.plot(IV_Driftfusion[:, 1], IV_Driftfusion[:, 2], label = "Driftfusion", linewidth = 3, markersize="7", marker= "o", linestyle="--", color = "green")
        PyPlot.plot(IV_Driftfusion_interface[:, 1], IV_Driftfusion_interface[:, 2], label = "Driftfusion (surface reco)", linewidth = 3, markersize="7", marker= "o", linestyle="--", color = "blue")
        PyPlot.ylim(0.0, 0.006)
        Plotter.legend()
        Plotter.title("Forward;\$ E_a =\$$(textEa)eV;  \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$  ")
        Plotter.xlabel("Applied Voltage [V]")
        Plotter.ylabel("current density [A \$ cm^{-2}\$ ]")
    end

    res = [biasValuesForward IVForward];
    writedlm("jl-IV-forward-without-internal-BC-Na-$(textNa)-Ea-$(textEa).dat",res)

end #  main

function test()
    #testval=-3.9737356405650717
    main(test = true, unknown_storage=:dense) ≈ testval #&& main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

