#=
# 104: 1D PSC p-i-n device with ions (linear I-V scan protocol).
([source code](SOURCE_URL))

Simulating a three layer PSC device SiO2| MAPI | SiO2 with mobile ions 
where the ion vacancy accumulation is limited by the Fermi-Dirac integral of order -1.
The simulations are performed out of equilibrium, time-dependent and with
abrupt interfaces. An linear I-V measurement protocol is included and the corresponding
solution vectors after the scan can be depicted.

This simulation coincides with the one made in Section 4.3
of Calado et al. (https://arxiv.org/abs/2009.04384).
The paramters can be found here:
https://github.com/barnesgroupICL/Driftfusion/blob/Methods-IonMonger-Comparison/Input_files/IonMonger_default_bulk.csv.
=#

module Example104_PSC_withIons_IVMeasurement

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize

function main(;n = 8, Plotter = nothing, plotting = false, verbose = false, test = false, unknown_storage=:dense)
    
    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # region numbers
    regionDonor           = 1                           # n doped region
    regionIntrinsic       = 2                           # intrinsic region
    regionAcceptor        = 3                           # p doped region
    regions               = [regionDonor, regionIntrinsic, regionAcceptor]

    # boundary region numbers
    bregionDonor          = 1
    bregionAcceptor       = 2
    bregions              = [bregionDonor, bregionAcceptor]

    # grid
    h_ndoping             = 9.90e-6 * cm 
    h_intrinsic           = 4.00e-5 * cm + 2.0e-7 * cm
    h_pdoping             = 1.99e-5 * cm
    heightLayers          = [h_ndoping,
                             h_ndoping + h_intrinsic,
                             h_ndoping + h_intrinsic + h_pdoping]

    x0              = 0.0 * cm 
    δ               = 4*n        # the larger, the finer the mesh
    t               = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k               = 1.5        # the closer to 1, the closer to the boundary geomspace

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
    
    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],        [heightLayers[1]], regionDonor, tol = 1.0e-12)       # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionIntrinsic, tol = 1.0e-12)   # intrinsic region = 2  
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionAcceptor, tol = 1.0e-12)    # p-doped region   = 3  

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
    
    # indices
    iphin, iphip, iphia, ipsi      = 1:4
    species                        = [iphin, iphip, iphia, ipsi]
    
    # number of (boundary) regions and carriers
    numberOfRegions         = length(regions)
    numberOfBoundaryRegions = length(bregions) 
    numberOfCarriers        = length(species) - 1
    
    # temperature
    T               = 300.0                 *  K
    
    # band edge energies    
    Ec_d            = -4.0                  *  eV 
    Ev_d            = -5.8                  *  eV 
    
    Ec_i            = -3.7                  *  eV 
    Ev_i            = -5.4                  *  eV 
    
    Ec_a            = -3.4                  *  eV 
    Ev_a            = -5.1                  *  eV 

    ###################### adjust Na, Ea here #####################
    Nanion          = 1.0e21                / (cm^3)
    Ea_i            = -4.45                 *  eV 
    # for the labels in the figures
    textEa          = Ea_i./eV
    textNa          = Nanion.*cm^3
    ###################### adjust Na, Ea here #####################
        
    EC              = [Ec_d, Ec_i, Ec_a] 
    EV              = [Ev_d, Ev_i, Ev_a] 
    EA              = [0.0,  Ea_i,  0.0]
    
    # effective densities of state
    Nc_d            = 5.0e19                / (cm^3)
    Nv_d            = 5.0e19                / (cm^3)
    
    Nc_i            = 8.1e18                / (cm^3)
    Nv_i            = 5.8e18                / (cm^3)

    Nc_a            = 5.0e19                / (cm^3)
    Nv_a            = 5.0e19                / (cm^3)

    NC              = [Nc_d, Nc_i,  Nc_a]
    NV              = [Nv_d, Nv_i,  Nv_a]
    NAnion          = [0.0,  Nanion, 0.0]

    # mobilities 
    μn_d            = 3.89                  * (cm^2) / (V * s)  
    μp_d            = 3.89                  * (cm^2) / (V * s)  

    μn_i            = 6.62e1                * (cm^2) / (V * s)  
    μp_i            = 6.62e1                * (cm^2) / (V * s)

    μa_i            = 3.93e-12              * (cm^2) / (V * s)

    μn_a            = 3.89e-1               * (cm^2) / (V * s) 
    μp_a            = 3.89e-1               * (cm^2) / (V * s)
    
    μn              = [μn_d, μn_i, μn_a] 
    μp              = [μp_d, μp_i, μp_a] 
    μa              = [0.0,  μa_i, 0.0 ] 

    # relative dielectric permittivity  
    ε_d             = 10.0                  *  1.0  
    ε_i             = 24.1                  *  1.0 
    ε_a             = 3.0                   *  1.0 

    ε               = [ε_d, ε_i, ε_a] 

    # recombination model
    recombinationOn = true

    # radiative recombination
    r0_d            = 0.0e+0               * cm^3 / s 
    r0_i            = 1.0e-12              * cm^3 / s  
    r0_a            = 0.0e+0               * cm^3 / s
        
    r0              = [r0_d, r0_i, r0_a]
        
    # life times and trap densities 
    τn_d            = 1.0e100              * s 
    τp_d            = 1.0e100              * s
        
    τn_i            = 3.0e-10              * s
    τp_i            = 3.0e-8               * s
    τn_a            = τn_d
    τp_a            = τp_d

    τn              = [τn_d, τn_i, τn_a]
    τp              = [τp_d, τp_i, τp_a]
        
    # SRH trap energies (needed for calculation of recombinationSRHTrapDensity)
    Ei_d            = -5.0                 * eV   
    Ei_i            = -4.55                * eV   
    Ei_a            = -4.1                 * eV   

    EI              = [Ei_d, Ei_i, Ei_a]
        
    # Auger recombination
    Auger           = 0.0

    # doping (doping values are from Driftfusion)
    Nd              =   1.03e18             / (cm^3) 
    Na              =   1.03e18             / (cm^3) 
    C0              =   1.6e19              / (cm^3) 

    # contact voltages: we impose an applied voltage only on one boundary.
    # At the other boundary the applied voltage is zero.
    voltageAcceptor =  1.2                  * V 

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define ChargeTransport data and fill in previously defined data")
    end
    ################################################################################

    # initialize ChargeTransport instance
    data      = ChargeTransportInSolids.ChargeTransportData(numberOfNodes,
                                                            numberOfRegions,
                                                            numberOfBoundaryRegions,
                                                            ;numberOfSpecies = numberOfCarriers + 1)

    # region independent data
    data.F                               = [Boltzmann, Boltzmann,FermiDiracMinusOne] # Boltzmann, FermiDiracOneHalf, FermiDiracMinusOne, Blakemore
    data.temperature                     = T
    data.UT                              = (kB * data.temperature) / q
    data.chargeNumbers[iphin]            = -1
    data.chargeNumbers[iphip]            =  1
    data.chargeNumbers[iphia]            =  1
    data.recombinationOn                 = recombinationOn

    # interior region data
    for ireg in 1:numberOfRegions

        data.dielectricConstant[ireg]                 = ε[ireg]

        # dos, band edge energy and mobilities
        data.densityOfStates[iphin, ireg]             = NC[ireg]
        data.densityOfStates[iphip, ireg]             = NV[ireg]
        data.densityOfStates[iphia, ireg]             = NAnion[ireg]

        data.bandEdgeEnergy[iphin, ireg]              = EC[ireg]
        data.bandEdgeEnergy[iphip, ireg]              = EV[ireg]
        data.bandEdgeEnergy[iphia, ireg]              = EA[ireg]

        data.mobility[iphin, ireg]                    = μn[ireg]
        data.mobility[iphip, ireg]                    = μp[ireg]
        data.mobility[iphia, ireg]                    = μa[ireg]

        # recombination parameters
        data.recombinationRadiative[ireg]             = r0[ireg]
        data.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        data.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        data.recombinationSRHTrapDensity[iphin, ireg] = ChargeTransportInSolids.trapDensity(iphin, ireg, data, EI[ireg])
        data.recombinationSRHTrapDensity[iphip, ireg] = ChargeTransportInSolids.trapDensity(iphip, ireg, data, EI[ireg])
        data.recombinationAuger[iphin, ireg]          = Auger
        data.recombinationAuger[iphip, ireg]          = Auger
    end

    # boundary region data
    data.bDensityOfStates[iphin, bregionDonor]        = Nc_d
    data.bDensityOfStates[iphip, bregionDonor]        = Nv_d

    data.bDensityOfStates[iphin, bregionAcceptor]     = Nc_a
    data.bDensityOfStates[iphip, bregionAcceptor]     = Nv_a

    data.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_d
    data.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_d

    data.bBandEdgeEnergy[iphin, bregionAcceptor]      = Ec_a
    data.bBandEdgeEnergy[iphip, bregionAcceptor]      = Ev_a

    # interior doping
    data.doping[iphin, regionDonor]                   = Nd
    data.doping[iphia, regionIntrinsic]               = C0
    data.doping[iphip, regionAcceptor]                = Na  

    # boundary doping
    data.bDoping[iphip, bregionAcceptor]              = Na      
    data.bDoping[iphin, bregionDonor]                 = Nd      

    if test == false
        println("*** done\n")
    end

    if plotting == true
        ################################################################################
        println("Plot electroneutral potential, band-edge energies and doping")
        ################################################################################
        ChargeTransportInSolids.plotEnergies(Plotter, grid, data)
        Plotter.figure()
        ChargeTransportInSolids.plotDoping(Plotter, grid, data)
        Plotter.figure()
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define physics and system")
    end
    ################################################################################

    ## initializing physics environment ##
    physics = VoronoiFVM.Physics(
    data        = data,
    num_species = numberOfCarriers + 1,
    flux        = ChargeTransportInSolids.Sedan!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
    reaction    = ChargeTransportInSolids.reaction!,
    breaction   = ChargeTransportInSolids.breactionOhmic!,
    storage     = ChargeTransportInSolids.storage!
    )

    sys         = VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)

    # enable all three species in all regions
    enable_species!(sys, ipsi,  regions)
    enable_species!(sys, iphin, regions)
    enable_species!(sys, iphip, regions)
    enable_species!(sys, iphia, [regionIntrinsic])

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
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    data.inEquilibrium             = true

    # initialize solution and starting vectors
    initialGuess                   = unknowns(sys)
    solution                       = unknowns(sys)
    @views initialGuess[ipsi,  :] .= 0.0
    @views initialGuess[iphin, :] .= 0.0
    @views initialGuess[iphip, :] .= 0.0
    @views initialGuess[iphia, :] .= 0.0

    control.damp_initial      = 0.1
    control.damp_growth       = 1.61 # >= 1
    control.max_round         = 5

    # set Dirichlet boundary conditions (Ohmic contacts), in Equilibrium we impose homogeneous Dirichlet conditions,
    # i.e. the boundary values at outer boundaries are zero.
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
        initialGuess .= solution
    end

    if plotting 
        ChargeTransportInSolids.plotEnergies(Plotter, grid, data, solution, "Equilibrium; \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$")
        Plotter.figure()
        ChargeTransportInSolids.plotDensities(Plotter, grid, data, solution,"Equilibrium; \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$")
        Plotter.figure()
        ChargeTransportInSolids.plotSolution(Plotter, coord, solution, 0.0, "Equilibrium; \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$")
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    data.inEquilibrium        = false

    control.damp_initial      = 0.5
    control.damp_growth       = 1.61 # >= 1
    control.max_round         = 7
 
    # there are different way to control timestepping
    # Here we assume these primary data
    scanrate                  = 1.0 * V/s
    ntsteps                   = 101
    vend                      = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
    v0                        = 0.0
    IV                        = zeros(0) # for IV values
    biasValues                = zeros(0) # for bias values

    # The end time then is calculated here:
    tend                      = vend/scanrate

    # with fixed timestep sizes we can calculate the times
    # a priori
    tvalues                   = range(0, stop = tend, length = ntsteps)

    for istep = 2:ntsteps
        
        t                     = tvalues[istep] # Actual time
        Δu                    = v0 + t*scanrate # Applied voltage 
        Δt                    = t - tvalues[istep-1] # Time step size
        
        # Apply new voltage
        sys.boundary_values[iphin, bregionAcceptor]      = Δu
        sys.boundary_values[iphip, bregionAcceptor]      = Δu
        
        if test == false
            println("time value: t = $(t)")
        end

        # Solve time step problems with timestep Δt. initialGuess plays the role of the solution
        # from last timestep
        solve!(solution, initialGuess, sys, control = control, tstep = Δt)
        
        # get IV curve
        factory = VoronoiFVM.TestFunctionFactory(sys)

        # testfunction zero in bregionAcceptor and one in bregionDonor
        tf     = testfunction(factory, [bregionDonor], [bregionAcceptor])
        I      = integrate(sys, tf, solution, initialGuess, Δt)

        current = I[ipsi] + I[iphin] + I[iphip] + I[iphia]

        push!(IV,  current)
        push!(biasValues, Δu)

        initialGuess .= solution

    end # time loop

    # here in res the biasValues and the corresponding current are stored.
    #res = [biasValues IV];

    if plotting 
        ChargeTransportInSolids.plotEnergies(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(v0+vend); \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$")
        Plotter.figure()
        ChargeTransportInSolids.plotDensities(Plotter, grid, data, solution,"bias \$\\Delta u\$ = $(v0+vend); \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$")
        Plotter.figure()
        ChargeTransportInSolids.plotSolution(Plotter, coord, solution, 0.0, "bias \$\\Delta u\$ = $(v0+vend); \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$")
    end

    testval = solution[ipsi, 15]
    return testval
    
    println("*** done\n")

end #  main

function test()
    testval = -4.100369963410221
    main(test = true, unknown_storage=:dense) ≈ testval  #&& main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
