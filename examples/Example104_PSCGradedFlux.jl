"""
Simulating a three layer PSC device without mobile ions.
The simulations are performed out of equilibrium and with
two junctions between perovskite layer and transport layers, to 
which we refer as graded interfaces.
Hence, a graded flux discretizations with space dependent
band-edge energies and density of states are tested here.

This simulation coincides with the one made in Section 4.3
of Calado et al. (https://arxiv.org/abs/2009.04384).
The paramters can be found here:
https://github.com/barnesgroupICL/Driftfusion/blob/Methods-IonMonger-comparison/Input_files/IonMonger_default_noIR.csv.

"""

module Example104_PSCGradedFlux

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using Printf

function gradingParameter(physicalParameter, coord, regionTransportLayers, regionJunctions, h, heightLayers, lengthLayers, values)
    for ireg in regionTransportLayers

        xcoord                     = lengthLayers[ireg]:lengthLayers[ireg+1]
        physicalParameter[xcoord] .= values[ireg]

    end

    for ireg in regionJunctions 

        xcoord   = lengthLayers[ireg]:lengthLayers[ireg+1] 
        left     = lengthLayers[ireg]-3  
        junction = h[ireg]
        right    = lengthLayers[ireg+2]-3

        gradient = ( physicalParameter[right] - physicalParameter[left] )/ junction

        for index in xcoord
            physicalParameter[index] = physicalParameter[left] + (coord[index] - heightLayers[ireg-1]) * gradient
        end

    end

    return physicalParameter
end

function main(;n = 4, Plotter = nothing, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # region numbers
    regionDonor           = 1                           # n doped region
    regionJunction1       = 2
    regionIntrinsic       = 3                           # intrinsic region
    regionJunction2       = 4
    regionAcceptor        = 5                           # p doped region
    regions               = [regionDonor, regionJunction1, regionIntrinsic, regionJunction2, regionAcceptor]
    regionTransportLayers = [regionDonor, regionIntrinsic, regionAcceptor]
    regionJunctions       = [regionJunction1, regionJunction2]

    # boundary region numbers
    bregionDonor          = 1
    bregionAcceptor       = 2
    bregions              = [bregionDonor, bregionAcceptor]

    # grid
    h_ndoping             = 9.90e-6 * cm 
    h_junction1           = 1.0e-7  * cm
    h_intrinsic           = 4.00e-5 * cm
    h_junction2           = 1.0e-7  * cm
    h_pdoping             = 1.99e-5 * cm
    h                     = [h_ndoping, h_junction1, h_intrinsic, h_junction2, h_pdoping]
    heightLayers          = [h_ndoping,
                             h_ndoping + h_junction1,
                             h_ndoping + h_junction1 + h_intrinsic,
                             h_ndoping + h_junction1 + h_intrinsic + h_junction2,
                             h_ndoping + h_junction1 + h_intrinsic + h_junction2 + h_pdoping]
    refinementfactor      = 2^(n-1)

    coord_ndoping         = collect(range(0.0, stop = h_ndoping, length = 4 * refinementfactor))
    length_n              = length(coord_ndoping)
    coord_junction1       = collect(range(h_ndoping,
                                        stop = h_ndoping + h_junction1,
                                        length = 3 * refinementfactor))
    coord_intrinsic       = collect(range(h_ndoping + h_junction1,
                                        stop = (h_ndoping + h_junction1 + h_intrinsic),
                                        length = 10 * refinementfactor))
    coord_junction2       = collect(range(h_ndoping + h_junction1 + h_intrinsic,
                                        stop = (h_ndoping + h_junction1 + h_intrinsic + h_junction2),
                                        length = 3 * refinementfactor))
    coord_pdoping         = collect(range((h_ndoping + h_junction1 + h_intrinsic + h_junction2),
                                        stop = (h_ndoping + h_junction1 + h_intrinsic + h_junction2 + h_pdoping),
                                        length = 4 * refinementfactor))

    coord                 = glue(coord_ndoping, coord_junction1)
    length_j1             = length(coord)
    coord                 = glue(coord, coord_intrinsic)
    length_i              = length(coord)
    coord                 = glue(coord, coord_junction2)
    length_j2             = length(coord)
    coord                 = glue(coord, coord_pdoping)

    grid                  = ExtendableGrids.simplexgrid(coord)
    numberOfNodes         = length(coord)
    lengthLayers          = [1, length_n, length_j1, length_i, length_j2, numberOfNodes]

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],        [heightLayers[1]], regionDonor)       # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionJunction1)   # first junction   = 2  
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionIntrinsic)   # intrinsic region = 3  
    cellmask!(grid, [heightLayers[3]], [heightLayers[4]], regionJunction2)   # sec. junction    = 4
    cellmask!(grid, [heightLayers[4]], [heightLayers[5]], regionAcceptor)    # p-doped region   = 5  

    if plotting
        ExtendableGrids.plot(grid, Plotter = Plotter, p = Plotter.plot()) 
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
    iphin, iphip, ipsi      = 1:3
    species                 = [iphin, iphip, ipsi]

    # number of (boundary) regions and carriers
    numberOfRegions         = length(regions)
    numberOfBoundaryRegions = length(bregions) 
    numberOfCarriers        = length(species) - 1

    # temperature
    T                       = 300.0                 *  K

    # band edge energies    
    Ec_d                    = -4.0                  *  eV 
    Ev_d                    = -6.0                  *  eV 
        
    Ec_i                    = -3.7                  *  eV 
    Ev_i                    = -5.4                  *  eV 
        
    Ec_a                    = -3.1                  *  eV 
    Ev_a                    = -5.1                  *  eV 
    
    # these parameters at the junctions for E_\alpha and N_\alpha will be overwritten.
    Ec_j1                   = Ec_d;     Ec_j2     = Ec_i
    Ev_j1                   = Ev_d;     Ev_j2     = Ev_i
          
    EC                      = [Ec_d, Ec_j1, Ec_i, Ec_j2, Ec_a] 
    EV                      = [Ev_d, Ev_j1, Ev_i, Ev_j2, Ev_a] 

    # effective densities of state
    Nc_d                   = 5.0e19                / (cm^3)
    Nv_d                   = 5.0e19                / (cm^3)

    Nc_i                   = 8.1e18                / (cm^3)
    Nv_i                   = 5.8e18                / (cm^3)

    Nc_a                   = 5.0e19                / (cm^3)
    Nv_a                   = 5.0e19                / (cm^3)

    Nc_j1                  = Nc_d;     Nc_j2      = Nc_i
    Nv_j1                  = Nv_d;     Nv_j2      = Nv_i

    NC                     = [Nc_d, Nc_j1, Nc_i, Nc_j2, Nc_a]
    NV                     = [Nv_d, Nv_j1, Nv_i, Nv_j2, Nv_a]
 
    # mobilities 
    μn_d                   = 3.89                  * (cm^2) / (V * s)  
    μp_d                   = 3.89                  * (cm^2) / (V * s)  

    μn_i                   = 6.62e1                * (cm^2) / (V * s)  
    μp_i                   = 6.62e1                * (cm^2) / (V * s)

    μn_a                   = 3.89e-1               * (cm^2) / (V * s) 
    μp_a                   = 3.89e-1               * (cm^2) / (V * s)
    
    μn_j1                  = μn_d;     μn_j2      = μn_i
    μp_j1                  = μp_d;     μp_j2      = μp_i

    μn                     = [μn_d, μn_j1, μn_i, μn_j2, μn_a] 
    μp                     = [μp_d, μp_j1, μp_i, μp_j2, μp_a] 

    # relative dielectric permittivity  
    ε_d                    = 10.0                  *  1.0  
    ε_i                    = 24.1                  *  1.0 
    ε_a                    = 3.0                   *  1.0 

    ε_j1                   = ε_d;       ε_j2      = ε_a

    ε                      = [ε_d, ε_j1, ε_i, ε_j2, ε_a] 

    # recombination model
    recombinationOn        = true

    # radiative recombination
    r0_d                   = 0.0e+0               * cm^3 / s 
    r0_i                   = 1.0e-12              * cm^3 / s  
    r0_a                   = 0.0e+0               * cm^3 / s

    r0_j1                  = r0_i;      r0_j2     = r0_i
        
    r0                     = [r0_d, r0_j1, r0_i, r0_j2, r0_a]
        
    # life times and trap densities 
    τn_d                   = 1.0e100              * s 
    τp_d                   = 1.0e100              * s
        
    τn_i                   = 3.0e-10              * s
    τp_i                   = 3.0e-8               * s
    τn_a                   = τn_d
    τp_a                   = τp_d

    τn_j1                  = τn_i;     τn_j2      = τn_a
    τp_j1                  = τp_i;     τp_j2      = τp_a
        
    τn                     = [τn_d, τn_j1, τn_i, τn_j2, τn_a]
    τp                     = [τp_d, τp_j1, τp_i, τp_j2, τp_a]
        
    # SRH trap energies (needed for calculation of recombinationSRHTrapDensity)
    Ei_d                   = -5.0                 * eV   
    Ei_i                   = -4.55                * eV   
    Ei_a                   = -4.1                 * eV   

    Ei_j1                  = Ei_i;      Ei_j2     = Ei_a

    EI                     = [Ei_d, Ei_j1, Ei_i, Ei_j2, Ei_a]
        
    # Auger recombination
    Auger                  = 0.0

    # doping (doping values are from Driftfusion)
    Nd                     =   1.03e18             / (cm^3) 
    Na                     =   1.03e18             / (cm^3) 
    Ni_acceptor            =   8.32e7              / (cm^3) 

    # contact voltages
    voltageAcceptor        =  1.2                 * V 
    voltageDonor           =  0.0                 * V 

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
    data.F                              .= Boltzmann # Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA, FermiDiracMinusOne, Blakemore
    data.temperature                     = T
    data.UT                              = (kB * data.temperature) / q
    data.contactVoltage[bregionAcceptor] = voltageAcceptor
    data.contactVoltage[bregionDonor]    = voltageDonor
    data.chargeNumbers[iphin]            = -1
    data.chargeNumbers[iphip]            =  1
    data.recombinationOn                 = recombinationOn

    # band-edge energies
    data.bandEdgeEnergyNode[iphin, :]    = gradingParameter(data.bandEdgeEnergyNode[iphin, :],
                                                            coord, regionTransportLayers, regionJunctions, h,
                                                            heightLayers, lengthLayers, EC)
    data.bandEdgeEnergyNode[iphip, :]    = gradingParameter(data.bandEdgeEnergyNode[iphip, :],
                                                            coord, regionTransportLayers, regionJunctions, h,
                                                            heightLayers, lengthLayers, EV)
    # # density of states
    data.densityOfStatesNode[iphin, :]   = gradingParameter(data.densityOfStatesNode[iphin, :],
                                                            coord, regionTransportLayers, regionJunctions, h,
                                                            heightLayers, lengthLayers, NC)
    data.densityOfStatesNode[iphip, :]   = gradingParameter(data.densityOfStatesNode[iphip, :],
                                                            coord, regionTransportLayers, regionJunctions, h,
                                                            heightLayers, lengthLayers, NV)
    # region dependent data
    for ireg in 1:numberOfRegions

        # mobility
        data.mobility[iphin, ireg]                    = μn[ireg]
        data.mobility[iphip, ireg]                    = μp[ireg]

        data.dielectricConstant[ireg]                 = ε[ireg]
        # recombination parameters
        data.recombinationRadiative[ireg]             = r0[ireg]
        data.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        data.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        data.recombinationSRHTrapDensity[iphin, ireg] = ChargeTransportInSolids.trapDensity(iphin, ireg, data, EI[ireg])
        data.recombinationSRHTrapDensity[iphip, ireg] = ChargeTransportInSolids.trapDensity(iphip, ireg, data, EI[ireg])
        data.recombinationAuger[iphin, ireg]          = Auger
        data.recombinationAuger[iphip, ireg]          = Auger
    end           
    
    # interior doping
    data.doping[iphin, regionDonor]               = Nd
    data.doping[iphip, regionIntrinsic]           = Ni_acceptor    
    data.doping[iphip, regionAcceptor]            = Na     
                             
    # boundary doping
    data.bDoping[iphip, bregionAcceptor]          = Na        # data.bDoping  = [Na  0.0;
    data.bDoping[iphin, bregionDonor]             = Nd        #                  0.0  Nd]

    if test == false
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
    flux        = ChargeTransportInSolids.ScharfetterGummelGraded!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
    reaction    = ChargeTransportInSolids.reaction!,
    breaction   = ChargeTransportInSolids.breactionOhmic!
    )

    sys         = VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)

    # enable all three species in all regions
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

    # initialize solution and starting vectors
    initialGuess                   = unknowns(sys)
    solution                       = unknowns(sys)
    @views initialGuess[ipsi,  :] .= 0.0
    @views initialGuess[iphin, :] .= 0.0
    @views initialGuess[iphip, :] .= 0.0

    control.damp_initial           = 0.5
    control.damp_growth            = 1.61 # >= 1
    control.max_round              = 5

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
        ChargeTransportInSolids.plotEnergies(Plotter, grid, data, solution, "EQULIBRIUM (NO illumination)")
        Plotter.figure()
        ChargeTransportInSolids.plotDensities(Plotter, grid, data, solution, "EQULIBRIUM (NO illumination)")
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

    control.damp_initial                             = 0.5
    control.damp_growth                              = 1.21 # >= 1
    control.max_round                                = 7

    # set non equilibrium boundary conditions
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
        ChargeTransportInSolids.plotEnergies(Plotter, grid, data, solution, maxBias)
        Plotter.figure()
        ChargeTransportInSolids.plotDensities(Plotter, grid, data, solution, maxBias)
    end

    if test == false
        println("*** done\n")
    end

    testval = solution[ipsi, 20]
    return testval

end #  main

function test()
    testval=-4.100368066229441
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
