#=
# 103: 1D PSC p-i-n device with graded interfaces.
([source code](SOURCE_URL))

Simulating a three layer PSC device SiO2| MAPI | SiO2 without mobile ions.
The simulations are performed out of equilibrium, stationary and with
two junctions between perovskite layer and transport layers, to 
which we refer as graded interfaces.
Hence, a graded flux discretization with space dependent
band-edge energies and density of states is tested here.

This simulation coincides with the one made in Section 4.3
of Calado et al. (https://arxiv.org/abs/2009.04384).
The paramters can be found in Table S.13 or slightly modified than the one
in the publication here:
https://github.com/barnesgroupICL/Driftfusion/blob/Methods-IonMonger-Comparison/Input_files/IonMonger_default_bulk.csv
=#

module Example104_PSC_gradedFlux

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize


# function for grading the physical parameters 
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
    regionDonor             = 1                           # n doped region
    regionJunction1         = 2
    regionIntrinsic         = 3                           # intrinsic region
    regionJunction2         = 4
    regionAcceptor          = 5                           # p doped region
    regions                 = [regionDonor, regionJunction1, regionIntrinsic, regionJunction2, regionAcceptor]
    regionTransportLayers   = [regionDonor, regionIntrinsic, regionAcceptor]
    regionJunctions         = [regionJunction1, regionJunction2]
    numberOfRegions         = length(regions)

    # boundary region numbers
    bregionDonor            = 1
    bregionAcceptor         = 2
    bregions                = [bregionDonor, bregionAcceptor]
    numberOfBoundaryRegions = length(bregions)

    # grid
    h_ndoping               = 9.90e-6 * cm 
    h_junction1             = 1.0e-7  * cm
    h_intrinsic             = 4.00e-5 * cm
    h_junction2             = 1.0e-7  * cm
    h_pdoping               = 1.99e-5 * cm
    h                       = [h_ndoping, h_junction1, h_intrinsic, h_junction2, h_pdoping]
    heightLayers            = [h_ndoping,
                               h_ndoping + h_junction1,
                               h_ndoping + h_junction1 + h_intrinsic,
                               h_ndoping + h_junction1 + h_intrinsic + h_junction2,
                               h_ndoping + h_junction1 + h_intrinsic + h_junction2 + h_pdoping]
    refinementfactor        = 2^(n-1)

    coord_ndoping           = collect(range(0.0, stop = h_ndoping, length = 4 * refinementfactor))
    length_n                = length(coord_ndoping)
    coord_junction1         = collect(range(h_ndoping,
                                           stop = h_ndoping + h_junction1,
                                           length = 3 * refinementfactor))
    coord_intrinsic         = collect(range(h_ndoping + h_junction1,
                                           stop = (h_ndoping + h_junction1 + h_intrinsic),
                                           length = 10 * refinementfactor))
    coord_junction2         = collect(range(h_ndoping + h_junction1 + h_intrinsic,
                                           stop = (h_ndoping + h_junction1 + h_intrinsic + h_junction2),
                                           length = 3 * refinementfactor))
    coord_pdoping           = collect(range((h_ndoping + h_junction1 + h_intrinsic + h_junction2),
                                            stop = (h_ndoping + h_junction1 + h_intrinsic + h_junction2 + h_pdoping),
                                            length = 4 * refinementfactor))

    coord                   = glue(coord_ndoping, coord_junction1)
    length_j1               = length(coord)
    coord                   = glue(coord, coord_intrinsic)
    length_i                = length(coord)
    coord                   = glue(coord, coord_junction2)
    length_j2               = length(coord)
    coord                   = glue(coord, coord_pdoping)

    grid                    = simplexgrid(coord)
    numberOfNodes           = length(coord)
    lengthLayers            = [1, length_n, length_j1, length_i, length_j2, numberOfNodes]

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],        [heightLayers[1]], regionDonor)       # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionJunction1)   # first junction   = 2  
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionIntrinsic)   # intrinsic region = 3  
    cellmask!(grid, [heightLayers[3]], [heightLayers[4]], regionJunction2)   # sec. junction    = 4
    cellmask!(grid, [heightLayers[4]], [heightLayers[5]], regionAcceptor)    # p-doped region   = 5  

    if plotting
        gridplot(grid, Plotter = Plotter)
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

    ##########      physical data      ##########
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
    Nc_d                    = 5.0e19                / (cm^3)
    Nv_d                    = 5.0e19                / (cm^3)

    Nc_i                    = 8.1e18                / (cm^3)
    Nv_i                    = 5.8e18                / (cm^3)

    Nc_a                    = 5.0e19                / (cm^3)
    Nv_a                    = 5.0e19                / (cm^3)

    Nc_j1                   = Nc_d;     Nc_j2      = Nc_i
    Nv_j1                   = Nv_d;     Nv_j2      = Nv_i

    NC                      = [Nc_d, Nc_j1, Nc_i, Nc_j2, Nc_a]
    NV                      = [Nv_d, Nv_j1, Nv_i, Nv_j2, Nv_a]
 
    # mobilities 
    μn_d                    = 3.89                  * (cm^2) / (V * s)  
    μp_d                    = 3.89                  * (cm^2) / (V * s)  

    μn_i                    = 6.62e1                * (cm^2) / (V * s)  
    μp_i                    = 6.62e1                * (cm^2) / (V * s)

    μn_a                    = 3.89e-1               * (cm^2) / (V * s) 
    μp_a                    = 3.89e-1               * (cm^2) / (V * s)
    
    μn_j1                   = μn_d;     μn_j2      = μn_i
    μp_j1                   = μp_d;     μp_j2      = μp_i

    μn                      = [μn_d, μn_j1, μn_i, μn_j2, μn_a] 
    μp                      = [μp_d, μp_j1, μp_i, μp_j2, μp_a] 

    # relative dielectric permittivity  
    ε_d                     = 10.0                  *  1.0  
    ε_i                     = 24.1                  *  1.0 
    ε_a                     = 3.0                   *  1.0 

    ε_j1                    = ε_d;       ε_j2      = ε_a
 
    ε                       = [ε_d, ε_j1, ε_i, ε_j2, ε_a] 

    # recombination model
    bulk_recombination = bulk_recombination_full

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
        
    # SRH trap energies (needed for calculation of trap_density! (SRH))
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

    # contact voltages: we impose an applied voltage only on one boundary.
    # At the other boundary the applied voltage is zero.
    voltageAcceptor        =  1.2                 * V 

    # interface model (this is needed for giving the user the correct index set)
    interface_reaction = interface_model_none

    # set the correct indices for each species (this is needed for giving the user the correct index set)
    # but likewise it is possible to define one owns index set, i.e. iphin, iphin, iphip = 1:3
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
    data.boundary_type[bregionDonor]    = ohmic_contact  
    data.boundary_type[bregionAcceptor] = ohmic_contact                       
     
    # Following choices are possible for the flux_discretization scheme: ScharfetterGummel, ScharfetterGummel_Graded,
    # excessChemicalPotential, excessChemicalPotential_Graded, diffusionEnhanced, generalized_SG
    data.flux_approximation             = ScharfetterGummel_Graded
    
    ################################################################################
    if test == false
        println("Define ChargeTransportParams and fill in physical parameters")
    end
    ################################################################################

    # for region dependent parameters
    params                                              = ChargeTransportParams(grid, numberOfCarriers)
    # for space dependent parameters
    paramsnodal                                         = ChargeTransportParamsNodal(grid, numberOfCarriers)


    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1

    # nodal band-edge energies
    paramsnodal.bandEdgeEnergy[iphin, :]  = gradingParameter(paramsnodal.bandEdgeEnergy[iphin, :],
                                                             coord, regionTransportLayers, regionJunctions, h,
                                                             heightLayers, lengthLayers, EC)
    paramsnodal.bandEdgeEnergy[iphip, :]  = gradingParameter(paramsnodal.bandEdgeEnergy[iphip, :],
                                                             coord, regionTransportLayers, regionJunctions, h,
                                                             heightLayers, lengthLayers, EV)
    # nodal effective density of states
    paramsnodal.densityOfStates[iphin, :] = gradingParameter(paramsnodal.densityOfStates[iphin, :],
                                                             coord, regionTransportLayers, regionJunctions, h,
                                                             heightLayers, lengthLayers, NC)
    paramsnodal.densityOfStates[iphip, :] = gradingParameter(paramsnodal.densityOfStates[iphip, :],
                                                             coord, regionTransportLayers, regionJunctions, h,
                                                             heightLayers, lengthLayers, NV)
    # region dependent data
    for ireg in 1:numberOfRegions

        # mobility
        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]

        params.dielectricConstant[ireg]                 = ε[ireg]
        # recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, data, EI[ireg])
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

    end           
    
    # interior doping
    params.doping[iphin, regionDonor]               = Nd
    params.doping[iphip, regionIntrinsic]           = Ni_acceptor    
    params.doping[iphip, regionAcceptor]            = Na     
                             
    # boundary doping
    params.bDoping[iphip, bregionAcceptor]          = Na        # data.bDoping  = [Na  0.0;
    params.bDoping[iphin, bregionDonor]             = Nd        #                  0.0  Nd]

    # Region dependent params is now a substruct of data which is again a substruct of the system and will be parsed 
    # in next step.
    data.params                                     = params

    # same holds true for space dependent params
    data.paramsnodal                                = paramsnodal

    # in the last step, we initialize our system with previous data which is likewise dependent on the parameters. 
    # important that this is in the end, otherwise our VoronoiFVMSys is not dependent on the data we initialized
    # but rather on default data.
    ctsys                                           = ChargeTransportSystem(grid, data, unknown_storage=unknown_storage)
    

    ########### It is also possible to print the nodal dependent data, but for the sake of readibility
    ########### we neglect this here.
    # print data
    if test == false
        # show region dependent physical parameters. show_params() only supports region dependent parameters, but, if one wishes to
        # print nodal dependent parameters, currently this is possible with println(ctsys.data.paramsnodal). We neglected here, since
        # in most applications where the numberOfNodes is >> 10 this would results in a large output in the terminal.
        show_params(ctsys)
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define boundary conditions and enabled layers")
    end
    ################################################################################

    # set ohmic contacts for each charge carrier at all outerior boundaries. First, 
    # we compute equilibrium solutions. Hence the boundary values at the ohmic contacts
    # are zero.
    set_ohmic_contact!(ctsys, iphin, bregionDonor, 0.0)
    set_ohmic_contact!(ctsys, iphip, bregionDonor, 0.0)
    set_ohmic_contact!(ctsys, iphin, bregionAcceptor, 0.0)
    set_ohmic_contact!(ctsys, iphip, bregionAcceptor, 0.0)

    # enable all three species in all regions
    enable_species!(ctsys, ipsi,  regions)
    enable_species!(ctsys, iphin, regions)
    enable_species!(ctsys, iphip, regions)

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

    control.damp_initial           = 0.5
    control.damp_growth            = 1.61 # >= 1
    control.max_round              = 5

    initialGuess          = unknowns(ctsys)
    solution              = unknowns(ctsys)

    solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

    initialGuess         .= solution 

    if plotting
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

    ctsys.data.calculation_type                      = outOfEquilibrium

    control.damp_initial                             = 0.9
    control.damp_growth                              = 1.61 # >= 1
    control.max_round                                = 7

    maxBias    = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 13)

    for Δu in biasValues
        if verbose
            println("Bias value: Δu = $(Δu) (no illumination)")
        end

        # set non equilibrium boundary conditions
        set_ohmic_contact!(ctsys, iphin, bregionAcceptor, Δu)
        set_ohmic_contact!(ctsys, iphip, bregionAcceptor, Δu)

        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        initialGuess .= solution

    end # bias loop

    #plotting
    if plotting
        plot_energies(Plotter, grid, data, solution, "Applied voltage Δu = $maxBias")
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution, "Applied voltage Δu = $maxBias")
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "Applied voltage Δu = $maxBias")
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
