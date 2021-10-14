#=
# 102: 1D PSC p-i-n device without mobile ions.
([source code](SOURCE_URL))

Simulating a three layer PSC device SiO2| MAPI | SiO2 without mobile ions
and in stationary state (i.e. no time-dependency). This means we consider
hetero interfaces.
The simulations are performed out of equilibrium and with
abrupt interfaces. For simplicity, the generation is off.

This simulation coincides with the one made in Section 4.3
of Calado et al. (https://arxiv.org/abs/2009.04384).
The paramters can be found in Table S.13 or slightly modified than the one
in the publication here:
https://github.com/barnesgroupICL/Driftfusion/blob/Methods-IonMonger-Comparison/Input_files/IonMonger_default_bulk.csv
=#

module Example103_PSC

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize

function main(;n = 8, Plotter = nothing, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # region numbers
    regionDonor             = 1                           # n doped region
    regionIntrinsic         = 2                           # intrinsic region
    regionAcceptor          = 3                           # p doped region
    regions                 = [regionDonor, regionIntrinsic, regionAcceptor]
    numberOfRegions         = length(regions)

    # boundary region numbers
    bregionDonor            = 1
    bregionAcceptor         = 2
    bregions                = [bregionAcceptor, bregionDonor]
    numberOfBoundaryRegions = length(bregions)

    # grid
    # NB: Using geomspace to create uniform mesh is not a good idea. It may create virtual duplicates at boundaries.
    h_ndoping               = 9.90e-6 * cm 
    h_intrinsic             = 4.00e-5 * cm + 2.0e-7 * cm # add 2.e-7 cm to this layer for agreement with grid of Driftfusion
    h_pdoping               = 1.99e-5 * cm

    x0                      = 0.0 * cm 
    δ                       = 2*n        # the larger, the finer the mesh
    t                       = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k                       = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_n_u               = collect(range(x0, h_ndoping/2, step=h_ndoping/(0.8*δ)))
    coord_n_g               = geomspace(h_ndoping/2, 
                                        h_ndoping, 
                                        h_ndoping/(1.1*δ), 
                                        h_ndoping/(1.1*δ), 
                                        tol=t)
    coord_i_g1              = geomspace(h_ndoping, 
                                        h_ndoping+h_intrinsic/k, 
                                        h_intrinsic/(2.8*δ), 
                                        h_intrinsic/(2.8*δ), 
                                        tol=t)
    coord_i_g2              = geomspace(h_ndoping+h_intrinsic/k, 
                                        h_ndoping+h_intrinsic,               
                                        h_intrinsic/(2.8*δ),    
                                        h_intrinsic/(2.8*δ), 
                                        tol=t)
    coord_p_g               = geomspace(h_ndoping+h_intrinsic,               
                                        h_ndoping+h_intrinsic+h_pdoping/2, 
                                        h_pdoping/(1.6*δ),   
                                        h_pdoping/(1.6*δ),      
                                        tol=t)
    coord_p_u               = collect(range(h_ndoping+h_intrinsic+h_pdoping/2, h_ndoping+h_intrinsic+h_pdoping, step=h_pdoping/(1.3*δ)))

    coord                   = glue(coord_n_u, coord_n_g,  tol=10*t) 
    coord                   = glue(coord,     coord_i_g1, tol=10*t)
    coord                   = glue(coord,     coord_i_g2, tol=10*t) 
    coord                   = glue(coord,     coord_p_g,  tol=10*t)
    coord                   = glue(coord,     coord_p_u,  tol=10*t)
    grid                    = simplexgrid(coord)

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],                [h_ndoping],                           regionDonor)     # n-doped region   = 1
    cellmask!(grid, [h_ndoping],               [h_ndoping + h_intrinsic],             regionIntrinsic) # intrinsic region = 2
    cellmask!(grid, [h_ndoping + h_intrinsic], [h_ndoping + h_intrinsic + h_pdoping], regionAcceptor)  # p-doped region   = 3

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

    # set indices of the quasi Fermi potentials
    iphin             = 2 # electron quasi Fermi potential
    iphip             = 1 # hole quasi Fermi potential
    numberOfCarriers  = 2 

    ##########      physical data      ##########
    # temperature
    T                 = 300.0                 *  K

    # band edge energies    
    Ec_d              = -4.0                  *  eV 
    Ev_d              = -6.0                  *  eV 

    Ec_i              = -3.7                  *  eV 
    Ev_i              = -5.4                  *  eV 

    Ec_a              = -3.1                  *  eV 
    Ev_a              = -5.1                  *  eV 

    EC                = [Ec_d, Ec_i, Ec_a] 
    EV                = [Ev_d, Ev_i, Ev_a] 

    # effective densities of state
    Nc_d               = 5.0e19                / (cm^3)
    Nv_d               = 5.0e19                / (cm^3)

    Nc_i               = 8.1e18                / (cm^3)
    Nv_i               = 5.8e18                / (cm^3)

    Nc_a               = 5.0e19                / (cm^3)
    Nv_a               = 5.0e19                / (cm^3)

    NC                 = [Nc_d, Nc_i, Nc_a]
    NV                 = [Nv_d, Nv_i, Nv_a]

    # mobilities 
    μn_d               = 3.89                  * (cm^2) / (V * s)  
    μp_d               = 3.89                  * (cm^2) / (V * s)  

    μn_i               = 6.62e1                * (cm^2) / (V * s)  
    μp_i               = 6.62e1                * (cm^2) / (V * s)

    μn_a               = 3.89e-1               * (cm^2) / (V * s) 
    μp_a               = 3.89e-1               * (cm^2) / (V * s) 

    μn                 = [μn_d, μn_i, μn_a] 
    μp                 = [μp_d, μp_i, μp_a] 

    # relative dielectric permittivity  

    ε_d                = 10.0                  *  1.0  
    ε_i                = 24.1                  *  1.0 
    ε_a                = 3.0                   *  1.0 

    ε                  = [ε_d, ε_i, ε_a] 

    # recombination model
    bulk_recombination = bulk_recomb_model_full

    # radiative recombination
    r0_d               = 0.0e+0               * cm^3 / s 
    r0_i               = 1.0e-12              * cm^3 / s  
    r0_a               = 0.0e+0               * cm^3 / s
        
    r0                 = [r0_d, r0_i, r0_a]
        
    # life times and trap densities 
    τn_d               = 1.0e100              * s 
    τp_d               = 1.0e100              * s
        
    τn_i               = 3.0e-10              * s
    τp_i               = 3.0e-8               * s
    τn_a               = τn_d
    τp_a               = τp_d
        
    τn                 = [τn_d, τn_i, τn_a]
    τp                 = [τp_d, τp_i, τp_a]
        
    # SRH trap energies (needed for calculation of trap_density! (SRH))
    Ei_d               = -5.0                 * eV   
    Ei_i               = -4.55                * eV   
    Ei_a               = -4.1                 * eV   

    EI                 = [Ei_d, Ei_i, Ei_a]
        
    # Auger recombination
    Auger              = 0.0

    # doping 
    Nd                  =   1.03e18             / (cm^3) 
    Na                  =   1.03e18             / (cm^3) 
    Ni_acceptor         =   8.32e7              / (cm^3) 

    # contact voltages: we impose an applied voltage only on one boundary.
    # At the other boundary the applied voltage is zero.
    voltageAcceptor     =  1.2                  * V 

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

    # Here the user can specify, if they assume continuous or discontinuous charge carriers.
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
    data.flux_approximation             = excessChemicalPotential

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define ChargeTransportParams and fill in physical parameters")
    end
    ################################################################################

    params                                              = ChargeTransportParams(grid, numberOfCarriers)

    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1

    # boundary region data
    params.bDensityOfStates[iphin, bregionDonor]    = Nc_d
    params.bDensityOfStates[iphip, bregionDonor]    = Nv_d

    params.bDensityOfStates[iphin, bregionAcceptor] = Nc_a
    params.bDensityOfStates[iphip, bregionAcceptor] = Nv_a

    params.bBandEdgeEnergy[iphin, bregionDonor]     = Ec_d
    params.bBandEdgeEnergy[iphip, bregionDonor]     = Ev_d 

    params.bBandEdgeEnergy[iphin, bregionAcceptor]  = Ec_a 
    params.bBandEdgeEnergy[iphip, bregionAcceptor]  = Ev_a 

    # interior region data
    for ireg in 1:numberOfRegions

        params.dielectricConstant[ireg]                 = ε[ireg]

        # effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = NC[ireg]
        params.densityOfStates[iphip, ireg]             = NV[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = EC[ireg] 
        params.bandEdgeEnergy[iphip, ireg]              = EV[ireg] 

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]

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

    # in the last step, we initialize our system with previous data which is likewise dependent on the parameters. 
    # important that this is in the end, otherwise our VoronoiFVMSys is not dependent on the data we initialized
    # but rather on default data.
    ctsys                                           = ChargeTransportSystem(grid, data, unknown_storage=unknown_storage)

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
        println("Define outerior boundary conditions")
    end
    ################################################################################

    # set ohmic contacts for each charge carrier at all outerior boundaries. First, 
    # we compute equilibrium solutions. Hence the boundary values at the ohmic contacts
    # are zero.
    set_ohmic_contact!(ctsys, bregionAcceptor, 0.0)
    set_ohmic_contact!(ctsys, bregionDonor, 0.0)

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


    control.damp_initial      = 0.8
    control.damp_growth       = 1.61 # >= 1
    control.max_round         = 5

    # initialize solution and starting vectors
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

    control.damp_initial                             = 0.6
    control.damp_growth                              = 1.21 # >= 1
    control.max_round                                = 7

    maxBias    = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 13)

    for Δu in biasValues
        if verbose
            println("Bias value: Δu = $(Δu) (no illumination)")
        end

        # set non equilibrium boundary conditions
        set_ohmic_contact!(ctsys, bregionAcceptor, Δu)

        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        initialGuess .= solution

    end # bias loop

    # plotting
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

    testval = solution[3, 20] # the 3 corresponds to the index of ipsi
    return testval

end #  main

function test()
    testval = -4.053036050777252
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
