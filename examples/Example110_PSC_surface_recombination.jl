#=
# 110: 1D PSC p-i-n device with surface recombination.
([source code](SOURCE_URL))

Simulating a three layer PSC device Pedot| MAPI | PCBM.
The simulations are performed out of equilibrium, time-dependent
and with abrupt interfaces. 
A linear I-V measurement protocol is included and the corresponding
solution vectors after the scan protocol can be depicted.

The paramters are from Calado et al. and can be found here:
https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/pedotpss_mapi_pcbm.csv.
(with adjustments on layer lengths)
=#

module Example110_PSC_surface_recombination

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize
using PyPlot

function main(;n = 6, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:dense)

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
    bregionJunction1        = 3
    bregionJunction2        = 4
    bregions                = [bregionAcceptor, bregionDonor, bregionJunction1, bregionJunction2]
    numberOfBoundaryRegions = length(bregions)

    # grid (the nearer to interface, the finer)
    h_pdoping               = 3.00e-6 * cm + 1.0e-7 * cm
    h_intrinsic             = 3.00e-5 * cm 
    h_ndoping               = 8.50e-6 * cm + 1.0e-7 * cm

    x0                      = 0.0 * cm 
    δ                       = 4*n        # the larger, the finer the mesh
    t                       = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k                       = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_p_u               = collect(range(x0, h_pdoping/2, step=h_pdoping/(0.9*δ)))
    coord_p_g               = geomspace(h_pdoping/2, 
                                        h_pdoping, 
                                        h_pdoping/(1.0*δ), 
                                        h_pdoping/(1.5*δ), 
                                        tol=t)
    coord_i_g1              = geomspace(h_pdoping, 
                                        h_pdoping+h_intrinsic/k, 
                                        h_intrinsic/(6.1*δ), 
                                        h_intrinsic/(3.1*δ), 
                                        tol=t)
    coord_i_g2              = geomspace(h_pdoping+h_intrinsic/k, 
                                        h_pdoping+h_intrinsic,               
                                        h_intrinsic/(3.1*δ),    
                                        h_intrinsic/(6.1*δ), 
                                        tol=t)
    coord_n_g               = geomspace(h_pdoping+h_intrinsic,               
                                        h_pdoping+h_intrinsic+h_ndoping/2, 
                                        h_ndoping/(3.0*δ),   
                                        h_ndoping/(1.0*δ),      
                                        tol=t)
    coord_n_u               = collect(range(h_pdoping+h_intrinsic+h_ndoping/2, h_pdoping+h_intrinsic+h_ndoping, step=h_pdoping/(0.8*δ)))

    coord                   = glue(coord_p_u,coord_p_g,  tol=10*t) 
    icoord_p                = length(coord)
    coord                   = glue(coord,    coord_i_g1, tol=10*t)
    coord                   = glue(coord,    coord_i_g2, tol=10*t) 
    icoord_pi               = length(coord)
    coord                   = glue(coord,    coord_n_g,  tol=10*t)
    coord                   = glue(coord,    coord_n_u,  tol=10*t)
    grid                    = ExtendableGrids.simplexgrid(coord)
    numberOfNodes           = length(coord)

    # # cellmask! for defining the subregions and assigning region number (doping profiles do not intersect)
    cellmask!(grid, [0.0 * μm],                 [h_pdoping],                           regionAcceptor, tol = 1.0e-18)   # p-doped region   = 1
    cellmask!(grid, [h_pdoping],                [h_pdoping + h_intrinsic],             regionIntrinsic, tol = 1.0e-18)  # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic],  [h_pdoping + h_intrinsic + h_ndoping], regionDonor, tol = 1.0e-18)      # n-doped region   = 3

    # bfacemask! for ``active'' boundary regions, i.e. internal interfaces. On the outer boudary regions, the 
    # conditions will be formulated later
    bfacemask!(grid, [h_pdoping],               [h_pdoping],                           bregionJunction1)  # first  inner interface
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic],             bregionJunction2)  # second inner interface

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

    # charge carriers (by construction the last index is automatically assigned to the electric potential, i.e.
    # ipsi = length(ichargeCarriers) + 1)

    iphin               = 1 # electron quasi Fermi potential
    iphip               = 2 # hole quasi Fermi potential
    iphia               = 3 # anion vacancy quasi Fermi potential

    ichargeCarriers     = [iphin, iphip, iphia]   # this is an Array of indices of the charge carriers 
    numberOfCarriers    = length(ichargeCarriers) # electrons, holes and anion vacancies

    # temperature
    T                   =  300.0                *  K

    # band edge energies    
    Ec_a                = -3.0                  *  eV 
    Ev_a                = -5.1                  *  eV 

    Ec_i                = -3.8                  *  eV 
    Ev_i                = -5.4                  *  eV 

    Ec_d                = -3.8                  *  eV 
    Ev_d                = -6.2                  *  eV 

    EC                  = [Ec_a, Ec_i, Ec_d] 
    EV                  = [Ev_a, Ev_i, Ev_d]
    
    # effective densities of state
    Nc_a                = 1.0e20                / (cm^3)
    Nv_a                = 1.0e20                / (cm^3)

    Nc_i                = 1.0e19                / (cm^3)
    Nv_i                = 1.0e19                / (cm^3)

    ###################### adjust Na, Ea here #####################
    Nanion              = 1.21e22               / (cm^3)
    Ea_i                = -5.175                *  eV 

    # for the labels in the figures
    textEa              = Ea_i                 ./  eV
    textNa              = Nanion               .* (cm^3)
    ###################### adjust Na, Ea here #####################

    EA                  = [0.0,  Ea_i,  0.0]

    Nc_d                = 1.0e19                / (cm^3)
    Nv_d                = 1.0e19                / (cm^3)

    NC                  = [Nc_a, Nc_i, Nc_d]
    NV                  = [Nv_a, Nv_i, Nv_d]
    NAnion              = [0.0,  Nanion, 0.0]

    # mobilities 
    μn_a                = 0.1                   * (cm^2) / (V * s)  
    μp_a                = 0.1                   * (cm^2) / (V * s)  

    μn_i                = 2.00e1                * (cm^2) / (V * s)  
    μp_i                = 2.00e1                * (cm^2) / (V * s)
    μa_i                = 1.00e-10              * (cm^2) / (V * s)

    μn_d                = 1.0e-3                * (cm^2) / (V * s) 
    μp_d                = 1.0e-3                * (cm^2) / (V * s) 

    μn                  = [μn_a, μn_i, μn_d] 
    μp                  = [μp_a, μp_i, μp_d] 
    μa                  = [0.0,  μa_i, 0.0 ] 

    # relative dielectric permittivity  
    ε_a                 = 4.0                   *  1.0  
    ε_i                 = 23.0                  *  1.0 
    ε_d                 = 3.0                   *  1.0 

    ε                   = [ε_a, ε_i, ε_d] 

    # radiative recombination
    r0_a                = 6.3e-11               * cm^3 / s 
    r0_i                = 3.6e-12               * cm^3 / s  
    r0_d                = 6.8e-11               * cm^3 / s
        
    r0                  = [r0_a, r0_i, r0_d]
        
    # life times and trap densities 
    τn_a                = 1.0e-6                * s 
    τp_a                = 1.0e-6                * s
        
    τn_i                = 1.0e-7                * s
    τp_i                = 1.0e-7                * s
    τn_d                = τn_a
    τp_d                = τp_a
        
    τn                  = [τn_a, τn_i, τn_d]
    τp                  = [τp_a, τp_i, τp_d]
        
    # SRH trap energies (needed for calculation of trap_density! (SRH))
    Ei_a                = -4.05                 * eV   
    Ei_i                = -4.60                 * eV   
    Ei_d                = -5.00                 * eV   

    EI                  = [Ei_a, Ei_i, Ei_d]
        
    # doping
    Nd                  = 2.089649130192123e17  / (cm^3) 
    Na                  = 4.529587947185444e18  / (cm^3) 
    C0                  = 1.0e18                / (cm^3) 

    # contact voltages: we impose an applied voltage only on one outer boundary.
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
    data                                 = ChargeTransportData(grid, numberOfCarriers)

    ##############################################################################
    ####     declare here all necessary information concerning the model      ####
    ##############################################################################

    # Following variable declares, if we want to solve a stationary or transient problem.
    # Choose between: model_transient, model_stationary
    data.model_type                      = model_transient

    # Following choices are possible for F: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                               = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    # Here, we need to specify which numbers are associated with electron and hole quasi Fermi potential. Further, the desired recombination 
    # processes can be chosen here. Note that, if you choose a SRH recombination you can further specify a transient SRH recombination by 
    # the method enable_traps! and adjusting the model_type. Otherwise, by default we use the stationary model for this type of recombination.
    data.bulk_recombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip, 
                                                                  bulk_recomb_Auger = false,
                                                                  bulk_recomb_radiative = true,
                                                                  bulk_recomb_SRH = true)

    # Following choices are possible for boundary model: For contacts currently only ohmic_contact and schottky_contact are possible.
    # For inner boundaries we have interface_model_none, interface_model_surface_recombination.
    data.boundary_type[bregionAcceptor]  = ohmic_contact  
    data.boundary_type[bregionJunction1] = interface_model_surface_recombination
    data.boundary_type[bregionJunction2] = interface_model_surface_recombination                   
    data.boundary_type[bregionDonor]     = ohmic_contact   

    # Here, the user gives information on which indices belong to ionic charge carriers and in which regions these charge carriers are present.
    # In this application ion vacancies only live in active perovskite layer
    data.enable_ionic_carriers            = enable_ionic_carriers(ionic_carriers = [iphia], regions = [regionIntrinsic])
    
    # Following choices are possible for the flux_discretization scheme: scharfetter_gummel, scharfetter_gummel_graded,
    # excess_chemical_potential, excess_chemical_potential_graded, diffusion_enhanced, generalized_sg
    data.flux_approximation              = excess_chemical_potential
    
    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define ChargeTransportParams and fill in physical parameters")
    end
    ################################################################################

    # Params is a struct which contains all necessary physical parameters. If one wants to simulate
    # space-dependent variables, one additionally needs to generate a ParamsNodal struct, see Example102 or Example104.
    params                                                       = ChargeTransportParams(grid, numberOfCarriers)

    params.temperature                                           = T
    params.UT                                                    = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                                  = -1
    params.chargeNumbers[iphip]                                  =  1
    params.chargeNumbers[iphia]                                  =  1

    # interior region data
    for ireg in 1:numberOfRegions

        params.dielectricConstant[ireg]                          = ε[ireg]

        # effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]                      = NC[ireg]
        params.densityOfStates[iphip, ireg]                      = NV[ireg]
        params.densityOfStates[iphia, ireg]                      = NAnion[ireg]

        params.bandEdgeEnergy[iphin, ireg]                       = EC[ireg] 
        params.bandEdgeEnergy[iphip, ireg]                       = EV[ireg] 
        params.bandEdgeEnergy[iphia, ireg]                       = EA[ireg] 

        params.mobility[iphin, ireg]                             = μn[ireg]
        params.mobility[iphip, ireg]                             = μp[ireg]
        params.mobility[iphia, ireg]                             = μa[ireg]

        # recombination parameters
        params.recombinationRadiative[ireg]                      = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]             = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]             = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg]          = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg]          = trap_density!(iphip, ireg, data, EI[ireg])
    end

    ## outer boundary region data
    params.bDensityOfStates[iphin, bregionAcceptor]              = Nc_a
    params.bDensityOfStates[iphip, bregionAcceptor]              = Nv_a

    params.bDensityOfStates[iphin, bregionDonor]                 = Nc_d
    params.bDensityOfStates[iphip, bregionDonor]                 = Nv_d

    params.bBandEdgeEnergy[iphin, bregionAcceptor]               = Ec_a
    params.bBandEdgeEnergy[iphip, bregionAcceptor]               = Ev_a

    params.bBandEdgeEnergy[iphin, bregionDonor]                  = Ec_d
    params.bBandEdgeEnergy[iphip, bregionDonor]                  = Ev_d

    ##############################################################
    ## inner boundary region data
    params.bDensityOfStates[iphin, bregionJunction1]             = Nc_i
    params.bDensityOfStates[iphip, bregionJunction1]             = Nv_i

    params.bDensityOfStates[iphin, bregionJunction2]             = Nc_i
    params.bDensityOfStates[iphip, bregionJunction2]             = Nv_i

    params.bBandEdgeEnergy[iphin, bregionJunction1]              = Ec_i 
    params.bBandEdgeEnergy[iphip, bregionJunction1]              = Ev_i 

    params.bBandEdgeEnergy[iphin, bregionJunction2]              = Ec_i 
    params.bBandEdgeEnergy[iphip, bregionJunction2]              = Ev_i 

    # for surface recombination
    params.recombinationSRHvelocity[iphin, bregionJunction1]     = 1.0e1  * cm / s
    params.recombinationSRHvelocity[iphip, bregionJunction1]     = 1.0e5  * cm / s

    params.recombinationSRHvelocity[iphin, bregionJunction2]     = 1.0e7  * cm / s
    params.recombinationSRHvelocity[iphip, bregionJunction2]     = 1.0e1  * cm / s

    params.bRecombinationSRHTrapDensity[iphin, bregionJunction1] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJunction1] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    params.bRecombinationSRHTrapDensity[iphin, bregionJunction2] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJunction2] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    ##############################################################
    
    # interior doping
    params.doping[iphin,  regionDonor]                           = Nd 
    params.doping[iphip,  regionAcceptor]                        = Na      
    params.doping[iphia,  regionIntrinsic]                       = C0                 
    # boundary doping
    params.bDoping[iphin, bregionDonor]                          = Nd 
    params.bDoping[iphip, bregionAcceptor]                       = Na    

    # Region dependent params is now a substruct of data which is again a substruct of the system and will be parsed 
    # in a next step.
    data.params                                                  = params

    # in the last step, we initialize our system with previous data which is likewise dependent on the parameters. 
    # important that this is in the end, otherwise our VoronoiFVMSys is not dependent on the data we initialized
    # but rather on default data.
    ctsys                                                        = ChargeTransportSystem(grid, data, unknown_storage=unknown_storage)

    # print all params stored in ctsys.data.params
    if test == false
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

    # set ohmic contacts for electrons and holes at all outerior boundaries. Here, 
    # we compute equilibrium solutions, i.e. the boundary values at the ohmic contacts
    # are zero.
    set_ohmic_contact!(ctsys, bregionAcceptor, 0.0)
    set_ohmic_contact!(ctsys, bregionDonor,    0.0)

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
    control.tol_absolute      = 1.0e-10
    control.tol_relative      = 1.0e-10
    control.handle_exceptions = true
    control.tol_round         = 1.0e-10
    control.max_round         = 5

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    control.damp_initial      = 0.05
    control.damp_growth       = 1.21 # >= 1
    control.max_round         = 5

    # initialize solution and starting vectors
    initialGuess              = unknowns(ctsys)
    solution                  = unknowns(ctsys)

    solution                  = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

    initialGuess             .= solution 

    if test == false
        println("*** done\n")
    end


    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("I-V Measurement Loop")
    end
    ################################################################################

    ctsys.data.calculation_type   = outOfEquilibrium

    control.damp_initial          = 0.5
    control.damp_growth           = 1.21
    control.max_round             = 5

    # there are different way to control timestepping
    # Here we assume these primary data
    scanrate                      = 0.04 * V/s
    ntsteps                       = 101
    vend                          = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
    v0                            = 0.0

    # The end time then is calculated here:
    tend                          = vend/scanrate

    # with fixed timestep sizes we can calculate the times
    # a priori
    tvalues                       = range(0, stop = tend, length = ntsteps)

    # for saving I-V data
    IV                            = zeros(0) # for IV values
    biasValues                    = zeros(0) # for bias values

    for istep = 2:ntsteps
        
        t                         = tvalues[istep]       # Actual time
        Δu                        = v0 + t*scanrate      # Applied voltage 
        Δt                        = t - tvalues[istep-1] # Time step size
        
        # Apply new voltage (set non-equilibrium values)
        set_ohmic_contact!(ctsys, bregionAcceptor, Δu)
        
        if test == false
            println("time value: Δt = $(t)")
        end

        # Solve time step problems with timestep Δt. initialGuess plays the role of the solution
        # from last timestep
        solve!(solution, initialGuess, ctsys, control = control, tstep = Δt)
        
        initialGuess .= solution

        # get I-V data
        current = get_current_val(ctsys, solution, initialGuess, Δt)

        push!(IV, current)
        push!(biasValues, Δu)

        if plotting
            plot_solution(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(Δu); \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$")
        end

    end # time loop

    #res = [biasValues, IV]

    if test == false
        println("*** done\n")
    end

    testval = VoronoiFVM.norm(ctsys.fvmsys, solution, 2)
    return testval


end #  main

function test()
    testval = 54.913949735570306
    main(test = true, unknown_storage=:dense) ≈ testval #&& main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
