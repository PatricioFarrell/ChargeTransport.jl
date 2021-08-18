#=
# 302: 1D PSC p-i-n device with discontinuous quasi Fermi potentials.
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

module Example333_PSC_pedot_mapi_pcbm_contqF

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize
using DelimitedFiles
using PyPlot


function main(;n = 14, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:dense)

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


    # grid
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
                                        h_pdoping/(1.2*δ), 
                                        h_pdoping/(1.2*δ), 
                                        tol=t)
    coord_i_g1              = geomspace(h_pdoping, 
                                        h_pdoping+h_intrinsic/k, 
                                        h_intrinsic/(7.1*δ), 
                                        h_intrinsic/(7.1*δ), 
                                        tol=t)
    coord_i_g2              = geomspace(h_pdoping+h_intrinsic/k, 
                                        h_pdoping+h_intrinsic,               
                                        h_intrinsic/(7.1*δ),    
                                        h_intrinsic/(7.1*δ), 
                                        tol=t)
    coord_n_g               = geomspace(h_pdoping+h_intrinsic,               
                                        h_pdoping+h_intrinsic+h_ndoping/4, 
                                        h_ndoping/(6.1*δ),   
                                        h_ndoping/(6.1*δ),      
                                        tol=t)
    coord_n_u               = collect(range(h_pdoping+h_intrinsic+h_ndoping/4, h_pdoping+h_intrinsic+h_ndoping, step=h_pdoping/(1.0*δ)))

    coord                   = glue(coord_p_u,coord_p_g,  tol=10*t) 
    coord                   = glue(coord,    coord_i_g1, tol=10*t)
    coord                   = glue(coord,    coord_i_g2, tol=10*t) 
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

    # recombination model
    bulk_recombination  = bulk_recomb_model_full

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
        
    # Auger recombination
    Auger               = 0.0

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
    data.F                               = [Boltzmann, Boltzmann, FermiDiracMinusOne]

    # Here the user can specify, if they assume continuous or discontinuous charge carriers. We note that for a surface recombination model,
    # we need to use discontinuous electron and hole quasi Fermi potentials.
    data.isContinuous[iphin]             = true
    data.isContinuous[iphip]             = true
    data.isContinuous[iphia]             = true

    # Following choices are possible for bulk_recombination_model:bulk_recomb_model_none, bulk_recomb_model_trap_assisted, bulk_recomb_radiative, bulk_recomb_full <: bulk_recombination_model 
    # The input iphin, iphip refers to the indices set by the user and needs to be specified
    data.bulk_recombination              = set_bulk_recombination(iphin = iphin, iphip = iphip, bulk_recombination_model = bulk_recombination)

    # Following choices are possible for boundary model: For contacts currently only ohmic_contact and schottky_contact are possible.
    # For inner boundaries we have interface_model_none, interface_model_surface_recombination, interface_model_ion_charge
    # (distinguish between left and right).
    data.boundary_type[bregionAcceptor]  = ohmic_contact  
    data.boundary_type[bregionJunction1] = interface_model_surface_reco_Cont
    data.boundary_type[bregionJunction2] = interface_model_surface_reco_Cont                   
    data.boundary_type[bregionDonor]     = ohmic_contact   

    # Here, the user gives information on which indices belong to ionic charge carriers and in which regions these charge carriers are present.
    # In this application ion vacancies only live in active perovskite layer
    data.enable_ion_vacancies            = enable_ion_vacancies(ionic_vacancies = [iphia], regions = [regionIntrinsic])
    
    # Following choices are possible for the flux_discretization scheme: ScharfetterGummel, ScharfetterGummel_Graded,
    # excessChemicalPotential, excessChemicalPotential_Graded, diffusionEnhanced, generalized_SG
    data.flux_approximation              = excessChemicalPotential
    
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
    params                                              = ChargeTransportParams(grid, numberOfCarriers)

    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1
    params.chargeNumbers[iphia]                         =  1

    # boundary region data
    params.bDensityOfStates[iphin, bregionDonor]        = Nc_d
    params.bDensityOfStates[iphip, bregionDonor]        = Nv_d

    params.bDensityOfStates[iphin, bregionAcceptor]     = Nc_a
    params.bDensityOfStates[iphip, bregionAcceptor]     = Nv_a

    params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_d
    params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_d

    params.bBandEdgeEnergy[iphin, bregionAcceptor]      = Ec_a
    params.bBandEdgeEnergy[iphip, bregionAcceptor]      = Ev_a

    # interior region data
    for ireg in 1:numberOfRegions

        params.dielectricConstant[ireg]                 = ε[ireg]

        # effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = NC[ireg]
        params.densityOfStates[iphip, ireg]             = NV[ireg]
        params.densityOfStates[iphia, ireg]             = NAnion[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = EC[ireg] 
        params.bandEdgeEnergy[iphip, ireg]              = EV[ireg] 
        params.bandEdgeEnergy[iphia, ireg]              = EA[ireg] 

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]
        params.mobility[iphia, ireg]                    = μa[ireg]

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
    params.doping[iphin,  regionDonor]                  = Nd 
    params.doping[iphip,  regionAcceptor]               = Na      
    params.doping[iphia,  regionIntrinsic]              = C0                 
    # boundary doping
    params.bDoping[iphin, bregionDonor]                 = Nd 
    params.bDoping[iphip, bregionAcceptor]              = Na    

    # Region dependent params is now a substruct of data which is again a substruct of the system and will be parsed 
    # in a next step.
    data.params                                         = params

    # in the last step, we initialize our system with previous data which is likewise dependent on the parameters. 
    # important that this is in the end, otherwise our VoronoiFVMSys is not dependent on the data we initialized
    # but rather on default data.
    ctsys                                               = ChargeTransportSystem(grid, data, unknown_storage=unknown_storage)

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
    set_ohmic_contact!(ctsys, bregionDonor, 0.0)

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
        
        t                         = tvalues[istep] # Actual time
        Δu                        = v0 + t*scanrate # Applied voltage 
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
            dfusion_grid            = 1.0e-2.*readdlm("data/Driftfusion-pedotpss-grid.dat")
            dfusion_psi             = readdlm("data/Driftfusion-pedotpss-Na-1p21e22-psi-t-end.dat")
            dfusion_phin            = readdlm("data/Driftfusion-pedotpss-Na-1p21e22-Efn-t-end.dat")
            dfusion_phip            = readdlm("data/Driftfusion-pedotpss-Na-1p21e22-Efp-t-end.dat")
        
            dfusion_psi_interface   = readdlm("data/Driftfusion-pedotpss-Na-1p21e22-interface-reco-psi-t-end.dat")
            dfusion_phin_interface  = readdlm("data/Driftfusion-pedotpss-Na-1p21e22-interface-reco-Efn-t-end.dat")
            dfusion_phip_interface  = readdlm("data/Driftfusion-pedotpss-Na-1p21e22-interface-reco-Efp-t-end.dat")
            ##########################
            plot_solution(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(Δu); \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$")
            PyPlot.plot(dfusion_grid', (dfusion_psi[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="black")
            PyPlot.plot(dfusion_grid', (-dfusion_phin[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="black")
            PyPlot.plot(dfusion_grid', (-dfusion_phip[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="black")
            #####
            PyPlot.plot(dfusion_grid', (dfusion_psi_interface[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="grey")
            PyPlot.plot(dfusion_grid', (-dfusion_phin_interface[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="grey")
            PyPlot.plot(dfusion_grid', (-dfusion_phip_interface[1,:]- 3.82*(ones(length(dfusion_grid)))), linewidth = 3, linestyle= ":", color="grey")
            PyPlot.axvline(h_pdoping, color="black", linestyle="solid")
            PyPlot.axvline(h_pdoping + h_intrinsic, color="black", linestyle="solid")
            #PyPlot.xlim(h_pdoping-1.2e-8, h_pdoping + h_intrinsic+1.2e-8)
            #PyPlot.ylim(-0.1, 1.25)
            
            ##########
            
        end
        #savefig("cont-qF.eps")

    end # time loop


    Plotter.figure()
    IV_measured         = readdlm("data/Driftfusion-IV-measurement-pcb-forward.dat")
    IV_Driftfusion_reco = readdlm("data/Driftfusion-pedotpss-Na-1p21e22-interface-reco-IV-forward.dat")
    IV_Driftfusion      = readdlm("data/Driftfusion-pedotpss-Na-1p21e22-IV-forward.dat")
        
    PyPlot.plot(biasValues, abs.(IV.*(cm)^2*1.0e3), label = "simulation",  linewidth= 3, linestyle="--", color="red")
    PyPlot.plot(IV_measured[:, 1], IV_measured[:, 2], label = "measurement",  linestyle="--", color = "black")
    PyPlot.plot(IV_Driftfusion[:, 1], IV_Driftfusion[:, 2].*1.0e3, label = "Driftfusion (without interface reco)", linewidth = 3, markersize="7", marker= "o", linestyle="--", color = "green")
    PyPlot.plot(IV_Driftfusion_reco[:, 1], IV_Driftfusion_reco[:, 2].*1.0e3, label = "Driftfusion (with interface reco)", markersize=7,marker= "x",  linestyle=":", color="blue")

    PyPlot.legend()
    PyPlot.title("Forward; \$ E_a =\$$(textEa)eV;  \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$ ")
    Plotter.ylabel("total current [mA]") # 
    Plotter.xlabel("Applied Voltage [V]")
    PyPlot.ylim(0.0, 0.006*1.0e3)
    PyPlot.tight_layout()
    #savefig("cont-qF-IV.eps")



    if test == false
        println("*** done\n")
    end


    testval = VoronoiFVM.norm(ctsys.fvmsys, solution, 2)
    return testval


end #  main

function test()
    testval = 97.57205176140376
    main(test = true, unknown_storage=:dense) ≈ testval #&& main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
