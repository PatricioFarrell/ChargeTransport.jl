
ENV["LC_NUMERIC"]="C"

module L_shape_linear_IV_without_ions

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using GridVisualize
using DelimitedFiles
using PyPlot
using SimplexGridFactory
using Triangulate


function main(;Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:dense)

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
    bregionNoFlux           = 5
    bregions                = [bregionAcceptor, bregionDonor, bregionJunction1, bregionJunction2]
    numberOfBoundaryRegions = length(bregions)


    # grid
    h_pdoping               = 3.00e-6 * cm + 1.0e-7 * cm
    h_intrinsic             = 3.00e-5 * cm 
    h_ndoping               = 8.50e-6 * cm + 1.0e-7 * cm
    height                  = 1.00e-5 * cm
     
    function unsuitable(x1,y1,x2,y2,x3,y3,area)
        bary_x=(x1+x2+x3)/3.0
        bary_y=(y1+y2+y3)/3.0
        dx=bary_x-refinement_center[1]
        dy=bary_y-refinement_center[2]
        qdist=dx^2+dy^2
        area>0.1*max(9.0e-18,qdist)
        # min(area1, area2) > 0.1*max(4.e-17,qdist)
    end

    # function unsuitable2(x1,y1,x2,y2,x3,y3,area)

    #     bary_x=(x1+x2+x3)/3.0
    #     bary_y=(y1+y2+y3)/3.0
    #     dx=bary_x-refinement_center2[1]
    #     dy=bary_y-refinement_center2[2]
    #     qdist=dx^2+dy^2
    #     area>0.1*max(5.e-18,qdist)
    # end
    
    b                = SimplexGridBuilder(Generator=Triangulate)

    # specify boundary nodes
    length_0    = point!(b, 0.0, 0.0)
    length_p    = point!(b, h_pdoping, 0.0)
    length_pi   = point!(b, h_pdoping + h_intrinsic, 0.0)
    length_pin  = point!(b, h_pdoping + h_intrinsic + h_ndoping, 0.0)
    height_0    = point!(b, 0.0, height)
    height_p    = point!(b, h_pdoping, height)
    
    # for L shape
    height_pi12 = point!(b, h_pdoping + h_intrinsic/2, height)
    height_pi2  = point!(b, h_pdoping + h_intrinsic/2, height/2)
    height_pi   = point!(b, h_pdoping + h_intrinsic, height/2)
    height_pin  = point!(b, h_pdoping + h_intrinsic + h_ndoping, height/2)
    
    ## specify boundary regions
    # metal interface
    facetregion!(b, bregionAcceptor)
    facet!(b, height_0, height_p)
    facetregion!(b, bregionDonor)
    facet!(b, length_pin, height_pin) 
          
    # no flux
    facetregion!(b, bregionNoFlux) # this is the old facet where previously the metal interface was defined.
    facet!(b, length_0, height_0)

    facetregion!(b, bregionNoFlux)
    facet!(b, length_0, length_pin)    

    facetregion!(b, bregionNoFlux)
    facet!(b, height_0, height_pi12)


    facetregion!(b, bregionNoFlux)
    facet!(b, height_pi12, height_pi2)
    facetregion!(b, bregionNoFlux)
    facet!(b, height_pi2, height_pin)
  
    # inner interface
    facetregion!(b, bregionJunction1)
    facet!(b, length_p, height_p)
    facetregion!(b, bregionJunction2)
    facet!(b, length_pi, height_pi)

    refinement_center = [h_pdoping + h_intrinsic/2, height/2]
    #refinement_center2 = [h_pdoping/2, height]
    # Activate unsuitable callback
    options!(b,unsuitable=unsuitable)

    #options!(b,unsuitable=unsuitable2)
    
    # cell regions
    cellregion!(b, regionAcceptor)
    regionpoint!(b, h_pdoping-1.0e-6*cm, height/2-1.0e-6*cm) 
    cellregion!(b,regionIntrinsic)
    regionpoint!(b, h_pdoping + h_intrinsic -1.0e-6*cm, height/2-1.0e-6*cm) 
    cellregion!(b,regionDonor)
    regionpoint!(b, h_pdoping + h_intrinsic + h_ndoping -1.0e-6*cm, height/2-1.0e-6*cm) 

    options!(b,maxvolume=9.0e-18)

    grid            = simplexgrid(b)
    numberOfNodes   = size(grid[Coordinates])[2]
    X               = grid[Coordinates][1,:]
    Y               = grid[Coordinates][2,:]

    if plotting
        PyPlot.figure()
        GridVisualize.gridplot(grid, Plotter= Plotter, resolution=(600,400),linewidth=0.6, legend=:lt) #, legend=:lt
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.title("Grid")
        Plotter.tight_layout()
        #savefig("grid-l-shape.eps")
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

    ichargeCarriers     = [iphin, iphip]   # this is an Array of indices of the charge carriers 
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

    Nc_d                = 1.0e19                / (cm^3)
    Nv_d                = 1.0e19                / (cm^3)

    NC                  = [Nc_a, Nc_i, Nc_d]
    NV                  = [Nv_a, Nv_i, Nv_d]

    # mobilities 
    μn_a                = 0.1                   * (cm^2) / (V * s)  
    μp_a                = 0.1                   * (cm^2) / (V * s)  

    μn_i                = 2.00e1                * (cm^2) / (V * s)  
    μp_i                = 2.00e1                * (cm^2) / (V * s)

    μn_d                = 1.0e-3                * (cm^2) / (V * s) 
    μp_d                = 1.0e-3                * (cm^2) / (V * s) 

    μn                  = [μn_a, μn_i, μn_d] 
    μp                  = [μp_a, μp_i, μp_d] 

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
    data.F                               = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA]

    # Following choices are possible for bulk_recombination_model:bulk_recomb_model_none, bulk_recomb_model_trap_assisted, bulk_recomb_radiative, bulk_recomb_full <: bulk_recombination_model 
    # The input iphin, iphip refers to the indices set by the user and needs to be specified
    data.bulk_recombination              = set_bulk_recombination(iphin = iphin, iphip = iphip, bulk_recombination_model = bulk_recombination)

    # Following choices are possible for boundary model: For contacts currently only ohmic_contact and schottky_contact are possible.
    # For inner boundaries we have interface_model_none, interface_model_surface_recombination, interface_model_ion_charge
    # (distinguish between left and right).
    data.boundary_type[bregionAcceptor]  = ohmic_contact  
    data.boundary_type[bregionJunction1] = interface_model_surface_recombination_and_tangential_flux#interface_model_surface_recombination#
    data.boundary_type[bregionJunction2] = interface_model_surface_recombination_and_tangential_flux##interface_model_surface_recombination#                  
    data.boundary_type[bregionDonor]     = ohmic_contact   
    
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

    # interior region data
    for ireg in 1:numberOfRegions

        params.dielectricConstant[ireg]                 = ε[ireg]

        # effective dos, band edge energy and mobilities
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

    ## outer boundary region data
    params.bDensityOfStates[iphin, bregionAcceptor]     = Nc_a
    params.bDensityOfStates[iphip, bregionAcceptor]     = Nv_a

    params.bDensityOfStates[iphin, bregionDonor]        = Nc_d
    params.bDensityOfStates[iphip, bregionDonor]        = Nv_d

    params.bBandEdgeEnergy[iphin, bregionAcceptor]      = Ec_a
    params.bBandEdgeEnergy[iphip, bregionAcceptor]      = Ev_a

    params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_d
    params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_d

    ##############################################################
    ## inner boundary region data

    # mobility needed for bflux! 
    prefactor_mobility                                           = 1.0
    params.bMobility[iphin, bregionJunction1]                    = prefactor_mobility * μn_i
    params.bMobility[iphip, bregionJunction1]                    = prefactor_mobility * μp_i

    params.bMobility[iphin, bregionJunction2]                    = prefactor_mobility * μn_i
    params.bMobility[iphip, bregionJunction2]                    = prefactor_mobility * μp_i

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
    
    # interior doping
    params.doping[iphin,  regionDonor]                  = Nd 
    params.doping[iphip,  regionAcceptor]               = Na      
              
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
    initialGuess              = unknowns(ctsys)
    solution                  = unknowns(ctsys)

    solution                  = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

    initialGuess             .= solution 
                  
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
    ntsteps                       = 51
    vend                          = 1.0#voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
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
       
    end # time loop

    #writedlm("sol-2D-tangential-flux-off-surface-reco-off.dat", [X Y solution'])

    res = [biasValues IV]
    #writedlm("IV-2D-tangential-flux-off-surface-reco-off.dat", res)
     
    function tran32!(a,b)
        a[1] = b[2]
    end

    p                 = GridVisualizer(;Plotter=Plotter,layout=(2,2),clear=true,resolution=(800,500))

    bgrid1            = subgrid(grid, [1], boundary = true)
    bgrid2            = subgrid(grid, [2], boundary = true, transform = tran32!)
    bgrid3            = subgrid(grid, [3], boundary = true, transform = tran32!)
    bgrid4            = subgrid(grid, [4], boundary = true, transform = tran32!)

    U_iphin_left      = view(solution[iphin,:], bgrid1)
    U_iphin_right     = view(solution[iphin,:], bgrid2)
    U_iphin_junction1 = view(solution[iphin,:], bgrid3)
    U_iphin_junction2 = view(solution[iphin,:], bgrid4)
    
    scalarplot!(p[1,1],bgrid1, U_iphin_left,      show=true, cellwise = true, xlabel = "x-axis [m]", ylabel = "voltage [V]", title = "Left contact")
    scalarplot!(p[2,1],bgrid3, U_iphin_junction1, show=true, cellwise = true, xlabel = "y-axis [m]", ylabel = "voltage [V]", title = "Left junction")
    scalarplot!(p[2,2],bgrid4, U_iphin_junction2, show=true, cellwise = true, xlabel = "y-axis [m]", ylabel = "voltage [V]", title = "Right junction")
    scalarplot!(p[1,2],bgrid2, U_iphin_right,     show=true, cellwise = true, xlabel = "x-axis [m]", ylabel = "voltage [V]", title = "Right contact")
    Plotter.tight_layout()
    #savefig("boundary-iphin-tangential-flux-on-prefactor-mobility-$(prefactor_mobility)-surface-reco-on.eps")

    if plotting

        lengthNaN = 61
        xx        = collect( range(1.83e-5 * cm, stop = h_pdoping + h_intrinsic + h_ndoping, length = lengthNaN) )
        yy        = collect( range(0.63e-5 * cm, stop = height,                              length = lengthNaN) )     
        ipsi      = 3

        phin      = append!(solution[iphin, :], NaN*ones( lengthNaN ))
        psi       = append!(solution[ipsi, :],  NaN*ones( lengthNaN ))
        Xsurf     = X
        Ysurf     = Y
        Xsurf     = append!(Xsurf, xx)
        Ysurf     = append!(Ysurf, yy)


        Plotter.figure()
        Plotter.surf(X[:], Y[:], psi)
        Plotter.title("Electrostatic potential \$ \\psi \$ at end time")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        #savefig("ipsi-2D-tangential-flux-on-prefactor-mobility-$(prefactor_mobility)-surface-reco-on.eps")
        ##########################
        Plotter.figure()
        Plotter.surf(X[:], Y[:], phin )
        Plotter.title("quasi Fermi potential \$ \\varphi_n \$ at end time")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        #savefig("iphin-2D-tangential-flux-on-prefactor-mobility-$(prefactor_mobility)-surface-reco-on.eps")

        #########################
        Plotter.figure()
        ########
        Plotter.plot(biasValues, abs.(IV.*(cm)^2), label = "simulation (prefactor mobility = $(prefactor_mobility) [\$ \\frac{\\mathrm{cm}^2}{\\mathrm{Vs}}\$] )",  linewidth= 3, linestyle="--", color="red")
        PyPlot.legend()
        Plotter.title("Forward I-V scan protocol")
        Plotter.ylabel("total current [A]") # 
        Plotter.xlabel("Applied Voltage [V]")
        #PyPlot.ylim(0.0, 0.006*1.0e3)
    end


    if test == false
        println("*** done\n")
    end

    testval = VoronoiFVM.norm(ctsys.fvmsys, solution, 2)
    return testval


end #  main

function test()
    # voltageAcceptor = 1.0 and ntsteps: 51, mesh: maxvolume=9.0e-18
    ################################################################################################################
    #
    #  norm values: (mobility_interface = prefactor_mobility * mobility_bulk)
    # 1. without tangential flux and without surface recombination:                 215.04468531184676
    #
    # 2. without tangential flux and with surface recombination:                    215.04457300209455
    #
    # 3. with tangential flux and without surface recombination:
    #    prefactor_mobility = 1.0:                                                  214.8638294731805
    #    prefactor_mobility = 1.0e2:                                                214.96762540843613
    #
    #
    # 4. with tangential flux and with surface recombination: 
    #    prefactor_mobility = 1.0:                                                  214.86382567944307
    #    prefactor_mobility = 1.0e2:                                                214.96762920024418
    ################################################################################################################
    main(test = true, unknown_storage=:dense) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
