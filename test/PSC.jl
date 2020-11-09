"""
Simulating charge transport in a perovskite solar cell (PSC)
without interfacial recombination.

"""
# Simulation values are from
# https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/Evidence_for_ion_migration_SRH.csv

module PSC

using VoronoiFVM
using ChargeTransportInSolids
using ExtendableGrids
using PyPlot; PyPlot.pygui(true)
using Printf
using DelimitedFiles

function main(;n = 7, pyplot = false, verbose = false, dense = true)

    # close all windows
    PyPlot.close("all")

    ################################################################################
    println("Set up grid and regions")
    ################################################################################

    # region numbers
    regionAcceptor  = 1                           # p doped region
    regionIntrinsic = 2                           # intrinsic region
    regionDonor     = 3                           # n doped region
    regions         = [regionAcceptor, regionIntrinsic, regionDonor]

    # boundary region numbers
    bregionAcceptor = 1
    bregionDonor    = 2
    bregions        = [bregionAcceptor, bregionDonor]

    # #inner regions
    # iregionAcceptor = 3
    # iregionDonor    = 4

    # iregions        = [iregionAcceptor, iregionDonor]

    # grid
    # NB: Using geomspace to create uniform mesh is not a good idea. It may create virtual duplicates at boundaries.
    h_pdoping       = 0.2 * μm
    h_intrinsic     = 0.4 * μm
    h_ndoping       = 0.2 * μm
    x0              = 0.0 * μm 
    δ               = 4*n        # the larger, the finer the mesh
    t               = 0.5*μm/δ   # tolerance for geomspace and glue (with factor 10)
    k               = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_p_u       = collect(range(x0, h_pdoping/2, step=h_pdoping/(0.5*δ)))
    coord_p_g       = geomspace(h_pdoping/2, 
                                h_pdoping, 
                                h_pdoping/(δ), 
                                h_pdoping/(1.5*δ), 
                                tol=t)
    coord_i_g1      = geomspace(h_pdoping, 
                                h_pdoping+h_intrinsic/k, 
                                h_intrinsic/(3.1*δ), 
                                h_intrinsic/(1.1*δ), 
                                tol=t)
    coord_i_g2      = geomspace(h_pdoping+h_intrinsic/k, 
                                h_pdoping+h_intrinsic,               
                                h_intrinsic/(1.1*δ),    
                                h_intrinsic/(3.0*δ), 
                                tol=t)
    coord_n_g       = geomspace(h_pdoping+h_intrinsic,               
                                h_pdoping+h_intrinsic+h_ndoping/2, 
                                h_ndoping/(1.5*δ),   
                                h_ndoping/(1δ),      
                                tol=t)
    coord_n_u       = collect(range(h_pdoping+h_intrinsic+h_ndoping/2, h_pdoping+h_intrinsic+h_ndoping, step=h_ndoping/(0.5*δ)))

    coord           = glue(coord_p_u,coord_p_g,  tol=10*t) 
    coord           = glue(coord,    coord_i_g1, tol=10*t)
    coord           = glue(coord,    coord_i_g2, tol=10*t) 
    coord           = glue(coord,    coord_n_g,  tol=10*t)
    coord           = glue(coord,    coord_n_u,  tol=10*t)
    grid            = ExtendableGrids.simplexgrid(coord)
    numberOfNodes   = length(coord)

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],                [h_pdoping],                           regionAcceptor)  # p-doped region   = 1
    cellmask!(grid, [h_pdoping],               [h_pdoping + h_intrinsic],             regionIntrinsic) # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor)     # n-doped region   = 3

    # # set interior boundary regions
    # bfacemask!(grid, [h_pdoping],           [h_pdoping],           iregionAcceptor)
    # bfacemask!(grid, [h_pdoping+h_intrinsic], [h_pdoping+h_intrinsic], iregionDonor)

    if pyplot
        ExtendableGrids.plot(grid, Plotter = PyPlot, p = PyPlot.plot()) 
        PyPlot.title("Grid")
    end
    println("*** done\n")

    ################################################################################
    println("Define physical parameters and model")
    ################################################################################

    # indices
    iphin, iphip, iphia, ipsi = 1:4
    species         = [iphin, iphip, iphia, ipsi]

    # number of (boundary) regions and carriers
    numberOfRegions         = length(regions)
    numberOfBoundaryRegions = length(bregions) #+ length(iregions)
    numberOfCarriers        = length(species) - 1

    # temperature
    T    = 300.0                *  K

    # band edge energies    
    Eref =  0.0#5.4        *  eV # reference energy 
    Ec_a = -3.8        *  eV 
    Ev_a = -5.4        *  eV 

    Ec_i = -3.8        *  eV 
    Ev_i = -5.4        *  eV 
    Ea_i = -4.8        *  eV 

    Ec_d = -3.8        *  eV 
    Ev_d = -5.4        *  eV 

    EC   = [Ec_a, Ec_i, Ec_d] 
    EV   = [Ev_a, Ev_i, Ev_d] 
    EA   = [0.0,  Ea_i, 0.0] 

    # effective densities of state
    Nc       = 1.0e20               / (cm^3)
    Nv       = 1.0e20               / (cm^3)
    Nanion   = 1.0e17               / (cm^3)

    NC   = [Nc, Nc, Nc]
    NV   = [Nv, Nv, Nv]
    NA   = [0.0, Nanion, 0.0]

    # mobilities 
    μn_a = 20.0                 * (cm^2) / (V * s)  
    μp_a = 20.0                 * (cm^2) / (V * s)  

    μn_i = 20.0                 * (cm^2) / (V * s)  
    μp_i = 20.0                 * (cm^2) / (V * s)  
    μa_i = 1.0e-10              * (cm^2) / (V * s)

    μn_d = 20.0                 * (cm^2) / (V * s) 
    μp_d = 20.0                 * (cm^2) / (V * s) 
 
    μn   = [μn_a, μn_i, μn_d] 
    μp   = [μp_a, μp_i, μp_d] 
    μa   = [0.0,  μa_i, 0.0] 

    # relative dielectric permittivity  

    ε_a  = 20.0                 *  1.0  
    ε_d  = 20.0                 *  1.0 
    ε_i  = 20.0                 *  1.0 

    ε   = [ε_a, ε_i, ε_d] 

    # recombination model
    recombinationOn = true

    # radiative recombination
    r0_a = 1.0e-10              * cm^3 / s 
    r0_i = 1.0e-10              * cm^3 / s  
    r0_d = 1.0e-10              * cm^3 / s

    r0   = [r0_a, r0_i, r0_d]

    # life times and trap densities 
    τn_a = 2.0e-15            * s 
    τp_a = 2.0e-15            * s

    τn_i = 1.0e6              * s
    τp_i = 1.0e6              * s
    τn_d = τn_a
    τp_d = τp_a

    τn   = [τn_a, τn_i, τn_d]
    τp   = [τp_a, τp_i, τp_d]

    # SRH trap energies (needed for calculation of recombinationSRHTrapDensity)
    Ei_a = -5.2                 * eV   + Eref
    Ei_i = -4.6                 * eV   + Eref
    Ei_d = -4.0                 * eV   + Eref

    EI   = [Ei_a, Ei_i, Ei_d]

    # Auger recombination
    Auger = 0.0

    # generation (only intrinsically present)
    generation_a = 0.0
    generation_i = 2.5e21 / (cm^3 * s)
    generation_d = 0.0
    generationEmittedLight = [generation_a, generation_i, generation_d]

    # doping (doping values are from Phils paper, not stated in the parameter list online)
    Nd             =   3.0e17   / (cm^3) 
    Na             =   3.0e17   / (cm^3) 
    C0             =   1.0e17   / (cm^3) 
    

    # intrinsic concentration (not doping!)
    ni             =   sqrt(Nc * Nv) * exp(-(Ec_i - Ev_i) / (2 * kB * T))

    # contact voltages
    voltageAcceptor     = 1.1*V #0.079 * V
    voltageDonor        = 0.0 * V 

    println("*** done\n")

    ################################################################################
    println("Define ChargeTransport data and fill in previously defined data")
    ################################################################################

    # initialize ChargeTransport instance
    data      = ChargeTransportInSolids.ChargeTransportData(numberOfNodes,
                                                            numberOfRegions,
                                                            numberOfBoundaryRegions,
                                                            numberOfSpecies = numberOfCarriers + 1)

    # region independent data
    data.F                              .= Boltzmann # Boltzmann, FermiDiracOneHalf, Blakemore
    data.temperature                     = T
    data.UT                              = (kB * data.temperature) / q
    data.contactVoltage[bregionAcceptor] = voltageAcceptor
    data.contactVoltage[bregionDonor]    = voltageDonor
    data.chargeNumbers[iphin]            = -1
    data.chargeNumbers[iphip]            =  1
    data.chargeNumbers[iphia]            =  1
    data.Eref                            =  Eref

    data.recombinationOn                 = recombinationOn


    # boundary region data
    for ibreg in bregions
        data.bDensityOfStates[iphin, ibreg] = Nc
        data.bDensityOfStates[iphip, ibreg] = Nv
        data.bDensityOfStates[iphia, ibreg] = 0.0
    end

    data.bBandEdgeEnergy[iphin, bregionAcceptor]  = Ec_a + data.Eref
    data.bBandEdgeEnergy[iphip, bregionAcceptor]  = Ev_a + data.Eref

    # data.bBandEdgeEnergy[iphin, iregionAcceptor]        =  Ec_a + data.Eref
    # data.bBandEdgeEnergy[iphip, iregionAcceptor]        =  Ev_a + data.Eref

    data.bBandEdgeEnergy[iphin, bregionDonor]     = Ec_d + data.Eref
    data.bBandEdgeEnergy[iphip, bregionDonor]     = Ev_d + data.Eref

    # data.bBandEdgeEnergy[iphin, iregionDonor]         =  Ec_a + data.Eref
    # data.bBandEdgeEnergy[iphip, iregionDonor]         =  Ev_a + data.Eref

    # interior region data
    for ireg in 1:numberOfRegions

        data.dielectricConstant[ireg]    = ε[ireg]

        # dos, band edge energy and mobilities
        data.densityOfStates[iphin,ireg] = NC[ireg]
        data.densityOfStates[iphip, ireg] = NV[ireg]
        data.densityOfStates[iphia, ireg] = NA[ireg]

        data.bandEdgeEnergy[iphin, ireg]  = EC[ireg] + data.Eref
        data.bandEdgeEnergy[iphip,ireg]  = EV[ireg] + data.Eref
        data.bandEdgeEnergy[iphia,ireg]  = EA[ireg] + data.Eref
        data.mobility[iphin,ireg]        = μn[ireg]
        data.mobility[iphip, ireg]        = μp[ireg]
        data.mobility[iphia, ireg]        = μa[ireg]

        # recombination parameters
        data.recombinationRadiative[ireg]            = r0[ireg]
        data.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        data.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        data.recombinationSRHTrapDensity[iphin, ireg] = ChargeTransportInSolids.trapDensity(iphin, ireg, data, EI[ireg])
        data.recombinationSRHTrapDensity[iphip, ireg] = ChargeTransportInSolids.trapDensity(iphip, ireg, data, EI[ireg])
        data.recombinationAuger[iphin, ireg]          = Auger
        data.recombinationAuger[iphip, ireg]          = Auger

        data.generationEmittedLight[ireg]            = generationEmittedLight[ireg]

    end

    # interior doping
    data.doping[iphip, regionAcceptor]   = Na 
    data.doping[iphia, regionAcceptor]   = 0.0        #                  ni   ni  C0; 
    data.doping[iphin, regionIntrinsic]  = 0.0#ni        #                  0.0  Na  0.0]
    data.doping[iphip, regionIntrinsic]  = 0.0#ni        
    data.doping[iphia, regionIntrinsic]  = C0        
    data.doping[iphin, regionDonor]      = Nd        # data.doping   = [Nd  0.0  0.0;                   
    data.doping[iphia,regionDonor]      = 0.0

    # boundary doping
    data.bDoping[iphip, bregionAcceptor] = Na        # data.bDoping  = [Nd  0.0;
    data.bDoping[iphin, bregionDonor]    = Nd        #                  0.0  Na]

    # print data
    println(data)
    println("*** done\n")

    ################################################################################
    println("Define physics and system")
    ################################################################################

    ## initializing physics environment ##
    physics = VoronoiFVM.Physics(
    data        = data,
    num_species = numberOfCarriers + 1,
    flux        = ChargeTransportInSolids.ScharfetterGummel!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
    reaction    = ChargeTransportInSolids.reaction!,
    breaction   = ChargeTransportInSolids.breaction!
    )

    if dense
        sys = VoronoiFVM.System(grid, physics, unknown_storage = :dense)
    else
        sys = VoronoiFVM.System(grid, physics, unknown_storage = :sparse)
    end

    # enable all three species in all regions
    enable_species!(sys, ipsi,  regions)
    enable_species!(sys, iphin, regions)
    enable_species!(sys, iphip, regions)
    enable_species!(sys, iphia, [regionIntrinsic])

    sys.boundary_values[iphin,  bregionAcceptor] = data.contactVoltage[bregionAcceptor]
    sys.boundary_factors[iphin, bregionAcceptor] = VoronoiFVM.Dirichlet

    sys.boundary_values[iphin,  bregionDonor]    = data.contactVoltage[bregionDonor]
    sys.boundary_factors[iphin, bregionDonor]    = VoronoiFVM.Dirichlet

    sys.boundary_values[iphip,  bregionAcceptor] = data.contactVoltage[bregionAcceptor]
    sys.boundary_factors[iphip, bregionAcceptor] = VoronoiFVM.Dirichlet

    sys.boundary_values[iphip,  bregionDonor]    = data.contactVoltage[bregionDonor]
    sys.boundary_factors[iphip, bregionDonor]    = VoronoiFVM.Dirichlet

    sys.boundary_values[iphia,  bregionAcceptor] = 0.0
    sys.boundary_factors[iphia, bregionAcceptor] = VoronoiFVM.Dirichlet

    sys.boundary_values[iphia,  bregionDonor]    = 0.0
    sys.boundary_factors[iphia, bregionDonor]    = VoronoiFVM.Dirichlet


    println("*** done\n")


    ################################################################################
    println("Define control parameters for Newton solver")
    ################################################################################

    control                   = VoronoiFVM.NewtonControl()
    control.verbose           = verbose
    control.max_iterations    = 200
    control.tol_absolute      = 1.0e-13
    control.tol_relative      = 1.0e-13
    control.handle_exceptions = true
    control.tol_round         = 1.0e-13
    control.max_round         = 5

    println("*** done\n")

    ################################################################################
    println("Compute solution in thermodynamic equilibrium for Boltzmann")
    ################################################################################

    data.inEquilibrium = true

    # initialize solution and starting vectors
    initialGuess                   = unknowns(sys)
    solution                       = unknowns(sys)
    @views initialGuess[ipsi,  :] .= 0.0
    @views initialGuess[iphin, :] .= 0.0
    @views initialGuess[iphip, :] .= 0.0
    @views initialGuess[iphia, :] .= 0.0

    control.damp_initial      = 0.9
    control.damp_growth       = 1.61 # >= 1
    control.max_round         = 5

    sys.boundary_values[iphin, bregionAcceptor] = 0.0 * V
    sys.boundary_values[iphip, bregionAcceptor] = 0.0 * V
    sys.physics.data.contactVoltage             = 0.0 * sys.physics.data.contactVoltage

    I = collect(1.0:-1:0.0)
    LAMBDA = 10 .^ (-I) 
    prepend!(LAMBDA,0.0)
    #node = 10
    for i in 1:length(LAMBDA)
        println("λ1 = $(LAMBDA[i])")
        sys.physics.data.λ1 = LAMBDA[i]
        solve!(solution, initialGuess, sys, control = control, tstep=Inf)
        initialGuess .= solution
        # if LAMBDA[i] == 1.0
        #     ChargeTransportInSolids.printJacobi(node, sys)
        # end
    end

    println("*** done\n")

    ################################################################################
    println("Bias loop")
    ################################################################################

    data.inEquilibrium = false

    control.damp_initial      = 0.5
    control.damp_growth       = 1.2 # >= 1
    control.max_round         = 7

    # set non equilibrium boundary conditions
    sys.physics.data.contactVoltage[bregionDonor]    = voltageDonor
    sys.physics.data.contactVoltage[bregionAcceptor] = voltageAcceptor
    sys.boundary_values[iphin, bregionAcceptor]      = data.contactVoltage[bregionAcceptor]
    sys.boundary_values[iphip, bregionAcceptor]      = data.contactVoltage[bregionAcceptor]

    maxBias    = data.contactVoltage[bregionAcceptor]
    biasValues = range(0, stop = maxBias, length = 11)
    IV         = zeros(0)

    w_device = 1.0 #0.5 * μm     # width of device
    z_device = 1.0 #1.0e-4 * cm  # depth of device


    for Δu in biasValues
        println("Bias value: Δu = $(Δu) (no illumination)")

        data.contactVoltage[bregionAcceptor] = Δu
        sys.boundary_values[iphin, bregionAcceptor] = Δu
        sys.boundary_values[iphip, bregionAcceptor] = Δu

        solve!(solution, initialGuess, sys, control = control, tstep = Inf)
        initialGuess .= solution

        # get IV curve
        factory = VoronoiFVM.TestFunctionFactory(sys)

        # testfunction zero in bregionDonor and one in bregionAcceptor
        tf     = testfunction(factory, [bregionDonor], [bregionAcceptor])
        I      = integrate(sys, tf, solution)

        push!(IV,  abs.(w_device * z_device * (I[iphin] + I[iphip])))

    #     # plotting
    #     if pyplot
    #         PyPlot.figure()
    #         ChargeTransportInSolids.plotDensities(grid, data, solution, "$Δu (no illumination)")
    #         PyPlot.figure()
    #         ChargeTransportInSolids.plotEnergies(grid, data, solution, "$Δu (no illumination)")
    #         #PyPlot.figure()
    #         #ChargeTransportInSolids.plotIV(biasValues,IV, Δu)
    # end


    end # bias loop

    # return IV
    println("*** done\n")
    ################################################################################
    println("Illumination loop")
    ################################################################################ 
    control.damp_initial      = 0.005
    control.damp_growth       = 1.2 # >= 1
    control.max_round         = 7
    I = collect(15.0:-1:0.0)
    LAMBDA = 10 .^ (-I) 
    prepend!(LAMBDA,0.0)

    for i in 1:length(LAMBDA)
        println("λ2 = $(LAMBDA[i])")
        sys.physics.data.λ2 = LAMBDA[i]
        solve!(solution, initialGuess, sys, control = control, tstep=Inf)
        initialGuess = solution
    end

        sol = [coord solution']
        writedlm("PSC-illuminated-sol.dat", sol)
    
    # if pyplot
    #     PyPlot.figure()
    #     ChargeTransportInSolids.plotDensities(grid, data, solution, "$(maxBias) (illumination)")
    #     PyPlot.figure()
    #     ChargeTransportInSolids.plotEnergies(grid, data, solution, "$(maxBias) (illumination)")
    #     #PyPlot.figure()
    #     #ChargeTransportInSolids.plotSolution(coord, solution, data.Eref, "$(maxBias) (illuminated)")
    # end
    
    println("*** done\n")
     ################################################################################
     println("Reverse Bias loop")
     ################################################################################
 
     data.inEquilibrium = false
 
     control.damp_initial      = 0.05
     control.damp_growth       = 1.2 # >= 1
     control.max_round         = 7

     IV         = zeros(0)
 
     w_device = 1.0 #0.5 * μm     # width of device
     z_device = 1.0 #1.0e-4 * cm  # depth of device

 
     for Δu in reverse(biasValues)

         println("Bias value: Δu = $(Δu) (illumination)")
 
         data.contactVoltage[bregionAcceptor] = Δu
         sys.boundary_values[iphin, bregionAcceptor] = Δu
         sys.boundary_values[iphip, bregionAcceptor] = Δu
 
         solve!(solution, initialGuess, sys, control = control, tstep = Inf)
         # if Δu == 0.0
         #     ChargeTransportInSolids.printJacobi(node, sys)
         # end
         # return
         initialGuess .= solution
 
         # get IV curve
         factory = VoronoiFVM.TestFunctionFactory(sys)
 
         # testfunction zero in bregionDonor and one in bregionAcceptor
         tf     = testfunction(factory, [bregionDonor], [bregionAcceptor])
         I      = integrate(sys, tf, solution)
 
         push!(IV,  abs.(w_device * z_device * (I[iphin] + I[iphip])))
 
         # if Δu==maxBias
         #     sol = [coord solution']
         #     writedlm("jl-PSC-sol-Boltzmann-SG.dat", sol)
         # end
 
         # plotting
         if pyplot
             if Δu == maxBias #>0.2 # == maxBias
             #PyPlot.figure()
             ChargeTransportInSolids.plotDensities(grid, data, solution, "$Δu (illumination)")
             #savefig("PSC-sol-$Δu-dens.eps")
             PyPlot.figure()
             ChargeTransportInSolids.plotEnergies(grid, data, solution, "$Δu (illumination)")
             PyPlot.figure()
             #ChargeTransportInSolids.plotSolution(coord, solution, data.Eref, "$Δu (illumination)")
             #ChargeTransportInSolids.plotIV(biasValues,IV, Δu)
             end
         #end
     end
 
 
     end # bias loop
     println("*** done\n")

    # ################################################################################
    # println("Transient solution")
    # ################################################################################

    # tstep = 0.5*1e-16
    # tstep_max = 0.5*1e0
    # dV = 0.25


    #  # Solve the stationary state system
    # control=VoronoiFVM.NewtonControl()
    # control.Δt=tstep
    # control.Δt_max=tstep_max
    # control.Δt_grow=1.5
    # control.Δu_opt=dV
    # control.verbose=false
    # control.max_lureuse=0
    # control.edge_cutoff=1.0e-16

    # # inival.=initial_solution
    # initial_solution[1:2,2:end-1].=0.0
    # println(initial_solution)
    # sampling_times=[t for t in 0.0:1e-1:1e1]

    # solution_transient = unknowns(sys)

    # evolve!(solution_transient,initial_solution,sys,sampling_times, control=control)

    # ChargeTransportInSolids.plotDensities(grid, data, solution_transient, "FINAL")

    # PyPlot.figure()
    # ChargeTransportInSolids.plotEnergies(grid, data, solution_transient, "FINAL")

    # println("*** done\n")

    # return solution_transient



    

end #  main

println("This message should show when the PSC module is successfully recompiled.")

end # module


#     # Number of periods to run
#     nper=1

#     # Number of samples per period
#     nsamp=4
#     ###############################################################

#     # BV kinetic constants
#     phi_min=0*V
#     phi_max=-1.6*V

#     # Scan rate
#     scan=scan_rate*mV/s

#     per=2*abs(phi_max-phi_min)/scan

#     sampling_times=[t for t in 0.0:0.5*per/nsamp:per*nper]


#             function pre(sol,time) # vielleich t nicht nötig bei uns
#             #theta=Nernst_const(time)
#             #omega=freq*2.0*π
#             eplus,eminus=BV_rate_constants(time)
#         end
#         I_disk=[0.0,0.0]
#         I_disk_old=[0.0,0.0]
#         I_ring=[]
#         di=0.0
#         function delta(solution, oldsolution, time,tstep) # discrete Energie
#             I_disk=VoronoiFVM.integrate(rdcell,tfc_disk,solution,oldsolution,tstep)
#             I_ring=VoronoiFVM.integrate(rdcell,tfc_ring,solution,oldsolution,tstep)
#             di=FaradayConstant*abs(I_disk_old[specB]-I_disk[specB])/mA
#         end

        
#         function post(solution, oldsolution, time,tstep) # vielleich t nicht nötig bei uns, plotten
#             push!(vdisk,phi_cv(time))
#             push!(iring,-I_ring[specB]*FaradayConstant)
#             push!(idisk,I_disk[specB]*FaradayConstant)
#             push!(time_discretization,time)
#             if verbose
#                 ProgressMeter.next!(pmeter,showvalues=[
#                     (:Δϕ,phi_cv(time)),
#                     (:dI,di),
#                     (:t,time),
#                 ],valuecolor=:yellow)
#             end
#             I_disk_old=I_disk
#         end

# evolve!(solution,inival,sys,sampling_times, control=control, pre=pre,post=post,delta=delta)