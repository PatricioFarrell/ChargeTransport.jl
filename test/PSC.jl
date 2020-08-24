"""
Simulating charge transport in a perovskite solar cell (PSC).
"""

module PSC

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using PyPlot; PyPlot.pygui(true)
using Printf

function main(;n = 4, pyplot = false, verbose = false, dense = true)

    # close all windows
    PyPlot.close("all")

    ################################################################################
    println("Set up grid and regions")
    ################################################################################

    # region numbers
    regionDonor     = 1                           # n doped region
    regionIntrinsic = 2                           # intrinsic region
    regionAcceptor  = 3                           # p doped region
    regions         = [regionDonor, regionIntrinsic, regionAcceptor]

    # boundary region numbers
    bregionDonor    = 1
    bregionAcceptor = 2
    bregions        = [bregionDonor, bregionAcceptor]

    # inner regions
    # iregionDonor    = 3
    # iregionAcceptor = 4
    # iregions        = [iregionDonor, iregionAcceptor]

    # grid
    # NB: Using geomspace to create uniform mesh is not a good idea. It may create virtual duplicates at boundaries.
    h_ndoping       = 0.2 * μm
    h_intrinsic     = 0.4 * μm
    h_pdoping       = 0.2 * μm
    x0              = 0 * μm 
    δ               = 15        # the larger, the finer the mesh
    t               = 0.2*μm/δ  # tolerance for geomspace and glue (with factor 10)
    k               = 1.5       # the closer to 1, the closer to the boundary geomspace works

    coord_n_u       = collect(range(x0, h_ndoping/2, step=h_ndoping/δ))
    coord_n_g       = geomspace(h_ndoping/2, 
                                h_ndoping, 
                                h_ndoping/δ, 
                                h_ndoping/(2δ), 
                                tol=t)
    coord_i_g1      = geomspace(h_ndoping, 
                                h_ndoping+h_intrinsic/k, 
                                h_intrinsic/(4δ), 
                                h_intrinsic/δ, 
                                tol=t)
    coord_i_g2      = geomspace(h_ndoping+h_intrinsic/k, 
                                h_ndoping+h_intrinsic,               
                                h_intrinsic/δ,    
                                h_intrinsic/(4δ), 
                                tol=t)
    coord_p_g       = geomspace(h_ndoping+h_intrinsic,               
                                h_ndoping+h_intrinsic+h_pdoping/2, 
                                h_pdoping/(2δ),   
                                h_pdoping/δ,      
                                tol=t)
    coord_p_u       = collect(range(h_ndoping+h_intrinsic+h_pdoping/2, h_ndoping+h_intrinsic+h_pdoping, step=h_pdoping/δ))

    coord           = glue(coord_n_u,coord_n_g,  tol=10*t)
    coord           = glue(coord,    coord_i_g1, tol=10*t)
    coord           = glue(coord,    coord_i_g2, tol=10*t)
    coord           = glue(coord,    coord_p_g,  tol=10*t)
    coord           = glue(coord,    coord_p_u,  tol=10*t)
    grid            = ExtendableGrids.simplexgrid(coord)
    numberOfNodes   = length(coord)

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],                [h_ndoping],                           regionDonor)     # n-doped region   = 1
    cellmask!(grid, [h_ndoping],               [h_ndoping + h_intrinsic],             regionIntrinsic) # intrinsic region = 2
    cellmask!(grid, [h_ndoping + h_intrinsic], [h_ndoping + h_intrinsic + h_pdoping], regionAcceptor)  # p-doped region   = 3

    # set interior boundary regions
    # bfacemask!(grid, [h_ndoping],           [h_ndoping],           iregionDonor)
    # bfacemask!(grid, [h_ndoping+h_pdoping], [h_ndoping+h_pdoping], iregionIntrinsic)

    if pyplot
        ExtendableGrids.plot(VoronoiFVM.Grid(coord, collect(0.0 : 0.25*μm : 0.5*μm)), Plotter = PyPlot, p = PyPlot.plot()) 
        PyPlot.title("Grid")
    end
    println("*** done\n")

    ################################################################################
    println("Define physical parameters and model")
    ################################################################################

    # indices

    iphin, iphip, ipsi        = 1:3
    species                   = [iphin, iphip, ipsi]
    # iphin, iphip, iphia, ipsi = 1:4
    # species         = [iphin, iphip, iphia, ipsi]

    # number of (boundary) regions and carriers
    numberOfRegions         = length(regions)
    numberOfBoundaryRegions = length(bregions) 
    numberOfSpecies         = length(species)

    # temperature
    T    = 300.0                *  K

    # band edge energies
    Eref =  3.8              
    Ec_d = (-3.8 + Eref)        *  eV # -4.1 # 3.3
    Ev_d = (-5.4 + Eref)        *  eV # -7.4 # 0.0

    Ec_i = (-3.8 + Eref)        *  eV # -3.8 # 1.6
    Ev_i = (-5.4 + Eref)        *  eV # -5.4 # 0.0
    Ea_i = (-1.8 + Eref)        *  eV ## ????? # -1.8

    Ec_a = (-3.8 + Eref)        *  eV # -1.9 # 3.0
    Ev_a = (-5.4 + Eref)        *  eV # -4.9 # 0.0

    EC   = [Ec_d, Ec_i, Ec_a] 
    EV   = [Ev_d, Ev_i, Ev_a] 
    EA   = [0.0,  Ea_i, 0.0] 

    # effective densities of state
    Nc   = 1.0e20               / (cm^3)
    Nv   = 1.0e20               / (cm^3)
    Na   = 1.0e19               / (cm^3)

    NC   = [Nc,  Nc, Nc]
    NV   = [Nv,  Nv, Nv]
    NA   = [0.0, Na, 0.0]

    # mobilities 
    μn_d = 20.0                 * (cm^2) / (V * s)  # 20 in Phil's paper  
    μp_d = 20.0                 * (cm^2) / (V * s)  # 20 in Phil's paper  

    μn_i = 20.0                 * (cm^2) / (V * s)  
    μp_i = 20.0                 * (cm^2) / (V * s)  
    μa_i = 1.0e-10              * (cm^2) / (V * s)

    μn_a = 20.0                 * (cm^2) / (V * s)  # 20 in Phil's paper  
    μp_a = 20.0                 * (cm^2) / (V * s)  # 20 in Phil's paper  
 
    μn   = [μn_d, μn_i, μn_a] 
    μp   = [μp_d, μp_i, μp_a] 
    μa   = [0.0,  μa_i, 0.0] 

    # relative dielectric permittivity  
    ε_d  = 20.0                 *  1.0         # 20 in Phil's paper      
    ε_i  = 20.0                 *  1.0         # 20 in Phil's paper      
    ε_a  = 20.0                 *  1.0         # 20 in Phil's paper  
    
    ε   = [ε_d, ε_i, ε_a] 

    # radiative recombination
    r0_d = 1.0e-10              * cm^3 / s   # 1e-10 in Phil's paper 
    r0_i = 1.0e-10              * cm^3 / s   # 1e-10 in Phil's paper  
    r0_a = 1.0e-10              * cm^3 / s   # 1e-10 in Phil's paper 

    r0   = [r0_d, r0_i, r0_a]
    ra   = [0.0, 0.0, 0.0]

    # life times
    τn_d = 1e-14                * s
    τp_d = 1e-15                * s
    τn_i = 1e6                  * s   # 2.0e-15 *s
    τp_i = 1e6                  * s   # 2.0e-15 *s
    τn_a = τn_d
    τp_a = τp_d

    τn   = [τn_d, τn_i, τn_a]
    τp   = [τp_d, τp_i, τp_a]

    # intrinsic Fermi energies
  
    Ei_d = -5.2                 * eV   
    Ei_i = -4.6                 * eV  # paper: intrinsic: -0.2+eV (n), -1.4 *eV (p)
    Ei_a = -4.0                 * eV

    EI   = [Ei_d, Ei_i, Ei_a]

    # recombination and generation parameters
    # noch in dd_system.jl direkt gesetzt
    # G               = 2.5e21    / (cm^3 * s)        # uniform generation rate (but only in the intrinsic layer)

    # doping
    Nd             =   3.0e18   / (cm^3)
    Na             =   3.0e18   / (cm^3)
    # E0             =   -0.5 *eV
    # Nd             =   Nc*exp((E0-Ec_i)/(kB*T))
    # Na             =   Nv*exp((Ev_i-E0)/(kB*T))
    # println(Nd)
    # println(Na)
    # @assert 1==0
    C0             =   1.0e19   / (cm^3)        # average anion vacancy density
    

    # intrinsic concentration (not doping!)
    ni             =   sqrt(Nc * Nv) * exp(-(Ec_i - Ev_i) / (2 * kB * T)) #/ (cm^3)

    # contact voltages
    voltageDonor     = 0.0 * V
    voltageAcceptor  = 1.0 * V # 1.1

    println("*** done\n")

    # PIN values
    εr   = 12.9 
    Ec   = 1.424                *  eV
    Ev   = 0.0                  *  eV
    Nc   = 4.351959895879690e17 / (cm^3)
    Nv   = 9.139615903601645e18 / (cm^3)
    mun  = 8500.0               * (cm^2) / (V * s)
    mup  = 400.0                * (cm^2) / (V * s)
    ε   = [εr,εr, εr] 
    EC   = [Ec, Ec, Ec] 
    EV   = [Ev, Ev, Ev]
    Ec_d = Ec; Ec_a = Ec; Ev_d = Ev; Ev_a = Ev;
    NC   = [Nc,  Nc, Nc]
    NV   = [Nv,  Nv, Nv]
    μn   = [mun, mun, mun] 
    μp   = [mup, mup, mup]
    dopingFactorNd =   1.0
    dopingFactorNa =   0.46
    Nd             =   dopingFactorNd * Nc
    Na             =   dopingFactorNa * Nv
    Auger           = 1.0e-29   * cm^6 / s          # 1.0e-41
    SRH_TrapDensity = 1.0e10    / cm^3              # 1.0e16
    SRH_LifeTime    = 1.0       * ns                # 1.0e10
    Radiative       = 1.0e-10   * cm^3 / s          # 1.0e-16


    ################################################################################
    println("Define ChargeTransport data and fill in previously defined data")
    ################################################################################

    # initialize ChargeTransport instance
    data      = ChargeTransport.ChargeTransportData(numberOfNodes, numberOfRegions, numberOfBoundaryRegions, numberOfSpecies)

    # region independent data
    data.F                              .= Boltzmann # Boltzmann, FermiDiracOneHalf, Blakemore
    data.temperature                     = T
    data.UT                              = (kB * data.temperature) / q
    data.contactVoltage[bregionDonor]    = voltageDonor
    data.contactVoltage[bregionAcceptor] = voltageAcceptor
    data.chargeNumbers[iphin]            = -1
    data.chargeNumbers[iphip]            =  1
    # data.chargeNumbers[iphia] =  1

    # boundary region data
    for ibreg in bregions

        data.bDensityOfStates[ibreg,iphin] = Nc
        data.bDensityOfStates[ibreg,iphip] = Nv

    end

    data.bBandEdgeEnergy[bregionDonor,iphin]     = Ec_d
    data.bBandEdgeEnergy[bregionDonor,iphip]     = Ev_d
    data.bBandEdgeEnergy[bregionAcceptor,iphin]  = Ec_a
    data.bBandEdgeEnergy[bregionAcceptor,iphip]  = Ev_a
    # data.bBandEdgeEnergy[bregionAcceptor,iphia]  = 0.0 # should not be necessary
    # data.bBandEdgeEnergy[bregionDonor,iphia]     = 0.0 # should not be necessary

    # interior region data
    for ireg in 1:numberOfRegions

        data.dielectricConstant[ireg]    = ε[ireg]

        # dos, band edge energy and mobilities
        data.densityOfStates[ireg,iphin] = NC[ireg]
        data.densityOfStates[ireg,iphip] = NV[ireg]
        # data.densityOfStates[ireg,iphia] = NA[ireg]
        data.bandEdgeEnergy[ireg,iphin]  = EC[ireg]
        data.bandEdgeEnergy[ireg,iphip]  = EV[ireg]
        # data.bandEdgeEnergy[ireg,iphia]  = EA[ireg]
        data.mobility[ireg,iphin]        = μn[ireg]
        data.mobility[ireg,iphip]        = μp[ireg]
        # data.mobility[ireg,iphia]        = μa[ireg]

        # recombination parameters
        data.recombinationRadiative[ireg]            = r0[ireg]
        data.recombinationRadiative[ireg]            = r0[ireg]
        data.recombinationSRHLifetime[ireg,iphin]    = τn[ireg]
        data.recombinationSRHLifetime[ireg,iphip]    = τp[ireg]
        data.recombinationSRHTrapDensity[ireg,iphin] = ChargeTransport.trapDensity(iphin, ireg, data, EI[ireg])
        data.recombinationSRHTrapDensity[ireg,iphip] = ChargeTransport.trapDensity(iphin, ireg, data, EI[ireg])
        # data.recombinationAuger[ireg,iphin]          = Auger
        # data.recombinationAuger[ireg,iphip]          = Auger

        # PIN recombination parameters
        # data.recombinationRadiative[ireg]      = Radiative
        # data.recombinationRadiative[ireg]      = Radiative
        # data.recombinationSRHLifetime[ireg,iphin]    = SRH_LifeTime
        # data.recombinationSRHLifetime[ireg,iphip]    = SRH_LifeTime
        # data.recombinationSRHTrapDensity[ireg,iphin] = SRH_TrapDensity
        # data.recombinationSRHTrapDensity[ireg,iphip] = SRH_TrapDensity
        # data.recombinationAuger[ireg,iphin]          = Auger
        # data.recombinationAuger[ireg,iphip]          = Auger

    end

    # interior doping
    data.doping[regionDonor,iphin]      = Nd        # data.doping   = [Nd  0.0  0.0;                   
    # data.doping[regionDonor,iphia]      = 0.0       #                  ni   ni  C0; 
    data.doping[regionIntrinsic,iphin]  = ni        #                  0.0  Na  0.0]
    data.doping[regionIntrinsic,iphip]  = ni        
    # data.doping[regionIntrinsic,iphia]  = C0        
    data.doping[regionAcceptor,iphip]   = Na
    # data.doping[regionAcceptor,iphia]   = 0.0

    # boundary doping
    data.bDoping[bregionDonor,iphin]    = Nd        # data.bDoping  = [Nd  0.0;
    data.bDoping[bregionAcceptor,iphip] = Na        #                  0.0  Na]

    # print data
    println(data)

    println("*** done\n")

    # psi0 = ChargeTransport.electroNeutralSolutionBoltzmann(grid, data)
    # psi0 = ChargeTransport.electroNeutralSolution!(data, grid)
    # println(psi0)
    if pyplot
        ################################################################################
        println("Plot electroneutral potential and doping")
        ################################################################################
        # ChargeTransport.plotEnergies(grid, data)
        PyPlot.figure()
        ChargeTransport.plotDoping(grid, data)
        # ChargeTransport.plotElectroNeutralSolutionBoltzmann(grid, psi0)

        println("*** done\n")
    end

    ################################################################################
    println("Define physics and system")
    ################################################################################

    ## initializing physics environment ##
    physics = VoronoiFVM.Physics(
    data        = data,
    num_species = numberOfSpecies,
    flux        = ChargeTransport.Sedan!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
    reaction    = ChargeTransport.reaction!,
    breaction   = ChargeTransport.breaction!
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
    # enable_species!(sys, iphia, [regionIntrinsic])


    sys.boundary_values[iphin,  bregionDonor]    = data.contactVoltage[bregionDonor]
    sys.boundary_factors[iphin, bregionDonor]    = VoronoiFVM.Dirichlet

    sys.boundary_values[iphin,  bregionAcceptor] = data.contactVoltage[bregionAcceptor]
    sys.boundary_factors[iphin, bregionAcceptor] = VoronoiFVM.Dirichlet

    sys.boundary_values[iphip,  bregionDonor]    = data.contactVoltage[bregionDonor]
    sys.boundary_factors[iphip, bregionDonor]    = VoronoiFVM.Dirichlet

    sys.boundary_values[iphip,  bregionAcceptor] = data.contactVoltage[bregionAcceptor]
    sys.boundary_factors[iphip, bregionAcceptor] = VoronoiFVM.Dirichlet

    println("*** done\n")


    ################################################################################
    println("Define control parameters for Newton solver")
    ################################################################################

    control = VoronoiFVM.NewtonControl()
    control.verbose           = verbose
    control.damp_initial      = 0.001
    control.damp_growth       = 1.21
    control.max_iterations    = 250
    control.tol_absolute      = 1.0e-14
    control.tol_relative      = 1.0e-14
    control.handle_exceptions = true
    control.tol_round         = 1.0e-14
    control.max_round         = 5
    control.tol_linear        = 1.0e-14
    # control.verbose           = verbose
    # control.damp_initial      = 0.001
    # control.damp_growth       = 1.5
    # control.max_iterations    = 250
    # control.tol_absolute      = 1.0e-10
    # control.tol_relative      = 1.0e-10
    # control.handle_exceptions = true
    # control.max_round         = 5
    # control.Δp                = 1e-16
    # control.Δp_min            = 1e-15


    println("*** done\n")

    ################################################################################
    println("Compute solution in thermodynamic equilibrium for Boltzmann")
    ################################################################################

    data.inEquilibrium = true

    # initialize solution and starting vectors
    initialGuess                   = unknowns(sys)
    solution                       = unknowns(sys)
    @views initialGuess[ipsi,  :] .= 0.0 # psi0
    @views initialGuess[iphin, :] .= 0.0
    @views initialGuess[iphip, :] .= 0.0

    # ChargeTransport.solveEquilibriumBoltzmann!(solution, initialGuess, data, grid, control, dense)

    function pre(u,lambda)
        sys.physics.data.λ1 = lambda
        # sys.physics.data.contactVoltage[bregionAcceptor] = lambda * 3.0
        sys.boundary_values[iphin, bregionAcceptor] = 0.0
        sys.boundary_values[iphip, bregionAcceptor] = 0.0
    end

    control.damp_initial      = 0.05
    control.damp_growth       = 1.2 # >= 1
    control.max_round         = 3

    sys.boundary_values[iphin, bregionAcceptor] = 0.0*V
    sys.boundary_values[iphip, bregionAcceptor] = 0.0*V
    sys.physics.data.contactVoltage             = 0.0 * sys.physics.data.contactVoltage

    I = collect(30.0:-1:0.0)
    LAMBDA = 10 .^ (-I) 
    prepend!(LAMBDA,0.0)


    for i in 1:length(LAMBDA)
        println("λ1 = $(LAMBDA[i])")
        sys.physics.data.λ1 = LAMBDA[i]
        solve!(solution, initialGuess, sys, control = control, tstep=Inf)
        initialGuess = solution
    end


    if pyplot
        ChargeTransport.plotDensities(grid, data, solution, "EQUILIBRIUM")
        PyPlot.figure()
        ChargeTransport.plotEnergies(grid, data, solution, "EQUILIBRIUM")
        PyPlot.figure()
    end

    # @assert 1==0


    println("*** done\n")



    ################################################################################
    println("Bias loop")
    ################################################################################

    data.inEquilibrium = false

    # set non equilibrium boundary conditions
    sys.physics.data.contactVoltage[bregionDonor]    = voltageDonor
    sys.physics.data.contactVoltage[bregionAcceptor] = voltageAcceptor
    sys.boundary_values[iphin, bregionAcceptor]      = data.contactVoltage[bregionAcceptor]
    sys.boundary_values[iphip, bregionAcceptor]      = data.contactVoltage[bregionAcceptor]

    function pre(u,lambda)
        # sys.physics.data.contactVoltage[bregionAcceptor] = lambda * 3.0
        sys.boundary_values[iphin, bregionAcceptor] = lambda * data.contactVoltage[bregionAcceptor]
        sys.boundary_values[iphip, bregionAcceptor] = lambda * data.contactVoltage[bregionAcceptor]
    end

    maxBias    = data.contactVoltage[bregionAcceptor]
    biasValues = range(0, stop = maxBias, length = 41)
    IV         = zeros(0)

    w_device = 0.5 * μm     # width of device
    z_device = 1.0e-4 * cm  # depth of device


    for Δu in biasValues

        println("Bias value: Δu = $(Δu) (no illumination)")

        data.contactVoltage[bregionAcceptor] = Δu
        sys.boundary_values[iphin, bregionAcceptor] = Δu
        sys.boundary_values[iphip, bregionAcceptor] = Δu

        solve!(solution, initialGuess, sys, control = control, tstep = Inf)
        # embed!(solution,initialGuess,sys,pre=pre, control=control)

        initialGuess .= solution

        # get IV curve
        factory = VoronoiFVM.TestFunctionFactory(sys)

        # testfunction zero in bregionAcceptor and one in bregionDonor
        tf     = testfunction(factory, [bregionAcceptor], [bregionDonor])
        I      = integrate(sys, tf, solution)

        push!(IV,  abs.(w_device * z_device * (I[iphin] + I[iphip])))

        # plot solution and IV curve
        if pyplot
            #ChargeTransport.plotEnergies(grid, data, solution, Δu)
            ChargeTransport.plotDensities(grid, data, solution, Δu)
            # if Δu == 0.0 || Δu == 1.5 Δu == 3
            #     savefig("psc-densities-nref-$n-deltaU-$Δu.eps")
            # end
            #ChargeTransport.plotIV(biasValues,IV)
        end


    end # bias loop

    # return IV
    println("*** done\n")

    ################################################################################
    println("Illumination loop")
    ################################################################################ 

    I = collect(20.0:-1:0.0)
    LAMBDA = 10 .^ (-I) 
    prepend!(LAMBDA,0.0)

    for i in 1:length(LAMBDA)
        println("λ2 = $(LAMBDA[i])")
        sys.physics.data.λ2 = LAMBDA[i]
        solve!(solution, initialGuess, sys, control = control, tstep=Inf)
        initialGuess = solution
    end

    if pyplot
        PyPlot.figure()
        ChargeTransport.plotDensities(grid, data, solution, "$(maxBias) (illuminated)")
        PyPlot.figure()
        ChargeTransport.plotEnergies(grid, data, solution, "$(maxBias) (illuminated)")
    end

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

    # ChargeTransport.plotDensities(grid, data, solution_transient, "FINAL")

    # PyPlot.figure()
    # ChargeTransport.plotEnergies(grid, data, solution_transient, "FINAL")

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