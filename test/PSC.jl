"""
Simulating charge transport in a perovskite solar cell (PSC).
"""

module PSC

using VoronoiFVM
using DDFermi
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
    regionIntrinsic = 2
    regionAcceptor  = 3                           # p doped region
    regions         = [regionDonor, regionIntrinsic, regionAcceptor]

    # boundary region numbers
    bregionDonor    = 1
    bregionAcceptor = 2
    bregions        = [bregionDonor, bregionAcceptor]

    # inner regions
    iregionDonor    = 3
    iregionAcceptor = 4
    iregions        = [iregionDonor, iregionAcceptor]

    # grid
    # NB: Using geomspace to create uniform mesh is not a good idea. It may create virtual duplicates at boundaries.
    h_ndoping       = 2 * μm
    h_intrinsic     = 2 * μm
    h_pdoping       = 2 * μm
    x0              = 0 * μm 
    δ               = 5        # the larger, the finer the mesh
    t               = μm/δ     # tolerance for geomspace and glue (with factor 10)
    k               = 1.5      # the closer to 1, the closer to the boundary geomspace works

    coord_n_u       = collect(range(x0, h_ndoping/k, step=h_ndoping/δ))
    coord_n_g       = geomspace(h_ndoping/k, 
                                h_ndoping, 
                                h_ndoping/δ, 
                                h_ndoping/(2δ), 
                                tol=t)
    coord_i_g1      = geomspace(h_ndoping, 
                                h_ndoping+h_intrinsic/k, 
                                h_intrinsic/(2δ), 
                                h_intrinsic/δ, 
                                tol=t)
    coord_i_g2      = geomspace(h_ndoping+h_intrinsic/k, 
                                h_ndoping+h_intrinsic,               
                                h_intrinsic/δ,    
                                h_intrinsic/(2δ), 
                                tol=t)
    coord_p_g       = geomspace(h_ndoping+h_intrinsic,               
                                h_ndoping+h_intrinsic+h_pdoping/k, 
                                h_pdoping/(2δ),   
                                h_pdoping/δ,      
                                tol=t)
    coord_p_u       = collect(range(h_ndoping+h_intrinsic+h_pdoping/k, h_ndoping+h_intrinsic+h_pdoping, step=h_pdoping/δ))

    coord           = glue(coord_n_u,coord_n_g,  tol=10*t)
    coord           = glue(coord,    coord_i_g1, tol=10*t)
    coord           = glue(coord,    coord_i_g2, tol=10*t)
    coord           = glue(coord,    coord_p_g,  tol=10*t)
    coord           = glue(coord,    coord_p_u,  tol=10*t)
    grid            = VoronoiFVM.Grid(coord)
    numberOfNodes   = length(coord)

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],                [h_ndoping],                           regionDonor)     # n-doped region   = 1
    cellmask!(grid, [h_ndoping],               [h_ndoping + h_intrinsic],             regionIntrinsic) # intrinsic region = 2
    cellmask!(grid, [h_ndoping + h_intrinsic], [h_ndoping + h_intrinsic + h_pdoping], regionAcceptor)  # p-doped region   = 3

    # set interior boundary regions
    bfacemask!(grid, [h_ndoping],           [h_ndoping],           iregionDonor)
    bfacemask!(grid, [h_ndoping+h_pdoping], [h_ndoping+h_pdoping], iregionAcceptor)

    if pyplot
        ExtendableGrids.plot(VoronoiFVM.Grid(coord, collect(0.0 : 0.25*μm : 0.5*μm)), Plotter = PyPlot, p = PyPlot.plot()) 
    end
    println("*** done\n")

    ################################################################################
    println("Define physical parameters and model")
    ################################################################################

    # indices
    iphin           = 1
    iphip           = 2
    ipsi            = 3
    species         = [iphin, iphip, ipsi]

    # number of (boundary) regions and carriers
    numberOfRegions         = length(regions)
    numberOfBoundaryRegions = length(bregions) + length(iregions)
    numberOfSpecies         = length(species)

    # physical data
    Ec   = 1.424                *  eV
    Ev   = 0.0                  *  eV
    Nc   = 4.351959895879690e17 / (cm^3)
    Nv   = 9.139615903601645e18 / (cm^3)
    mun  = 8500.0               * (cm^2) / (V * s)
    mup  = 400.0                * (cm^2) / (V * s)
    εr   = 12.9                 *  1.0              # relative dielectric permittivity of GAs
    T    = 300.0                *  K


    # recombination parameters
    Auger           = 1.0e-29   * cm^6 / s          # 1.0e-41
    SRH_TrapDensity = 1.0e10    / cm^3              # 1.0e16
    SRH_LifeTime    = 1.0       * ns                # 1.0e10
    Radiative       = 1.0e-10   * cm^3 / s          # 1.0e-16

    # doping
    dopingFactorNd =   1.0
    dopingFactorNa =   0.46
    Nd             =   dopingFactorNd * Nc
    Na             =   dopingFactorNa * Nv

    # intrinsic concentration (not doping!)
    ni             =   sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * kB * T)) / (cm^3)

    # contact voltages
    voltageDonor     = 0.0 * V
    voltageAcceptor  = 3.0 * V

    println("*** done\n")


    ################################################################################
    println("Define ddfermi data and fill in previously defined data")
    ################################################################################

    # initialize ddfermi instance
    data      = DDFermi.DDFermiData(numberOfNodes, numberOfRegions, numberOfBoundaryRegions, numberOfSpecies)

    # region independent data
    data.F                   .= Blakemore # Boltzmann, FermiDiracOneHalf, Blakemore
    data.temperature          = T
    data.UT                   = (kB * data.temperature) / q
    data.contactVoltage       = [voltageDonor, voltageAcceptor]
    data.chargeNumbers[iphin] = -1
    data.chargeNumbers[iphip] =  1

    # boundary region data
    for ibreg in bregions

        data.bDensityOfStates[ibreg,iphin] = Nc
        data.bDensityOfStates[ibreg,iphip] = Nv
        data.bBandEdgeEnergy[ibreg,iphin]  = Ec
        data.bBandEdgeEnergy[ibreg,iphip]  = Ev

    end

    # interior region data
    for ireg in 1:numberOfRegions

        data.dielectricConstant[ireg]    = εr

        # dos, band edge energy and mobilities
        data.densityOfStates[ireg,iphin] = Nc
        data.densityOfStates[ireg,iphip] = Nv
        data.bandEdgeEnergy[ireg,iphin]  = Ec
        data.bandEdgeEnergy[ireg,iphip]  = Ev
        data.mobility[ireg,iphin]        = mun
        data.mobility[ireg,iphip]        = mup

        # recombination parameters
        data.recombinationRadiative[ireg]            = Radiative
        data.recombinationSRHLifetime[ireg,iphin]    = SRH_LifeTime
        data.recombinationSRHLifetime[ireg,iphip]    = SRH_LifeTime
        data.recombinationSRHTrapDensity[ireg,iphin] = SRH_TrapDensity
        data.recombinationSRHTrapDensity[ireg,iphip] = SRH_TrapDensity
        data.recombinationAuger[ireg,iphin]          = Auger
        data.recombinationAuger[ireg,iphip]          = Auger

    end

    # interior doping
    data.doping[regionDonor,iphin]      = Nd        # data.doping   = [Nd  0.0;
    data.doping[regionIntrinsic,iphin]  = ni        #                  ni   ni;
    data.doping[regionIntrinsic,iphip]  = ni        #                  0.0  Na]
    data.doping[regionAcceptor,iphip]   = Na

    # boundary doping
    data.bDoping[bregionDonor,iphin]    = Nd        # data.bDoping  = [Nd  0.0;
    data.bDoping[bregionAcceptor,iphip] = Na        #                  0.0  Na]

    # print data
    println(data)

    println("*** done\n")

    psi0 = DDFermi.electroNeutralSolutionBoltzmann(grid, data)
    if pyplot
        ################################################################################
        println("Plot electroneutral potential and doping")
        ################################################################################
        #DDFermi.plotEnergies(grid, data)
        DDFermi.plotDoping(grid, data)
        #DDFermi.plotElectroNeutralSolutionBoltzmann(grid, psi0)

        println("*** done\n")
    end

    ################################################################################
    println("Define physics and system")
    ################################################################################

    ## initializing physics environment ##
    physics = VoronoiFVM.Physics(
    data        = data,
    num_species = numberOfSpecies,
    flux        = DDFermi.Sedan!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
    reaction    = DDFermi.reaction!,
    breaction   = DDFermi.breaction!
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
    control.tol_round         = 1.0e-8
    control.max_round         = 5


    println("*** done\n")

    ################################################################################
    println("Compute solution in thermodynamic equilibrium for Boltzmann")
    ################################################################################

    # initialize solution and starting vectors
    initialGuess                   = unknowns(sys)
    solution                       = unknowns(sys)
    @views initialGuess[ipsi,  :] .= psi0 #0.0
    @views initialGuess[iphin, :] .= 0.0
    @views initialGuess[iphip, :] .= 0.0

    DDFermi.solveEquilibriumBoltzmann!(solution, initialGuess, data, grid, control, dense)

    println("*** done\n")

    ################################################################################
    println("Bias loop")
    ################################################################################
    if !(data.F == DDFermi.Boltzmann) # adjust control, when not using Boltzmann
        control.damp_initial      = 0.5
        control.damp_growth       = 1.2
        control.max_iterations    = 30
    end

    maxBias    = data.contactVoltage[bregionAcceptor]
    biasValues = range(0, stop = maxBias, length = 41)
    IV         = zeros(0)

    w_device = 0.5 * μm     # width of device
    z_device = 1.0e-4 * cm  # depth of device


    # put the below values in comments, if using Boltzmann statistics.
    control.damp_initial      = 0.5
    control.damp_growth       = 1.2
    control.max_iterations    = 30

    for Δu in biasValues
        data.contactVoltage[bregionAcceptor] = Δu

        sys.boundary_values[iphin, bregionAcceptor] = Δu
        sys.boundary_values[iphip, bregionAcceptor] = Δu

        # solve!(solution, initialGuess, sys, control = control, tstep = Inf)
        solve!(solution, initialGuess, sys, control = control)

        initialGuess .= solution

        # get IV curve
        factory = VoronoiFVM.TestFunctionFactory(sys)

        # testfunction zero in bregionAcceptor and one in bregionDonor
        tf     = testfunction(factory, [bregionAcceptor], [bregionDonor])
        I      = integrate(sys, tf, solution)

        push!(IV,  abs.(w_device * z_device * (I[iphin] + I[iphip])))

        # plot solution and IV curve
        if pyplot
            #DDFermi.plotDensities(grid, sys, solution, Δu)
            #PyPlot.figure()
            DDFermi.plotDensities(grid, sys, solution , Δu)
            if Δu == 0.0 || Δu == 1.5 Δu == 3
                savefig("psc-densities-nref-$n-deltaU-$Δu.eps")
            end
            #DDFermi.plotIV(biasValues,IV)
        end

    end # bias loop

    # return IV
    println("*** done\n")

end #  main

println("This message should show when the PSC module is successfully recompiled.")

end # module
