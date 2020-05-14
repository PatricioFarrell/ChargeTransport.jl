"""
Simulating charge transport in a GAs pin diode in the non-stationary case.
"""

module PINtimedependent

using VoronoiFVM
using DDFermi
using ExtendableGrids
using PyPlot; PyPlot.pygui(true)
using Printf


function main(;n = 3, pyplot = false, verbose = false, dense = true)

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

    # grid
    refinementfactor = 2^(n-1)
    coord_pdoping    = collect(range(0, stop = 2 * μm, length = 3 * refinementfactor))
    coord_intrinsic  = collect(range(2* μm, stop = 4* μm, length = 3 * refinementfactor))
    coord_intrinsic  = filter!(x->x≠2.0e-6, coord_intrinsic)
    coord_ndoping    = collect(range(4* μm, stop = 6* μm, length = 3 * refinementfactor))
    coord_ndoping    = filter!(x->x≠4.0e-6, coord_ndoping)
    coord            = vcat(coord_pdoping, coord_intrinsic)
    coord            = vcat(coord, coord_ndoping)
    grid             = VoronoiFVM.Grid(coord)

    numberOfNodes = length(coord)

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm], [2.0 * μm], regionDonor)        # n-doped region = 1
    cellmask!(grid, [2.0 * μm], [4.0 * μm], regionIntrinsic)    # intrinsic region = 2
    cellmask!(grid, [4.0 * μm], [6.0 * μm], regionAcceptor)     # p-doped region = 3

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
    numberOfBoundaryRegions = length(bregions)
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
    data.F                    = Blakemore # Boltzmann, FermiDiracOneHalf, Blakemore
    data.temperature          = T
    data.UT                   = (kB * data.temperature) / q
    data.contactVoltage       = [voltageDonor, voltageAcceptor]
    data.chargeNumbers[iphin] = -1
    data.chargeNumbers[iphip] =  1

    # boundary region data
    for ibreg in 1:numberOfBoundaryRegions

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


    ################################################################################
    println("Define physics and system")
    ################################################################################

    ## initializing physics environment ##
    physics = VoronoiFVM.Physics(
    data        = data,
    num_species = numberOfSpecies,
    flux        = DDFermi.diffusionEnhanced!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
    reaction    = DDFermi.reaction!,
    breaction   = DDFermi.breaction!,
    storage     = DDFermi.storage!
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
    control.damp_initial      = 0.5
    control.damp_growth       = 1.2
    control.max_iterations    = 250
    control.tol_absolute      = 1.0e-14
    control.tol_relative      = 1.0e-14
    control.handle_exceptions = true
    control.tol_round         = 1.0e-8
    control.max_round         = 5


    println("*** done\n")


    ################################################################################
    println("Declare initial guess")
    ################################################################################

    psi0 = DDFermi.electroNeutralSolutionBoltzmann(grid, data)

    # initialize solution and starting vectors
    initialGuessBiasLoop                   = unknowns(sys)
    initialGuessTimeLoop                   = unknowns(sys)
    solution                               = unknowns(sys)
    @views initialGuessBiasLoop[ipsi,  :] .= psi0
    @views initialGuessBiasLoop[iphin, :] .= 0.0
    @views initialGuessBiasLoop[iphip, :] .= 0.0

    println("*** done\n")

    ################################################################################
    println("Bias and time loop")
    ################################################################################

    maxBias    = data.contactVoltage[bregionAcceptor]
    biasValues = range(0, stop = maxBias, length = 41)
    IV         = zeros(0)
    timeCoord  = zeros(0)


    w_device = 0.5 * μm# width of device
    z_device = 1.0e-4 * cm  # depth of device

    timeEnd  = 5.0
    timeStep = 0.5

    for Δu in biasValues

        time = 0.0 # for next bias loop

        data.contactVoltage[bregionAcceptor] = Δu

        sys.boundary_values[iphin, bregionAcceptor] = Δu
        sys.boundary_values[iphip, bregionAcceptor] = Δu

        solve!(solution, initialGuessBiasLoop, sys, control = control, tstep = Inf)
        initialGuessTimeLoop .= solution
        initialGuessBiasLoop .= solution

        while time < timeEnd
            push!(timeCoord, time)

            if verbose
                @printf("time = %g\n",time)
            end

            solve!(solution, initialGuessTimeLoop, sys, control = control, tstep = timeStep)
            initialGuessTimeLoop .= solution

            # plot solution and IV curve
            if pyplot
                DDFermi.plotSolution(grid, sys, solution, Δu, time)
            end

            time = time + timeStep
            #timeStep*= 1.4
        end # time loop

    end # bias loop


    println("*** done\n")

end #  main

println("This message should show when the time-dependent PIN module is successfully recompiled.")

end # module
