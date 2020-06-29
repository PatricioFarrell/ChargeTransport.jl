"""
Simulating charge transport in a GAs pin diode with (psi, n, p) as set of unknowns.
"""

module PINdensities

using VoronoiFVM
using DDFermi
using ExtendableGrids
using PyPlot; PyPlot.pygui(true)
using Printf

# function for initializing the grid for a possble extension to other
# p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_intrinsic, h_pdoping)
    # without filter! we would have 2.0e-6 and 4.0e-6 twice
    #  -> yields to singularexception
    coord_ndoping    = collect(range(0, stop = h_ndoping, length = 3 * refinementfactor))
    coord_intrinsic  = collect(range(h_ndoping, stop = (h_ndoping + h_intrinsic), length = 3 * refinementfactor))
    coord_intrinsic  = filter!(x->x≠h_pdoping, coord_intrinsic)
    coord_pdoping    = collect(range((h_ndoping + h_intrinsic), stop = (h_ndoping + h_intrinsic + h_pdoping), length = 3 * refinementfactor))
    coord_pdoping    = filter!(x->x≠(h_ndoping+h_intrinsic), coord_pdoping)
    coord            = vcat(coord_ndoping, coord_intrinsic)
    coord            = vcat(coord, coord_pdoping)
    return coord
end


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

    # grid
    refinementfactor = 2^(n-1)
    h_ndoping        = 2 * μm
    h_intrinsic      = 2 * μm
    h_pdoping        = 2 * μm
    coord            = initialize_pin_grid(refinementfactor,
    h_ndoping,
    h_intrinsic,
    h_pdoping)

    grid             = VoronoiFVM.Grid(coord)
    numberOfNodes    = length(coord)
    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm], [h_ndoping], regionDonor)        # n-doped region = 1
    cellmask!(grid, [h_ndoping], [h_ndoping + h_intrinsic], regionIntrinsic)    # intrinsic region = 2
    cellmask!(grid, [h_ndoping + h_intrinsic], [h_ndoping + h_intrinsic + h_pdoping], regionAcceptor)     # p-doped region = 3

    #ExtendableGrids.plot(VoronoiFVM.Grid(coord, collect(0.0:0.25 * μm:0.5 * μm)), Plotter = PyPlot, p = PyPlot.plot()) # Plot grid
    println("*** done\n")

    ################################################################################
    println("Define physical parameters and model")
    ################################################################################

    # indices
    in           = 1
    ip           = 2
    ipsi            = 3
    species         = [in, ip, ipsi]

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
    data.F                    = Boltzmann # Boltzmann, FermiDiracOneHalf, Blakemore
    data.temperature          = T
    data.UT                   = (kB * data.temperature) / q
    data.contactVoltage       = [voltageDonor, voltageAcceptor]
    data.chargeNumbers[in] = -1
    data.chargeNumbers[ip] =  1

    # boundary region data
    for ibreg in 1:numberOfBoundaryRegions

        data.bDensityOfStates[ibreg,in] = Nc
        data.bDensityOfStates[ibreg,ip] = Nv
        data.bBandEdgeEnergy[ibreg,in]  = Ec
        data.bBandEdgeEnergy[ibreg,ip]  = Ev

    end

    # interior region data
    for ireg in 1:numberOfRegions

        data.dielectricConstant[ireg]    = εr

        # dos, band edge energy and mobilities
        data.densityOfStates[ireg,in] = Nc
        data.densityOfStates[ireg,ip] = Nv
        data.bandEdgeEnergy[ireg,in]  = Ec
        data.bandEdgeEnergy[ireg,ip]  = Ev
        data.mobility[ireg,in]        = mun
        data.mobility[ireg,ip]        = mup

        # recombination parameters
        data.recombinationRadiative[ireg]            = Radiative
        data.recombinationSRHLifetime[ireg,in]    = SRH_LifeTime
        data.recombinationSRHLifetime[ireg,ip]    = SRH_LifeTime
        data.recombinationSRHTrapDensity[ireg,in] = SRH_TrapDensity
        data.recombinationSRHTrapDensity[ireg,ip] = SRH_TrapDensity
        data.recombinationAuger[ireg,in]          = Auger
        data.recombinationAuger[ireg,ip]          = Auger

    end

    # interior doping
    data.doping[regionDonor,in]      = Nd        # data.doping   = [Nd  0.0;
    data.doping[regionIntrinsic,in]  = ni        #                  ni   ni;
    data.doping[regionIntrinsic,ip]  = ni        #                  0.0  Na]
    data.doping[regionAcceptor,ip]   = Na

    # boundary doping
    data.bDoping[bregionDonor,in]    = Nd        # data.bDoping  = [Nd  0.0;
    data.bDoping[bregionAcceptor,ip] = Na        #                  0.0  Na]

    # print data
    println(data)

    println("*** done\n")
    psi0 = zeros(length(coord))
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
    flux        = DDFermi.ScharfetterGummelDensities!, #Sedan!, ScharfetterGummel!, diffusionEnhanced!, KopruckiGaertner!
    reaction    = DDFermi.reactionDensities!,
    breaction   = DDFermi.breactionDensities!
    )

    if dense
        sys = VoronoiFVM.System(grid, physics, unknown_storage = :dense)
    else
        sys = VoronoiFVM.System(grid, physics, unknown_storage = :sparse)
    end

    # enable all three species in all regions
    enable_species!(sys, ipsi,  regions)
    enable_species!(sys, in, regions)
    enable_species!(sys, ip, regions)

    bfaceregions  = grid[BFaceRegions]
    bfacenodes    = grid[BFaceNodes]

    for icc = 1:data.numberOfSpecies-1
        EDonor = data.bBandEdgeEnergy[bfaceregions[bregionDonor],icc] + data.bandEdgeEnergyNode[bfacenodes[bregionDonor],icc]
        EAcceptor = data.bBandEdgeEnergy[bfaceregions[bregionAcceptor],icc] + data.bandEdgeEnergyNode[bfacenodes[bregionAcceptor],icc]
        etaDonor = data.chargeNumbers[icc] / data.UT * ( (data.contactVoltage[bregionDonor]- psi0[1]) + EDonor / q )
        etaAcceptor = data.chargeNumbers[icc] / data.UT * ( (data.contactVoltage[bregionAcceptor]- psi0[length(coord)]) + EAcceptor / q )

        sys.boundary_values[icc,  bregionDonor]    = data.bDensityOfStates[bregionDonor, icc] * data.F(etaDonor)
        sys.boundary_factors[icc, bregionDonor]    = VoronoiFVM.Dirichlet

        sys.boundary_values[icc,  bregionAcceptor] = data.bDensityOfStates[bregionAcceptor, icc] * data.F(etaAcceptor)
        sys.boundary_factors[icc, bregionAcceptor] = VoronoiFVM.Dirichlet
    end


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
    @views initialGuess[ipsi,  :]  = psi0
    @views initialGuess[in, :] .= 0.0
    @views initialGuess[ip, :] .= 0.0

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
    biasValues = range(0, stop = maxBias, length = 51)
    IV         = zeros(0)

    w_device = 0.5 * μm     # width of device
    z_device = 1.0e-4 * cm  # depth of device


    # put the below values in comments, if using Boltzmann statistics.
    control.damp_initial      = 0.5
    control.damp_growth       = 1.2
    control.max_iterations    = 30



    for Δu in biasValues
        data.contactVoltage[bregionAcceptor] = Δu

        for icc = 1:data.numberOfSpecies-1
            EAcceptor = data.bBandEdgeEnergy[bfaceregions[bregionAcceptor],icc] + data.bandEdgeEnergyNode[bfacenodes[bregionAcceptor],icc]
            etaAcceptor = data.chargeNumbers[icc] / data.UT * ( (data.contactVoltage[bregionAcceptor]- psi0[length(coord)]) + EAcceptor / q )

            sys.boundary_values[icc,  bregionAcceptor] = data.bDensityOfStates[bregionAcceptor, icc] * data.F(etaAcceptor)
        end

        solve!(solution, initialGuess, sys, control = control, tstep = Inf)

        initialGuess .= solution

        # get IV curve
        factory = VoronoiFVM.TestFunctionFactory(sys)

        # testfunction zero in bregionAcceptor and one in bregionDonor
        tf     = testfunction(factory, [bregionAcceptor], [bregionDonor])
        I      = integrate(sys, tf, solution)

        push!(IV,  abs.(w_device * z_device * (I[in] + I[ip])))

        # plot solution and IV curve
        if pyplot
            DDFermi.plotIV(biasValues,IV)
        end

    end # bias loop

    # return IV
    println("*** done\n")

end #  main

println("This message should show when the PIN module is successfully recompiled.")

end # module
