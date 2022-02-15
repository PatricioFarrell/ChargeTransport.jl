#=
# GaAs diode (1D).
([source code](SOURCE_URL))

We simulate charge transport in a GaAs pin diode, where use the van Roosbroeck
system of equations as charge transport model. The unknowns are given by the quasi Fermi
potentials of electrons and holes $\varphi_n$, $\varphi_p$ and the electric potential $\psi$.
The simulations are performed out of equilibrium and for the
stationary problem.
=#

module TestAllocations

using VoronoiFVM       # PDE solver with a FVM spatial discretization
using ChargeTransport  # drift-diffusion solver
using ExtendableGrids  # grid initializer
using GridVisualize    # grid visualizer
using PyPlot           # solution visualizer
using BenchmarkTools


## This function is used to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_intrinsic, h_pdoping)
    coord_ndoping    = collect(range(0.0, stop = h_ndoping, length = 3 * refinementfactor))
    coord_intrinsic  = collect(range(h_ndoping, stop = (h_ndoping + h_intrinsic), length = 3 * refinementfactor))
    coord_pdoping    = collect(range((h_ndoping + h_intrinsic), stop = (h_ndoping + h_intrinsic + h_pdoping), length = 3 * refinementfactor))
    coord            = glue(coord_ndoping, coord_intrinsic)
    coord            = glue(coord, coord_pdoping)

    return coord
end


function main(;n = 4, Plotter = PyPlot, plotting = false, verbose = false, test = true, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionAcceptor          = 1          # p doped region
    regionIntrinsic         = 2          # intrinsic region
    regionDonor             = 3          # n doped region
    regions                 = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions         = length(regions)

    ## boundary region numbers
    bregionAcceptor         = 1
    bregionDonor            = 2
    bregions                = [bregionAcceptor, bregionDonor]
    numberOfBoundaryRegions = length(bregions)

    ## grid
    refinementfactor        = 2^(n-1)
    h_pdoping               = 2 * μm
    h_intrinsic             = 2 * μm
    h_ndoping               = 2 * μm
    coord                   = initialize_pin_grid(refinementfactor,
                                                  h_pdoping,
                                                  h_intrinsic,
                                                  h_ndoping)

    grid                    = simplexgrid(coord)

    ## cellmask! for defining the subregions and assigning region number
    cellmask!(grid, [0.0 * μm],                [h_pdoping],                           regionAcceptor)  # p-doped region = 1
    cellmask!(grid, [h_pdoping],               [h_pdoping + h_intrinsic],             regionIntrinsic) # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor)     # n-doped region = 3

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

    ## set indices of the quasi Fermi potentials
    iphin              = 1 # electron quasi Fermi potential
    iphip              = 2 # hole quasi Fermi potential
    numberOfCarriers   = 2

    # We define the physical data.
    Ec                 = 1.424                *  eV
    Ev                 = 0.0                  *  eV
    Nc                 = 4.351959895879690e17 / (cm^3)
    Nv                 = 9.139615903601645e18 / (cm^3)
    mun                = 8500.0               * (cm^2) / (V * s)
    mup                = 400.0                * (cm^2) / (V * s)
    εr                 = 12.9                 *  1.0              # relative dielectric permittivity of GAs
    T                  = 300.0                *  K

    ## recombination parameters
    Auger             = 1.0e-29              * cm^6 / s
    SRH_TrapDensity   = 1.0e10               / cm^3
    SRH_LifeTime      = 1.0                  * ns
    Radiative         = 1.0e-10              * cm^3 / s

    ## doping
    dopingFactorNd    = 1.0
    dopingFactorNa    = 0.46
    Nd                = dopingFactorNd * Nc
    Na                = dopingFactorNa * Nv

    ## intrinsic concentration
    ni                = sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * kB * T))

    ## contact voltages: we impose an applied voltage only on one boundary.
    ## At the other boundary the applied voltage is zero.
    voltageAcceptor   = 1.5                  * V

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    text = ["StandardFuncSet", "Functions"]
    i = 0
    for functionset in [StandardFuncSet, Function]
        i = i + 1
        println(" ")
        println(text[i], ": ")

        if functionset==Function
            data                            = Data(grid, numberOfCarriers, statfunctions=Function)
        else
            data                            = Data(grid, numberOfCarriers)
        end


        data.model_type                     = Stationary
        data.F                             .= Boltzmann

        data.bulk_recombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                     bulk_recomb_Auger = true,
                                                                     bulk_recomb_radiative = true,
                                                                     bulk_recomb_SRH = true)

        data.boundary_type[bregionAcceptor] = OhmicContact
        data.boundary_type[bregionDonor]    = OhmicContact

        data.flux_approximation             = ExcessChemicalPotential

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Define Params and fill in physical parameters")
        end
        ################################################################################

        params                                              = Params(grid, numberOfCarriers)

        params.temperature                                  = T
        params.UT                                           = (kB * params.temperature) / q
        params.chargeNumbers[iphin]                         = -1
        params.chargeNumbers[iphip]                         =  1

        for ibreg in 1:numberOfBoundaryRegions   # boundary region data
            params.bDensityOfStates[iphin, ibreg]           = Nc
            params.bDensityOfStates[iphip, ibreg]           = Nv
            params.bBandEdgeEnergy[iphin, ibreg]            = Ec
            params.bBandEdgeEnergy[iphip, ibreg]            = Ev
        end

        for ireg in 1:numberOfRegions           # interior region data

            params.dielectricConstant[ireg]                 = εr

            ## effective DOS, band-edge energy and mobilities
            params.densityOfStates[iphin, ireg]             = Nc
            params.densityOfStates[iphip, ireg]             = Nv
            params.bandEdgeEnergy[iphin, ireg]              = Ec
            params.bandEdgeEnergy[iphip, ireg]              = Ev
            params.mobility[iphin, ireg]                    = mun
            params.mobility[iphip, ireg]                    = mup

            ## recombination parameters
            params.recombinationRadiative[ireg]             = Radiative
            params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
            params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
            params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
            params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity
            params.recombinationAuger[iphin, ireg]          = Auger
            params.recombinationAuger[iphip, ireg]          = Auger

        end

        ## interior doping
        params.doping[iphin, regionDonor]                   = Nd     # data.doping   = [0.0  Na;
        params.doping[iphin, regionIntrinsic]               = ni     #                  ni   0.0;
        params.doping[iphip, regionIntrinsic]               = 0.0    #                  Nd  0.0]
        params.doping[iphip, regionAcceptor]                = Na

        ## boundary doping
        params.bDoping[iphin, bregionDonor]                 = Nd     # data.bDoping  = [0.0  Na;
        params.bDoping[iphip, bregionAcceptor]              = Na     #                  Nd  0.0]


        data.params                                         = params
        ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

        if test == false
            show_params(ctsys)
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Define outerior boundary conditions")
        end
        ################################################################################
        set_contact!(ctsys, bregionAcceptor, Δu = 0.0)
        set_contact!(ctsys, bregionDonor, Δu = 0.0)

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
        control.damp_initial      = 0.5
        control.damp_growth       = 1.21
        control.max_iterations    = 250
        control.tol_absolute      = 1.0e-14
        control.tol_relative      = 1.0e-14
        control.handle_exceptions = true
        control.tol_round         = 1.0e-8
        control.max_round         = 3

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Compute solution in thermodynamic equilibrium")
        end
        ################################################################################

        ## initialize solution and starting vectors
        initialGuess          = unknowns(ctsys)
        solution              = unknowns(ctsys)

        solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

        initialGuess         .= solution


        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Bias loop")
        end
        ################################################################################

        # Set calculation type to OutOfEquilibrium for starting with respective simulation.
        ctsys.data.calculation_type      = OutOfEquilibrium

        maxBias    = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
        biasValues = range(0, stop = maxBias, length = 4)

        for Δu in biasValues

            println("Δu  = ", Δu )

            ## set non equilibrium boundary conditions
            set_contact!(ctsys, bregionAcceptor, Δu = Δu)

            @btime solve!($solution, $initialGuess, $ctsys, control = $control, tstep = $Inf)

            initialGuess .= solution


        end # bias loop

    end


end #  main


end # module
