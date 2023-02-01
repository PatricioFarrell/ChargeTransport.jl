#=
# CIGS pn junction: stationary with traps and Schottky contacts.
([source code](SOURCE_URL))
Simulating stationary charge transport in a pn junction with hole traps and mixed Schottky/Ohmic contact conditions.
=#

module Ex112_CIGS_WithTraps

using ChargeTransport
using ExtendableGrids
using PyPlot

## function to initialize the grid for a possble extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_pdoping_left, h_pdoping_trap, h_pdoing_right)
    coord_ndoping    = collect(range(0.0, stop = h_ndoping, length = 2 * refinementfactor))
    coord_pdoping_left  = collect(range(h_ndoping, stop = (h_ndoping + h_pdoping_left), length = 3 * refinementfactor))
    coord_pdoping_plus  = collect(range((h_ndoping + h_pdoping_left),
                                        stop = (h_ndoping + h_pdoping_left + h_pdoping_trap),
                                        length =  refinementfactor))
    coord_pdoping_right = collect(range((h_ndoping + h_pdoping_left + h_pdoping_trap),
                                        stop = (h_ndoping + h_pdoping_left + h_pdoping_trap + h_pdoing_right),
                                        length = 3 * refinementfactor))
    coord            = glue(coord_ndoping, coord_pdoping_left)
    coord            = glue(coord, coord_pdoping_plus)
    coord            = glue(coord, coord_pdoping_right)

    return coord
end

function main(;n = 3, Plotter = PyPlot, plotting = false, verbose = false, test = false)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionDonor             = 1                           # n doped region
    regionAcceptorLeft      = 2                           # p doped region
    regionAcceptorTrap      = 3                           # p doped region with trap
    regionAcceptorRight     = 4                           # p doped region
    regions                 = [regionDonor, regionAcceptorLeft, regionAcceptorTrap, regionAcceptorRight]
    numberOfRegions         = length(regions)

    ## boundary region numbers
    bregionDonor            = 1
    bregionAcceptor         = 2
    bregions                = [bregionDonor, bregionAcceptor]
    numberOfBoundaryRegions = length(bregions)

    ## grid
    refinementfactor        = 2^(n-1)
    h_ndoping               = 0.5    * μm
    h_pdoping_left          = 1.0    * μm
    h_pdoping_trap          = 0.01   * μm
    h_pdoing_right          = 1.0    * μm
    w_device                = 0.5    * μm  # width of device
    z_device                = 1.0e-4 * cm  # depth of device
    h_total                 = h_ndoping + h_pdoping_left + h_pdoping_trap + h_pdoing_right
    coord                   = initialize_pin_grid(refinementfactor,
                                                  h_ndoping,
                                                  h_pdoping_left,
                                                  h_pdoping_trap,
                                                  h_pdoing_right)

    grid                    = simplexgrid(coord)

    ## set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm], [h_ndoping], regionDonor) # n doped
    cellmask!(grid, [h_ndoping], [h_ndoping + h_pdoping_left], regionAcceptorLeft) # p doped
    cellmask!(grid, [h_ndoping + h_pdoping_left], [h_ndoping + h_pdoping_left + h_pdoping_trap], regionAcceptorTrap) # p doped with traps
    cellmask!(grid, [h_ndoping + h_pdoping_left + h_pdoping_trap], [h_total], regionAcceptorRight) # p doped

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

    iphin             = 1 # index electron quasi Fermi potential
    iphip             = 2 # index hole quasi Fermi potential
    numberOfCarriers  = 2 # electrons and holes

    ## physical data
    T                 = 300.0                *  K

    ## band edge energies
    Ec_ZnO            = 3.4                  *  eV
    Ev_ZnO            = 0.0                  *  eV

    Ec_CIGS           = 3.4                  *  eV
    Ev_CIGS           = 2.3                  *  eV

    EC                = [Ec_ZnO, Ec_CIGS, Ec_CIGS, Ec_CIGS]
    EV                = [Ev_ZnO, Ev_CIGS, Ev_CIGS, Ev_CIGS]

    ## effective densities of states
    Nc                = 4.351959895879690e17 / (cm^3)
    Nv                = 9.139615903601645e18 / (cm^3)

    NC                = [Nc, Nc, Nc, Nc]
    NV                = [Nv, Nv, Nv, Nv]

    ## mobilities
    mun_CIGS          = 100.0                * (cm^2) / (V * s)
    mup_CIGS          = 25                   * (cm^2) / (V * s)
    mun_ZnO           = 100                  * (cm^2) / (V * s)
    mup_ZnO           = 25                   * (cm^2) / (V * s)

    μn                = [mun_ZnO, mun_CIGS, mun_CIGS, mun_CIGS]
    μp                = [mup_ZnO, mup_CIGS, mup_CIGS, mup_CIGS]

    ## relative dielectric permittivity
    εr_CIGS           = 13.6                 *  1.0
    εr_ZnO            = 9                    *  1.0

    ε                 = [εr_ZnO, εr_CIGS, εr_CIGS, εr_CIGS]

    ## trap information
    zt                = 1 # hole traps
    Et                = 2.8                  *  eV
    Nt                = 5e14                 / (cm^3)
    NT                = [Nt, Nt, Nt, Nt]

    ## recombination information parameters
    ni_CIGS           = sqrt(Nc * Nv) * exp(-(Ec_CIGS - Ev_CIGS) / (2 * kB * T)) # intrinsic concentration
    n0_CIGS           = Nc * Boltzmann( (Et-Ec_CIGS) / (kB*T) )                  # Boltzmann equilibrium concentration
    p0_CIGS           = ni_CIGS^2 / n0_CIGS                                      # Boltzmann equilibrium concentration
    ni_ZnO            = sqrt(Nc * Nv) * exp(-(Ec_ZnO - Ev_ZnO) / (2 * kB * T))   # intrinsic concentration
    n0_ZnO            = Nc * Boltzmann( (Et-Ec_ZnO) / (kB*T) )                   # Boltzmann equilibrium concentration
    p0_ZnO            = ni_ZnO^2 / n0_ZnO                                        # Boltzmann equilibrium concentration

    p0                = [p0_ZnO, p0_CIGS, p0_CIGS, p0_CIGS]
    n0                = [n0_ZnO, n0_CIGS, n0_CIGS, n0_CIGS]

    Auger             = 1.0e-29  * cm^6 / s
    SRH_LifeTime      = 1.0e-3   * ns
    Radiative         = 1.0e-10  * cm^3 / s

    ## Schottky contact information
    An                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    Ap                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    vn                = An * T^2 / (q*Nc)
    vp                = Ap * T^2 / (q*Nv)
    barrier           = 0.7 * eV

    ## doping information
    Nd                = 1.0e18 / (cm^3)
    Na                = 5.5e15 / (cm^3)

    ## we will impose this applied voltage on one boundary
    voltageAcceptor   = 1.0 * V

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)
    data.modelType                      = Stationary
    data.F                             .= FermiDiracOneHalfTeSCA
    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = true,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)
    data.boundaryType[bregionAcceptor]  = SchottkyContact
    data.boundaryType[bregionDonor]     = OhmicContact
    data.fluxApproximation             .= ExcessChemicalPotential

    ## pass trap data in stationary setting since there is no separate trap species
    add_trap_density!(data=data, zt = zt, Nt = NT)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    ## physical parameters
    params                                              = Params(grid, numberOfCarriers)
    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1

    for ibreg in 1:numberOfBoundaryRegions   # boundary region data
        params.bDensityOfStates[iphin, ibreg]           = Nc
        params.bDensityOfStates[iphip, ibreg]           = Nv
    end

    params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_ZnO
    params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_ZnO
    params.bBandEdgeEnergy[iphin, bregionAcceptor]      = Ec_CIGS
    params.bBandEdgeEnergy[iphip, bregionAcceptor]      = Ev_CIGS

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = NC[ireg]
        params.densityOfStates[iphip, ireg]             = NV[ireg]
        params.bandEdgeEnergy[iphin, ireg]              = EC[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = EV[ireg]
        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg]             = Radiative
        params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = n0[ireg]
        params.recombinationSRHTrapDensity[iphip, ireg] = p0[ireg]
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

    end

    ## doping -- since we do not set any doping for the traps it is automatically zero
    params.doping[iphin, regionDonor]                   = Nd
    params.doping[iphip, regionAcceptorLeft]            = Na
    params.doping[iphip, regionAcceptorTrap]            = Na
    params.doping[iphip, regionAcceptorRight]           = Na

    ## boundary doping
    params.bDoping[iphin, bregionDonor]                 = Nd
    params.bDoping[iphip, bregionAcceptor]              = Na

    ## values for the schottky contacts
    params.SchottkyBarrier[bregionAcceptor]             = barrier
    params.bVelocity[iphin, bregionAcceptor]            = vn
    params.bVelocity[iphip, bregionAcceptor]            = vp

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=:sparse)

    if test == false
        show_params(ctsys)
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control                = NewtonControl()
    control.verbose        = verbose
    control.tol_absolute   = 1.0e-10
    control.tol_relative   = 1.0e-10
    control.tol_round      = 1.0e-7
    control.damp_initial   = 0.5
    control.damp_growth    = 1.2
    control.max_iterations = 30
    control.max_round      = 3

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    ## initialize starting vector
    initialGuess          = unknowns(ctsys)

    ## solve thermodynamic equilibrium and update initial guess
    solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    initialGuess         .= solution

    if plotting
        label_solution, label_density, label_energy = set_plotting_labels(data)

        ## ##### set legend for plotting routines #####
        plot_energies(Plotter, ctsys, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "Equilibrium", label_solution)
        Plotter.figure()
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Stationary bias loop")
    end
    ################################################################################

    ## set calculationType to OutOfEquilibrium for starting with respective simulation.
    data.calculationType = OutOfEquilibrium      # Rn = Rp = R, since the model type is stationary
    endVoltage           = voltageAcceptor       # final bias value
    biasValues           = collect(range(0, stop = endVoltage, length = 52))

    IV                   = zeros(0)
    chargeDensities      = zeros(0)

    for i in eachindex(biasValues)

        Δu = biasValues[i] # bias

        ## Apply new voltage: set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("bias: Δu = $(Δu)")
        end

        ## solve time step problems with timestep Δt
        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)
        initialGuess .= solution

        ## save IV data
        current = get_current_val(ctsys, solution)
        push!(IV, w_device * z_device * current)

        ## store charge density in donor region (ZnO)
        push!(chargeDensities, charge_density(ctsys, solution)[regionDonor])


    end # bias loop

    ## compute static capacitance: check this is correctly computed
    staticCapacitance = diff(chargeDensities) ./ diff(biasValues)

    ## plot solution and IV curve
    if plotting
        plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage) V", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"bias \$\\Delta u\$ = $(endVoltage) V", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage) V", label_solution)
        Plotter.figure()
        plot_IV(Plotter, biasValues,IV, "bias \$\\Delta u\$ = $(biasValues[end]) V", plotGridpoints = true)
        Plotter.figure()
        plot_IV(Plotter, biasValues,chargeDensities, "bias \$\\Delta u\$ = $(biasValues[end]) V", plotGridpoints = true)
        Plotter.title("Charge density in donor region")
        Plotter.ylabel("Charge density [C]")
        Plotter.figure()
        plot_IV(Plotter, biasValues, staticCapacitance, "bias \$\\Delta u\$ = $(biasValues[end]) V", plotGridpoints = true)
        Plotter.title("Static capacitance in donor region")
        Plotter.ylabel("Static capacitance [C/V]")

    end

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solution))/length(solution) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = 1.2612532579953044
    main(test = true) ≈ testval
end

if test == false
    println("This message should show when this module has successfully recompiled.")
end


end # module

