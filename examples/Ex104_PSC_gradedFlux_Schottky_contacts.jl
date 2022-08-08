#=
# PSC device with graded interfaces & Schottky contacts (1D).
([source code](SOURCE_URL))

Simulating a three layer PSC device SiO2| MAPI | SiO2 without mobile ions. The simulations are
performed out of equilibrium, stationary and with two junctions between perovskite layer and
transport layers, to which we refer as graded interfaces. Hence, a graded flux discretization
with space dependent band-edge energies and density of states is tested here.
Additionally to that instead of Ohmic contacts we have Schottky contacts.

The parameters can be found in Table S.13, https://arxiv.org/abs/2009.04384. Or here:
https://github.com/barnesgroupICL/Driftfusion/blob/Methods-IonMonger-Comparison/Input_files/IonMonger_default_bulk.csv
=#

module Ex104_PSC_gradedFlux_Schottky_contacts

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot

## function for grading the physical parameters
function grading_parameter!(physicalParameter, coord, regionTransportLayers, regionJunctions, h, heightLayers, lengthLayers, values)
    for ireg in regionTransportLayers

        xcoord                     = lengthLayers[ireg]:lengthLayers[ireg+1]
        physicalParameter[xcoord] .= values[ireg]

    end

    for ireg in regionJunctions

        xcoord   = lengthLayers[ireg]:lengthLayers[ireg+1]
        left     = lengthLayers[ireg]-3
        junction = h[ireg]
        right    = lengthLayers[ireg+2]-3

        gradient = ( physicalParameter[right] - physicalParameter[left] ) / junction

        for index in xcoord
            physicalParameter[index] = physicalParameter[left] + (coord[index] - heightLayers[ireg-1]) * gradient
        end

    end

    return physicalParameter
end

function main(;n = 2, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionDonor           = 1          # n doped region
    regionJunction1       = 2
    regionIntrinsic       = 3          # intrinsic region
    regionJunction2       = 4
    regionAcceptor        = 5          # p doped region
    regions               = [regionDonor, regionJunction1, regionIntrinsic, regionJunction2, regionAcceptor]
    regionTransportLayers = [regionDonor, regionIntrinsic, regionAcceptor]
    regionJunctions       = [regionJunction1, regionJunction2]
    numberOfRegions       = length(regions)

    ## boundary region numbers
    bregionDonor          = 1
    bregionAcceptor       = 2

    ## grid
    h_ndoping             = 9.90e-6 * cm
    h_junction1           = 1.0e-7  * cm
    h_intrinsic           = 4.00e-5 * cm
    h_junction2           = 1.0e-7  * cm
    h_pdoping             = 1.99e-5 * cm
    h                     = [h_ndoping, h_junction1, h_intrinsic, h_junction2, h_pdoping]
    heightLayers          = [h_ndoping,
                             h_ndoping + h_junction1,
                             h_ndoping + h_junction1 + h_intrinsic,
                             h_ndoping + h_junction1 + h_intrinsic + h_junction2,
                             h_ndoping + h_junction1 + h_intrinsic + h_junction2 + h_pdoping]
    refinementfactor      = 2^(n-1)

    coord_ndoping         = collect(range(0.0, stop = h_ndoping, length = 4 * refinementfactor))
    length_n              = length(coord_ndoping)
    coord_junction1       = collect(range(h_ndoping,
                                         stop = h_ndoping + h_junction1,
                                         length = 3 * refinementfactor))
    coord_intrinsic       = collect(range(h_ndoping + h_junction1,
                                         stop = (h_ndoping + h_junction1 + h_intrinsic),
                                         length = 10 * refinementfactor))
    coord_junction2       = collect(range(h_ndoping + h_junction1 + h_intrinsic,
                                         stop = (h_ndoping + h_junction1 + h_intrinsic + h_junction2),
                                         length = 3 * refinementfactor))
    coord_pdoping         = collect(range((h_ndoping + h_junction1 + h_intrinsic + h_junction2),
                                          stop = (h_ndoping + h_junction1 + h_intrinsic + h_junction2 + h_pdoping),
                                          length = 4 * refinementfactor))

    coord                 = glue(coord_ndoping, coord_junction1)
    length_j1             = length(coord)
    coord                 = glue(coord, coord_intrinsic)
    length_i              = length(coord)
    coord                 = glue(coord, coord_junction2)
    length_j2             = length(coord)
    coord                 = glue(coord, coord_pdoping)

    grid                  = simplexgrid(coord)
    numberOfNodes         = length(coord)
    lengthLayers          = [1, length_n, length_j1, length_i, length_j2, numberOfNodes]

    ## set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],        [heightLayers[1]], regionDonor)     # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionJunction1) # first junction   = 2
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionIntrinsic) # intrinsic region = 3
    cellmask!(grid, [heightLayers[3]], [heightLayers[4]], regionJunction2) # sec. junction    = 4
    cellmask!(grid, [heightLayers[4]], [heightLayers[5]], regionAcceptor)  # p-doped region   = 5

    if plotting
        gridplot(grid, Plotter = Plotter)
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
    iphin            = 1 # electron quasi Fermi potential
    iphip            = 2 # hole quasi Fermi potential
    numberOfCarriers = 2

    ## temperature
    T                = 300.0                 *  K

    Eref             = 4.1                   *  eV  # reference energy
    ## band edge energies
    Ec_d             = -4.0                  *  eV  + Eref
    Ev_d             = -6.0                  *  eV  + Eref

    Ec_i             = -3.7                  *  eV  + Eref
    Ev_i             = -5.4                  *  eV  + Eref

    Ec_a             = -3.1                  *  eV  + Eref
    Ev_a             = -5.1                  *  eV  + Eref

    ## these parameters at the junctions for E_\alpha and N_\alpha will be overwritten.
    Ec_j1            = Ec_d;     Ec_j2     = Ec_i
    Ev_j1            = Ev_d;     Ev_j2     = Ev_i

    EC               = [Ec_d, Ec_j1, Ec_i, Ec_j2, Ec_a]
    EV               = [Ev_d, Ev_j1, Ev_i, Ev_j2, Ev_a]

    ## effective densities of state
    Nc_d             = 5.0e19                / (cm^3)
    Nv_d             = 5.0e19                / (cm^3)

    Nc_i             = 8.1e18                / (cm^3)
    Nv_i             = 5.8e18                / (cm^3)

    Nc_a             = 5.0e19                / (cm^3)
    Nv_a             = 5.0e19                / (cm^3)

    Nc_j1            = Nc_d;     Nc_j2      = Nc_i
    Nv_j1            = Nv_d;     Nv_j2      = Nv_i

    NC               = [Nc_d, Nc_j1, Nc_i, Nc_j2, Nc_a]
    NV               = [Nv_d, Nv_j1, Nv_i, Nv_j2, Nv_a]

    ## mobilities
    μn_d             = 3.89                  * (cm^2) / (V * s)
    μp_d             = 3.89                  * (cm^2) / (V * s)

    μn_i             = 6.62e1                * (cm^2) / (V * s)
    μp_i             = 6.62e1                * (cm^2) / (V * s)

    μn_a             = 3.89e-1               * (cm^2) / (V * s)
    μp_a             = 3.89e-1               * (cm^2) / (V * s)

    μn_j1            = μn_d;     μn_j2      = μn_i
    μp_j1            = μp_d;     μp_j2      = μp_i

    μn               = [μn_d, μn_j1, μn_i, μn_j2, μn_a]
    μp               = [μp_d, μp_j1, μp_i, μp_j2, μp_a]

    ## relative dielectric permittivity
    ε_d              = 10.0                  *  1.0
    ε_i              = 24.1                  *  1.0
    ε_a              = 3.0                   *  1.0

    ε_j1             = ε_d;       ε_j2      = ε_a

    ε                = [ε_d, ε_j1, ε_i, ε_j2, ε_a]

    ## radiative recombination
    r0_d             = 0.0e+0               * cm^3 / s
    r0_i             = 1.0e-12              * cm^3 / s
    r0_a             = 0.0e+0               * cm^3 / s

    r0_j1            = r0_i;      r0_j2     = r0_i

    r0               = [r0_d, r0_j1, r0_i, r0_j2, r0_a]

    ## life times and trap densities
    τn_d             = 1.0e100              * s
    τp_d             = 1.0e100              * s

    τn_i             = 3.0e-10              * s
    τp_i             = 3.0e-8               * s
    τn_a             = τn_d
    τp_a             = τp_d

    τn_j1            = τn_i;     τn_j2      = τn_a
    τp_j1            = τp_i;     τp_j2      = τp_a

    τn               = [τn_d, τn_j1, τn_i, τn_j2, τn_a]
    τp               = [τp_d, τp_j1, τp_i, τp_j2, τp_a]

    ## SRH trap energies
    Ei_d             = -5.0                 * eV
    Ei_i             = -4.55                * eV
    Ei_a             = -4.1                 * eV

    Ei_j1            = Ei_i;      Ei_j2     = Ei_a

    EI               = [Ei_d, Ei_j1, Ei_i, Ei_j2, Ei_a]

    ## Auger recombination
    Auger            = 0.0

    ## doping (doping values are from Driftfusion)
    Nd               = 1.03e18             / (cm^3)
    Na               = 1.03e18             / (cm^3)

    ## contact voltage
    voltageAcceptor  =  1.2                 * V

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in data
    data                               = Data(grid, numberOfCarriers)

    ## Possible choices: Stationary, Transient
    data.modelType                     = Stationary

    ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    ## FermiDiracMinusOne, Blakemore
    data.F                            .= Boltzmann

    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = true,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceModelNone,
    ## InterfaceModelSurfaceReco (inner boundary).
    data.boundaryType[bregionDonor]    = SchottkyContact
    data.boundaryType[bregionAcceptor] = SchottkyContact

    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation            .= ScharfetterGummelGraded

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    ## for region dependent parameters
    params                                = Params(grid, numberOfCarriers)
    ## for space dependent parameters
    paramsnodal                           = ParamsNodal(grid, numberOfCarriers)


    params.temperature                    = T
    params.UT                             = (kB * params.temperature) / q
    params.chargeNumbers[iphin]           = -1
    params.chargeNumbers[iphip]           =  1

    ## nodal band-edge energies
    paramsnodal.bandEdgeEnergy[iphin, :]  = grading_parameter!(paramsnodal.bandEdgeEnergy[iphin, :],
                                                              coord, regionTransportLayers, regionJunctions, h,
                                                              heightLayers, lengthLayers, EC)
    paramsnodal.bandEdgeEnergy[iphip, :]  = grading_parameter!(paramsnodal.bandEdgeEnergy[iphip, :],
                                                              coord, regionTransportLayers, regionJunctions, h,
                                                              heightLayers, lengthLayers, EV)
    ## nodal effective density of states
    paramsnodal.densityOfStates[iphin, :] = grading_parameter!(paramsnodal.densityOfStates[iphin, :],
                                                              coord, regionTransportLayers, regionJunctions, h,
                                                              heightLayers, lengthLayers, NC)
    paramsnodal.densityOfStates[iphip, :] = grading_parameter!(paramsnodal.densityOfStates[iphip, :],
                                                              coord, regionTransportLayers, regionJunctions, h,
                                                              heightLayers, lengthLayers, NV)
    ## region dependent data
    for ireg in 1:numberOfRegions

        ## mobility
        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0
        ## recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, data, EI[ireg])
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

    end

    ## interior doping
    params.doping[iphin, regionDonor]        = Nd
    params.doping[iphip, regionAcceptor]     = Na

    ## boundary doping
    params.bDoping[iphip, bregionAcceptor]   = Na
    params.bDoping[iphin, bregionDonor]      = Nd

    ## values for the schottky contacts
    params.SchottkyBarrier[bregionDonor]     =  0.0
    params.SchottkyBarrier[bregionAcceptor]  = -5.0  * eV  -  (-4.1  * eV) # difference between boundary Fermi level
    params.bVelocity                         = [1.0e7 * cm/s    0.0   * cm/s;
                                                0.0   * cm/s    1.0e7 * cm/s]

    data.params                              = params
    data.paramsnodal                         = paramsnodal
    ctsys                                    = System(grid, data, unknown_storage=unknown_storage)

    if test == false
        show_params(ctsys)
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define boundary conditions")
    end
    ################################################################################

    ## set Schottky contacts.
    set_contact!(ctsys, bregionAcceptor, Δu = 0.0)
    set_contact!(ctsys, bregionDonor,    Δu = 0.0)

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
    control.max_iterations    = 300
    control.tol_absolute      = 1.0e-13
    control.tol_relative      = 1.0e-13
    control.handle_exceptions = true
    control.tol_round         = 1.0e-13

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    control.damp_initial  = 0.5
    control.damp_growth   = 1.61 # >= 1
    control.max_round     = 5

    initialGuess          = unknowns(ctsys)
    solution              = unknowns(ctsys)

    solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    initialGuess         .= solution

    if plotting
        label_solution, label_density, label_energy = set_plotting_labels(data)

        plot_energies(Plotter,  grid, data, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution, "Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter,  grid, data, solution, "Equilibrium", label_solution)
        Plotter.figure()
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    data.calculationType = OutOfEquilibrium

    control.damp_initial = 0.9
    control.damp_growth  = 1.61 # >= 1
    control.max_round    = 7

    maxBias    = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 25)

    for Δu in biasValues
        if test == false
            println("Bias value: Δu = $(Δu)")
        end

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        initialGuess .= solution

    end # bias loop

    ## plotting
    if plotting
        plot_energies(Plotter,  grid, data, solution, "Applied voltage Δu = $maxBias", label_energy)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution, "Applied voltage Δu = $maxBias", label_density)
        Plotter.figure()
        plot_solution(Plotter,  grid, data, solution, "Applied voltage Δu = $maxBias", label_solution)
    end

    if test == false
        println("*** done\n")
    end

    testval = solution[data.index_psi, 20]
    return testval

end #  main

function test()
    testval=0.11725154137943199
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
