#=
# PSC device with uniform generation rate (1D).
([source code](SOURCE_URL))

Simulating a three layer PSC device Pedot| MAPI | PCBM. The simulations are performed out of
equilibrium, time-dependent, with abrupt interfaces and with a uniform generation.
A linear I-V measurement protocol is included and the corresponding solution vectors after
the scan protocol can be depicted.

The paramters are from Calado et al.:
https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/pedotpss_mapi_pcbm.csv.
(with adjustments on layer lengths)
=#

module Ex108_PSC_uniform_Generation

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot

function main(;n = 4, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:dense)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionAcceptor  = 1                           # p doped region
    regionIntrinsic = 2                           # intrinsic region
    regionDonor     = 3                           # n doped region
    regions         = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions = length(regions)

    ## boundary region numbers
    bregionAcceptor = 1
    bregionDonor    = 2

    ## grid
    h_pdoping       = 3.00e-6 * cm + 1.0e-7 *cm # add 1.e-7 cm to this layer for agreement with grid of Driftfusion
    h_intrinsic     = 3.00e-5 * cm
    h_ndoping       = 8.50e-6 * cm + 1.0e-7 *cm # add 1.e-7 cm to this layer for agreement with grid of Driftfusion

    x0              = 0.0 * cm
    δ               = 4*n        # the larger, the finer the mesh
    t               = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k               = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_p_u       = collect(range(x0, h_pdoping/2, step=h_pdoping/(0.6*δ)))
    coord_p_g       = geomspace(h_pdoping/2,
                                h_pdoping,
                                h_pdoping/(0.8*δ),
                                h_pdoping/(0.6*δ),
                                tol=t)
    coord_i_g1      = geomspace(h_pdoping,
                                h_pdoping+h_intrinsic/k,
                                h_intrinsic/(6.1*δ),
                                h_intrinsic/(2.1*δ),
                                tol=t)
    coord_i_g2      = geomspace(h_pdoping+h_intrinsic/k,
                                h_pdoping+h_intrinsic,
                                h_intrinsic/(2.1*δ),
                                h_intrinsic/(6.1*δ),
                                tol=t)
    coord_n_g       = geomspace(h_pdoping+h_intrinsic,
                                h_pdoping+h_intrinsic+h_ndoping/2,
                                h_ndoping/(1.5*δ),
                                h_ndoping/(0.9*δ),
                                tol=t)
    coord_n_u       = collect(range(h_pdoping+h_intrinsic+h_ndoping/2, h_pdoping+h_intrinsic+h_ndoping, step=h_pdoping/(0.5*δ)))

    coord           = glue(coord_p_u, coord_p_g,  tol=10*t)
    coord           = glue(coord,     coord_i_g1, tol=10*t)
    coord           = glue(coord,     coord_i_g2, tol=10*t)
    coord           = glue(coord,     coord_n_g,  tol=10*t)
    coord          = glue(coord,     coord_n_u,  tol=10*t)
    grid           = simplexgrid(coord)

    ## set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm],                [h_pdoping],                           regionAcceptor, tol = 1.0e-18)  # n-doped region   = 1
    cellmask!(grid, [h_pdoping],               [h_pdoping + h_intrinsic],             regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor, tol = 1.0e-18)     # p-doped region   = 3

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
    iphin              = 1 # electron quasi Fermi potential
    iphip              = 2 # hole quasi Fermi potential
    iphia              = 3 # anion vacancy quasi Fermi potential
    numberOfCarriers   = 3 # electrons, holes and anion vacancies

    ## temperature
    T                  =  300.0                 *  K

    ## band edge energies
    Ec_a               = -3.0                  *  eV
    Ev_a               = -5.1                  *  eV

    Ec_i               = -3.8                  *  eV
    Ev_i               = -5.4                  *  eV

    Ec_d               = -3.8                  *  eV
    Ev_d               = -6.2                  *  eV

    EC                 = [Ec_a, Ec_i, Ec_d]
    EV                 = [Ev_a, Ev_i, Ev_d]

    ## effective densities of state
    Nc_a               = 1.0e20                / (cm^3)
    Nv_a               = 1.0e20                / (cm^3)

    Nc_i               = 1.0e19                / (cm^3)
    Nv_i               = 1.0e19                / (cm^3)

    ## ############ adjust Na, Ea for anion vacancies here ###########
    Nanion             = 1.21e22                / (cm^3)
    Ea_i               = -5.175                 *  eV
    ## for the labels in the figures
    textEa             = Ea_i./eV
    textNa             = Nanion.*cm^3
    ## ############ adjust Na, Ea for anion vacancies here ###########

    EA                 = [0.0,  Ea_i,  0.0]

    Nc_d               = 1.0e19                 / (cm^3)
    Nv_d               = 1.0e19                 / (cm^3)

    NC                 = [Nc_a, Nc_i, Nc_d]
    NV                 = [Nv_a, Nv_i, Nv_d]
    NAnion             = [0.0,  Nanion, 0.0]

    ## mobilities
    μn_a               = 0.1                    * (cm^2) / (V * s)
    μp_a               = 0.1                    * (cm^2) / (V * s)

    μn_i               = 2.00e1                 * (cm^2) / (V * s)
    μp_i               = 2.00e1                 * (cm^2) / (V * s)
    μa_i               = 1.00e-10               * (cm^2) / (V * s)

    μn_d               = 1.0e-3                 * (cm^2) / (V * s)
    μp_d               = 1.0e-3                 * (cm^2) / (V * s)

    μn                 = [μn_a, μn_i, μn_d]
    μp                 = [μp_a, μp_i, μp_d]
    μa                 = [0.0,  μa_i, 0.0 ]

    ## relative dielectric permittivity
    ε_a                = 4.0                    *  1.0
    ε_i                = 23.0                   *  1.0
    ε_d                = 3.0                    *  1.0

    ε                  = [ε_a, ε_i, ε_d]

    ## radiative recombination
    r0_a               = 6.3e-11                * cm^3 / s
    r0_i               = 3.6e-12                * cm^3 / s
    r0_d               = 6.8e-11                * cm^3 / s

    r0                 = [r0_a, r0_i, r0_d]

    ## life times and trap densities
    τn_a               = 1.0e-6                 * s
    τp_a               = 1.0e-6                 * s

    τn_i               = 1.0e-7                 * s
    τp_i               = 1.0e-7                 * s
    τn_d               = τn_a
    τp_d               = τp_a

    τn                 = [τn_a, τn_i, τn_d]
    τp                 = [τp_a, τp_i, τp_d]

    ## SRH trap energies
    Ei_a               = -4.05                 * eV
    Ei_i               = -4.60                 * eV
    Ei_d               = -5.00                 * eV

    EI                 = [Ei_a, Ei_i, Ei_d]

    ## Auger recombination
    Auger              = 0.0

    ## generation
    generation_a       = 0.0
    generation_i       = 2.64e21               / (cm^3 * s)
    generation_d       = 0.0

    generation_uniform = [generation_a, generation_i, generation_d]

    ## doping
    Nd                 = 2.089649130192123e17  / (cm^3)
    Na                 = 4.529587947185444e18  / (cm^3)
    C0                 = 1.0e18                / (cm^3)

    ## contact voltage
    voltageAcceptor    = 1.1                   * V

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in predefined data
    data                               = Data(grid, numberOfCarriers)

    ## Possible choices: Stationary, Transient
    data.modelType                     = Transient

    ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    ## FermiDiracMinusOne, Blakemore
    data.F                             = [Boltzmann, Boltzmann, FermiDiracMinusOne]

    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = true,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)

    ## Possible choices: GenerationNone, GenerationUniform
    data.generationModel               = GenerationUniform

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceModelNone,
    ## InterfaceModelSurfaceReco (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor]    = OhmicContact

    ## Present ionic vacancies in perovskite layer
    data.enableIonicCarriers           = enable_ionic_carriers(ionic_carriers = [iphia], regions = [regionIntrinsic])

    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation             = ExcessChemicalPotential

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
    params.chargeNumbers[iphia]                         =  1

    ## boundary region data
    params.bDensityOfStates[iphin, bregionDonor]        = Nc_d
    params.bDensityOfStates[iphip, bregionDonor]        = Nv_d

    params.bDensityOfStates[iphin, bregionAcceptor]     = Nc_a
    params.bDensityOfStates[iphip, bregionAcceptor]     = Nv_a

    params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_d
    params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_d

    params.bBandEdgeEnergy[iphin, bregionAcceptor]      = Ec_a
    params.bBandEdgeEnergy[iphip, bregionAcceptor]      = Ev_a


    for ireg in 1:numberOfRegions # interior region data

        params.dielectricConstant[ireg]                 = ε[ireg]

        ## effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = NC[ireg]
        params.densityOfStates[iphip, ireg]             = NV[ireg]
        params.densityOfStates[iphia, ireg]             = NAnion[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = EC[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = EV[ireg]
        params.bandEdgeEnergy[iphia, ireg]              = EA[ireg]

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]
        params.mobility[iphia, ireg]                    = μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, data, EI[ireg])
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

        params.generationUniform[ireg]                  = generation_uniform[ireg]
    end

    ## interior doping
    params.doping[iphin, regionDonor]                   = Nd
    params.doping[iphia, regionIntrinsic]               = C0
    params.doping[iphip, regionAcceptor]                = Na

    ## boundary doping
    params.bDoping[iphip, bregionAcceptor]              = Na
    params.bDoping[iphin, bregionDonor]                 = Nd

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    if test == false
        show_params(ctsys)
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define outer boundary conditions")
    end
    ################################################################################

    ## set zero voltage ohmic contacts for electrons and holes at all outer boundaries.
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

    control                   = VoronoiFVM.NewtonControl()
    control.verbose           = verbose
    control.max_iterations    = 300
    control.tol_absolute      = 1.0e-10
    control.tol_relative      = 1.0e-10
    control.handle_exceptions = true
    control.tol_round         = 1.0e-10
    control.max_round         = 5
    control.damp_initial      = 0.5
    control.damp_growth       = 1.21 # >= 1

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium for Boltzmann")
    end
    ################################################################################

    ## initialize solution and starting vectors
    initialGuess  = unknowns(ctsys)
    solution      = unknowns(ctsys)

    solution      = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

    initialGuess .= solution

    if plotting
        label_solution, label_density, label_energy = set_plotting_labels(data)

        ## add labels for anion vacancy
        label_energy[1, iphia] = "\$E_a-q\\psi\$"; label_energy[2, iphia] = "\$ - q \\varphi_a\$"
        label_density[iphia]   = "a";              label_solution[iphia]  = "\$ \\varphi_a\$"

        plot_energies(Plotter, grid, data, solution, "Equilibrium; \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$", label_energy)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution,"Equilibrium; \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$", label_density)
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "Equilibrium; \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$", label_solution)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Loop for putting generation on")
    end
    ################################################################################

    data.calculationType = OutOfEquilibrium

    ## primary data for I-V scan protocol
    scanrate             = 0.04 * V/s
    number_tsteps        = 31
    endVoltage           = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend                 = endVoltage/scanrate

    ## with fixed timestep sizes we can calculate the times priori
    tvalues              = range(0, stop = tend, length = number_tsteps)

    ## these values are needed for putting the generation slightly on
    I      = collect(length(tvalues):-1:0.0)
    LAMBDA = 10 .^ (-I)

    for istep = 2:number_tsteps

        t  = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep-1] # Time step size

        ## Apply new voltage by setting non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        ## turn slowly generation on
        data.λ2   = LAMBDA[istep + 1]

        if test == false
            println("increase generation with λ2 = $(data.λ2)")
            println("time value: t = $(t)")
        end

        solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)

        initialGuess .= solution

    end # time loop

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Reverse IV Scan Protocol")
    end
    ################################################################################

    ## for saving I-V data
    IVReverse          = zeros(0) # for IV values
    biasValuesReverse  = zeros(0) # for bias values

    for istep = number_tsteps:-1:2

        t  = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep-1] # Time step size

        ## Apply new voltage
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: t = $(t)")
        end

        solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)

        initialGuess .= solution

        ## get I-V data
        current = get_current_val(ctsys, solution, initialGuess, Δt)

        push!(IVReverse, current)
        push!(biasValuesReverse, Δu)

    end # time loop

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Forward IV Scan Protocol")
    end
    ################################################################################

    ## for saving I-V data
    IVForward          = zeros(0) # for IV values
    biasValuesForward  = zeros(0) # for bias values

    for istep = 2:number_tsteps

        t  = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep-1] # Time step size

        ## Apply new voltage
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: t = $(t)")
        end

        solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)

        ## get I-V data
        current = get_current_val(ctsys, solution, initialGuess, Δt)

        push!(IVForward, current)
        push!(biasValuesForward, Δu)

        initialGuess .= solution


    end # time loop

    if test == false
        println("*** done\n")
    end

    ##resForward = [biasValuesForward IVForward]

    if plotting
        Plotter.figure()

        plot_densities(Plotter, grid, data, solution, "\$ \\Delta u = $(biasValuesForward[end])\$; \$ E_a =\$$(textEa)eV;  \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$", label_density)
        ## ###############
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "\$ \\Delta u = $(biasValuesForward[end])\$; \$ E_a =\$$(textEa)eV;  \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$", label_solution)
        ## ###############
        Plotter.figure()

        Plotter.plot(biasValuesForward, IVForward.*(cm^2), label = "forward",  linewidth= 3, linestyle="--", color="red")
        Plotter.plot(biasValuesReverse, IVReverse.*(cm^2), label = "reverse",  linewidth= 3, linestyle="--", color="blue")
        Plotter.legend()
        Plotter.xlabel("Applied Voltage [V]")
        Plotter.ylabel("current density [A \$ cm^{-2}\$ ]")
        Plotter.title("\$ E_a =\$$(textEa)eV;  \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$ (without internal BC)")
        Plotter.tight_layout()
    end

    testval = VoronoiFVM.norm(ctsys.fvmsys, solution, 2)
    return testval

end #  main

function test()
    testval = 39.067907030939914
    main(test = true, unknown_storage=:dense) ≈ testval #&& main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
