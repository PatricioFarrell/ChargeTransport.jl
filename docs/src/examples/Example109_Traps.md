# GaAs diode: transient with traps (1D).
([source code](https://github.com/PatricioFarrell/ChargeTransport.jl/tree/master/examplesExample109_Traps.jl))

Simulating transient charge transport in a GaAs pin diode with an electron trap.

````julia
module Example109_Traps

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot

# function to initialize the grid for a possble extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_intrinsic, h_pdoping)
    coord_ndoping    = collect(range(0.0, stop = h_ndoping, length = 3 * refinementfactor))
    coord_intrinsic  = collect(range(h_ndoping, stop = (h_ndoping + h_intrinsic), length = 3 * refinementfactor))
    coord_pdoping    = collect(range((h_ndoping + h_intrinsic),
                                        stop = (h_ndoping + h_intrinsic + h_pdoping),
                                        length = 3 * refinementfactor))
    coord            = glue(coord_ndoping, coord_intrinsic)
    coord            = glue(coord, coord_pdoping)

    return coord
end

function main(;n = 3, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

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
    bregions                = [bregionAcceptor, bregionDonor]
    numberOfBoundaryRegions = length(bregions)

    # grid
    refinementfactor        = 2^(n-1)
    h_pdoping               = 2.0    * μm
    h_intrinsic             = 2.0    * μm
    h_ndoping               = 2.0    * μm
    w_device                = 0.5    * μm  # width of device
    z_device                = 1.0e-4 * cm  # depth of device

    coord                   = initialize_pin_grid(refinementfactor,
                                                  h_pdoping,
                                                  h_intrinsic,
                                                  h_ndoping)

    grid                    = simplexgrid(coord)

    # set different regions in grid, doping profiles do not intersect
    cellmask!(grid, [0.0 * μm], [h_pdoping], regionAcceptor)                                        # p-doped
    cellmask!(grid, [h_pdoping], [h_pdoping + h_intrinsic], regionIntrinsic)                        # intrinsic
    cellmask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic + h_ndoping], regionDonor)  # n-doped

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

    iphin              = 1 # index electron quasi Fermi potential
    iphip              = 2 # index hole quasi Fermi potential
    iphit              = 3 # index trap quasi Fermi potential
    numberOfCarriers   = 3 # electrons, holes and traps

    # physical data
    Ec                  = 1.424                             *  eV
    Ev                  = 0.0                               *  eV
    Et                  = 0.6                               *  eV
    Nc                  = 4.351959895879690e17              / (cm^3)
    Nv                  = 9.139615903601645e18              / (cm^3)
    Nt                  = 1e16                              / (cm^3)
    mun                 = 8500.0                            * (cm^2) / (V * s)
    mup                 = 400.0                             * (cm^2) / (V * s)
    mut                 = 0.0                               * (cm^2) / (V * s)  # such that there is no flux
    εr                  = 12.9                              *  1.0              # relative dielectric permittivity of GAs
    T                   = 300.0                             *  K

    # recombination parameters
    ni                  = sqrt(Nc * Nv) * exp(-(Ec - Ev) / (2 * kB * T))        # intrinsic concentration
    n0                  = Nc * Boltzmann( (Et-Ec) / (kB*T) )                    # Boltzmann equilibrium concentration
    p0                  = ni^2 / n0                                             # Boltzmann equilibrium concentration
    Auger               = 1.0e-29                           * cm^6 / s          # 1.0e-41
    SRH_LifeTime        = 1.0e-3                            * ns
    Radiative           = 1.0e-10                           * cm^3 / s          # 1.0e-16
    G                   = 1.0e25                            / (cm^3 * s)

    # doping -- trap doping will not be set and thus automatically zero
    dopingFactorNd      = 1.0
    dopingFactorNa      = 0.46
    Nd                  = dopingFactorNd * Nc
    Na                  = dopingFactorNa * Nv

    # contact voltages
    voltageAcceptor     = 1.5                               * V

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    # initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)

    # possible choices: model_stationary, model_transient
    data.model_type                     = model_transient

    # possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                             .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    data.bulk_recombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                  bulk_recomb_Auger = true,
                                                                  bulk_recomb_radiative = true,
                                                                  bulk_recomb_SRH = true)

    # Here we enable the traps and parse the respective index and the regions where the trap is defined.
    enable_traps!(data = data, traps = iphit, regions = regions)

    # possible choices: generation_none, generation_uniform
    data.generation_model               = generation_uniform

    # possible choices: ohmic_contact, schottky_contact (outer boundary) and interface_model_none,
    # interface_model_surface_recombination (inner boundary).
    data.boundary_type[bregionAcceptor] = ohmic_contact
    data.boundary_type[bregionDonor]    = ohmic_contact

    # possible choices: scharfetter_gummel, scharfetter_gummel_graded, excess_chemical_potential,
    # excess_chemical_potential_graded, diffusion_enhanced, generalized_sg
    data.flux_approximation             = excess_chemical_potential

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
    params.chargeNumbers[iphit]                         = -1 # trap charge number determines whether hole or electron trap is used

    for ibreg in 1:numberOfBoundaryRegions   # boundary region data

        params.bDensityOfStates[iphin, ibreg]           = Nc
        params.bDensityOfStates[iphip, ibreg]           = Nv
        params.bDensityOfStates[iphit, ibreg]           = Nt
        params.bBandEdgeEnergy[iphin, ibreg]            = Ec
        params.bBandEdgeEnergy[iphip, ibreg]            = Ev
        params.bBandEdgeEnergy[iphit, ibreg]            = Et
    end

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg]                 = εr

        # effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nc
        params.densityOfStates[iphip, ireg]             = Nv
        params.densityOfStates[iphit, ireg]             = Nt
        params.bandEdgeEnergy[iphin, ireg]              = Ec
        params.bandEdgeEnergy[iphip, ireg]              = Ev
        params.bandEdgeEnergy[iphit, ireg]              = Et
        params.mobility[iphin, ireg]                    = mun
        params.mobility[iphip, ireg]                    = mup
        params.mobility[iphit, ireg]                    = mut

        # recombination parameters
        params.recombinationRadiative[ireg]             = Radiative
        params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = n0
        params.recombinationSRHTrapDensity[iphip, ireg] = p0
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger
        params.generationUniform[ireg]                  = G

    end

    # doping -- since we do not set any doping for the traps it is automatically zero
    params.doping[iphin, regionDonor]                   = Nd
    params.doping[iphin, regionIntrinsic]               = ni
    params.doping[iphip, regionIntrinsic]               = 0.0
    params.doping[iphip, regionAcceptor]                = Na

    # boundary doping
    params.bDoping[iphin, bregionDonor]                 = Nd
    params.bDoping[iphip, bregionAcceptor]              = Na

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    if test == false
        show_params(ctsys)
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define outerior boundary conditions and enabled layers")
    end
    ################################################################################

    # set zero voltage ohmic contacts for each charge carrier at all outerior boundaries.
    set_ohmic_contact!(ctsys, bregionAcceptor, 0.0)
    set_ohmic_contact!(ctsys, bregionDonor,    0.0)

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
    control.damp_growth       = 1.21    #>= 1
    control.max_iterations    = 250
    control.tol_absolute      = 1.0e-14
    control.tol_relative      = 1.0e-14
    control.handle_exceptions = true
    control.tol_round         = 1.0e-8
    control.max_round         = 5

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    data.calculation_type = inEquilibrium

    # initialize solution and starting vectors
    initialGuess          = unknowns(ctsys)
    solution              = unknowns(ctsys)

    # solve thermodynamic equilibrium and update initial guess
    solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    initialGuess         .= solution

    if test == false
        println("*** done\n")
    end

    if plotting
        # ##### set legend for plotting routines #####
        label_energy   = Array{String, 2}(undef, 2, numberOfCarriers) # band-edge energies and potential
        label_density  = Array{String, 1}(undef, numberOfCarriers)
        label_solution = Array{String, 1}(undef, numberOfCarriers)

        # for electrons
        label_energy[1, iphin] = "\$E_c-q\\psi\$";       label_energy[2, iphin] = "\$ - q \\varphi_n\$"
        label_density[iphin]   = "n";                    label_solution[iphin]  = "\$ \\varphi_n\$"

        # for holes
        label_energy[1, iphip] = "\$E_v-q\\psi\$";       label_energy[2, iphip] = "\$ - q \\varphi_p\$"
        label_density[iphip]   = "p";                    label_solution[iphip]  = "\$ \\varphi_p\$"

        # for traps
        label_energy[1, iphit] = "\$E_{\\tau}-q\\psi\$"; label_energy[2, iphit] = "\$ - q \\varphi_{\\tau}\$"
        label_density[iphit]   = "\$n_{\\tau}\$";        label_solution[iphit]  = "\$ \\varphi_{\\tau}\$"
        # ##### set legend for plotting routines #####
        plot_energies(Plotter, grid, data, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution,"Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "Equilibrium", label_solution)
    end

    ################################################################################
    if test == false
        println("Adjust Newton parameters")
    end
    ################################################################################

    control.tol_absolute          = 1.0e-10
    control.tol_relative          = 1.0e-10
    control.tol_round             = 1.0e-4
    control.damp_initial          = 0.5
    control.damp_growth           = 1.2
    control.max_iterations        = 30
    control.max_round             = 3

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Embed life times and increase generation")
    end
    ################################################################################

    # set calculation type to outOfEquilibrium for starting with respective simulation.
    ctsys.data.calculation_type   = outOfEquilibrium

    # Scan rate and time steps
    scanrate                      = 1.0 * V/s
    number_tsteps                 = 81
    endVoltage                    = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary

    IV                            = zeros(0) # for IV values
    biasValues                    = zeros(0) # for bias values

    # The end time then is calculated here:
    tend                          = endVoltage/scanrate

    # with fixed timestep sizes we can calculate the times
    # a priori
    tvalues                       = range(0.0, stop = tend, length = number_tsteps)

    steps                         = 35
    I                             = collect(steps:-1:0.0)
    LAMBDA                        = 10 .^ (I)
    Δt                            = tvalues[2] - tvalues[1]

    for i in 1:length(LAMBDA)
        if control.verbose
            println("λ = $(LAMBDA[i]), Nt = $(Nt), life time = $(SRH_LifeTime), n0 = $(n0), p0 = $(p0), G = $(G*1 / (LAMBDA[i] * SRH_LifeTime * Nt))")
        end
        ctsys.data.params.recombinationSRHLifetime[iphin,regions] .= LAMBDA[i] * SRH_LifeTime
        ctsys.data.params.recombinationSRHLifetime[iphip,regions] .= LAMBDA[i] * SRH_LifeTime

        data.λ2 = 1 / (LAMBDA[i] )
        VoronoiFVM.solve!(solution, initialGuess, ctsys, control = control, tstep=Δt)
        initialGuess = solution
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    for istep = 2:number_tsteps

        t                     = tvalues[istep]          # Actual time
        Δu                    = t*scanrate              # Applied voltage
        Δt                    = t - tvalues[istep-1]    # Time step size

        # Apply new voltage: set non equilibrium boundary conditions
        set_ohmic_contact!(ctsys, bregionAcceptor, Δu)

        if verbose
            println("time value: t = $(t)")
        end

        # Solve time step problems with timestep Δt. initialGuess plays the role of the solution
        # from last timestep
        solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)

        # get I-V data
        current = get_current_val(ctsys, solution, initialGuess, Δt)

        push!(IV, w_device * z_device * current)
        push!(biasValues, Δu)

        initialGuess .= solution

        ######## CHECK IF TRAPS CONTRIBUTE TO CURRENT!

    end # bias loop

    if test == false
        println("*** done\n")
    end

    # plot solution and IV curve
    if plotting
        plot_energies(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tvalues[number_tsteps])\$", label_energy)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution,"bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tvalues[number_tsteps])\$", label_density)
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tvalues[number_tsteps])\$", label_solution)
        Plotter.figure()
        plot_IV(Plotter, biasValues,IV, biasValues[end], plotGridpoints = true)
    end

    testval = solution[iphit, 17]
    return testval

    if test == false
        println("*** done\n")
    end

end #  main

function test()
    testval = 1.0245795906936692
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when the PIN module has successfully recompiled.")
end

end # module
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

