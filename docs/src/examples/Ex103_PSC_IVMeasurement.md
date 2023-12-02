# PSC device with ions and different I-V scan protocols (1D).
([source code](https://github.com/PatricioFarrell/ChargeTransport.jl/tree/master/examples/Ex103_PSC_IVMeasurement.jl))

Simulating a three layer PSC device Ti02| MAPI | spiro-OMeTAD with mobile ions where
the ion vacancy accumulation is limited by the Fermi-Dirac integral of order -1.

The time-dependent simulations are performed with abrupt interfaces.
Two different I-V measurement protocols are included and the corresponding solution vectors
and the I-V curve after the scan can be depicted.

````julia
module Ex103_PSC_IVMeasurement

using ChargeTransport
using ExtendableGrids
using PyPlot

function main(;n = 3, Plotter = PyPlot, plotting = false, verbose = false, test = false,
    parameter_file = "../parameter_files/Params_PSC_TiO2_MAPI_spiro.jl", # choose the parameter file
    otherScanProtocol = false) # you can choose between two scan protocols

    PyPlot.close("all")

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    include(parameter_file) # include the parameter file we specified

    # contact voltage
    voltageAcceptor = 1.2 * V

    # primary data for I-V scan protocol
    scanrate        = 1.0 * V/s
    number_tsteps   = 31
    endVoltage      = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend            = endVoltage/scanrate

    # Define scan protocol function
    function linearScanProtocol(t)
        if t == Inf
            0.0
        else
            scanrate * t
        end
    end

    # Apply zero voltage on left boundary and a linear scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, linearScanProtocol]

    # Instead of a linear scan protocol, we can also apply other scan protocols which we
    # define by our own and parse to the model generator via the struct Data
    if otherScanProtocol
        # scan protocol parameter
        number_tsteps = 40
        frequence     = 10.0 * Hz
        amplitude     = 0.2  * V
        tend          = 1/frequence

        # Define sinusoidal applied voltage
        function sinusoidalScanProtocol(t)
            if t == Inf
                0.0
            else
                amplitude * sin(2.0 * pi * frequence * t)
            end
        end

        # Apply zero voltage on left boundary and a linear scan protocol on right boundary
        contactVoltageFunction = [zeroVoltage, sinusoidalScanProtocol]
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    δ          = 4*n        # the larger, the finer the mesh
    t          = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k          = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_n_u  = collect(range(0.0, h_ndoping/2, step=h_ndoping/(0.8*δ)))
    coord_n_g  = geomspace(h_ndoping/2, h_ndoping,
                           h_ndoping/(0.7*δ), h_ndoping/(1.1*δ),
                           tol=t)
    coord_i_g1 = geomspace(h_ndoping, h_ndoping+h_intrinsic/k,
                           h_intrinsic/(2.8*δ), h_intrinsic/(2.1*δ),
                           tol=t)
    coord_i_g2 = geomspace(h_ndoping+h_intrinsic/k, h_ndoping+h_intrinsic,
                            h_intrinsic/(2.1*δ), h_intrinsic/(2.8*δ),
                           tol=t)
    coord_p_g  = geomspace(h_ndoping+h_intrinsic, h_ndoping+h_intrinsic+h_pdoping/2,
                           h_pdoping/(1.6*δ), h_pdoping/(1.6*δ),
                           tol=t)
    coord_p_u  = collect(range(h_ndoping+h_intrinsic+h_pdoping/2, h_ndoping+h_intrinsic+h_pdoping, step=h_pdoping/(1.3*δ)))

    coord      = glue(coord_n_u, coord_n_g,  tol=10*t)
    coord      = glue(coord,     coord_i_g1, tol=10*t)
    coord      = glue(coord,     coord_i_g2, tol=10*t)
    coord      = glue(coord,     coord_p_g,  tol=10*t)
    coord      = glue(coord,     coord_p_u,  tol=10*t)
    grid       = ExtendableGrids.simplexgrid(coord)

    # set different regions in grid
    cellmask!(grid,  [0.0 * μm],        [heightLayers[1]], regionDonor, tol = 1.0e-18)     # n-doped region   = 1
    cellmask!(grid,  [heightLayers[1]], [heightLayers[2]], regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid,  [heightLayers[2]], [heightLayers[3]], regionAcceptor, tol = 1.0e-18)  # p-doped region   = 3

    # bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0],             [0.0],             bregionDonor, tol = 1.0e-18)     # outer left boundary
    bfacemask!(grid, [h_total],         [h_total],         bregionAcceptor, tol = 1.0e-18)  # outer right boundary
    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1, tol = 1.0e-18) # first  inner interface
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJ2, tol = 1.0e-18) # second inner interface

    if plotting
        gridplot(grid, Plotter = Plotter, legend=:lt)
        Plotter.title("Grid")
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    # Initialize Data instance and fill in predefined data
    # Currently, the way to go is to pass a contact voltage function exactly here.
    data                               = Data(grid, numberOfCarriers)

    # Possible choices: Stationary, Transient
    data.modelType                     = Transient

    # Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    # FermiDiracMinusOne, Blakemore
    data.F                             = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = false,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)

    # Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    # InterfaceRecombination (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor]    = OhmicContact

    # With this method, the user enable the ionic carrier parsed to ionicCarrier and gives
    # gives the information on which regions this ionic carrier is defined.
    # In this application ion vacancies only live in active perovskite layer.
    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

    # Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    # ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation            .= ExcessChemicalPotential

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
    params.chargeNumbers[iphin]                         = zn
    params.chargeNumbers[iphip]                         = zp
    params.chargeNumbers[iphia]                         = za

    for ireg in 1:numberOfRegions # region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

        # effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nn[ireg]
        params.densityOfStates[iphip, ireg]             = Np[ireg]
        params.densityOfStates[iphia, ireg]             = Na[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = En[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = Ep[ireg]
        params.bandEdgeEnergy[iphia, ireg]              = Ea[ireg]

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]
        params.mobility[iphia, ireg]                    = μa[ireg]

        # recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, params, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, params, EI[ireg])
    end

    # doping
    params.doping[iphin, regionDonor]                   = Cn
    params.doping[iphia, regionIntrinsic]               = Ca
    params.doping[iphip, regionAcceptor]                = Cp

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=:sparse)

    if test == false
        show_params(ctsys)
        println("*** done\n")
    end

    if plotting == true
        ################################################################################
        println("Plot electroneutral potential, band-edge energies and doping")
        ################################################################################
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        # add labels for anion vacancy
        label_energy[1, iphia] = "\$E_a-q\\psi\$"; label_energy[2, iphia] = "\$ - q \\varphi_a\$"; label_BEE[iphia] = "\$E_a\$"
        label_density[iphia]   = "\$ n_a \$";      label_solution[iphia]  = "\$ \\varphi_a\$"

    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control              = SolverControl()
    control.verbose      = verbose
    control.max_round    = 5
    control.damp_initial = 0.1
    control.damp_growth  = 1.21 # >= 1

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    # calculate equilibrium solution and as initial guess
    solution = equilibrium_solve!(ctsys, control = control)
    inival   = solution

    if plotting
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "Equilibrium", label_solution)
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    # with fixed timestep sizes we can calculate the times a priori
    tvalues    = range(0, stop = tend, length = number_tsteps)

    # for saving I-V data
    IV         = zeros(0) # for IV values

    for istep = 2:number_tsteps

        t  = tvalues[istep]                                    # Actual time
        Δu = contactVoltageFunction[bregionAcceptor](t) # Applied voltage
        Δt = t - tvalues[istep-1]                              # Time step size

        # Apply new voltage by setting non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: t = $(t) s")
        end

        # Solve time step problems with timestep Δt. inival plays the role of the solution
        # from last timestep
        solution = solve(ctsys; inival = inival, control = control, tstep = Δt)
        # get I-V data
        current  = get_current_val(ctsys, solution, inival, Δt)

        push!(IV, current)

        inival = solution

    end # time loop

    if test == false
        println("*** done\n")
    end

    # here in res the biasValues and the corresponding current are stored.
    # res = [biasValues IV];

    if plotting
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage)", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution,"bias \$\\Delta u\$ = $(endVoltage)", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(endVoltage)", label_solution)
    end

    biasValues = contactVoltageFunction[bregionAcceptor].(tvalues)

    if plotting
        Plotter.figure()
        Plotter.plot(tvalues, biasValues, marker ="o")
        Plotter.xlabel("time [s]")
        Plotter.ylabel("bias [V]")
        Plotter.figure()
        plot_IV(Plotter, biasValues[2:end], IV, "bias \$\\Delta u\$ = $(endVoltage)")
    end

    testval = sum(filter(!isnan, solution))/length(solution) # when using sparse storage, we get NaN values in solution
    return testval


end #  main

function test()
    testval = -0.6302819608784171; testvalOther = -1.123710261723505
    main(test = true, otherScanProtocol = false) ≈ testval && main(test = true, otherScanProtocol = true) ≈ testvalOther
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

