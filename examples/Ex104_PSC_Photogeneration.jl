
#=
# PSC device with photogeneration rate (1D).
([source code](SOURCE_URL))

Simulating a three layer PSC device TiO2 | MAPI | Pedot with mobile ions where the
ion vacancy accumulation is limited by the Fermi-Dirac integral of order -1.

We perform a linear scan protocol and try out different photogeneration rates.
=#

module Ex104_PSC_Photogeneration

using ChargeTransport
using ExtendableGrids
using PyPlot

function main(;n = 5, Plotter = PyPlot, plotting = false, verbose = false, test = false,
              ########################
              parameter_file = "../parameter_files/Params_PSC_TiO2_MAPI_spiro.jl", # choose the parameter file
              ########################
              userdefinedGeneration = false) # you can choose between predefined and user-defined generation profiles

    Plotter.close("all")

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    include(parameter_file) # include the parameter file we specified

    ## contact voltage
    voltageAcceptor = 1.2 * V

    ## primary data for I-V scan protocol
    scanrate        = 0.04 * V/s
    number_tsteps   = 31
    endVoltage      = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend            = endVoltage/scanrate

    ## Define scan protocol function
    function scanProtocol(t)

        if    0.0 <= t  && t <= tend
            biasVal = 0.0 + scanrate * t
        elseif  t > tend  && t <= 2*tend
            biasVal = scanrate * tend .+ scanrate * (tend - t)
        else
            biasVal = 0.0
        end

        return biasVal

    end

    ## Apply zero voltage on left boundary and a linear scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, scanProtocol]

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    δ                = 4*n        # the larger, the finer the mesh
    t                = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k                = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_n_u        = collect(range(0.0, h_ndoping/2, step=h_ndoping/(0.8*δ)))
    coord_n_g        = geomspace(h_ndoping/2,
                                 h_ndoping,
                                 h_ndoping/(0.7*δ),
                                 h_ndoping/(1.1*δ),
                                 tol=t)
    coord_i_g1       = geomspace(h_ndoping,
                                 h_ndoping+h_intrinsic/k,
                                 h_intrinsic/(2.8*δ),
                                 h_intrinsic/(2.1*δ),
                                 tol=t)
    coord_i_g2       = geomspace(h_ndoping+h_intrinsic/k,
                                 h_ndoping+h_intrinsic,
                                 h_intrinsic/(2.1*δ),
                                 h_intrinsic/(2.8*δ),
                                 tol=t)
    coord_p_g        = geomspace(h_ndoping+h_intrinsic,
                                 h_ndoping+h_intrinsic+h_pdoping/2,
                                 h_pdoping/(1.6*δ),
                                 h_pdoping/(1.6*δ),
                                 tol=t)
    coord_p_u        = collect(range(h_ndoping+h_intrinsic+h_pdoping/2, h_ndoping+h_intrinsic+h_pdoping, step=h_pdoping/(1.3*δ)))

    coord            = glue(coord_n_u, coord_n_g,  tol=10*t)
    coord            = glue(coord,     coord_i_g1, tol=10*t)
    coord            = glue(coord,     coord_i_g2, tol=10*t)
    coord            = glue(coord,     coord_p_g,  tol=10*t)
    coord            = glue(coord,     coord_p_u,  tol=10*t)
    grid             = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm],        [heightLayers[1]], regionDonor,     tol = 1.0e-18) # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionAcceptor,  tol = 1.0e-18) # p-doped region   = 3

    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1,      tol = 1.0e-18)
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJ2,      tol = 1.0e-18)

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

    ## Initialize Data instance and fill in predefined data
    if userdefinedGeneration

        subg1          = subgrid(grid, [regionDonor]); subg2 = subgrid(grid, [regionIntrinsic]); subg3 = subgrid(grid, [regionAcceptor])

        gen1           = zeros(length(subg1[Coordinates])-1); gen3 = zeros(length(subg3[Coordinates])-1)
        gen2           = incidentPhotonFlux[regionIntrinsic] .* absorption[regionIntrinsic] .* exp.( - absorption[regionIntrinsic] .* (subg2[Coordinates] .- generationPeak))

        weight1        = (subg2[Coordinates][1] - subg1[Coordinates][end-1]) / (subg2[Coordinates][2]-subg1[Coordinates][end-1])
        weight2        = (subg2[Coordinates][end] - subg2[Coordinates][end-1]) / (subg3[Coordinates][2]-subg2[Coordinates][end-1])

        gen2[1]        = weight1 * gen2[1]; gen2[end] = weight2 * gen2[end]

        generationData = [gen1; gen2'; gen3]

        data                          = Data(grid, numberOfCarriers,
                                             contactVoltageFunction = contactVoltageFunction,
                                             generationData = generationData)
    else

        data                          = Data(grid, numberOfCarriers,
                                             contactVoltageFunction = contactVoltageFunction)

    end

    data.modelType                     = Transient
    data.F                             = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = false,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor]    = OhmicContact
    data.fluxApproximation            .= ExcessChemicalPotential

    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

    if userdefinedGeneration
        data.generationModel           = GenerationUserDefined
    else
        data.generationModel           = GenerationBeerLambert
    end

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

    for ireg in 1:numberOfRegions # interior region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

        ## effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nn[ireg]
        params.densityOfStates[iphip, ireg]             = Np[ireg]
        params.densityOfStates[iphia, ireg]             = Na[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = En[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = Ep[ireg]
        params.bandEdgeEnergy[iphia, ireg]              = Ea[ireg]

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]
        params.mobility[iphia, ireg]                    = μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, params, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, params, EI[ireg])

        ## generation parameters
        params.generationIncidentPhotonFlux[ireg]       = incidentPhotonFlux[ireg]
        params.generationAbsorption[ireg]               = absorption[ireg]
        params.generationUniform[ireg]                  = generation_uniform[ireg]
    end

    # parameter which passes the shift information in the Beer-Lambert generation
    params.generationPeak                               = generationPeak

    ## interior doping
    params.doping[iphin, regionDonor]                   = Cn
    params.doping[iphia, regionIntrinsic]               = Ca
    params.doping[iphip, regionAcceptor]                = Cp

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=:sparse)

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control              = SolverControl()
    if verbose == true
        control.verbose  = verbose
    else
        control.verbose  = "eda" # still print the time values
    end
    if test == true
        control.verbose  = false # do not show time values in testing case
    end
    control.maxiters     = 300
    control.max_round    = 5
    control.damp_initial = 0.5
    control.damp_growth  = 1.21 # >= 1
    control.Δt_max       = 5.0e-1

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    ## calculate equilibrium solution and as initial guess
    solution = equilibrium_solve!(ctsys, control = control)
    inival   = solution

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Loop for generation")
    end
    ################################################################################

    # these values are needed for putting the generation slightly on
    I      = collect(20:-1:0.0)
    LAMBDA = 10 .^ (-I)

    ## since the constant which represents the constant quasi Fermi potential of anion vacancies is undetermined, we need
    ## to fix it in the bias loop, since we have no applied bias. Otherwise we get convergence errors
    ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 1.0e30
    ctsys.fvmsys.boundary_values[iphia, bregionJ2]  = 0.0

    for istep = 1:length(I)-1

        ## turn slowly generation on
        ctsys.data.λ2   = LAMBDA[istep + 1]

        if test == false
            println("increase generation with λ2 = $(data.λ2)")
        end

        solution = solve(ctsys, inival = inival, control = control)
        inival   = solution

    end # generation loop

    solutionEQ = inival

    if plotting
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        ## add labels for anion vacancy
        label_energy[1, iphia] = "\$E_a-q\\psi\$"; label_energy[2, iphia] = "\$ - q \\varphi_a\$"; label_BEE[iphia] = "\$E_a\$"
        label_density[iphia]   = "\$ n_a \$";      label_solution[iphia]  = "\$ \\varphi_a\$"

        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Initial condition", label_density)
        Plotter.figure()
        plot_solution(Plotter, ctsys, solution, "Initial condition", label_solution)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    ## put here back the homogenous Neumann boundary conditions.
    ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 0.0
    ctsys.fvmsys.boundary_values[iphia, bregionJ2]  = 0.0

    sol = solve(ctsys, inival = inival, times=(0.0, tend), control = control)

    if plotting
        tsol = sol(tend)
        Plotter.figure()
        plot_densities(Plotter, ctsys, tsol, "Densities at end time", label_density)
        Plotter.tight_layout()
        Plotter.figure()
        plot_solution(Plotter, ctsys, tsol, "Solution at end time", label_solution)
        Plotter.tight_layout()
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Reverse scan protocol")
    end
    ################################################################################

    inivalReverse = sol(tend)
    solReverse    = solve(ctsys, inival = inivalReverse, times=(tend, 2 * tend), control = control)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("IV Curve calculation")
    end
    ################################################################################

    factory       = TestFunctionFactory(ctsys)
    tf            = testfunction(factory, [bregionDonor], [bregionAcceptor])

    tvalues       = sol.t
    number_tsteps = length(tvalues)
    biasValues    = scanProtocol.(tvalues)
    IV            = zeros(0)

    for istep = 2:number_tsteps
        Δt       = tvalues[istep] - tvalues[istep-1] # Time step size
        inival   = sol[istep-1]
        solution = sol[istep]

        I        = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii = 1:numberOfCarriers+1
            current = current + I[ii]
        end

        push!(IV, current)

    end

    tvaluesReverse       = solReverse.t
    number_tstepsReverse = length(tvaluesReverse)
    biasValuesReverse    = scanProtocol.(tvaluesReverse)
    IVReverse            = zeros(0)

    for istep = 2:number_tstepsReverse
        Δt       = tvaluesReverse[istep] - tvaluesReverse[istep-1] # Time step size
        inival   = solReverse[istep-1]
        solution = solReverse[istep]

        I        = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii = 1:numberOfCarriers+1
            current = current + I[ii]
        end

        push!(IVReverse, current)

    end

    if plotting
        Plotter.figure()
        Plotter.plot([tvalues tvaluesReverse], [biasValues biasValuesReverse], marker = "x")
        Plotter.xlabel("time [s]")
        Plotter.ylabel("voltage [V]")
        Plotter.grid()

        Plotter.figure()
        Plotter.plot(biasValues[2:end], -IV, linewidth = 5, label = "forward")
        Plotter.plot(biasValuesReverse[2:end], -IVReverse, linewidth = 5, label = "reverse")
        Plotter.grid()
        Plotter.legend()
        Plotter.xlabel("applied bias [V]")
        Plotter.ylabel("total current [A]")

        Plotter.figure()
        if userdefinedGeneration
            Plotter.plot(coord, data.generationData)
        else
            for ireg = 1:numberOfRegions
                subg = subgrid(grid, [ireg])
                Plotter.plot(subg[Coordinates]', BeerLambert(ctsys, ireg, subg[Coordinates])', label = "region $ireg")
            end

        end
        Plotter.legend()
        Plotter.grid()
        Plotter.xlabel("space [\$m\$]")
        Plotter.ylabel("photogeneration [\$\\frac{1}{cm^3s}\$]")
        Plotter.tight_layout()
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute fill factor and efficiency")
    end
    ################################################################################

    bias                      = biasValues[2:end]
    IV                        = -IV

    powerDensity              = bias .* (IV)           # power density function
    MaxPD, indexPD            = findmax(powerDensity)

    open_circuit              = compute_open_circuit_voltage(bias, IV)

    IncidentLightPowerDensity = 1000.0 * W/m^2

    efficiency                =  bias[indexPD] * IV[indexPD]  / IncidentLightPowerDensity
    fillfactor                = (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

    if test == false
        println("The fill factor is $fillfactor %.")
        println("The efficiency  is $efficiency %.")
        println("The open circuit voltage  is $open_circuit.")
    end

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solutionEQ))/length(solutionEQ) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = -1.052813874410313; testvalUserdefined = -1.0528971495353738
    main(test = true, userdefinedGeneration = false) ≈ testval && main(test = true, userdefinedGeneration = true) ≈ testvalUserdefined
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
