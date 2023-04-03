# PSC device without mobile ions (1D).
([source code](https://github.com/PatricioFarrell/ChargeTransport.jl/tree/master/examplesEx103_PSC.jl))

Simulating a three layer PSC device SiO2| MAPI | SiO2 without mobile ions and in stationary
state. We consider heterojunctions. The simulations are performed out of equilibrium and with
abrupt interfaces. For simplicity, the generation is off.

This simulation coincides with the one made in Section 4.3
of Calado et al. (https://arxiv.org/abs/2009.04384) with the parameters in Table S.13. Or here:
https://github.com/barnesgroupICL/Driftfusion/blob/Methods-IonMonger-Comparison/Input_files/IonMonger_default_bulk.csv

````julia
module Ex103_PSC

using ChargeTransport
using ExtendableGrids
using PyPlot

function main(;n = 3, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # region numbers
    regionDonor      = 1                           # n doped region
    regionIntrinsic  = 2                           # intrinsic region
    regionAcceptor   = 3                           # p doped region
    regions          = [regionDonor, regionIntrinsic, regionAcceptor]
    numberOfRegions  = length(regions)

    # boundary region numbers
    bregionDonor     = 1
    bregionAcceptor  = 2
    bregionJunction1 = 3
    bregionJunction2 = 4

    h_ndoping        = 9.90e-6 * cm
    h_intrinsic      = 4.00e-5 * cm + 2.0e-7 * cm # add 2.e-7 cm to this layer for agreement with grid of Driftfusion
    h_pdoping        = 1.99e-5 * cm
    h_total          = h_ndoping + h_intrinsic + h_pdoping

    x0               = 0.0 * cm
    δ                = 2*n        # the larger, the finer the mesh
    t                = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k                = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_n_u        = collect(range(x0, h_ndoping/2, step=h_ndoping/(0.8*δ)))
    coord_n_g        = geomspace(h_ndoping/2,
                                 h_ndoping,
                                 h_ndoping/(1.0*δ),
                                 h_ndoping/(1.0*δ),
                                 tol=t)
    coord_i_g1       = geomspace(h_ndoping,
                                 h_ndoping+h_intrinsic/k,
                                 h_intrinsic/(2.8*δ),
                                 h_intrinsic/(2.0*δ),
                                 tol=t)
    coord_i_g2       = geomspace(h_ndoping+h_intrinsic/k,
                                 h_ndoping+h_intrinsic,
                                 h_intrinsic/(2.0*δ),
                                 h_intrinsic/(2.8*δ),
                                 tol=t)
    coord_p_g        = geomspace(h_ndoping+h_intrinsic,
                                 h_ndoping+h_intrinsic+h_pdoping/2,
                                 h_pdoping/(1.6*δ),
                                 h_pdoping/(1.6*δ),
                                 tol=t)
    coord_p_u        = collect(range(h_ndoping+h_intrinsic+h_pdoping/2, h_total, step=h_pdoping/(1.3*δ)))

    coord            = glue(coord_n_u, coord_n_g,  tol=10*t)
    coord            = glue(coord,     coord_i_g1, tol=10*t)
    coord            = glue(coord,     coord_i_g2, tol=10*t)
    coord            = glue(coord,     coord_p_g,  tol=10*t)
    coord            = glue(coord,     coord_p_u,  tol=10*t)
    grid             = simplexgrid(coord)


    # set different regions in grid
    cellmask!(grid, [0.0 * μm],                [h_ndoping],               regionDonor)       # n-doped region   = 1
    cellmask!(grid, [h_ndoping],               [h_ndoping + h_intrinsic], regionIntrinsic)   # intrinsic region = 2
    cellmask!(grid, [h_ndoping + h_intrinsic], [h_total],                 regionAcceptor)    # p-doped region   = 3

    # set different boundary regions
    bfacemask!(grid, [0.0],                     [0.0],                     bregionDonor)     # outer left boundary
    bfacemask!(grid, [h_total],                 [h_total],                 bregionAcceptor)  # outer right boundary
    bfacemask!(grid, [h_ndoping],               [h_ndoping],               bregionJunction1) # first  inner interface
    bfacemask!(grid, [h_ndoping + h_intrinsic], [h_ndoping + h_intrinsic], bregionJunction2) # second inner interface

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

    # set indices of the quasi Fermi potentials
    iphin            = 2 # electron quasi Fermi potential
    iphip            = 1 # hole quasi Fermi potential
    numberOfCarriers = 2
````

Define the physical data.

````julia
    # temperature
    T                = 300.0                 *  K

    # band edge energies
    Ec_d             = -4.0                  *  eV
    Ev_d             = -6.0                  *  eV

    Ec_i             = -3.7                  *  eV
    Ev_i             = -5.4                  *  eV

    Ec_a             = -3.1                  *  eV
    Ev_a             = -5.1                  *  eV

    EC               = [Ec_d, Ec_i, Ec_a]
    EV               = [Ev_d, Ev_i, Ev_a]

    # effective densities of state
    Nc_d             = 5.0e19                / (cm^3)
    Nv_d             = 5.0e19                / (cm^3)

    Nc_i             = 8.1e18                / (cm^3)
    Nv_i             = 5.8e18                / (cm^3)

    Nc_a             = 5.0e19                / (cm^3)
    Nv_a             = 5.0e19                / (cm^3)

    NC               = [Nc_d, Nc_i, Nc_a]
    NV               = [Nv_d, Nv_i, Nv_a]

    # mobilities
    μn_d             = 3.89                  * (cm^2) / (V * s)
    μp_d             = 3.89                  * (cm^2) / (V * s)

    μn_i             = 6.62e1                * (cm^2) / (V * s)
    μp_i             = 6.62e1                * (cm^2) / (V * s)

    μn_a             = 3.89e-1               * (cm^2) / (V * s)
    μp_a             = 3.89e-1               * (cm^2) / (V * s)

    μn               = [μn_d, μn_i, μn_a]
    μp               = [μp_d, μp_i, μp_a]

    # relative dielectric permittivity

    ε_d              = 10.0                  *  1.0
    ε_i              = 24.1                  *  1.0
    ε_a              = 3.0                   *  1.0

    ε                = [ε_d, ε_i, ε_a]

    # radiative recombination
    r0_d             = 0.0e+0               * cm^3 / s
    r0_i             = 1.0e-12              * cm^3 / s
    r0_a             = 0.0e+0               * cm^3 / s

    r0               = [r0_d, r0_i, r0_a]

    # life times and trap densities
    τn_d             = 1.0e100              * s
    τp_d             = 1.0e100              * s

    τn_i             = 3.0e-10              * s
    τp_i             = 3.0e-8               * s
    τn_a             = τn_d
    τp_a             = τp_d

    τn               = [τn_d, τn_i, τn_a]
    τp               = [τp_d, τp_i, τp_a]

    # SRH trap energies (needed for calculation of trap_density! (SRH))
    Ei_d             = -5.0                 * eV
    Ei_i             = -4.55                * eV
    Ei_a             = -4.1                 * eV

    EI               = [Ei_d, Ei_i, Ei_a]

    # Auger recombination
    Auger            = 0.0

    # doping
    Nd               =   1.03e18             / (cm^3)
    Na               =   1.03e18             / (cm^3)
    Ni_acceptor      =   8.32e7              / (cm^3)

    # contact voltage: we impose an applied voltage only on one boundary.
    # At the other boundary the applied voltage is zero.
    voltageAcceptor  =  1.2                  * V

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################
````

We initialize the Data instance and fill in predefined data.

````julia
    data                               = Data(grid, numberOfCarriers)

    # Possible choices: Stationary, Transient
    data.modelType                     = Stationary

    # Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk,
    # FermiDiracOneHalfTeSCA, FermiDiracMinusOne, Blakemore
    data.F                            .= FermiDiracOneHalfTeSCA

    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = true,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = false)

    # Possible choices: OhmicContact, SchottkyContact(outer boundary) and InterfaceNone,
    # InterfaceRecombination (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor]    = OhmicContact

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
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1

    for ireg in 1:numberOfRegions ## region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

        # effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = NC[ireg]
        params.densityOfStates[iphip, ireg]             = NV[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = EC[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = EV[ireg]

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]

        # recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, data, EI[ireg])
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger
    end

    #  doping
    params.doping[iphin, regionDonor]                   = Nd
    params.doping[iphip, regionIntrinsic]               = Ni_acceptor
    params.doping[iphip, regionAcceptor]                = Na

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    # print data
    if test == false
        show_params(ctsys)
    end

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control              = SolverControl()
    control.verbose      = verbose
    control.abstol       = 1.0e-13
    control.reltol       = 1.0e-13
    control.tol_round    = 1.0e-13
    control.damp_initial = 0.6
    control.damp_growth  = 1.21 # >= 1
    control.max_round    = 7

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    solution = equilibrium_solve!(ctsys, control = control)
    inival   = solution

    if plotting
        label_solution, label_density, label_energy = set_plotting_labels(data)

        plot_energies(Plotter,  ctsys, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter,  ctsys, solution, "Equilibrium", label_solution)
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

    maxBias    = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 13)

    for Δu in biasValues

        if test == false
            println("bias value: Δu = $(Δu)", " V")
        end

        # set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        solution  = solve(ctsys; inival = inival, control = control)
        inival   .= solution

    end # bias loop

    # plotting
    if plotting
        plot_energies(Plotter,  ctsys, solution, "Applied voltage Δu = $maxBias", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Applied voltage Δu = $maxBias", label_density)
        Plotter.figure()
        plot_solution(Plotter,  ctsys, solution, "Applied voltage Δu = $maxBias", label_solution)
    end

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solution))/length(solution) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = -0.9114885947340186
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

