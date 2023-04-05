#=
# PSC device with surface recombination (1D).
([source code](SOURCE_URL))

Simulating a three layer PSC device Pedot| MAPI | PCBM. The simulations are performed out of
equilibrium, time-dependent, with abrupt interfaces and with surface recombination at the
internal boundaries.

The parameters are from Calado et al.:
https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/pedotpss_mapi_pcbm.csv.
(with adjustments on layer lengths)
=#

module Ex110_PSC_surface_recombination

using ChargeTransport
using ExtendableGrids
using PyPlot

function main(;n = 6, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionAcceptor   = 1                           # p doped region
    regionIntrinsic  = 2                           # intrinsic region
    regionDonor      = 3                           # n doped region
    regions          = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions  = length(regions)

    ## boundary region numbers
    ## Note that by convention we have 1 for the left boundary and 2 for the right boundary. If
    ## adding additional interior boundaries, continue with 3, 4, ...
    bregionAcceptor  = 1
    bregionDonor     = 2
    bregionJunction1 = 3
    bregionJunction2 = 4

    ## grid (the nearer to interface, the finer)
    h_pdoping        = 3.00e-6 * cm + 1.0e-7 * cm
    h_intrinsic      = 3.00e-5 * cm
    h_ndoping        = 8.50e-6 * cm + 1.0e-7 * cm

    x0               = 0.0 * cm
    δ                = 4*n        # the larger, the finer the mesh
    t                = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k                = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_p_u        = collect(range(x0, h_pdoping/2, step=h_pdoping/(0.5*δ)))
    coord_p_g        = geomspace(h_pdoping/2,
                                 h_pdoping,
                                 h_pdoping/(0.8*δ),
                                 h_pdoping/(1.5*δ),
                                 tol=t)
    coord_i_g1       = geomspace(h_pdoping,
                                 h_pdoping+h_intrinsic/k,
                                 h_intrinsic/(6.1*δ),
                                 h_intrinsic/(2.1*δ),
                                 tol=t)
    coord_i_g2       = geomspace(h_pdoping+h_intrinsic/k,
                                 h_pdoping+h_intrinsic,
                                 h_intrinsic/(2.1*δ),
                                 h_intrinsic/(6.1*δ),
                                 tol=t)
    coord_n_g        = geomspace(h_pdoping+h_intrinsic,
                                 h_pdoping+h_intrinsic+h_ndoping/2,
                                 h_ndoping/(3.0*δ),
                                 h_ndoping/(1.0*δ),
                                 tol=t)
    coord_n_u        = collect(range(h_pdoping+h_intrinsic+h_ndoping/2, h_pdoping+h_intrinsic+h_ndoping, step=h_pdoping/(0.8*δ)))

    coord            = glue(coord_p_u,coord_p_g,  tol=10*t)
    coord            = glue(coord,    coord_i_g1, tol=10*t)
    coord            = glue(coord,    coord_i_g2, tol=10*t)
    coord            = glue(coord,    coord_n_g,  tol=10*t)
    coord            = glue(coord,    coord_n_u,  tol=10*t)
    grid             = ExtendableGrids.simplexgrid(coord)

    ## cellmask! for defining the subregions and assigning region number
    cellmask!(grid, [0.0 * μm],                 [h_pdoping],                           regionAcceptor, tol = 1.0e-18)   # p-doped region   = 1
    cellmask!(grid, [h_pdoping],                [h_pdoping + h_intrinsic],             regionIntrinsic, tol = 1.0e-18)  # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic],  [h_pdoping + h_intrinsic + h_ndoping], regionDonor, tol = 1.0e-18)      # n-doped region   = 3

    ## bfacemask! for ``active'' boundary regions, i.e. internal interfaces. On the outer boundary regions, the
    ## conditions will be formulated later
    bfacemask!(grid, [h_pdoping],               [h_pdoping],                           bregionJunction1)  # first  inner interface
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic],             bregionJunction2)  # second inner interface

    if plotting
        GridVisualize.gridplot(grid, Plotter = Plotter, legend=:lt)
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

    iphin            = 1 # electron quasi Fermi potential
    iphip            = 2 # hole quasi Fermi potential
    iphia            = 3 # anion vacancy quasi Fermi potential
    numberOfCarriers = 3

    ## temperature
    T                =  300.0                *  K

    ## band edge energies
    Ec_a             = -3.0                  *  eV
    Ev_a             = -5.1                  *  eV

    Ec_i             = -3.8                  *  eV
    Ev_i             = -5.4                  *  eV

    Ec_d             = -3.8                  *  eV
    Ev_d             = -6.2                  *  eV

    EC               = [Ec_a, Ec_i, Ec_d]
    EV               = [Ev_a, Ev_i, Ev_d]

    ## effective densities of state
    Nc_a             = 1.0e20                / (cm^3)
    Nv_a             = 1.0e20                / (cm^3)

    Nc_i             = 1.0e19                / (cm^3)
    Nv_i             = 1.0e19                / (cm^3)

    ## ###################### adjust Na, Ea here #####################
    Nanion           = 1.21e22               / (cm^3)
    Ea_i             = -5.175                *  eV

    ## for the labels in the figures
    textEa           = Ea_i                 ./  eV
    textNa           = Nanion               .* (cm^3)
    ## ###################### adjust Na, Ea here #####################

    EA               = [0.0,  Ea_i,  0.0]

    Nc_d             = 1.0e19                / (cm^3)
    Nv_d             = 1.0e19                / (cm^3)

    NC               = [Nc_a, Nc_i, Nc_d]
    NV               = [Nv_a, Nv_i, Nv_d]
    NAnion           = [0.0,  Nanion, 0.0]

    ## mobilities
    μn_a             = 0.1                   * (cm^2) / (V * s)
    μp_a             = 0.1                   * (cm^2) / (V * s)

    μn_i             = 2.00e1                * (cm^2) / (V * s)
    μp_i             = 2.00e1                * (cm^2) / (V * s)
    μa_i             = 1.00e-10              * (cm^2) / (V * s)

    μn_d             = 1.0e-3                * (cm^2) / (V * s)
    μp_d             = 1.0e-3                * (cm^2) / (V * s)

    μn               = [μn_a, μn_i, μn_d]
    μp               = [μp_a, μp_i, μp_d]
    μa               = [0.0,  μa_i, 0.0 ]

    ## relative dielectric permittivity
    ε_a              = 4.0                   *  1.0
    ε_i              = 23.0                  *  1.0
    ε_d              = 3.0                   *  1.0

    ε                = [ε_a, ε_i, ε_d]

    ## radiative recombination
    r0_a             = 6.3e-11               * cm^3 / s
    r0_i             = 3.6e-12               * cm^3 / s
    r0_d             = 6.8e-11               * cm^3 / s

    r0               = [r0_a, r0_i, r0_d]

    ## life times and trap densities
    τn_a             = 1.0e-6                * s
    τp_a             = 1.0e-6                * s

    τn_i             = 1.0e-7                * s
    τp_i             = 1.0e-7                * s
    τn_d             = τn_a
    τp_d             = τp_a

    τn               = [τn_a, τn_i, τn_d]
    τp               = [τp_a, τp_i, τp_d]

    ## SRH trap energies (needed for calculation of trap_density! (SRH))
    Ei_a             = -4.05                 * eV
    Ei_i             = -4.60                 * eV
    Ei_d             = -5.00                 * eV

    EI               = [Ei_a, Ei_i, Ei_d]

    ## doping
    Nd               = 2.089649130192123e17  / (cm^3)
    Na               = 4.529587947185444e18  / (cm^3)
    C0               = 1.0e18                / (cm^3)

    ## contact voltage
    voltageAcceptor  = 1.2                   * V

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)

    ## Possible choices: Stationary, Transient
    data.modelType                      = Transient

    ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    ## FermiDiracMinusOne, Blakemore
    data.F                              = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = false,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[bregionAcceptor]  = OhmicContact
    data.boundaryType[bregionJunction1] = InterfaceRecombination
    data.boundaryType[bregionJunction2] = InterfaceRecombination
    data.boundaryType[bregionDonor]     = OhmicContact

    ## Present ionic vacancies in perovskite layer
    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
    data.fluxApproximation             .= ExcessChemicalPotential

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                                       = Params(grid, numberOfCarriers)

    params.temperature                                           = T
    params.UT                                                    = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                                  = -1
    params.chargeNumbers[iphip]                                  =  1
    params.chargeNumbers[iphia]                                  =  1

    for ireg in 1:numberOfRegions ## interior region data

        params.dielectricConstant[ireg]                          = ε[ireg] * ε0

        ## effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]                      = NC[ireg]
        params.densityOfStates[iphip, ireg]                      = NV[ireg]
        params.densityOfStates[iphia, ireg]                      = NAnion[ireg]

        params.bandEdgeEnergy[iphin, ireg]                       = EC[ireg]
        params.bandEdgeEnergy[iphip, ireg]                       = EV[ireg]
        params.bandEdgeEnergy[iphia, ireg]                       = EA[ireg]

        params.mobility[iphin, ireg]                             = μn[ireg]
        params.mobility[iphip, ireg]                             = μp[ireg]
        params.mobility[iphia, ireg]                             = μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg]                      = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]             = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]             = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg]          = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg]          = trap_density!(iphip, ireg, data, EI[ireg])
    end

    ##############################################################
    ## inner boundary region data (we choose the intrinsic values)
    params.bDensityOfStates[iphin, bregionJunction1]             = Nc_i
    params.bDensityOfStates[iphip, bregionJunction1]             = Nv_i

    params.bDensityOfStates[iphin, bregionJunction2]             = Nc_i
    params.bDensityOfStates[iphip, bregionJunction2]             = Nv_i

    params.bBandEdgeEnergy[iphin, bregionJunction1]              = Ec_i
    params.bBandEdgeEnergy[iphip, bregionJunction1]              = Ev_i

    params.bBandEdgeEnergy[iphin, bregionJunction2]              = Ec_i
    params.bBandEdgeEnergy[iphip, bregionJunction2]              = Ev_i

    ## for surface recombination
    params.recombinationSRHvelocity[iphin, bregionJunction1]     = 1.0e1  * cm / s
    params.recombinationSRHvelocity[iphip, bregionJunction1]     = 1.0e5  * cm / s

    params.bRecombinationSRHTrapDensity[iphin, bregionJunction1] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJunction1] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    params.recombinationSRHvelocity[iphin, bregionJunction2]     = 1.0e7  * cm / s
    params.recombinationSRHvelocity[iphip, bregionJunction2]     = 1.0e1  * cm / s

    params.bRecombinationSRHTrapDensity[iphin, bregionJunction2] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJunction2] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    ##############################################################

    ## interior doping
    params.doping[iphin,  regionDonor]                           = Nd
    params.doping[iphip,  regionAcceptor]                        = Na
    params.doping[iphia,  regionIntrinsic]                       = C0

    data.params                                                  = params
    ctsys                                                        = System(grid, data, unknown_storage=unknown_storage)

    ## print all params stored in ctsys.data.params
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
    control.damp_initial = 0.9
    control.damp_growth  = 1.61 # >= 1
    control.max_round    = 5

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

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("I-V Measurement Loop")
    end
    ################################################################################

    ## there are different ways to control time stepping. Here we assume these primary data
    scanrate   = 1.0 * V/s
    ntsteps    = 31
    vend       = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend       = vend/scanrate

    ## with fixed timestep sizes we can calculate the times a priori
    tvalues    = range(0, stop = tend, length = ntsteps)

    ## for saving I-V data
    IV         = zeros(0) # for IV values
    biasValues = zeros(0) # for bias values

    for istep = 2:ntsteps

        t  = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep-1] # Time step size

        ## Apply new voltage (set non-equilibrium values)
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: Δt = $(t)")
        end

        solution = solve(ctsys, inival = inival, control = control, tstep = Δt)
        inival   = solution

        ## get I-V data
        current  = get_current_val(ctsys, solution, inival, Δt)

        push!(IV, current)
        push!(biasValues, Δu)

        if plotting
            label_solution = Array{String, 1}(undef, numberOfCarriers)
            label_solution[iphia]  = "\$ \\varphi_a\$"

            PyPlot.clf()
            plot_solution(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(Δu); \$E_a\$ =$(textEa)eV; \$N_a\$ =$textNa\$\\mathrm{cm}^{⁻3} \$", label_solution)
            PyPlot.pause(0.5)
        end

    end # time loop

    ##res = [biasValues, IV]

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solution))/length(solution) # when using sparse storage, we get NaN values in solution
    return testval


end # main

function test()
    testval = -0.5198379953833077
    main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
