#=
# PSC device on 2D domain (unstructured grid).
([source code](SOURCE_URL))

Simulating a three layer PSC device Pedot| MAPI | PCBM with mobile ions. The simulations are
performed in 2D on an unstructured grid, out of equilibrium and with abrupt interfaces.

The paramters are from Calado et al.:
https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/pedotpss_mapi_pcbm.csv.
(with adjustments on layer lengths)
=#

ENV["LC_NUMERIC"]="C" # put this in to work with Triangulate.jl, which is originally written in c++

module PSC_2D_unstructuredGrid

using ChargeTransport
using ExtendableGrids
using GridVisualize

## For using this example one additionally needs to add Triangulate.
## SimplexGridFactory is a wrapper for using this meshgenerator.
using SimplexGridFactory
using Triangulate

## problem with linux, when including PyPlot not until the end: "ERROR: LoadError: InitError: could not load library "/home/abdel/.julia/artifacts/8cc532f6a1ace8d1b756fc413f4ab340195ec3c3/lib/libgio-2.0.so"/home/abdel/.julia/artifacts/8cc532f6a1ace8d1b756fc413f4ab340195ec3c3/lib/libgobject-2.0.so.0: undefined symbol: g_uri_ref"
## It seems that this problem is common: https://discourse.julialang.org/t/could-not-load-library-librsvg-very-strange-error/21276
using PyPlot

function main(Plotter = PyPlot, ;plotting = false, verbose = false, test = true, unknown_storage=:sparse)

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
    bregionAcceptor  = 1
    bregionDonor     = 2
    bregionJunction1 = 3
    bregionJunction2 = 4
    bregionNoFlux    = 5

    ## grid
    h_pdoping        = 3.00e-6 * cm + 1.0e-7 *cm
    h_intrinsic      = 3.00e-5 * cm
    h_ndoping        = 8.50e-6 * cm + 1.0e-7 *cm
    height           = 1.00e-5 * cm

    b                = SimplexGridBuilder(Generator=Triangulate)

    ## specify boundary nodes
    length_0         = point!(b, 0.0, 0.0)
    length_p         = point!(b, h_pdoping, 0.0)
    length_pi        = point!(b, h_pdoping + h_intrinsic, 0.0)
    length_pin       = point!(b, h_pdoping + h_intrinsic + h_ndoping, 0.0)

    height_0         = point!(b, 0.0, height)
    height_p         = point!(b, h_pdoping, height)
    height_pi        = point!(b, h_pdoping + h_intrinsic, height)
    height_pin       = point!(b, h_pdoping + h_intrinsic + h_ndoping, height)

    ## specify boundary regions
    ## metal interface
    facetregion!(b, bregionAcceptor)
    facet!(b, length_0, height_0)
    facetregion!(b, bregionDonor)
    facet!(b, length_pin, height_pin)

    ## no flux
    facetregion!(b, bregionNoFlux)
    facet!(b, length_0, length_pin)
    facetregion!(b, bregionNoFlux)
    facet!(b, height_0, height_pin)

    ## inner interface
    facetregion!(b, bregionJunction1)
    facet!(b, length_p, height_p)
    facetregion!(b, bregionJunction2)
    facet!(b, length_pi, height_pi)

    ## cell regions
    cellregion!(b, regionAcceptor)
	regionpoint!(b, h_pdoping/2, height/2)
    cellregion!(b,regionIntrinsic)
	regionpoint!(b, (h_pdoping + h_intrinsic)/2, height/2)
    cellregion!(b,regionDonor)
	regionpoint!(b, h_pdoping + h_intrinsic + h_ndoping/2, height/2)

    options!(b,maxvolume=1.0e-16)

    grid           = simplexgrid(b)

    if plotting
        GridVisualize.gridplot(grid, Plotter= Plotter, resolution=(600,400),linewidth=0.5, legend=:lt)
        Plotter.title("Grid")
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
    iphia            = 3 # anion vacancy quasi Fermi potential
    numberOfCarriers = 3 # electrons, holes and anion vacancies

    ## temperature
    T                = 300.0                 *  K

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
    Nanion           = 1.0e18                / (cm^3)
    Ea_i             = -4.4                  *  eV
    ## for the labels in the figures
    textEa           = Ea_i./eV
    textNa           = Nanion.*cm^3
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

    ## SRH trap energies (needed for calculation of recombinationSRHTrapDensity)
    Ei_a             = -4.05                * eV
    Ei_i             = -4.60                * eV
    Ei_d             = -5.00                * eV

    EI               = [Ei_a, Ei_i, Ei_d]

    ## doping (doping values are from Phils paper, not stated in the parameter list online)
    Nd               = 2.089649130192123e17 / (cm^3)
    Na               = 4.529587947185444e18 / (cm^3)
    C0               = 1.0e18               / (cm^3)

    ## contact voltage
    voltageAcceptor  = 1.0                  * V

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
    data.modelType                     = Transient

    ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA, FermiDiracMinusOne, Blakemore
    data.F                             = [Boltzmann, Boltzmann, FermiDiracMinusOne]

    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = false,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor]    = OhmicContact

    ## Present ionic vacancies in perovskite layer
    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

    ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
    ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
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

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

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

    ## print data
    if test == false
        show_params(ctsys)
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
    control.tol_absolute      = 1.0e-10
    control.tol_relative      = 1.0e-10
    control.handle_exceptions = true
    control.tol_round         = 1.0e-10
    control.max_round         = 5

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

    if plotting # currently, plotting the solution was only tested with PyPlot.
        ipsi = data.index_psi
        X    = grid[Coordinates][1,:]
        Y    = grid[Coordinates][2,:]

        Plotter.figure()
        Plotter.surf(X[:], Y[:], solution[ipsi, :])
        Plotter.title("Electrostatic potential \$ \\psi \$ in Equilibrium")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        Plotter.tight_layout()
        ################
        Plotter.figure()
        Plotter.surf(X[:], Y[:], solution[iphin,:] )
        Plotter.title("quasi Fermi potential \$ \\varphi_n \$ in Equilibrium")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        Plotter.tight_layout()
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("I-V Measurement Loop")
    end
    ################################################################################

    data.calculationType = OutOfEquilibrium

    ## primary data for I-V scan protocol
    scanrate             = 0.04 * V/s
    number_tsteps        = 41
    endVoltage           = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend                 = endVoltage/scanrate

    ## with fixed timestep sizes we can calculate the times a priori
    tvalues              = range(0, stop = tend, length = number_tsteps)
    numbertsteps         = length(tvalues)

    ## for saving I-V data
    IV                   = zeros(0) # for IV values
    biasValues           = zeros(0) # for bias values

    for istep = 2:numbertsteps

        t  = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep-1] # Time step size

        ## Apply new voltage; set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: t = $(t)")
        end

        solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)

        ## get I-V data
        current = get_current_val(ctsys, solution, initialGuess, Δt)

        push!(IV, current)
        push!(biasValues, Δu)

        initialGuess .= solution
    end # time loop

    if test == false
        println("*** done\n")
    end

    if plotting
        Plotter.figure()
        Plotter.surf(X[:], Y[:], solution[ipsi, :])
        Plotter.title("Electrostatic potential \$ \\psi \$ at end time")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        ## ################
        Plotter.figure()
        Plotter.surf(X[:], Y[:], solution[iphin,:] )
        Plotter.title("quasi Fermi potential \$ \\varphi_n \$ at end time")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("height [m]")
        Plotter.zlabel("potential [V]")
        ## ################
        Plotter.figure()
        Plotter.plot(biasValues, IV.*(cm)^2/height, label = "\$ E_a =\$$(textEa)eV;  \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$ (without internal BC)",  linewidth= 3, linestyle="--", color="red")
        Plotter.title("Forward; \$ E_a =\$$(textEa)eV;  \$ N_a =\$ $textNa\$\\mathrm{cm}^{⁻3}\$ ")
        Plotter.ylabel("total current [A]") #
        Plotter.xlabel("Applied Voltage [V]")
    end

    testval = sum(filter(!isnan, solution))/length(solution) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = -0.6729879350578445
    main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
