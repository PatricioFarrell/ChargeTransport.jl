#=
# CIGS pn junction: stationary with trap and Schottky contacts.
([source code](SOURCE_URL))

Simulating stationary charge transport in a pn junction with hole traps and a Schottky boundary condition.
=#

module Example204_CIGS_stationary

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot

## function to initialize the grid for a possible extension to other p-i-n devices.
function initialize_pin_grid(refinementfactor, h_ndoping, h_pdoping_left, h_pdoping_trap, h_pdoping_right)
    coord_ndoping    = collect(range(0.0, stop = h_ndoping, length = 2 * refinementfactor))
    coord_pdoping_left  = collect(range(h_ndoping, stop = (h_ndoping + h_pdoping_left), length = 3 * refinementfactor))
    coord_pdoping_plus  = collect(range((h_ndoping + h_pdoping_left),
                                        stop = (h_ndoping + h_pdoping_left + h_pdoping_trap),
                                        length =  refinementfactor))
    coord_pdoping_right = collect(range((h_ndoping + h_pdoping_left + h_pdoping_trap),
                                        stop = (h_ndoping + h_pdoping_left + h_pdoping_trap + h_pdoping_right),
                                        length = 3 * refinementfactor))
    coord            = glue(coord_ndoping, coord_pdoping_left)
    coord            = glue(coord, coord_pdoping_plus)
    coord            = glue(coord, coord_pdoping_right)

    return coord
end

function main(;n = 3, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

    ################################################################################
    println("Set up grid and regions")
    ################################################################################

    ## region numbers
    regionDonor         = 1                           # n doped region
    regionAcceptorLeft  = 2                           # p doped region
    regionAcceptorTrap  = 3                           # p doped region with trap
    regionAcceptorRight = 4                           # p doped region
    regions             = [regionDonor, regionAcceptorLeft, regionAcceptorTrap, regionAcceptorRight]
    numberOfRegions     = length(regions)

    ## boundary region numbers
    bregionDonor            = 1
    bregionAcceptor         = 2
    bregions                = [bregionDonor, bregionAcceptor]
    numberOfBoundaryRegions = length(bregions)

    ## grid
    refinementfactor        = 2^(n-1)
    h_ndoping               = 0.5 * μm
    h_pdoping_left          = 1 * μm
    h_pdoping_trap          = 0.01 * μm#0.01 * μm
    h_pdoping_right         = 1 * μm
    coord                   = initialize_pin_grid(refinementfactor,
                                                  h_ndoping,
                                                  h_pdoping_left,
                                                  h_pdoping_trap,
                                                  h_pdoping_right)

    grid                    = simplexgrid(coord)

    ## set different regions in grid, doping profiles do not intersect
    ## n doped
    cellmask!(grid, [0.0 * μm], [h_ndoping], regionDonor)
    ## p doped
    cellmask!(grid, [h_ndoping], [h_ndoping + h_pdoping_left], regionAcceptorLeft)
    ## p doped with traps
    cellmask!(grid, [h_ndoping + h_pdoping_left], [h_ndoping + h_pdoping_left + h_pdoping_trap], regionAcceptorTrap)
    ## p doped
    cellmask!(grid, [h_ndoping + h_pdoping_left + h_pdoping_trap], [h_ndoping + h_pdoping_left + h_pdoping_trap + h_pdoping_right], regionAcceptorRight)

    if plotting
        gridplot(grid, Plotter = Plotter, legend=:lt)
        Plotter.title("Grid")
        Plotter.figure()
    end

    println("*** done\n")

    ################################################################################
    println("Define physical parameters and model")
    ################################################################################

    iphin             = 1 # index electron quasi Fermi potential
    iphip             = 2 # index hole quasi Fermi potential
    iphit             = 3 # index trap quasi Fermi potential
    numberOfCarriers  = 3 # electrons, holes and traps

    ## physical data
    # Ec                = 1.424                *  eV
    # Ev                = 0.0                  *  eV
    # Et                = 0.6                  *  eV

    Ec_CIGS           = 3.4                  *  eV
    Ev_CIGS           = 2.3                  *  eV
    Ec_ZnO            = 3.4                  *  eV
    Ev_ZnO            = 0.0                  *  eV
    Et                = 2.8                  *  eV

    Nc                = 4.351959895879690e17 / (cm^3)
    Nv                = 9.139615903601645e18 / (cm^3)
    Nt                = 5e14                 / (cm^3)
    Nt_low            = Nt#/1e3
    mun_CIGS          = 100.0                * (cm^2) / (V * s)
    mup_CIGS          = 25                   * (cm^2) / (V * s)
    mun_ZnO           = 100                  * (cm^2) / (V * s)
    mup_ZnO           = 25                   * (cm^2) / (V * s)
    mut               = 0                    * (cm^2) / (V * s)  # no flux for traps
    εr_CIGS           = 13.6                 *  1.0
    εr_ZnO            = 9                    *  1.0
    T                 = 300.0                *  K

    An                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    Ap                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
    vn                = An * T^2 / (q*Nc)
    vp                = Ap * T^2 / (q*Nv)
    barrier           = Ev_CIGS + 0.1 * eV

    ## recombination parameters
    ni_CIGS           = sqrt(Nc * Nv) * exp(-(Ec_CIGS - Ev_CIGS) / (2 * kB * T)) # intrinsic concentration
    n0_CIGS           = Nc * Boltzmann( (Et-Ec_CIGS) / (kB*T) )             # Boltzmann equilibrium concentration
    p0_CIGS           = ni_CIGS^2 / n0_CIGS                                      # Boltzmann equilibrium concentration
    ni_ZnO            = sqrt(Nc * Nv) * exp(-(Ec_ZnO - Ev_ZnO) / (2 * kB * T)) # intrinsic concentration
    n0_ZnO            = Nc * Boltzmann( (Et-Ec_ZnO) / (kB*T) )             # Boltzmann equilibrium concentration
    p0_ZnO            = ni_ZnO^2 / n0_ZnO                                      # Boltzmann equilibrium concentration
    Auger             = 1.0e-29  * cm^6 / s          # 1.0e-41 m^6 / s
    SRH_LifeTime      = 1.0e-3   * ns
    Radiative         = 1.0e-10  * cm^3 / s          # 1.0e-16 m^3 / s
    G                 = 1.0e20   / (cm^3 * s)
    A_CIGS            = 1.0e5    / cm
    A_ZnO             = 0.0      / cm
    N0                = 1e17     / cm^2/s

    ## doping -- trap doping will not be set and thus automatically zero
    Nd                = 1.0e18 / (cm^3)
    Na                = 5.5e15 / (cm^3)

    ## we will impose this applied voltage on one boundary
    voltageAcceptor   = 2.0 * V

    println("*** done\n")
    ################################################################################
    println("Define System and fill in information about model")
    ################################################################################

    ## initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)

    ## possible choices: model_stationary, model_transient
    data.model_type                     = model_stationary

    ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                             .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    data.bulk_recombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                  bulk_recomb_Auger = true,
                                                                  bulk_recomb_radiative = true,
                                                                  bulk_recomb_SRH = true)

    data.generation_model               = generation_beer_lambert #generation_uniform #generation_beer_lambert

    ## Here, the user gives information on which indices belong to ionic charge carriers and in which regions these charge carriers are present.
    ## In this application ion vacancies only live in active perovskite layer
    data.boundary_type[bregionAcceptor] = schottky_contact
    data.boundary_type[bregionDonor]    = ohmic_contact

    ## possible choices: scharfetter_gummel, scharfetter_gummel_graded, excess_chemical_potential,
    ## excess_chemical_potential_graded, diffusion_enhanced, generalized_sg
    data.flux_approximation             = excess_chemical_potential

    println("*** done\n")

    ################################################################################
    println("Define Params and fill in physical parameters")
    ################################################################################

    ## physical parameters
    params                                              = Params(grid, numberOfCarriers)
    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1
    params.chargeNumbers[iphit]                         =  1  # +1: hole trap is used

    for ibreg in 1:numberOfBoundaryRegions   # boundary region data
        params.bDensityOfStates[iphin, ibreg]           = Nc
        params.bDensityOfStates[iphip, ibreg]           = Nv
        params.bDensityOfStates[iphit, ibreg]           = Nt_low
        params.bBandEdgeEnergy[iphit, ibreg]            = Et
    end

    params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_ZnO
    params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_ZnO
    params.bBandEdgeEnergy[iphin, bregionAcceptor]      = Ec_CIGS
    params.bBandEdgeEnergy[iphip, bregionAcceptor]      = Ev_CIGS

    for ireg in 1:numberOfRegions           # interior region data

        params.dielectricConstant[ireg]                 = εr_CIGS

        ## effective DOS, band-edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nc
        params.densityOfStates[iphip, ireg]             = Nv
        params.bandEdgeEnergy[iphin, ireg]              = Ec_CIGS
        params.bandEdgeEnergy[iphip, ireg]              = Ev_CIGS
        params.bandEdgeEnergy[iphit, ireg]              = Et
        params.mobility[iphin, ireg]                    = mun_CIGS
        params.mobility[iphip, ireg]                    = mup_CIGS
        params.mobility[iphit, ireg]                    = mut

        ## recombination parameters
        params.recombinationRadiative[ireg]             = Radiative
        params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
        params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
        params.recombinationSRHTrapDensity[iphin, ireg] = n0_CIGS
        params.recombinationSRHTrapDensity[iphip, ireg] = p0_CIGS
        params.recombinationAuger[iphin, ireg]          = Auger
        params.recombinationAuger[iphip, ireg]          = Auger

        ## generation parameters
        params.generationAbsorption[ireg]               = A_CIGS
        params.generationIncidentPhotonFlux[ireg]       = N0
        params.generationUniform[ireg]                  = G

    end

    ## overwrite parameters in ZnO donor region
    params.generationUniform[regionDonor]                  = 0.0      # only used if for "generation_uniform"
    params.generationAbsorption[regionDonor]               = A_ZnO    # only used if for "generation_beer_lambert"
    params.generationIncidentPhotonFlux[regionDonor]       = N0
    params.recombinationSRHTrapDensity[iphin, regionDonor] = n0_ZnO
    params.recombinationSRHTrapDensity[iphip, regionDonor] = p0_ZnO
    params.bandEdgeEnergy[iphin, regionDonor]              = Ec_ZnO
    params.bandEdgeEnergy[iphip, regionDonor]              = Ev_ZnO
    params.dielectricConstant[regionDonor]                 = εr_ZnO
    params.mobility[iphin, regionDonor]                    = mun_ZnO
    params.mobility[iphip, regionDonor]                    = mup_ZnO

    ## hole trap density only high in grain
    params.densityOfStates[iphit, regionDonor]          = Nt_low
    params.densityOfStates[iphit, regionAcceptorLeft]   = Nt_low
    params.densityOfStates[iphit, regionAcceptorTrap]   = Nt
    params.densityOfStates[iphit, regionAcceptorRight]  = Nt_low

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
    params.bVelocity[iphin,bregionAcceptor]             = vn
    params.bVelocity[iphip,bregionAcceptor]             = vp

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    show_params(ctsys)
    println("*** done\n")

    ################################################################################
    println("Define outer boundary conditions and enabled layers")
    ################################################################################

    ## set ohmic contact in bregionDonor and schottky contact in bregionAcceptor
    set_ohmic_contact!(ctsys, bregionDonor, 0.0)
    set_schottky_contact!(ctsys, bregionAcceptor, appliedVoltage = 0.0)

    println("*** done\n")

    ################################################################################
    println("Define control parameters for Newton solver")
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

    println("*** done\n")

    ################################################################################
    println("Compute solution in thermodynamic equilibrium")
    ################################################################################

    ## initialize solution and starting vectors
    initialGuess          = unknowns(ctsys)
    solution              = unknowns(ctsys)

    data.calculation_type = inEquilibrium

    ## solve thermodynamic equilibrium and update initial guess
    solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    initialGuess         .= solution

    println("*** done\n")

    if plotting
        ## ##### set legend for plotting routines #####
        label_energy   = Array{String, 2}(undef, 2, numberOfCarriers) # band-edge energies and potential
        label_density  = Array{String, 1}(undef, numberOfCarriers)
        label_solution = Array{String, 1}(undef, numberOfCarriers)

        ## for electrons
        label_energy[1, iphin] = "\$E_c-q\\psi\$";       label_energy[2, iphin] = "\$ - q \\varphi_n\$"
        label_density[iphin]   = "n";                    label_solution[iphin]  = "\$ \\varphi_n\$"

        ## for holes
        label_energy[1, iphip] = "\$E_v-q\\psi\$";       label_energy[2, iphip] = "\$ - q \\varphi_p\$"
        label_density[iphip]   = "p";                    label_solution[iphip]  = "\$ \\varphi_p\$"

        ## for traps
        label_energy[1, iphit] = "\$E_{\\tau}-q\\psi\$"; label_energy[2, iphit] = "\$ - q \\varphi_{\\tau}\$"
        label_density[iphit]   = "\$n_{\\tau}\$";        label_solution[iphit]  = "\$ \\varphi_{\\tau}\$"
        ## ##### set legend for plotting routines #####
        plot_energies(Plotter, grid, data, solution, "Equilibrium", label_energy)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution,"Equilibrium", label_density)
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "Equilibrium", label_solution)
    end

    ################################################################################
    println("Stationary bias loop")
    ################################################################################

    ## set calculation type to outOfEquilibrium for starting with respective simulation.
    ctsys.data.calculation_type   = outOfEquilibrium      # Rn = Rp = R, since the model type is stationary
    endVoltage                    = voltageAcceptor       # final bias value

    IV         = zeros(0)
    maxBias    = voltageAcceptor
    biasSteps  = 52
    biasValues = collect(range(0, stop = maxBias, length = biasSteps))
    chargeDensities = zeros(0)

    w_device = 0.5    * μm  # width of device
    z_device = 1.0e-4 * cm  # depth of device

    ## adjust Newton parameters
    control.tol_absolute      = 1.0e-10
    control.tol_relative      = 1.0e-10
    control.tol_round         = 1.0e-7
    control.damp_initial      = 0.5
    control.damp_growth       = 1.2
    control.max_iterations    = 30
    control.max_round         = 3

    for i in eachindex(biasValues)

        Δu = biasValues[i] # bias

        ## set non equilibrium boundary condition
        set_schottky_contact!(ctsys, bregionAcceptor, appliedVoltage = Δu)

        ## increase generation rate with bias
        ctsys.data.λ2 = 10.0^(-biasSteps + i)

        if verbose
            println("bias: Δu = $(Δu)")
        end

        ## solve time step problems with timestep Δt
        solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

        ## save IV data
        current = get_current_val(ctsys, solution)
        push!(IV, w_device * z_device * current)

        ## store charge density in donor region (ZnO)
        push!(chargeDensities,chargeDensity(ctsys,solution)[regionDonor])

        initialGuess .= solution

    end # bias loop

    println("*** done\n")

    ## compute static capacitance: check this is correctly computed
    staticCapacitance = diff(chargeDensities) ./ diff(biasValues)

    ## plot solution and IV curve
    if plotting
        plot_energies(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(0)\$", label_energy)
        Plotter.figure()
        plot_densities(Plotter, grid, data, solution,"bias \$\\Delta u\$ = $(endVoltage), \$ t=$(0)\$", label_density)
        Plotter.figure()
        plot_solution(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(0)\$", label_solution)
        Plotter.figure()
        plot_IV(Plotter, biasValues,IV, biasValues[end], plotGridpoints = true)
        Plotter.figure()
        plot_IV(Plotter, biasValues,chargeDensities, biasValues[end], plotGridpoints = true)
        Plotter.title("Charge density in donor region")
        Plotter.ylabel("Charge density [C]")
        Plotter.figure()
        plot_IV(Plotter, biasValues,staticCapacitance, biasValues[end-1], plotGridpoints = true)
        Plotter.title("Static capacitance in donor region")
        Plotter.ylabel("Static capacitance [C/V]")

        ## Plotter.figure()
        ## Plotter.yscale("symlog")
        ## dens = compute_densities!(grid, data, solution)
        ## n = dens[iphin,:] #.*1e6 # cm^(-3)
        ## p = dens[iphip,:] #.*1e6 # cm^(-3)
        ## t = dens[iphit,:]
        ## # p_tr = N_t - n_tr
        ## plot(coord, 1e-6*(Nt .- (p0*Nt .+ n.*Nt) ./ ((p0 .+p) .+ (n0 .+n))) )
        ## Plotter.title("Check p_traps agrees with computed traps" )

        ## Plotter.figure()
        ## Plotter.yscale("symlog")
        ## plot(coord, 1e-6*t .- 1e-6*(Nt .- (p0*Nt .+ n.*Nt) ./ ((p0 .+p) .+ (n0 .+n))) )
        ## Plotter.title("Error" )

    end


    ## ipsi                          = data.index_psi
    ## number_tsteps                 = 41
    ## tend                          = 1*s
    ## tvalues                       = range(0,stop=tend,length=number_tsteps)

    ## time loop
    ## for istep = 2:number_tsteps

    ##     t                     = tvalues[istep]          # Actual time
    ##     Δt                    = t - tvalues[istep-1]    # Time step size

    ##     if verbose
    ##         println("time:  = $(t)")
    ##     end

    ##      # Solve time step problems with timestep Δt
    ##     solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)

    ##     # get I-V data
    ##     current = get_current_val(ctsys, solution)

    ##     push!(IV, w_device * z_device * current)
    ##     push!(biasValues,endVoltage)

    ##     initialGuess .= solution

    ## end # time loop

    ## if test == false
    ##     println("*** done\n")
    ## end

    ## plot solution and IV curve
    ## if plotting
    ##     plot_energies(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tend)\$", label_energy)
    ##     Plotter.figure()
    ##     plot_densities(Plotter, grid, data, solution,"bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tend)\$", label_density)
    ##     Plotter.figure()
    ##     plot_solution(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tend)\$", label_solution)
    ##     Plotter.figure()
    ##     plot_IV(Plotter, biasValues,IV, biasValues[end], plotGridpoints = true)
    ## end

    ## println("Max error")
    ## @show max(abs.(solution_stationary_bias[iphin,:].-solution[iphin,:])...)
    ## @show max(abs.(solution_stationary_bias[iphip,:].-solution[iphip,:])...)
    ## @show max(abs.(solution_stationary_bias[iphit,:].-solution[iphit,:])...)
    ## @show max(abs.(solution_stationary_bias[ipsi,:].-solution[ipsi,:])...)

    testval = solution[data.index_psi, 10]
    return testval

    println("*** done\n")

end #  main

function test()
    testval = 1.3214196490674017
    main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
end

println("This message should show when this module has successfully recompiled.")


end # module




#######################################



#=
# CIGS pn junction: stationary with trap and Schottky contacts.
([source code](SOURCE_URL))

Simulating stationary charge transport in a pn junction with hole traps and a Schottky boundary condition.
=#

# module Example204_CIGS_stationary

# using VoronoiFVM
# using ChargeTransport
# using ExtendableGrids
# using GridVisualize
# using PyPlot

# ## function to initialize the grid for a possible extension to other p-i-n devices.
# function initialize_pin_grid(refinementfactor, h_ndoping, h_pdoping_left, h_pdoping_trap, h_pdoping_right)
#     coord_ndoping    = collect(range(0.0, stop = h_ndoping, length = 2 * refinementfactor))
#     coord_pdoping_left  = collect(range(h_ndoping, stop = (h_ndoping + h_pdoping_left), length = 3 * refinementfactor))
#     coord_pdoping_plus  = collect(range((h_ndoping + h_pdoping_left),
#                                         stop = (h_ndoping + h_pdoping_left + h_pdoping_trap),
#                                         length =  refinementfactor))
#     coord_pdoping_right = collect(range((h_ndoping + h_pdoping_left + h_pdoping_trap),
#                                         stop = (h_ndoping + h_pdoping_left + h_pdoping_trap + h_pdoping_right),
#                                         length = 3 * refinementfactor))
#     coord            = glue(coord_ndoping, coord_pdoping_left)
#     coord            = glue(coord, coord_pdoping_plus)
#     coord            = glue(coord, coord_pdoping_right)

#     return coord
# end

# function main(;n = 3, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:sparse)

#     ################################################################################
#     println("Set up grid and regions")
#     ################################################################################

#     ## region numbers
#     regionDonor         = 1                           # n doped region
#     regionAcceptorLeft  = 2                           # p doped region
#     regionAcceptorTrap  = 3                           # p doped region with trap
#     regionAcceptorRight = 4                           # p doped region
#     regions             = [regionDonor, regionAcceptorLeft, regionAcceptorTrap, regionAcceptorRight]
#     numberOfRegions     = length(regions)

#     ## boundary region numbers
#     bregionDonor            = 1
#     bregionAcceptor         = 2
#     bregions                = [bregionDonor, bregionAcceptor]
#     numberOfBoundaryRegions = length(bregions)

#     ## grid
#     refinementfactor        = 2^(n-1)
#     h_ndoping               = 0.5 * μm
#     h_pdoping_left          = 1 * μm
#     h_pdoping_trap          = 0.01 * μm#0.01 * μm
#     h_pdoping_right          = 1 * μm
#     coord                   = initialize_pin_grid(refinementfactor,
#                                                   h_ndoping,
#                                                   h_pdoping_left,
#                                                   h_pdoping_trap,
#                                                   h_pdoping_right)

#     grid                    = simplexgrid(coord)

#     ## set different regions in grid, doping profiles do not intersect
#     ## n doped
#     cellmask!(grid, [0.0 * μm], [h_ndoping], regionDonor)
#     ## p doped
#     cellmask!(grid, [h_ndoping], [h_ndoping + h_pdoping_left], regionAcceptorLeft)
#     ## p doped with traps
#     cellmask!(grid, [h_ndoping + h_pdoping_left], [h_ndoping + h_pdoping_left + h_pdoping_trap], regionAcceptorTrap)
#     ## p doped
#     cellmask!(grid, [h_ndoping + h_pdoping_left + h_pdoping_trap], [h_ndoping + h_pdoping_left + h_pdoping_trap + h_pdoping_right], regionAcceptorRight)

#     if plotting
#         gridplot(grid, Plotter = Plotter, legend=:lt)
#         Plotter.title("Grid")
#         Plotter.figure()
#     end

#     println("*** done\n")

#     ################################################################################
#     println("Define physical parameters and model")
#     ################################################################################

#     iphin             = 1 # index electron quasi Fermi potential
#     iphip             = 2 # index hole quasi Fermi potential
#     iphit             = 3 # index trap quasi Fermi potential
#     numberOfCarriers  = 3 # electrons, holes and traps

#     ## physical data
#     Ec_CIGS           = 3.4                  *  eV
#     Ev_CIGS           = 2.3                  *  eV
#     Ec_ZnO            = 3.4                  *  eV
#     Ev_ZnO            = 0.0                  *  eV
#     Et                = 0.6                  *  eV
#     Nc                = 4.351959895879690e17 / (cm^3)
#     Nv                = 9.139615903601645e18 / (cm^3)
#     Nt                = 5e14                / (cm^3)
#     Nt_low            = Nt#/1e3
#     mun_CIGS          = 100.0                * (cm^2) / (V * s)
#     mup_CIGS          = 25                   * (cm^2) / (V * s)
#     mun_ZnO           = 100                  * (cm^2) / (V * s)
#     mup_ZnO           = 25                   * (cm^2) / (V * s)
#     mut               = 0                    * (cm^2) / (V * s)  # no flux for traps
#     εr_CIGS           = 13.6                 *  1.0
#     εr_ZnO            = 9                    *  1.0
#     T                 = 300.0                *  K
#     A_CIGS            = 1.0e5                / cm
#     A_ZnO             = 0.0                  / cm
#     N0                = 1e17                 / cm^2/s

#     An                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
#     Ap                = 4 * pi * q * mₑ * kB^2 / Planck_constant^3
#     vn                = An * T^2 / (q*Nc)
#     vp                = Ap * T^2 / (q*Nv)
#     barrier           = 0.1 * eV

#     ## recombination parameters
#     ni_CIGS           = sqrt(Nc * Nv) * exp(-(Ec_CIGS - Ev_CIGS) / (2 * kB * T)) # intrinsic concentration
#     n0_CIGS           = Nc * Boltzmann( (Et-Ec_CIGS) / (kB*T) )             # Boltzmann equilibrium concentration
#     p0_CIGS           = ni_CIGS^2 / n0_CIGS                                      # Boltzmann equilibrium concentration
#     ni_ZnO            = sqrt(Nc * Nv) * exp(-(Ec_ZnO - Ev_ZnO) / (2 * kB * T)) # intrinsic concentration
#     n0_ZnO            = Nc * Boltzmann( (Et-Ec_ZnO) / (kB*T) )             # Boltzmann equilibrium concentration
#     p0_ZnO            = ni_ZnO^2 / n0_ZnO                                      # Boltzmann equilibrium concentration
#     Auger             = 1.0e-29  * cm^6 / s          # 1.0e-41 m^6 / s
#     SRH_LifeTime      = 1.0e-3   * ns
#     Radiative         = 1.0e-10  * cm^3 / s          # 1.0e-16 m^3 / s
#     G                 = 1.0e20  / (cm^3 * s)

#     ## doping -- trap doping will not be set and thus automatically zero
#     Nd                = 1.0e18 / (cm^3)
#     Na                = 5.5e15 / (cm^3)

#     ## we will impose this applied voltage on one boundary
#     voltageAcceptor   = 2.0 * V

#     println("*** done\n")
#     ################################################################################
#     println("Define System and fill in information about model")
#     ################################################################################

#     ## initialize Data instance and fill in data
#     data                                = Data(grid, numberOfCarriers)

#     ## possible choices: model_stationary, model_transient
#     data.model_type                     = model_stationary

#     ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
#     data.F                             .= [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

#     data.bulk_recombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
#                                                                   bulk_recomb_Auger = true,
#                                                                   bulk_recomb_radiative = true,
#                                                                   bulk_recomb_SRH = true)

#     ## possible choices: ohmic_contact, schottky_contact (outer boundary) and interface_model_none,
#     ## interface_model_surface_recombination (inner boundary).
#     data.generation_model               = generation_beer_lambert # generation_uniform

#     ## Here, the user gives information on which indices belong to ionic charge carriers and in which regions these charge carriers are present.
#     ## In this application ion vacancies only live in active perovskite layer
#     data.boundary_type[bregionAcceptor] = schottky_contact
#     data.boundary_type[bregionDonor]    = ohmic_contact

#     ## possible choices: scharfetter_gummel, scharfetter_gummel_graded, excess_chemical_potential,
#     ## excess_chemical_potential_graded, diffusion_enhanced, generalized_sg
#     data.flux_approximation             = excess_chemical_potential

#     println("*** done\n")

#     ################################################################################
#     println("Define Params and fill in physical parameters")
#     ################################################################################

#     ## physical parameters
#     params                                              = Params(grid, numberOfCarriers)
#     params.temperature                                  = T
#     params.UT                                           = (kB * params.temperature) / q
#     params.chargeNumbers[iphin]                         = -1
#     params.chargeNumbers[iphip]                         =  1
#     params.chargeNumbers[iphit]                         =  1  # +1: hole trap is used

#     # boundary region data
#     for ibreg in 1:numberOfBoundaryRegions   # boundary region data

#         params.bDensityOfStates[iphin, ibreg]           = Nc
#         params.bDensityOfStates[iphip, ibreg]           = Nv
#         params.bDensityOfStates[iphit, ibreg]           = Nt_low
#         params.bBandEdgeEnergy[iphit, ibreg]            = Et
#     end

#     params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_ZnO
#     params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_ZnO
#     params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_CIGS
#     params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_CIGS

#     # only CIGS regions
#     for ireg in [regionAcceptorLeft, regionAcceptorTrap, regionAcceptorRight]

#         params.dielectricConstant[ireg] = εr_CIGS

#         ## effective DOS, band-edge energy and mobilities
#         params.densityOfStates[iphin, ireg] = Nc
#         params.densityOfStates[iphip, ireg] = Nv
#         params.bandEdgeEnergy[iphin, ireg] = Ec_CIGS
#         params.bandEdgeEnergy[iphip, ireg] = Ev_CIGS
#         params.bandEdgeEnergy[iphit, ireg] = Et
#         params.mobility[iphin, ireg] = mun_CIGS
#         params.mobility[iphip, ireg] = mup_CIGS
#         params.mobility[iphit, ireg] = mut

#         ## recombination parameters
#         params.recombinationRadiative[ireg] = Radiative
#         params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
#         params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
#         params.recombinationSRHTrapDensity[iphin, ireg] = n0_CIGS
#         params.recombinationSRHTrapDensity[iphip, ireg] = p0_CIGS
#         params.recombinationAuger[iphin, ireg] = Auger
#         params.recombinationAuger[iphip, ireg] = Auger

#         ## generation parameters
#         params.generationAbsorption[ireg]   = A_CIGS
#         params.generationIncidentPhotonFlux[ireg]  = N0

#     end

#     # ZnO region (donor region)
#     params.mobility[iphin, regionDonor] = mun_ZnO
#     params.mobility[iphip, regionDonor] = mup_ZnO
#     params.mobility[iphit, regionDonor] = mut
#     params.dielectricConstant[regionDonor] = εr_ZnO
#     params.bandEdgeEnergy[iphin, regionDonor] = Ec_ZnO
#     params.bandEdgeEnergy[iphip, regionDonor] = Ev_ZnO
#     params.recombinationSRHTrapDensity[iphin, regionDonor] = n0_ZnO
#     params.recombinationSRHTrapDensity[iphip, regionDonor] = p0_ZnO

#     ## (uniform) generation only in donor region: can be improved later with exp decay
#     # params.generationUniform[regionDonor]          = 0.0
#     # params.generationUniform[regionAcceptorLeft]   = G
#     # params.generationUniform[regionAcceptorTrap]   = G
#     # params.generationUniform[regionAcceptorRight]  = G
#     params.generationAbsorption[regionDonor]          = A_ZnO
#     params.generationIncidentPhotonFlux[regionDonor]  = N0

#     ## hole trap density only high in grain
#     params.densityOfStates[iphit, regionDonor]          = Nt_low
#     params.densityOfStates[iphit, regionAcceptorLeft]   = Nt_low
#     params.densityOfStates[iphit, regionAcceptorTrap]   = Nt
#     params.densityOfStates[iphit, regionAcceptorRight]  = Nt_low

#     ## doping -- since we do not set any doping for the traps it is automatically zero
#     params.doping[iphin, regionDonor]                   = Nd
#     params.doping[iphip, regionAcceptorLeft]            = Na
#     params.doping[iphip, regionAcceptorTrap]            = Na
#     params.doping[iphip, regionAcceptorRight]           = Na

#     ## boundary doping
#     params.bDoping[iphin, bregionDonor]                 = Nd
#     params.bDoping[iphip, bregionAcceptor]              = Na

#     ## values for the schottky contacts
#     params.SchottkyBarrier[bregionAcceptor]             = barrier
#     params.bVelocity[iphin,bregionAcceptor]             = vn
#     params.bVelocity[iphip,bregionAcceptor]             = vp

#     data.params                                         = params
#     ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

#     show_params(ctsys)
#     println("*** done\n")

#     ################################################################################
#     println("Define outer boundary conditions and enabled layers")
#     ################################################################################

#     ## set ohmic contact in bregionDonor and schottky contact in bregionAcceptor
#     set_ohmic_contact!(ctsys, bregionDonor, 0.0)
#     set_schottky_contact!(ctsys, bregionAcceptor, appliedVoltage = 0.0)

#     println("*** done\n")

#     ################################################################################
#     println("Define control parameters for Newton solver")
#     ################################################################################

#     control                   = NewtonControl()
#     control.verbose           = verbose
#     control.damp_initial      = 0.5
#     control.damp_growth       = 1.21    #>= 1
#     control.max_iterations    = 250
#     control.tol_absolute      = 1.0e-14
#     control.tol_relative      = 1.0e-14
#     control.handle_exceptions = true
#     control.tol_round         = 1.0e-8
#     control.max_round         = 5

#     println("*** done\n")

#     ################################################################################
#     println("Compute solution in thermodynamic equilibrium")
#     ################################################################################

#     ## initialize solution and starting vectors
#     initialGuess          = unknowns(ctsys)
#     solution              = unknowns(ctsys)

#     data.calculation_type = inEquilibrium

#     ## solve thermodynamic equilibrium and update initial guess
#     solution              = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
#     initialGuess         .= solution

#     println("*** done\n")

#     if plotting
#         ## ##### set legend for plotting routines #####
#         label_energy   = Array{String, 2}(undef, 2, numberOfCarriers) # band-edge energies and potential
#         label_density  = Array{String, 1}(undef, numberOfCarriers)
#         label_solution = Array{String, 1}(undef, numberOfCarriers)

#         ## for electrons
#         label_energy[1, iphin] = "\$E_c-q\\psi\$";       label_energy[2, iphin] = "\$ - q \\varphi_n\$"
#         label_density[iphin]   = "n";                    label_solution[iphin]  = "\$ \\varphi_n\$"

#         ## for holes
#         label_energy[1, iphip] = "\$E_v-q\\psi\$";       label_energy[2, iphip] = "\$ - q \\varphi_p\$"
#         label_density[iphip]   = "p";                    label_solution[iphip]  = "\$ \\varphi_p\$"

#         ## for traps
#         label_energy[1, iphit] = "\$E_{\\tau}-q\\psi\$"; label_energy[2, iphit] = "\$ - q \\varphi_{\\tau}\$"
#         label_density[iphit]   = "\$n_{\\tau}\$";        label_solution[iphit]  = "\$ \\varphi_{\\tau}\$"
#         ## ##### set legend for plotting routines #####
#         plot_energies(Plotter, grid, data, solution, "Equilibrium", label_energy)
#         Plotter.figure()
#         plot_densities(Plotter, grid, data, solution,"Equilibrium", label_density)
#         Plotter.figure()
#         plot_solution(Plotter, grid, data, solution, "Equilibrium", label_solution)
#     end

#     ################################################################################
#     println("Stationary bias loop")
#     ################################################################################

#     ## set calculation type to outOfEquilibrium for starting with respective simulation.
#     ctsys.data.calculation_type   = outOfEquilibrium      # Rn = Rp = R, since the model type is stationary
#     endVoltage                    = voltageAcceptor       # final bias value

#     IV         = zeros(0)
#     maxBias    = voltageAcceptor
#     biasValues = collect(range(0, stop = maxBias, length = 52))
#     chargeDensities = zeros(0)

#     w_device = 0.5    * μm  # width of device
#     z_device = 1.0e-4 * cm  # depth of device

#     ## adjust Newton parameters
#     control.tol_absolute      = 1.0e-10
#     control.tol_relative      = 1.0e-10
#     control.tol_round         = 1.0e-7
#     control.damp_initial      = 0.5
#     control.damp_growth       = 1.2
#     control.max_iterations    = 30
#     control.max_round         = 3

#     for i in eachindex(biasValues)

#         Δu = biasValues[i] # bias

#         ## set non equilibrium boundary condition
#         set_schottky_contact!(ctsys, bregionAcceptor, appliedVoltage = Δu)

#         if verbose
#             println("bias: Δu = $(Δu)")
#         end

#         ## solve time step problems with timestep Δt
#         solve!(solution, initialGuess, ctsys, control  = control, tstep = Inf)

#         ## save IV data
#         current = get_current_val(ctsys, solution)
#         push!(IV, w_device * z_device * current)

#         ## story charge density in donor region (ZnO)
#         push!(chargeDensities,chargeDensity(ctsys,solution)[regionDonor])

#         initialGuess .= solution

#     end # bias loop

#     println("*** done\n")

#     ## compute static capacitance: check this is correctly computed
#     staticCapacitance = diff(chargeDensities) ./ diff(biasValues)

#     ## plot solution and IV curve
#     if plotting
#         plot_energies(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(0)\$", label_energy)
#         Plotter.figure()
#         plot_densities(Plotter, grid, data, solution,"bias \$\\Delta u\$ = $(endVoltage), \$ t=$(0)\$", label_density)
#         Plotter.figure()
#         plot_solution(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(0)\$", label_solution)
#         Plotter.figure()
#         plot_IV(Plotter, biasValues,IV, biasValues[end], plotGridpoints = true)
#         Plotter.figure()
#         plot_IV(Plotter, biasValues,chargeDensities, biasValues[end], plotGridpoints = true)
#         Plotter.title("Charge density in donor region")
#         Plotter.ylabel("Charge density [C]")
#         Plotter.figure()
#         plot_IV(Plotter, biasValues,staticCapacitance, biasValues[end-1], plotGridpoints = true)
#         Plotter.title("Static capacitance in donor region")
#         Plotter.ylabel("Static capacitance [C/V]")

#         ## Plotter.figure()
#         ## Plotter.yscale("symlog")
#         ## dens = compute_densities!(grid, data, solution)
#         ## n = dens[iphin,:] #.*1e6 # cm^(-3)
#         ## p = dens[iphip,:] #.*1e6 # cm^(-3)
#         ## t = dens[iphit,:]
#         ## # p_tr = N_t - n_tr
#         ## plot(coord, 1e-6*(Nt .- (p0*Nt .+ n.*Nt) ./ ((p0 .+p) .+ (n0 .+n))) )
#         ## Plotter.title("Check p_traps agrees with computed traps" )

#         ## Plotter.figure()
#         ## Plotter.yscale("symlog")
#         ## plot(coord, 1e-6*t .- 1e-6*(Nt .- (p0*Nt .+ n.*Nt) ./ ((p0 .+p) .+ (n0 .+n))) )
#         ## Plotter.title("Error" )

#     end


#     ## ipsi                          = data.index_psi
#     ## number_tsteps                 = 41
#     ## tend                          = 1*s
#     ## tvalues                       = range(0,stop=tend,length=number_tsteps)

#     ## time loop
#     ## for istep = 2:number_tsteps

#     ##     t                     = tvalues[istep]          # Actual time
#     ##     Δt                    = t - tvalues[istep-1]    # Time step size

#     ##     if verbose
#     ##         println("time:  = $(t)")
#     ##     end

#     ##      # Solve time step problems with timestep Δt
#     ##     solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)

#     ##     # get I-V data
#     ##     current = get_current_val(ctsys, solution)

#     ##     push!(IV, w_device * z_device * current)
#     ##     push!(biasValues,endVoltage)

#     ##     initialGuess .= solution

#     ## end # time loop

#     ## if test == false
#     ##     println("*** done\n")
#     ## end

#     ## plot solution and IV curve
#     ## if plotting
#     ##     plot_energies(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tend)\$", label_energy)
#     ##     Plotter.figure()
#     ##     plot_densities(Plotter, grid, data, solution,"bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tend)\$", label_density)
#     ##     Plotter.figure()
#     ##     plot_solution(Plotter, grid, data, solution, "bias \$\\Delta u\$ = $(endVoltage), \$ t=$(tend)\$", label_solution)
#     ##     Plotter.figure()
#     ##     plot_IV(Plotter, biasValues,IV, biasValues[end], plotGridpoints = true)
#     ## end

#     ## println("Max error")
#     ## @show max(abs.(solution_stationary_bias[iphin,:].-solution[iphin,:])...)
#     ## @show max(abs.(solution_stationary_bias[iphip,:].-solution[iphip,:])...)
#     ## @show max(abs.(solution_stationary_bias[iphit,:].-solution[iphit,:])...)
#     ## @show max(abs.(solution_stationary_bias[ipsi,:].-solution[ipsi,:])...)

#     testval = solution[data.index_psi, 10]
#     return testval

#     println("*** done\n")

# end #  main

# function test()
#     testval = 1.3214196490674017
#     main(test = true, unknown_storage=:dense) ≈ testval && main(test = true, unknown_storage=:sparse) ≈ testval
# end

# println("This message should show when this module has successfully recompiled.")


# end # module
