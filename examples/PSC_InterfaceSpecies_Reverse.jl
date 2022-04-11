
module PSC_InterfaceSpecies_Reverse

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot
using DelimitedFiles

function main(;n = 19, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:dense)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionAcceptor          = 1                           # p doped region
    regionIntrinsic         = 2                           # intrinsic region
    regionDonor             = 3                           # n doped region
    regions                 = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions         = length(regions)

    ## boundary region numbers
    bregionAcceptor         = 1
    bregionDonor            = 2
    bregionJunction1        = 3
    bregionJunction2        = 4
    bregions                = [bregionAcceptor, bregionDonor, bregionJunction1, bregionJunction2]
    numberOfBoundaryRegions = length(bregions)

    ## grid (the nearer to interface, the finer)
    h_pdoping               = 3.00e-6 * cm + 1.0e-7 * cm
    h_intrinsic             = 3.00e-5 * cm
    h_ndoping               = 8.50e-6 * cm + 1.0e-7 * cm

    x0                      = 0.0 * cm
    δ                       = 4*n        # the larger, the finer the mesh
    t                       = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k                       = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_p_u               = collect(range(x0, 2/3 * h_pdoping, step=h_pdoping/(0.3*δ)))
    coord_p_g               = geomspace(2/3 * h_pdoping,
                                        h_pdoping,
                                        h_pdoping/(0.4*δ),
                                        h_pdoping/(1.2*δ),
                                        tol=t)
    coord_i_g1              = geomspace(h_pdoping,
                                        h_pdoping+h_intrinsic/k,
                                        h_intrinsic/(7.1*δ),
                                        h_intrinsic/(0.4*δ),
                                        tol=t)
    coord_i_g2              = geomspace(h_pdoping+h_intrinsic/k,
                                        h_pdoping+h_intrinsic,
                                        h_intrinsic/(0.4*δ),
                                        h_intrinsic/(7.1*δ),
                                        tol=t)
    coord_n_g               = geomspace(h_pdoping+h_intrinsic,
                                        h_pdoping+h_intrinsic+1/3 * h_ndoping,
                                        h_ndoping/(2.0*δ),
                                        h_ndoping/(0.4*δ),
                                        tol=t)
    coord_n_u               = collect(range(h_pdoping+h_intrinsic+1/3 * h_ndoping, h_pdoping+h_intrinsic+h_ndoping, step=h_pdoping/(0.3*δ)))

    coord                   = glue(coord_p_u,coord_p_g,  tol=10*t)
    icoord_p                = length(coord)
    coord                   = glue(coord,    coord_i_g1, tol=10*t)
    coord                   = glue(coord,    coord_i_g2, tol=10*t)
    icoord_pi               = length(coord)
    coord                   = glue(coord,    coord_n_g,  tol=10*t)
    coord                   = glue(coord,    coord_n_u,  tol=10*t)
    grid                    = simplexgrid(coord)

    ## cellmask! for defining the subregions and assigning region number
    cellmask!(grid, [0.0 * μm],                 [h_pdoping],                           regionAcceptor, tol = 1.0e-18)   # p-doped region   = 1
    cellmask!(grid, [h_pdoping],                [h_pdoping + h_intrinsic],             regionIntrinsic, tol = 1.0e-18)  # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic],  [h_pdoping + h_intrinsic + h_ndoping], regionDonor, tol = 1.0e-18)      # n-doped region   = 3

    bfacemask!(grid, [h_pdoping],               [h_pdoping],                           bregionJunction1)  # first  inner interface
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic],             bregionJunction2)  # second inner interface

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    iphin               = 1 # electron quasi Fermi potential
    iphip               = 2 # hole quasi Fermi potential
    iphin_b2            = 3
    iphip_b2            = 4
    numberOfCarriers    = 4

    ## temperature
    T                   =  300.0                *  K

    ## band edge energies
    Ec_a                = -3.0                  *  eV
    Ev_a                = -5.1                  *  eV

    Ec_i                = -3.8                  *  eV
    Ev_i                = -5.4                  *  eV

    Ec_d                = -3.8                  *  eV
    Ev_d                = -6.2                  *  eV

    EC                  = [Ec_a, Ec_i, Ec_d]
    EV                  = [Ev_a, Ev_i, Ev_d]

    ## effective densities of state
    Nc_a                = 1.0e20                / (cm^3)
    Nv_a                = 1.0e20                / (cm^3)

    Nc_i                = 1.0e19                / (cm^3)
    Nv_i                = 1.0e19                / (cm^3)

    Nc_d                = 1.0e19                / (cm^3)
    Nv_d                = 1.0e19                / (cm^3)

    NC                  = [Nc_a, Nc_i, Nc_d]
    NV                  = [Nv_a, Nv_i, Nv_d]

    ## mobilities
    μn_a                = 0.1                   * (cm^2) / (V * s)
    μp_a                = 0.1                   * (cm^2) / (V * s)

    μn_i                = 2.00e1                * (cm^2) / (V * s)
    μp_i                = 2.00e1                * (cm^2) / (V * s)

    μn_d                = 1.0e-3                * (cm^2) / (V * s)
    μp_d                = 1.0e-3                * (cm^2) / (V * s)

    μn                  = [μn_a, μn_i, μn_d]
    μp                  = [μp_a, μp_i, μp_d]

    ## relative dielectric permittivity
    ε_a                 = 4.0                   *  1.0
    ε_i                 = 23.0                  *  1.0
    ε_d                 = 3.0                   *  1.0

    ε                   = [ε_a, ε_i, ε_d]

    ## radiative recombination
    r0_a                = 6.3e-11               * cm^3 / s
    r0_i                = 3.6e-12               * cm^3 / s
    r0_d                = 6.8e-11               * cm^3 / s

    r0                  = [r0_a, r0_i, r0_d]

    ## life times and trap densities
    τn_a                = 1.0e-6                * s
    τp_a                = 1.0e-6                * s

    τn_i                = 1.0e-7                * s
    τp_i                = 1.0e-7                * s
    τn_d                = τn_a
    τp_d                = τp_a

    τn                  = [τn_a, τn_i, τn_d]
    τp                  = [τp_a, τp_i, τp_d]

    ## SRH trap energies
    Ei_a                = -4.05                 * eV
    Ei_i                = -4.60                 * eV
    Ei_d                = -5.00                 * eV

    EI                  = [Ei_a, Ei_i, Ei_d]

    ## doping
    Nd                  = 2.089649130192123e17  / (cm^3)
    Na                  = 4.529587947185444e18  / (cm^3)

    ## contact voltages
    voltageAcceptor     =  1.2                  * V

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

    ## possible choices: Stationary, Transient
    data.modelType                      = Stationary

    ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                              = [Boltzmann, Boltzmann]

    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                  bulk_recomb_Auger = false,
                                                                  bulk_recomb_radiative = true,
                                                                  bulk_recomb_SRH = true)

    data.boundaryType[bregionAcceptor]  = OhmicContact
    data.boundaryType[bregionJunction1] = InterfaceModelDiscontqFNoReaction
    data.boundaryType[bregionJunction2] = InterfaceModelDiscontqF
    data.boundaryType[bregionDonor]     = OhmicContact

    # wäre schöner, wenn pro iphin_b2 nur iphin
    enable_interface_carriers!(data, bulkSpecies = [iphin, iphip], interfaceSpecies = [iphin_b2, iphip_b2], boundaryRegion = bregionJunction2)

    data.isContinuous[iphin]            = false
    data.isContinuous[iphip]            = false

    data.fluxApproximation              = ScharfetterGummel

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

    for ireg in 1:numberOfRegions ## interior region data

        params.dielectricConstant[ireg]                 = ε[ireg]

        ## effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = NC[ireg]
        params.densityOfStates[iphip, ireg]             = NV[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = EC[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = EV[ireg]

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, data, EI[ireg])
    end

    ## outer boundary region data
    params.bDensityOfStates[iphin, bregionAcceptor]     = Nc_a
    params.bDensityOfStates[iphip, bregionAcceptor]     = Nv_a

    params.bDensityOfStates[iphin, bregionDonor]        = Nc_d
    params.bDensityOfStates[iphip, bregionDonor]        = Nv_d

    params.bBandEdgeEnergy[iphin, bregionAcceptor]      = Ec_a
    params.bBandEdgeEnergy[iphip, bregionAcceptor]      = Ev_a

    params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_d
    params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_d

    # ##############################################################

    ## inner boundary region data
    # delta with negative sign -> IV shifted to right.
    delta1                                              = +0.3 * eV
    delta2                                              = -0.0 * eV
    data.d                                              = 6.28 * 10e-8 * cm # lattice size perovskite
    params.bDensityOfStates[iphin_b2, bregionJunction2] = data.d * params.densityOfStates[iphin, regionDonor]
    params.bDensityOfStates[iphip_b2, bregionJunction2] = data.d * params.densityOfStates[iphip, regionDonor]

    params.bBandEdgeEnergy[iphin_b2, bregionJunction2]  = params.bandEdgeEnergy[iphin, regionDonor] + delta1
    params.bBandEdgeEnergy[iphip_b2, bregionJunction2]  = params.bandEdgeEnergy[iphip, regionDonor] + delta2

    ##############################################################

    ## interior doping
    params.doping[iphin,  regionDonor]                  = Nd
    params.doping[iphip,  regionAcceptor]               = Na
    ## boundary doping
    params.bDoping[iphin, bregionDonor]                 = Nd
    params.bDoping[iphip, bregionAcceptor]              = Na

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=unknown_storage)

    # set legend for plotting routines. Either you can use the predefined labels or write your own.
    label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
    plot_energies(Plotter, grid, data, label_BEE)
    PyPlot.ylim(-7.0, -2.0)
 
    #savefig("rev-Ec-M0p3-Ev-M0p3.eps")

    function compute_densities(icc, ireg,  phi, psi)
        eta = data.params.chargeNumbers[icc]/ data.params.UT .* ( (phi .- psi) .+ data.params.bandEdgeEnergy[icc, ireg]./q )

        return (data.params.densityOfStates[icc, ireg] .* data.F[icc].(eta))
    end


    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define outer boundary conditions")
    end
    ################################################################################

    ## We set zero voltage ohmic contacts for each charge carrier at all outerior boundaries
    ## for the equilibrium calculations.
    set_contact!(ctsys, bregionAcceptor, Δu = 0.0)
    set_contact!(ctsys, bregionDonor, Δu = 0.0)

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
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    ## initialize solution and starting vectors
    inival  = unknowns(ctsys)
    sol     = unknowns(ctsys)

    sol     = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    inival .= sol

    #writedlm("data/PSC-stationary-reference-sol-EQ.dat", [coord sol'])

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    # Set calculation type to outOfEquilibrium for starting with respective simulation.
    data.calculationType = OutOfEquilibrium

    maxBias    = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 51)
    IV         = zeros(0)

    for Δu in biasValues

        println("Δu  = ", Δu )

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        solve!(sol, inival, ctsys, control = control, tstep = Inf)

        inival .= sol

        factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
        tf      = testfunction(factory, [1], [2])
        I       = integrate(ctsys.fvmsys, tf, sol)

        val = 0.0
        for ii = 1:length(I)-1
            val = val + I[ii]
        end

        push!(IV, abs(val) )

    end # bias loop

    # writedlm("data/PSC-stationary-reference-sol.dat", [coord sol'])
    # res = [biasValues IV]
    # writedlm("data/PSC-stationary-reference-IV.dat", res)

    if test == false
        println("*** done\n")
    end

    if plotting == true

        vis = GridVisualizer(Plotter = PyPlot, layout=(2,1))

        subgrids    = VoronoiFVM.subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
        bgrid2      = subgrid(grid, [bregionJunction2], boundary=true)

        phin_sol    = VoronoiFVM.views(sol, data.chargeCarrierList[iphin], subgrids, ctsys.fvmsys)
        phip_sol    = VoronoiFVM.views(sol, data.chargeCarrierList[iphip], subgrids, ctsys.fvmsys)
        psi_sol     = VoronoiFVM.views(sol, data.index_psi, subgrids, ctsys.fvmsys)
        
        phin_b2_sol = view(sol[iphin_b2, :], bgrid2)
        phip_b2_sol = view(sol[iphip_b2, :], bgrid2)

        sol_ref     = readdlm("data/PSC-stationary-reference-sol.dat")

        for i = 1:length(phin_sol)
            scalarplot!(vis[1, 1], subgrids[i], phin_sol[i], clear = false, color=:green)
            scalarplot!(vis[1, 1], subgrids[i], phip_sol[i], clear = false, color=:red)
            if i == 3
                scalarplot!(vis[1, 1], subgrids[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                scalarplot!(vis[1, 1], subgrids[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$",  color=:red)
            end
        end
        PyPlot.plot(sol_ref[:, 1], sol_ref[:, 2], linestyle="--", color = "black")
        PyPlot.plot(sol_ref[:, 1], sol_ref[:, 3], linestyle="--", color = "black")
        PyPlot.plot(coord[icoord_pi], phin_b2_sol, marker = "x", markersize = 12, color = "darkgreen", label = "\$ \\bar{\\varphi}_n \$")
        PyPlot.plot(coord[icoord_pi], phip_b2_sol, marker = "x", markersize = 12, color = "darkred", label = "\$ \\bar{\\varphi}_p \$")
        println("value phin_b2 = ", phin_b2_sol)
        println("value phip_b2 = ", phip_b2_sol)

        Plotter.xlabel("space [m]")
        Plotter.ylabel("potential [V]")
        Plotter.legend(loc = "best", fontsize=11)
        Plotter.title("Solution with Bias")
        PyPlot.tight_layout()

        # for i = 1:length(phin_sol)
        #     scalarplot!(vis[2, 1], subgrids[i], psi_sol[i],  clear = false, color=:blue)
        #     if i == 3
        #         scalarplot!(vis[2, 1], subgrids[i], psi_sol[i],  clear = false, label = "\$ \\psi \$",color=:blue)
        #     end
        # end
        # PyPlot.plot(sol_ref[:, 1], sol_ref[:, 4], linestyle="--", color = "black")

        # Plotter.xlabel("space [m]")
        # Plotter.ylabel("potential [V]")
        # Plotter.legend(fancybox = true, loc = "best", fontsize=11)
        # Plotter.title("Solution with Bias")
        # PyPlot.tight_layout()
        xcoord = [1:icoord_p, icoord_p:icoord_pi, icoord_pi:length(coord)]

        for i = 1:length(phin_sol)
            scalarplot!(vis[2, 1], subgrids[i], compute_densities(iphin, subgrids[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, linestyle=:solid, yscale=:log)
            scalarplot!(vis[2, 1], subgrids[i], compute_densities(iphip, subgrids[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, color=:red)
            if i == 3
            scalarplot!(vis[2, 1], subgrids[i], compute_densities(iphin, subgrids[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, label ="n", yscale=:log)
            scalarplot!(vis[2, 1], subgrids[i], compute_densities(iphip, subgrids[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, label ="p", color=:red)
            end

        end

        for i = 1:length(phin_sol)
            scalarplot!(vis[2, 1], subgrids[i], compute_densities(iphin, subgrids[i][CellRegions][1], sol_ref[xcoord[i], 2], sol_ref[xcoord[i], 4]), clear = false, label=:none, linestyle=:dot, color=:black)
            scalarplot!(vis[2, 1], subgrids[i], compute_densities(iphip, subgrids[i][CellRegions][1], sol_ref[xcoord[i], 3], sol_ref[xcoord[i], 4]), clear = false, linestyle=:dot, color=:black)

        end

        eta_nb2 = -1/ data.params.UT * ( (phin_b2_sol[1] - psi_sol[2][end]) + data.params.bBandEdgeEnergy[iphin_b2, bregionJunction2]/q )
        eta_pb2 =  1/ data.params.UT * ( (phip_b2_sol[1] - psi_sol[2][end]) + data.params.bBandEdgeEnergy[iphip_b2, bregionJunction2]/q )

        nb2 = data.params.bDensityOfStates[iphin_b2, bregionJunction2] * data.F[iphin](eta_nb2)
        pb2 = data.params.bDensityOfStates[iphip_b2, bregionJunction2] * data.F[iphip](eta_pb2)

        println("value n_b2 = ", nb2)
        println("value p_b2 = ", pb2)

        PyPlot.semilogy(coord[icoord_pi], nb2, marker = "x", markersize = 12, color = "darkgreen", label = "\$ \\bar{n} \$")
        PyPlot.semilogy(coord[icoord_pi], pb2, marker = "x", markersize = 12, color = "darkred", label = "\$ \\bar{p} \$")

        Plotter.xlabel("space [\$m\$]")
        Plotter.ylabel("density [\$\\frac{1}{m^3}\$]")
        #Plotter.legend()
        Plotter.tight_layout()

        # ##########################################################
        # scalarplot!(vis[3, 1], biasValues, IV, clear = false, color=:green)
        # IV_ref         = readdlm("data/PSC-stationary-reference-IV.dat")
        # PyPlot.plot(IV_ref[:, 1], IV_ref[:, 2], linestyle="--", color = "black")

    end

    #savefig("rev-Ec-M0p3-Ev-M0p3-sol.eps")
    testval = VoronoiFVM.norm(ctsys.fvmsys, sol, 2)
    return testval


end #  main

function test()
    #testval = 49.92777921026983
    main(test = true, unknown_storage=:dense) ≈ testval #&& main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
