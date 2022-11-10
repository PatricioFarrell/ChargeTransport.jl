

module Ex402_PSC_InterfaceSpecies

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot
# using DelimitedFiles

function main(;n = 19, plotting = false, verbose = false, test = false,
              interfaceSpecies = true, leftInterface = true, interfaceReco = false)

    PyPlot.close("all")
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

    ## grid (the nearer to interface, the finer)
    h_pdoping               = 3.00e-6 * cm + 1.0e-7 * cm
    h_intrinsic             = 3.00e-5 * cm
    h_ndoping               = 8.50e-6 * cm + 1.0e-7 * cm

    δ                       = 4*n        # the larger, the finer the mesh
    t                       = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k                       = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_p_u               = collect(range(0.0, 2/3 * h_pdoping, step=h_pdoping/(0.3*δ)))
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

    iphin                = 1 # electron quasi Fermi potential
    iphip                = 2 # hole quasi Fermi potential

    if interfaceSpecies
        iphinb           = 3
        iphipb           = 4
        ipsi             = 5
        numberOfCarriers = 4
    else
        ipsi             = 3
        numberOfCarriers = 2
    end

    ## temperature
    T                    =  300.0                *  K

    ## band edge energies
    Ec_a                 = -3.0                  *  eV
    Ev_a                 = -5.1                  *  eV

    Ec_i                 = -3.8                  *  eV
    Ev_i                 = -5.4                  *  eV

    Ec_d                 = -3.8                  *  eV
    Ev_d                 = -6.2                  *  eV

    EC                   = [Ec_a, Ec_i, Ec_d]
    EV                   = [Ev_a, Ev_i, Ev_d]

    ## effective densities of state
    Nc_a                 = 1.0e20                / (cm^3)
    Nv_a                 = 1.0e20                / (cm^3)

    Nc_i                 = 1.0e19                / (cm^3)
    Nv_i                 = 1.0e19                / (cm^3)

    Nc_d                 = 1.0e19                / (cm^3)
    Nv_d                 = 1.0e19                / (cm^3)

    NC                   = [Nc_a, Nc_i, Nc_d]
    NV                   = [Nv_a, Nv_i, Nv_d]

    ## mobilities
    μn_a                 = 0.1                   * (cm^2) / (V * s)
    μp_a                 = 0.1                   * (cm^2) / (V * s)

    μn_i                 = 2.00e1                * (cm^2) / (V * s)
    μp_i                 = 2.00e1                * (cm^2) / (V * s)

    μn_d                 = 1.0e-3                * (cm^2) / (V * s)
    μp_d                 = 1.0e-3                * (cm^2) / (V * s)

    μn                   = [μn_a, μn_i, μn_d]
    μp                   = [μp_a, μp_i, μp_d]

    ## relative dielectric permittivity
    ε_a                  = 4.0                   *  1.0
    ε_i                  = 23.0                  *  1.0
    ε_d                  = 3.0                   *  1.0

    ε                    = [ε_a, ε_i, ε_d]

    ## radiative recombination
    r0_a                 = 6.3e-11               * cm^3 / s
    r0_i                 = 3.6e-12               * cm^3 / s
    r0_d                 = 6.8e-11               * cm^3 / s

    r0                   = [r0_a, r0_i, r0_d]

    ## life times and trap densities
    τn_a                 = 1.0e-6                * s
    τp_a                 = 1.0e-6                * s

    τn_i                 = 1.0e-7                * s
    τp_i                 = 1.0e-7                * s
    τn_d                 = τn_a
    τp_d                 = τp_a

    τn                   = [τn_a, τn_i, τn_d]
    τp                   = [τp_a, τp_i, τp_d]

    ## SRH trap energies
    Ei_a                 = -4.05                 * eV
    Ei_i                 = -4.60                 * eV
    Ei_d                 = -5.00                 * eV

    EI                   = [Ei_a, Ei_i, Ei_d]

    ## doping
    Nd                   = 2.089649130192123e17  / (cm^3)
    Na                   = 4.529587947185444e18  / (cm^3)

    ## contact voltages
    voltageAcceptor      =  1.2                  * V

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
    data.modelType                      = Stationary
    data.F                             .= Boltzmann
    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                  bulk_recomb_Auger = false,
                                                                  bulk_recomb_radiative = true,
                                                                  bulk_recomb_SRH = true)
    data.boundaryType[bregionAcceptor]  = OhmicContact
    data.boundaryType[bregionDonor]     = OhmicContact
    data.fluxApproximation             .= ScharfetterGummel

    if leftInterface == true
        bregActive = bregionJunction1
        bregDeact  = bregionJunction2
        icoordJ    = icoord_p

        regl       = regionAcceptor
        regr       = regionIntrinsic
    else
        bregActive = bregionJunction2
        bregDeact  = bregionJunction1
        icoordJ    = icoord_pi

        regl       = regionIntrinsic
        regr       = regionDonor
    end

    if interfaceSpecies
        enable_interface_carrier!(data, bulkCarrier = iphin, interfaceCarrier = iphinb, bregions = [bregActive])
        enable_interface_carrier!(data, bulkCarrier = iphip, interfaceCarrier = iphipb, bregions = [bregActive])
    end

    if interfaceReco

        data.boundaryType[bregActive] = InterfaceRecombination
        if interfaceSpecies

            data.interfaceRecombination = set_interface_recombination(;iphin = iphinb, iphip = iphipb,
                                                                       bregions = [bregActive])

        end
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
    params.chargeNumbers[iphin]                         = -1
    params.chargeNumbers[iphip]                         =  1

    for ireg in 1:numberOfRegions ## interior region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

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
    if interfaceSpecies

        data.d                                      = 6.28 * 10e-8 * cm # lattice size perovskite

        if leftInterface
            dopingN                                 = 0.0
            dopingP                                 = data.d * Na
        else
            dopingN                                 = data.d * Nd
            dopingP                                 = 0.0
        end

        δn                                          =    0.2 * eV
        δp                                          =  - 0.1 * eV

        params.bDoping[iphinb, bregActive]          = dopingN
        params.bDoping[iphipb, bregActive]          = dopingP

        params.bDensityOfStates[iphinb, bregActive] = data.d * params.densityOfStates[iphin, regl]
        params.bDensityOfStates[iphipb, bregActive] = data.d * params.densityOfStates[iphip, regl]

        params.bBandEdgeEnergy[iphinb, bregActive]  = params.bandEdgeEnergy[iphin, regl] + δn
        params.bBandEdgeEnergy[iphipb, bregActive]  = params.bandEdgeEnergy[iphip, regl] + δp

        zn  = params.chargeNumbers[iphin]
        zp  = params.chargeNumbers[iphip]
        UT  = params.UT

        Nc_l  = params.densityOfStates[iphin, regl]
        Nv_l  = params.densityOfStates[iphip, regl]
        Ec_l  = params.bandEdgeEnergy[iphin,  regl]
        Ev_l  = params.bandEdgeEnergy[iphip,  regl]

        Nc_r  = params.densityOfStates[iphin, regr]
        Nv_r  = params.densityOfStates[iphip, regr]
        Ec_r  = params.bandEdgeEnergy[iphin,  regr]
        Ev_r  = params.bandEdgeEnergy[iphip,  regr]

        # take values from intrinsic layer
        mun = params.mobility[iphin, regionIntrinsic]
        mup = params.mobility[iphip, regionIntrinsic]

        Nc_b = params.bDensityOfStates[iphinb, bregActive]
        Nv_b = params.bDensityOfStates[iphipb, bregActive]
        Ec_b = params.bBandEdgeEnergy[iphinb,  bregActive]
        Ev_b = params.bBandEdgeEnergy[iphipb,  bregActive]

        k0n = UT * mun * Nc_r/data.d
        k0p = UT * mup * Nv_r/data.d

        params.bReactionCoefficient[iphin, bregActive] = k0n
        params.bReactionCoefficient[iphip, bregActive] = k0p

        # If you decrease these values you will observe discontinuity in the qFs.
        params.bReactionCoefficient[iphin, bregDeact]        = 1.0e13
        params.bReactionCoefficient[iphip, bregDeact]        = 1.0e13

    end

    if interfaceReco

        # Caution! In case of interface species we have iphinb, iphipb and in bulk case iphin, iphip
        if interfaceSpecies

            params.bRecombinationSRHLifetime[iphinb, bregActive]    = 1/(1.0e1  * cm / s)
            params.bRecombinationSRHLifetime[iphipb, bregActive]    = 1/(1.0e5  * cm / s)

            params.bRecombinationSRHTrapDensity[iphinb, bregActive] = data.d * params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
            params.bRecombinationSRHTrapDensity[iphipb, bregActive] = data.d * params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

            params.bRecombinationSRHLifetime[iphinb, bregDeact]     = 1.0e100  * s#1.0e7  * cm / s
            params.bRecombinationSRHLifetime[iphipb, bregDeact]     = 1.0e100  * s#1.0e1  * cm / s

            params.bRecombinationSRHTrapDensity[iphinb, bregDeact]  = data.d * params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
            params.bRecombinationSRHTrapDensity[iphipb, bregDeact]  = data.d * params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

        else
            ## inner boundary region data (we choose the intrinsic values)
            params.bDensityOfStates[iphin, bregActive]             = Nc_i
            params.bDensityOfStates[iphip, bregActive]             = Nv_i

            params.bBandEdgeEnergy[iphin, bregActive]              = Ec_i
            params.bBandEdgeEnergy[iphip, bregActive]              = Ev_i

            params.recombinationSRHvelocity[iphin, bregActive]     = 1.0e1  * cm / s
            params.recombinationSRHvelocity[iphip, bregActive]     = 1.0e5  * cm / s

            params.bRecombinationSRHTrapDensity[iphin, bregActive] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
            params.bRecombinationSRHTrapDensity[iphip, bregActive] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

            # put these values small since the inverse enters
            params.recombinationSRHvelocity[iphin, bregDeact]      = 1.0e100  * cm / s#1.0e7  * cm / s
            params.recombinationSRHvelocity[iphip, bregDeact]      = 1.0e100  * cm / s#1.0e1  * cm / s

            params.bRecombinationSRHTrapDensity[iphin, bregDeact]  = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
            params.bRecombinationSRHTrapDensity[iphip, bregDeact]  = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]
        end

    end

    ##############################################################

    ## interior doping
    params.doping[iphin,  regionDonor]                  = Nd
    params.doping[iphip,  regionAcceptor]               = Na
    ## boundary doping
    params.bDoping[iphin, bregionDonor]                 = Nd
    params.bDoping[iphip, bregionAcceptor]              = Na

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=:sparse)

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
    control.tol_absolute      = 1.0e-8
    control.tol_relative      = 1.0e-8
    control.tol_round         = 1.0e-8
    control.damp_initial      = 0.9
    control.damp_growth       = 1.61 # >= 1

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

    if plotting

        function compute_densities(icc, ireg,  phi, psi)
            eta = data.params.chargeNumbers[icc]/ data.params.UT .* ( (phi .- psi) .+ data.params.bandEdgeEnergy[icc, ireg]./q )

            return (data.params.densityOfStates[icc, ireg] .* data.F[icc].(eta))
        end

        sol_ref  = readdlm("data/PSC-stationary-reference-sol-EQ.dat") # [coord sol_iphin sol_iphip sol_ipsi]
        vis1     = GridVisualizer(Plotter = PyPlot, layout=(2,1), size = (600,670), xlabel = "space [m]", ylabel = "potential [V]", fignumber=2)
        vis2     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "space [m]", ylabel = "density [m\$^{-3}\$]", fignumber=3)

        if interfaceSpecies
            subgridB = subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
            phin_sol = views(sol, data.chargeCarrierList[iphin], subgridB, ctsys.fvmsys)
            phip_sol = views(sol, data.chargeCarrierList[iphip], subgridB, ctsys.fvmsys)
            psi_sol  = views(sol, data.index_psi, subgridB, ctsys.fvmsys)

            # this is unfortunately not working soooo good ...
            bgrid    = subgrids(data.chargeCarrierList[iphinb], ctsys.fvmsys)

            phinb_sol = views(sol, data.chargeCarrierList[iphinb], bgrid, ctsys.fvmsys)
            phipb_sol = views(sol, data.chargeCarrierList[iphipb], bgrid, ctsys.fvmsys)
        else
            subgridp  = subgrid(grid, [1])
            subgridi  = subgrid(grid, [2])
            subgridn  = subgrid(grid, [3])
            # actually, we do not need this cases, since for this set-up we have continuity
            # in the densities with no effects on interfaces, but still for demonstrational
            # purpose.
            psi_solp  = view(sol[ipsi, :],  subgridp)
            psi_soli  = view(sol[ipsi, :],  subgridi)
            psi_soln  = view(sol[ipsi, :],  subgridn)

            phin_solp = view(sol[iphin, :], subgridp)
            phin_soli = view(sol[iphin, :], subgridi)
            phin_soln = view(sol[iphin, :], subgridn)

            phip_solp = view(sol[iphip, :], subgridp)
            phip_soli = view(sol[iphip, :], subgridi)
            phip_soln = view(sol[iphip, :], subgridn)


        end
        ###############################################################################
        ##########                         Potentials                        ##########
        ###############################################################################
        if interfaceSpecies
            for i in eachindex(phin_sol)
                scalarplot!(    vis1[1, 1], subgridB[i], psi_sol[i], clear = false, color=:blue, linewidth = 5)
                if i == 3
                    scalarplot!(vis1[1, 1], subgridB[i], psi_sol[i], clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")
                end
            end
            scalarplot!(vis1[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best)

            ##########################
            for i in eachindex(phin_sol)
                scalarplot!(vis1[2, 1], subgridB[i], phin_sol[i], clear = false, color=:green, linestyle=:solid, linewidth = 5)
                scalarplot!(vis1[2, 1], subgridB[i], phip_sol[i], clear = false, color=:red)

                if i == 3
                    scalarplot!(vis1[2, 1], subgridB[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                    scalarplot!(vis1[2, 1], subgridB[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$", color=:red)
                end
            end

            # DA: current way out, when waiting for changes within ExtendableGrids and GridVisualize
            PyPlot.figure(2)
            PyPlot.plot(coord[icoordJ], phinb_sol, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{\\varphi}_n \$")
            PyPlot.plot(coord[icoordJ], phipb_sol, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{\\varphi}_p \$")
            println("value phin_b = ", phinb_sol)
            println("value phip_b = ", phipb_sol)

            scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
            scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best, show = true)
        else
            scalarplot!(vis1[1, 1], subgridp, psi_solp,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis1[1, 1], subgridi, psi_soli,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis1[1, 1], subgridn, psi_soln,  clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")
            scalarplot!(vis1[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best)
            scalarplot!(vis1[2, 1], subgridp, phin_solp, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis1[2, 1], subgridp, phip_solp, clear = false, color=:red)
            scalarplot!(vis1[2, 1], subgridi, phin_soli, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis1[2, 1], subgridi, phip_soli, clear = false, color=:red)
            scalarplot!(vis1[2, 1], subgridn, phin_soln, clear = false, label = "\$ \\varphi_n \$", color=:green)
            scalarplot!(vis1[2, 1], subgridn, phip_soln, clear = false, label = "\$ \\varphi_p \$", color=:red)
            scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
            scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best, show = true)

        end

        ###############################################################################
        ##########                         Densities                         ##########
        ###############################################################################

        PyPlot.figure(3)
        if interfaceSpecies
            for i in eachindex(phin_sol)
                scalarplot!(vis2, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, linewidth = 5, yscale=:log)
                scalarplot!(vis2, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, color=:red)
                if i == 3
                    scalarplot!(vis2, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
                    scalarplot!(vis2, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, label ="\$ n_p \$", color=:red)
                end

            end

            if leftInterface
                psib         = psi_sol[1][end]
            else
                psib         = psi_sol[2][end]
            end

            eta_nb = -1/ data.params.UT * ( (phinb_sol[1] - psib) + data.params.bBandEdgeEnergy[iphinb, bregActive]/q )
            eta_pb =  1/ data.params.UT * ( (phipb_sol[1] - psib) + data.params.bBandEdgeEnergy[iphipb, bregActive]/q )

            nb     = data.params.bDensityOfStates[iphinb, bregActive] * data.F[iphin](eta_nb)./data.d
            pb     = data.params.bDensityOfStates[iphipb, bregActive] * data.F[iphip](eta_pb)./data.d

            println("value n_b (scaled by thickness) = ", nb)
            println("value p_b (scaled by thickness) = ", pb)

            PyPlot.semilogy(coord[icoordJ], nb, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{n}_n \$")
            PyPlot.semilogy(coord[icoordJ], pb, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{n}_p \$")

            n1 = compute_densities(iphin, regionAcceptor, sol_ref[1:icoord_p, 2], sol_ref[1:icoord_p, 4])
            p1 = compute_densities(iphip, regionAcceptor, sol_ref[1:icoord_p, 3], sol_ref[1:icoord_p, 4])

            n2 = compute_densities(iphin, regionIntrinsic, sol_ref[icoord_p:icoord_pi, 2], sol_ref[icoord_p:icoord_pi, 4])
            p2 = compute_densities(iphip, regionIntrinsic, sol_ref[icoord_p:icoord_pi, 3], sol_ref[icoord_p:icoord_pi, 4])

            n3 = compute_densities(iphin, regionDonor, sol_ref[icoord_pi:end, 2], sol_ref[icoord_pi:end, 4])
            p3 = compute_densities(iphip, regionDonor, sol_ref[icoord_pi:end, 3], sol_ref[icoord_pi:end, 4])

            scalarplot!(vis2, subgridB[1], n1, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
            scalarplot!(vis2, subgridB[1], p1, clear = false)
            scalarplot!(vis2, subgridB[2], n2, clear = false)
            scalarplot!(vis2, subgridB[2], p2, clear = false)
            scalarplot!(vis2, subgridB[3], n3, clear = false)
            scalarplot!(vis2, subgridB[3], p3, clear = false, label = "ref sol", legend =:best, show = true)

        else
            scalarplot!(vis2, subgridp, compute_densities(iphin, regionAcceptor,  phin_solp, psi_solp), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis2, subgridp, compute_densities(iphip, regionAcceptor,  phip_solp, psi_solp), clear = false, color=:red)
            scalarplot!(vis2, subgridi, compute_densities(iphin, regionIntrinsic, phin_soli, psi_soli), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis2, subgridi, compute_densities(iphip, regionIntrinsic, phip_soli, psi_soli), clear = false, color=:red)
            scalarplot!(vis2, subgridn, compute_densities(iphin, regionDonor,     phin_soln, psi_soln), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
            scalarplot!(vis2, subgridn, compute_densities(iphip, regionDonor,     phip_soln, psi_soln), clear = false, label ="\$ n_p \$", color=:red)

            n1 = compute_densities(iphin, regionAcceptor, sol_ref[1:icoord_p, 2], sol_ref[1:icoord_p, 4])
            p1 = compute_densities(iphip, regionAcceptor, sol_ref[1:icoord_p, 3], sol_ref[1:icoord_p, 4])

            n2 = compute_densities(iphin, regionIntrinsic, sol_ref[icoord_p:icoord_pi, 2], sol_ref[icoord_p:icoord_pi, 4])
            p2 = compute_densities(iphip, regionIntrinsic, sol_ref[icoord_p:icoord_pi, 3], sol_ref[icoord_p:icoord_pi, 4])

            n3 = compute_densities(iphin, regionDonor, sol_ref[icoord_pi:end, 2], sol_ref[icoord_pi:end, 4])
            p3 = compute_densities(iphip, regionDonor, sol_ref[icoord_pi:end, 3], sol_ref[icoord_pi:end, 4])

            scalarplot!(vis2, subgridp, n1, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
            scalarplot!(vis2, subgridp, p1, clear = false)
            scalarplot!(vis2, subgridi, n2, clear = false)
            scalarplot!(vis2, subgridi, p2, clear = false)
            scalarplot!(vis2, subgridn, n3, clear = false)
            scalarplot!(vis2, subgridn, p3, clear = false, label = "ref sol", legend =:best, show = true)

        end

    end

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
    biasValues           = range(0, stop = voltageAcceptor, length = 51)
    IV                   = zeros(0)

    control.damp_initial      = 0.5
    control.damp_growth       = 1.21 # >= 1

    for Δu in biasValues

        if test == false
            println("Δu  = ", Δu )
        end

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        solve!(sol, inival, ctsys, control = control, tstep = Inf)
        inival .= sol

        ## get I-V data
        val = get_current_val(ctsys, sol)

        push!(IV, val)

        if interfaceSpecies
            subgridB = subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
            bgrid    = subgrids(data.chargeCarrierList[iphinb], ctsys.fvmsys)

            phin_sol = views(sol, data.chargeCarrierList[iphin], subgridB, ctsys.fvmsys)
            phip_sol = views(sol, data.chargeCarrierList[iphip], subgridB, ctsys.fvmsys)
            psi_sol  = views(sol, data.index_psi,                subgridB, ctsys.fvmsys)

            phinb_sol = views(sol, data.chargeCarrierList[iphinb], bgrid, ctsys.fvmsys)
            phipb_sol = views(sol, data.chargeCarrierList[iphipb], bgrid, ctsys.fvmsys)

            # left values
            etanl = zn / UT * ( (phin_sol[regl][end] - psi_sol[regl][end]) + Ec_l / q ) # left
            etapl = zp / UT * ( (phip_sol[regr][end] - psi_sol[regr][end]) + Ev_l / q )  # left

            nl    = Nc_l * data.F[iphin](etanl)
            pl    = Nv_l * data.F[iphip](etapl)

            # interface values
            etan_b = zn / UT * ( (phinb_sol[1] - psi_sol[regl][end]) + Ec_b / q ) # interface
            etap_b = zp / UT * ( (phipb_sol[1] - psi_sol[regl][end]) + Ev_b / q ) # interface

            n_b    = Nc_b * data.F[iphin](etan_b)
            p_b    = Nv_b * data.F[iphip](etap_b)

            # right values
            etanr = zn / UT * ( (phin_sol[regr][1] - psi_sol[regr][1]) + Ec_r / q ) # left
            etapr = zp / UT * ( (phip_sol[regr][1] - psi_sol[regr][1]) + Ev_r / q )  # left

            nr    = Nc_r * data.F[iphin](etanr)
            pr    = Nv_r * data.F[iphip](etapr)

            Knleft  = exp(- zn/ (kB * data.params.temperature) * (Ec_l - Ec_b))
            Knright = exp(- zn/ (kB * data.params.temperature) * (Ec_r - Ec_b))
            Kpleft  = exp(- zp/ (kB * data.params.temperature) * (Ev_l - Ev_b))
            Kpright = exp(- zp/ (kB * data.params.temperature) * (Ev_r - Ev_b))

            #@show Knleft, Knright, Kpleft, Kpright

            reactnl  = - q * zn * (Knleft^(1/2)  * nl/Nc_l - Knleft^( - 1/2) * n_b/Nc_b)
            reactnr  = - q * zn * (Knright^(1/2) * nr/Nc_r - Knright^(- 1/2) * n_b/Nc_b)

            reactpl  = - q * zp * (Kpleft^(1/2)  * pl/Nv_l - Kpleft^( - 1/2) * p_b/Nv_b)
            reactpr  = - q * zp * (Kpright^(1/2) * pr/Nv_r - Kpright^(- 1/2) * p_b/Nv_b)

            # #println("main file: ")
            # @show reactnl, reactnr
            # @show reactpl, reactpr
        end

    end # bias loop

    # writedlm("data/PSC-stationary-reference-sol-surface-reco.dat", [coord sol'])
    # res = [biasValues IV]
    # writedlm("data/PSC-stationary-reference-IV-surface-reco.dat", res)

    # return

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Some plotting")
    end
    ################################################################################

    if plotting

        if interfaceReco
            sol_ref  = readdlm("data/PSC-stationary-reference-sol-surface-reco.dat") # [coord sol_iphin sol_iphip sol_ipsi]
            IV_ref   = readdlm("data/PSC-stationary-reference-IV-surface-reco.dat")
        else
            sol_ref  = readdlm("data/PSC-stationary-reference-sol.dat") # [coord sol_iphin sol_iphip sol_ipsi]
            IV_ref   = readdlm("data/PSC-stationary-reference-IV.dat")
        end

        vis3     = GridVisualizer(Plotter = PyPlot, layout=(2,1), size = (600,670), xlabel = "space [m]", ylabel = "potential [V]", fignumber=4)
        vis4     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "space [m]", ylabel = "density [m\$^{-3}\$]", fignumber=5)
        vis5     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "voltage [V]", ylabel = "current density [Am\$^{-2}\$]",fignumber=6)

        if interfaceSpecies
            subgridB = subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
            phin_sol = views(sol, data.chargeCarrierList[iphin], subgridB, ctsys.fvmsys)
            phip_sol = views(sol, data.chargeCarrierList[iphip], subgridB, ctsys.fvmsys)
            psi_sol  = views(sol, data.index_psi, subgridB, ctsys.fvmsys)

            # this is unfortunately not working soooo good ...
            bgrid    = subgrids(data.chargeCarrierList[iphinb], ctsys.fvmsys)

            phinb_sol = views(sol, data.chargeCarrierList[iphinb], bgrid, ctsys.fvmsys)
            phipb_sol = views(sol, data.chargeCarrierList[iphipb], bgrid, ctsys.fvmsys)
        else
            subgridp  = subgrid(grid, [1])
            subgridi  = subgrid(grid, [2])
            subgridn  = subgrid(grid, [3])
            # actually, we do not need this cases, since for this set-up we have continuity
            # in the densities with no effects on interfaces, but still for demonstrational
            # purpose.
            psi_solp  = view(sol[ipsi, :],  subgridp)
            psi_soli  = view(sol[ipsi, :],  subgridi)
            psi_soln  = view(sol[ipsi, :],  subgridn)

            phin_solp = view(sol[iphin, :], subgridp)
            phin_soli = view(sol[iphin, :], subgridi)
            phin_soln = view(sol[iphin, :], subgridn)

            phip_solp = view(sol[iphip, :], subgridp)
            phip_soli = view(sol[iphip, :], subgridi)
            phip_soln = view(sol[iphip, :], subgridn)


        end

        ###############################################################################
        ##########                         Potentials                        ##########
        ###############################################################################
        if interfaceSpecies
            for i in eachindex(phin_sol)
                scalarplot!(    vis3[1, 1], subgridB[i], psi_sol[i], clear = false, color=:blue, linewidth = 5)
                if i == 3
                    scalarplot!(vis3[1, 1], subgridB[i], psi_sol[i], clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")
                end
            end
            scalarplot!(vis3[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best)

            ##########################
            for i in eachindex(phin_sol)
                scalarplot!(vis3[2, 1], subgridB[i], phin_sol[i], clear = false, color=:green, linestyle=:solid, linewidth = 5)
                scalarplot!(vis3[2, 1], subgridB[i], phip_sol[i], clear = false, color=:red)

                if i == 3
                    scalarplot!(vis3[2, 1], subgridB[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                    scalarplot!(vis3[2, 1], subgridB[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$", color=:red)
                end
            end

            # DA: current way out, when waiting for changes within ExtendableGrids and GridVisualize
            PyPlot.figure(4)
            PyPlot.plot(coord[icoordJ], phinb_sol, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{\\varphi}_n \$")
            PyPlot.plot(coord[icoordJ], phipb_sol, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{\\varphi}_p \$")
            println("value phin_b = ", phinb_sol)
            println("value phip_b = ", phipb_sol)

            scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
            scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best, show = true)
        else
            scalarplot!(vis3[1, 1], subgridp, psi_solp,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis3[1, 1], subgridi, psi_soli,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis3[1, 1], subgridn, psi_soln,  clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")
            scalarplot!(vis3[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best)
            scalarplot!(vis3[2, 1], subgridp, phin_solp, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis3[2, 1], subgridp, phip_solp, clear = false, color=:red)
            scalarplot!(vis3[2, 1], subgridi, phin_soli, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis3[2, 1], subgridi, phip_soli, clear = false, color=:red)
            scalarplot!(vis3[2, 1], subgridn, phin_soln, clear = false, label = "\$ \\varphi_n \$", color=:green)
            scalarplot!(vis3[2, 1], subgridn, phip_soln, clear = false, label = "\$ \\varphi_p \$", color=:red)
            scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
            scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best, show = true)

        end

        ###############################################################################
        ##########                         Densities                         ##########
        ###############################################################################

        if interfaceSpecies
            for i in eachindex(phin_sol)
                scalarplot!(vis4, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, linewidth = 5, yscale=:log)
                scalarplot!(vis4, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, color=:red)
                if i == 3
                    scalarplot!(vis4, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
                    scalarplot!(vis4, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, label ="\$ n_p \$", color=:red)
                end

            end

            PyPlot.figure(5)


            if leftInterface
                psib         = psi_sol[1][end]
            else
                psib         = psi_sol[2][end]
            end

            eta_nb = -1/ data.params.UT * ( (phinb_sol[1] - psib) + data.params.bBandEdgeEnergy[iphinb, bregActive]/q )
            eta_pb =  1/ data.params.UT * ( (phipb_sol[1] - psib) + data.params.bBandEdgeEnergy[iphipb, bregActive]/q )

            # DA: divide by d such that it is three dimensional again?
            nb     = data.params.bDensityOfStates[iphinb, bregActive] * data.F[iphin](eta_nb)./data.d
            pb     = data.params.bDensityOfStates[iphipb, bregActive] * data.F[iphip](eta_pb)./data.d

            println("value n_b (scaled by thickness) = ", nb)
            println("value p_b (scaled by thickness) = ", pb)

            PyPlot.semilogy(coord[icoordJ], nb, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{n}_n \$")
            PyPlot.semilogy(coord[icoordJ], pb, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{n}_p \$")

            n1 = compute_densities(iphin, regionAcceptor, sol_ref[1:icoord_p, 2], sol_ref[1:icoord_p, 4])
            p1 = compute_densities(iphip, regionAcceptor, sol_ref[1:icoord_p, 3], sol_ref[1:icoord_p, 4])

            n2 = compute_densities(iphin, regionIntrinsic, sol_ref[icoord_p:icoord_pi, 2], sol_ref[icoord_p:icoord_pi, 4])
            p2 = compute_densities(iphip, regionIntrinsic, sol_ref[icoord_p:icoord_pi, 3], sol_ref[icoord_p:icoord_pi, 4])

            n3 = compute_densities(iphin, regionDonor, sol_ref[icoord_pi:end, 2], sol_ref[icoord_pi:end, 4])
            p3 = compute_densities(iphip, regionDonor, sol_ref[icoord_pi:end, 3], sol_ref[icoord_pi:end, 4])

            scalarplot!(vis4, subgridB[1], n1, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
            scalarplot!(vis4, subgridB[1], p1, clear = false)
            scalarplot!(vis4, subgridB[2], n2, clear = false)
            scalarplot!(vis4, subgridB[2], p2, clear = false)
            scalarplot!(vis4, subgridB[3], n3, clear = false)
            scalarplot!(vis4, subgridB[3], p3, clear = false, label = "ref sol", legend =:best, show = true)
        else
            scalarplot!(vis4, subgridp, compute_densities(iphin, regionAcceptor,  phin_solp, psi_solp), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis4, subgridp, compute_densities(iphip, regionAcceptor,  phip_solp, psi_solp), clear = false, color=:red)
            scalarplot!(vis4, subgridi, compute_densities(iphin, regionIntrinsic, phin_soli, psi_soli), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis4, subgridi, compute_densities(iphip, regionIntrinsic, phip_soli, psi_soli), clear = false, color=:red)
            scalarplot!(vis4, subgridn, compute_densities(iphin, regionDonor,     phin_soln, psi_soln), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
            scalarplot!(vis4, subgridn, compute_densities(iphip, regionDonor,     phip_soln, psi_soln), clear = false, label ="\$ n_p \$", color=:red)

            n1 = compute_densities(iphin, regionAcceptor, sol_ref[1:icoord_p, 2], sol_ref[1:icoord_p, 4])
            p1 = compute_densities(iphip, regionAcceptor, sol_ref[1:icoord_p, 3], sol_ref[1:icoord_p, 4])

            n2 = compute_densities(iphin, regionIntrinsic, sol_ref[icoord_p:icoord_pi, 2], sol_ref[icoord_p:icoord_pi, 4])
            p2 = compute_densities(iphip, regionIntrinsic, sol_ref[icoord_p:icoord_pi, 3], sol_ref[icoord_p:icoord_pi, 4])

            n3 = compute_densities(iphin, regionDonor, sol_ref[icoord_pi:end, 2], sol_ref[icoord_pi:end, 4])
            p3 = compute_densities(iphip, regionDonor, sol_ref[icoord_pi:end, 3], sol_ref[icoord_pi:end, 4])

            scalarplot!(vis4, subgridp, n1, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
            scalarplot!(vis4, subgridp, p1, clear = false)
            scalarplot!(vis4, subgridi, n2, clear = false)
            scalarplot!(vis4, subgridi, p2, clear = false)
            scalarplot!(vis4, subgridn, n3, clear = false)
            scalarplot!(vis4, subgridn, p3, clear = false, label = "ref sol", legend =:best, show = true)

        end

        ###############################################################################
        ##########                            IV                             ##########
        ###############################################################################
        scalarplot!(vis5, biasValues,   IV,           clear = false, color=:green, linewidth = 5)
        scalarplot!(vis5, IV_ref[:, 1], IV_ref[:, 2], clear = false, color=:black, linestyle=:dot)

    end

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, sol))/length(sol) # when using sparse storage, we get NaN values in solution
    return testval


end #  main



function test()
    testvalwithoutReco = -0.2988993117821689; testvalwithReco = -0.2988993117821689
    main(test = true, interfaceSpecies = true, leftInterface = true, interfaceReco=false) ≈ testvalwithoutReco && main(test = true, interfaceSpecies = true, leftInterface = true, interfaceReco=false) ≈ testvalwithReco

    # main(test = true, interfaceSpecies = false, leftInterface = true, interfaceReco=false) = -0.905305257878682
    # main(test = true, interfaceSpecies = false, leftInterface = true, interfaceReco=true)  = -0.892757917525781

end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
