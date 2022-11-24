

module Ex402_PSC_IM_InterfaceSpecies

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot
using DelimitedFiles

function main(;n = 19, plotting = false, verbose = false, test = false,
              plotCTWithoutIntSpec  = false,
              ionicSpecies      = true,
              interfaceSpecies  = true,
              interfaceReco     = false,
              leftInterface     = true)

    PyPlot.close("all")
    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionDonor             = 1                           # n doped region
    regionIntrinsic         = 2                           # intrinsic region
    regionAcceptor          = 3                           # p doped region
    regions                 = [regionDonor, regionIntrinsic, regionAcceptor]
    numberOfRegions         = length(regions)

    ## boundary region numbers
    bregionDonor            = 1
    bregionAcceptor         = 2
    bregionJunction1        = 3
    bregionJunction2        = 4
    bregions                = [bregionDonor, bregionAcceptor, bregionJunction1, bregionJunction2]

    ## grid (the nearer to interface, the finer)
    h_ndoping               = 1.00e-5 * cm
    h_intrinsic             = 4.00e-5 * cm
    h_pdoping               = 2.00e-5 * cm

    δ                       = 4*n        # the larger, the finer the mesh
    t                       = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k                       = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_n_u               = collect(range(0.0, 2/3 * h_ndoping, step=h_ndoping/(0.6*δ)))
    coord_n_g               = geomspace(2/3 * h_ndoping,
                                        h_ndoping,
                                        h_ndoping/(0.4*δ),
                                        h_ndoping/(4.2*δ),
                                        tol=t)
    coord_i_g1              = geomspace(h_ndoping,
                                        h_ndoping+h_intrinsic/k,
                                        h_intrinsic/(9.1*δ),
                                        h_intrinsic/(0.4*δ),
                                        tol=t)
    coord_i_g2              = geomspace(h_ndoping+h_intrinsic/k,
                                        h_ndoping+h_intrinsic,
                                        h_intrinsic/(0.4*δ),
                                        h_intrinsic/(9.1*δ),
                                        tol=t)
    coord_p_g               = geomspace(h_ndoping+h_intrinsic,
                                        h_ndoping+h_intrinsic+1/3 * h_pdoping,
                                        h_pdoping/(7.0*δ),
                                        h_pdoping/(0.4*δ),
                                        tol=t)
    coord_p_u               = collect(range(h_ndoping+h_intrinsic+1/3 * h_pdoping, h_ndoping+h_intrinsic+h_pdoping, step=h_ndoping/(0.6*δ)))

    coord                   = glue(coord_n_u,coord_n_g,  tol=10*t)
    icoord_n                = length(coord)
    coord                   = glue(coord,    coord_i_g1, tol=10*t)
    coord                   = glue(coord,    coord_i_g2, tol=10*t)
    icoord_ni               = length(coord)
    coord                   = glue(coord,    coord_p_g,  tol=10*t)
    coord                   = glue(coord,    coord_p_u,  tol=10*t)
    grid                    = simplexgrid(coord)

    ## cellmask! for defining the subregions and assigning region number
    cellmask!(grid, [0.0 * μm],                 [h_ndoping],                           regionDonor, tol = 1.0e-18)
    cellmask!(grid, [h_ndoping],                [h_ndoping + h_intrinsic],             regionIntrinsic, tol = 1.0e-18)
    cellmask!(grid, [h_ndoping + h_intrinsic],  [h_ndoping + h_intrinsic + h_pdoping], regionAcceptor, tol = 1.0e-18)

    bfacemask!(grid, [h_ndoping],               [h_ndoping],                           bregionJunction1)  # first  inner interface
    bfacemask!(grid, [h_ndoping + h_intrinsic], [h_ndoping + h_intrinsic],             bregionJunction2)  # second inner interface

    if plotting
        gridplot(grid, Plotter = PyPlot, legend=:lt, fignumber = 1)
        PyPlot.title("Grid")
    end

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

    if interfaceSpecies && ionicSpecies
        iphia            = 3
        iphinb           = 4
        iphipb           = 5
        numberOfCarriers = 5

        ipsi             = 6

    elseif interfaceSpecies && ionicSpecies == false
        iphinb           = 3
        iphipb           = 4
        numberOfCarriers = 4

        ipsi             = 5

    elseif interfaceSpecies == false && ionicSpecies
        iphia            = 3
        numberOfCarriers = 3

        ipsi             = 4

    else
        numberOfCarriers = 2

        ipsi             = 3

    end

    ## temperature
    T                  = 298.0                 *  K

    ## band edge energies
    Ec_d               = -4.0                  *  eV
    Ev_d               = -5.8                  *  eV

    Ec_i               = -3.7                  *  eV
    Ev_i               = -5.4                  *  eV
    Ea_i               = -4.45                 *  eV

    Ec_a               = -3.4                  *  eV
    Ev_a               = -5.1                  *  eV

    EC                 = [Ec_d, Ec_i, Ec_a]
    EV                 = [Ev_d, Ev_i, Ev_a]
    EA                 = [0.0,  Ea_i,  0.0]

    ## effective densities of state
    Nc_d               = 5.0e19                / (cm^3)
    Nv_d               = 5.0e19                / (cm^3)

    Nc_i               = 8.1e18                / (cm^3)
    Nv_i               = 5.8e18                / (cm^3)
    Nanion             = 1.6e19/0.01           / (cm^3)

    Nc_a               = 5.0e19                / (cm^3)
    Nv_a               = 5.0e19                / (cm^3)

    NC                 = [Nc_d, Nc_i,  Nc_a]
    NV                 = [Nv_d, Nv_i,  Nv_a]
    NAnion             = [0.0,  Nanion, 0.0]

    ## diffusivities (from Ionmonger)
    UT                 = (kB * T) / q

    Dn_e               = 1.0e-5                * m^2 / s
    Dp_e               = 1.0e-5                * m^2 / s

    Dn_i               = 1.7e-4                * m^2 / s
    Dp_i               = 1.7e-4                * m^2 / s
    Da_i               = 6.5e-8 * m^2/s * exp(-0.58*eV/(kB * T))

    Dn_a               = 1.0e-6                * m^2 / s
    Dp_a               = 1.0e-6                * m^2 / s

    ## mobilities (get from Einstein relation)
    μn_d               = Dn_e / UT
    μp_d               = Dp_e / UT

    μn_i               = Dn_i / UT
    μp_i               = Dp_i / UT
    μa_i               = Da_i / UT

    μn_a               = Dn_a / UT
    μp_a               = Dp_a / UT

    μn                 = [μn_d, μn_i, μn_a]
    μp                 = [μp_d, μp_i, μp_a]
    μa                 = [0.0,  μa_i, 0.0 ]

    ## relative dielectric permittivity
    ε_d                = 10.0                  *  1.0
    ε_i                = 24.1                  *  1.0
    ε_a                = 3.0                   *  1.0

    ε                  = [ε_d, ε_i, ε_a]

    ## radiative recombination
    r0_d               = 0.0e0                 * cm^3 / s
    r0_i               = 0.0e0                 * cm^3 / s
    r0_a               = 0.0e0                 * cm^3 / s

    r0                 = [r0_d, r0_i, r0_a]

    ## life times and trap densities
    τn_d               = 1.0e100               * s
    τp_d               = 1.0e100               * s

    τn_i               = 3.0e-9                * s
    τp_i               = 3.0e-7                * s
    τn_a               = τn_d
    τp_a               = τp_d

    τn                 = [τn_d, τn_i, τn_a]
    τp                 = [τp_d, τp_i, τp_a]

    ## SRH trap energies
    Ei_d               = -5.0                  * eV
    Ei_i               = -4.55                 * eV
    Ei_a               = -4.1                  * eV

    EI                 = [Ei_d, Ei_i, Ei_a]


    ## doping
    Nd                 = 1.00e18             / (cm^3)
    Na                 = 1.00e18             / (cm^3)
    C0                 = 1.6e19              / (cm^3)

    ## contact voltage
    voltageAcceptor    =  1.2                  * V

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
    data.modelType                      = Transient
    data.F                             .= Boltzmann
    if ionicSpecies # in case of present anion vacancies, put FermiDiracMinusOne.
        data.F[iphia]                   = FermiDiracMinusOne
    end
    data.bulkRecombination              = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                  bulk_recomb_Auger = false,
                                                                  bulk_recomb_radiative = false,
                                                                  bulk_recomb_SRH = true)
    data.boundaryType[bregionAcceptor]  = OhmicContact
    data.boundaryType[bregionDonor]     = OhmicContact
    data.fluxApproximation             .= ExcessChemicalPotential

    # declare left and right interface
    if leftInterface == true
        bregActive = bregionJunction1
        bregDeact  = bregionJunction2
        icoordJ    = icoord_n

        regl       = regionDonor
        regr       = regionIntrinsic
    else
        bregActive = bregionJunction2
        bregDeact  = bregionJunction1
        icoordJ    = icoord_ni

        regl       = regionIntrinsic
        regr       = regionAcceptor
    end

    if ionicSpecies     # enable ionic carrier
        enable_ionic_carrier!(data = data, ionicCarrier = iphia, regions = [regionIntrinsic])
    end

    if interfaceSpecies # present interface species
        enable_interface_carrier!(data = data, bulkCarrier = iphin, interfaceCarrier = iphinb, bregions = [bregActive])
        enable_interface_carrier!(data = data, bulkCarrier = iphip, interfaceCarrier = iphipb, bregions = [bregActive])
    end

    if interfaceReco
        data.boundaryType[bregActive] = InterfaceRecombination
        if interfaceSpecies
            set_interface_recombination!(data = data, iphin = iphinb, iphip = iphipb, bregions = [bregActive])
        else
            data.boundaryType[bregDeact] = InterfaceRecombination
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
    if ionicSpecies
        params.chargeNumbers[iphia]                     =  1
    end
    if interfaceSpecies
        params.chargeNumbers[iphinb]                    = -1
        params.chargeNumbers[iphipb]                    =  1
    end

    for ireg in 1:numberOfRegions ## interior region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

        ## effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = NC[ireg]
        params.densityOfStates[iphip, ireg]             = NV[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = EC[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = EV[ireg]

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]

        if ionicSpecies
            params.densityOfStates[iphia, ireg]         = NAnion[ireg]
            params.bandEdgeEnergy[iphia, ireg]          = EA[ireg]
            params.mobility[iphia, ireg]                = μa[ireg]
        end

        ## recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = params.densityOfStates[iphin, ireg] * exp( params.chargeNumbers[iphin] * (params.bandEdgeEnergy[iphin, ireg] - EI[ireg]) / (kB * params.temperature))
        params.recombinationSRHTrapDensity[iphip, ireg] = params.densityOfStates[iphip, ireg] * exp( params.chargeNumbers[iphip] * (params.bandEdgeEnergy[iphip, ireg] - EI[ireg]) / (kB * params.temperature))
    end


    ## interior doping
    params.doping[iphin,  regionDonor]                  = Nd
    params.doping[iphip,  regionAcceptor]               = Na
    if ionicSpecies
        params.doping[iphia, regionIntrinsic]           = C0
    end

    ##############################################################
    ## outer boundary region data
    params.bDoping[iphin, bregionDonor]                 = Nd
    params.bDoping[iphip, bregionAcceptor]              = Na


    params.bDensityOfStates[iphin, bregionDonor]        = Nc_d
    params.bDensityOfStates[iphip, bregionDonor]        = Nv_d

    params.bDensityOfStates[iphin, bregionAcceptor]     = Nc_a
    params.bDensityOfStates[iphip, bregionAcceptor]     = Nv_a

    params.bBandEdgeEnergy[iphin, bregionDonor]         = Ec_d
    params.bBandEdgeEnergy[iphip, bregionDonor]         = Ev_d

    params.bBandEdgeEnergy[iphin, bregionAcceptor]      = Ec_a
    params.bBandEdgeEnergy[iphip, bregionAcceptor]      = Ev_a

    # ##############################################################

    ## inner boundary region data
    if interfaceSpecies

        data.d                                      = 6.28 * 10e-8 * cm # lattice size perovskite

        if ionicSpecies
            params.bBandEdgeEnergy[iphia, bregionJunction1]  = Ea_i
            params.bBandEdgeEnergy[iphia, bregionJunction2]  = Ea_i

            params.bDensityOfStates[iphia, bregionJunction1] = data.d * Nanion
            params.bDensityOfStates[iphia, bregionJunction2] = data.d * Nanion

            params.bDoping[iphia, bregionJunction1]           = data.d * C0
            params.bDoping[iphia, bregionJunction2]           = data.d * C0
        end

        if leftInterface
            dopingN                                 = data.d * Nd
            dopingP                                 = 0.0

        else
            dopingN                                 = 0.0
            dopingP                                 = data.d * Na
        end

        δn                                          =  0.05 * eV#0.1 * eV
        δp                                          =  0.15 * eV#0.1 * eV

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

            params.bRecombinationSRHLifetime[iphinb, bregActive]    = 1/(1.0e5  * m / s)
            params.bRecombinationSRHLifetime[iphipb, bregActive]    = 1/(1.0e1  * m / s)

            params.bRecombinationSRHTrapDensity[iphinb, bregActive] = data.d * params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
            params.bRecombinationSRHTrapDensity[iphipb, bregActive] = data.d * params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

            params.bRecombinationSRHLifetime[iphinb, bregDeact]     = 1.0e100  * s#1.0e7  * cm / s
            params.bRecombinationSRHLifetime[iphipb, bregDeact]     = 1.0e100  * s#1.0e1  * cm / s

            params.bRecombinationSRHTrapDensity[iphinb, bregDeact]  = data.d * params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
            params.bRecombinationSRHTrapDensity[iphipb, bregDeact]  = data.d * params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

        else
            ## inner boundary region data (we choose the intrinsic values)
            params.bDensityOfStates[iphin, bregionJunction1]             = Nc_i
            params.bDensityOfStates[iphip, bregionJunction1]             = Nv_i

            params.bBandEdgeEnergy[iphin, bregionJunction1]              = Ec_i
            params.bBandEdgeEnergy[iphip, bregionJunction1]              = Ev_i

            params.recombinationSRHvelocity[iphin, bregionJunction1]     = 1.0e5  * m / s
            params.recombinationSRHvelocity[iphip, bregionJunction1]     = 1.0e1  * m / s

            params.bRecombinationSRHTrapDensity[iphin, bregionJunction1] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
            params.bRecombinationSRHTrapDensity[iphip, bregionJunction1] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

            ###################################################
            params.bDensityOfStates[iphin, bregionJunction2]             = Nc_a
            params.bDensityOfStates[iphip, bregionJunction2]             = Nv_i

            params.bBandEdgeEnergy[iphin, bregionJunction2]              = Ec_a
            params.bBandEdgeEnergy[iphip, bregionJunction2]              = Ev_i

            params.recombinationSRHvelocity[iphin, bregionJunction2]     = 1.0e-1 * m / s
            params.recombinationSRHvelocity[iphip, bregionJunction2]     = 1.0e5  * m / s

            params.bRecombinationSRHTrapDensity[iphin, bregionJunction2] = params.recombinationSRHTrapDensity[iphin, regionAcceptor]
            params.bRecombinationSRHTrapDensity[iphip, bregionJunction2] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]
        end

    end

    ##############################################################

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=:sparse)

    #return data.electricCarrierList, data.ionicCarrierList, data.interfaceCarrierList
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

        vis1     = GridVisualizer(Plotter = PyPlot, layout=(2,1), size = (600,670), xlabel = "space [m]", ylabel = "potential [V]", fignumber=2)
        vis2     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "space [m]", ylabel = "density [m\$^{-3}\$]", fignumber=3)

        if interfaceSpecies
            subgridB     = subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
            phin_sol     = views(sol, data.chargeCarrierList[iphin], subgridB, ctsys.fvmsys)
            phip_sol     = views(sol, data.chargeCarrierList[iphip], subgridB, ctsys.fvmsys)
            psi_sol      = views(sol, data.index_psi, subgridB, ctsys.fvmsys)
            if ionicSpecies
                phia_sol = views(sol, data.chargeCarrierList[iphia], subgridB, ctsys.fvmsys)
            end

            # this is unfortunately not working soooo good ...
            bgrid    = subgrids(data.chargeCarrierList[iphinb], ctsys.fvmsys)

            phinb_sol = views(sol, data.chargeCarrierList[iphinb], bgrid, ctsys.fvmsys)
            phipb_sol = views(sol, data.chargeCarrierList[iphipb], bgrid, ctsys.fvmsys)

        else
            subgridn  = subgrid(grid, [1])
            subgridi  = subgrid(grid, [2])
            subgridp  = subgrid(grid, [3])
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

            if ionicSpecies
                phia_soli = view(sol[iphia, :], subgridi)
            end

        end

        if plotCTWithoutIntSpec
            sol_ref  = readdlm("data/PSC-stationary-reference-sol-EQ.dat") # [coord sol_iphin sol_iphip sol_ipsi]
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

            if plotCTWithoutIntSpec
                scalarplot!(vis1[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol", legend =:best)
            end

            ##########################
            for i in eachindex(phin_sol)
                scalarplot!(vis1[2, 1], subgridB[i], phin_sol[i], clear = false, color=:green, linestyle=:solid, linewidth = 5)
                scalarplot!(vis1[2, 1], subgridB[i], phip_sol[i], clear = false, color=:red)

                if i == 3
                    scalarplot!(vis1[2, 1], subgridB[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                    scalarplot!(vis1[2, 1], subgridB[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$", color=:red)
                end
            end

            if ionicSpecies
                scalarplot!(vis1[2, 1], subgridB[regionIntrinsic], phia_sol[regionIntrinsic], clear = false, label = "\$ \\varphi_a \$", color=:gold)
            end

            # DA: current way out, when waiting for changes within ExtendableGrids and GridVisualize
            PyPlot.figure(2)
            PyPlot.plot(coord[icoordJ], phinb_sol, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{\\varphi}_n \$")
            PyPlot.plot(coord[icoordJ], phipb_sol, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{\\varphi}_p \$")
            println("value phin_b = ", phinb_sol)
            println("value phip_b = ", phipb_sol)

            if plotCTWithoutIntSpec
                scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
                scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                            legend =:best, show = true)
            end

        else
            scalarplot!(vis1[1, 1], subgridp, psi_solp,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis1[1, 1], subgridi, psi_soli,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis1[1, 1], subgridn, psi_soln,  clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")

            if plotCTWithoutIntSpec
                scalarplot!(vis1[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                            legend =:best)
            end

            scalarplot!(vis1[2, 1], subgridp, phin_solp, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis1[2, 1], subgridp, phip_solp, clear = false, color=:red)
            scalarplot!(vis1[2, 1], subgridi, phin_soli, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis1[2, 1], subgridi, phip_soli, clear = false, color=:red)
            if ionicSpecies
                scalarplot!(vis1[2, 1], subgridi, phia_soli, clear = false, label = "\$ \\varphi_a \$", color=:gold)
            end

            scalarplot!(vis1[2, 1], subgridn, phin_soln, clear = false, label = "\$ \\varphi_n \$", color=:green)
            scalarplot!(vis1[2, 1], subgridn, phip_soln, clear = false, label = "\$ \\varphi_p \$", color=:red,
            legend =:best, show = true)

            if plotCTWithoutIntSpec
                scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
                scalarplot!(vis1[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol")
            end

        end

        ##### Ionmonger solution

        coord_IM_ETL  = vec(1.0e-9 .* (readdlm("IM/IM-parameter-grid-ETL.dat")  .- readdlm("IM/IM-parameter-grid-ETL.dat")[1]))
        coord_IM_intr = vec(1.0e-9 .* (readdlm("IM/IM-parameter-grid-intr.dat") .- readdlm("IM/IM-parameter-grid-ETL.dat")[1]))
        coord_IM_HTL  = vec(1.0e-9 .* (readdlm("IM/IM-parameter-grid-HTL.dat")  .- readdlm("IM/IM-parameter-grid-ETL.dat")[1]))
        coordIM       = glue(coord_IM_ETL[:, 1], coord_IM_intr[:, 1])
        coordIM       = glue(coordIM[:, 1], coord_IM_HTL[:, 1])

        grid_ETL      = ExtendableGrids.simplexgrid(coord_IM_ETL)
        grid_intr     = ExtendableGrids.simplexgrid(coord_IM_intr)
        grid_HTL      = ExtendableGrids.simplexgrid(coord_IM_HTL)

        Vbi           = (Ec_d - Ev_a - kB * T * log(Nc_d * Nv_a/ (Nd * Na)))/q
        Nintr_d       = sqrt(Nc_d * Nv_d) * exp(-(Ec_d - Ev_d) / (2 * kB * T))
        psi0_d        = (Ec_d + Ev_d)/(2 * q) - 0.5 * UT * log(Nc_d/Nv_d) + UT * asinh(Nd/(2*Nintr_d)) # -4.100459358855549

        psi_ETL_imIR  = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-psi-ETL-t-0p0.dat")     .-Vbi/2 .+ psi0_d)
        psi_intr_imIR = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-psi-intr-t-0p0.dat")    .-Vbi/2 .+ psi0_d)
        psi_HTL_imIR  = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-psi-HTL-t-0p0.dat")     .-Vbi/2 .+ psi0_d)

        a_imIR        = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-a-intr-t-0p0.dat"))
        p_imIR        = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-p-intr-t-0p0.dat"))
        n_imIR        = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-n-intr-t-0p0.dat"))
        ##########################
        psi_ETL_im    = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-psi-ETL-t-0p0.dat")  .-Vbi/2 .+ psi0_d)
        psi_intr_im   = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-psi-intr-t-0p0.dat") .-Vbi/2 .+ psi0_d)
        psi_HTL_im    = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-psi-HTL-t-0p0.dat")  .-Vbi/2 .+ psi0_d)

        a_im          = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-a-intr-t-0p0.dat"))
        p_im          = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-p-intr-t-0p0.dat"))
        n_im          = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-n-intr-t-0p0.dat"))

        scalarplot!(vis1[1, 1], grid_ETL,  psi_ETL_im,  clear = false, marker ="", label = "", linestyle=:dash, color=:black, linewidth = 4)
        scalarplot!(vis1[1, 1], grid_intr, psi_intr_im, clear = false)
        scalarplot!(vis1[1, 1], grid_HTL,  psi_HTL_im,  clear = false, label = "IM (without reco)")

        scalarplot!(vis1[1, 1], grid_ETL,  psi_ETL_imIR,  clear = false, marker ="", label = "", linestyle=:dot, color=:gray, linewidth = 5)
        scalarplot!(vis1[1, 1], grid_intr, psi_intr_imIR, clear = false)
        scalarplot!(vis1[1, 1], grid_HTL,  psi_HTL_imIR,  clear = false, label = "IM (with reco)", legend =:best, show = true)

        ###############################################################################
        ##########                         Densities                         ##########
        ###############################################################################

        PyPlot.figure(3)
        if interfaceSpecies

            if ionicSpecies
                scalarplot!(vis2, subgridB[regionIntrinsic], compute_densities(iphia, subgridB[regionIntrinsic][CellRegions][1], phia_sol[regionIntrinsic], psi_sol[regionIntrinsic]), clear = false, color=:gold, label ="\$ n_a \$", linewidth = 5, yscale=:log)
            end

            for i in eachindex(phin_sol)
                scalarplot!(vis2, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), label = "", clear = false, color=:green, linewidth = 5, yscale=:log)
                scalarplot!(vis2, subgridB[i], compute_densities(iphip, subgridB[i][CellRegions][1], phip_sol[i], psi_sol[i]), clear = false, color=:red)
                if i == 3
                   scalarplot!(vis2, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), clear = false, color=:green, label ="\$ n_n \$", yscale=:log, show = true)
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

            if plotCTWithoutIntSpec

                PyPlot.semilogy(coord[icoordJ], nb, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{n}_n \$")
                PyPlot.semilogy(coord[icoordJ], pb, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{n}_p \$")

                n1 = compute_densities(iphin, regionDonor, sol_ref[1:icoord_n, 2], sol_ref[1:icoord_n, 4])
                p1 = compute_densities(iphip, regionDonor, sol_ref[1:icoord_n, 3], sol_ref[1:icoord_n, 4])

                n2 = compute_densities(iphin, regionIntrinsic, sol_ref[icoord_n:icoord_ni, 2], sol_ref[icoord_n:icoord_ni, 4])
                p2 = compute_densities(iphip, regionIntrinsic, sol_ref[icoord_n:icoord_ni, 3], sol_ref[icoord_n:icoord_ni, 4])

                n3 = compute_densities(iphin, regionAcceptor, sol_ref[icoord_ni:end, 2], sol_ref[icoord_ni:end, 4])
                p3 = compute_densities(iphip, regionAcceptor, sol_ref[icoord_ni:end, 3], sol_ref[icoord_ni:end, 4])

                scalarplot!(vis2, subgridB[1], n1, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
                scalarplot!(vis2, subgridB[1], p1, clear = false)
                scalarplot!(vis2, subgridB[2], n2, clear = false)
                scalarplot!(vis2, subgridB[2], p2, clear = false)
                scalarplot!(vis2, subgridB[3], n3, clear = false)
                scalarplot!(vis2, subgridB[3], p3, clear = false, label = "ref sol", legend =:best, show = true)
            end

        else


            scalarplot!(vis2, subgridp, compute_densities(iphin, regionAcceptor,  phin_solp, psi_solp), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis2, subgridp, compute_densities(iphip, regionAcceptor,  phip_solp, psi_solp), clear = false, color=:red)
            scalarplot!(vis2, subgridi, compute_densities(iphin, regionIntrinsic, phin_soli, psi_soli), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis2, subgridi, compute_densities(iphip, regionIntrinsic, phip_soli, psi_soli), clear = false, color=:red)
            scalarplot!(vis2, subgridn, compute_densities(iphin, regionDonor,     phin_soln, psi_soln), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
            if ionicSpecies
                scalarplot!(vis2, subgridi, compute_densities(iphia, regionIntrinsic, phia_soli, psi_soli), clear = false, color=:gold, label ="\$ n_a \$", linewidth = 5, yscale=:log)
            end
            scalarplot!(vis2, subgridn, compute_densities(iphip, regionDonor,     phip_soln, psi_soln), clear = false, label ="\$ n_p \$", color=:red)


            if plotCTWithoutIntSpec
                n1 = compute_densities(iphin, regionDonor, sol_ref[1:icoord_n, 2], sol_ref[1:icoord_n, 4])
                p1 = compute_densities(iphip, regionDonor, sol_ref[1:icoord_n, 3], sol_ref[1:icoord_n, 4])

                n2 = compute_densities(iphin, regionIntrinsic, sol_ref[icoord_n:icoord_ni, 2], sol_ref[icoord_n:icoord_ni, 4])
                p2 = compute_densities(iphip, regionIntrinsic, sol_ref[icoord_n:icoord_ni, 3], sol_ref[icoord_n:icoord_ni, 4])

                n3 = compute_densities(iphin, regionAcceptor, sol_ref[icoord_ni:end, 2], sol_ref[icoord_ni:end, 4])
                p3 = compute_densities(iphip, regionAcceptor, sol_ref[icoord_ni:end, 3], sol_ref[icoord_ni:end, 4])

                scalarplot!(vis2, subgridn, n1, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
                scalarplot!(vis2, subgridn, p1, clear = false)
                scalarplot!(vis2, subgridi, n2, clear = false)
                scalarplot!(vis2, subgridi, p2, clear = false)
                scalarplot!(vis2, subgridp, n3, clear = false)
                scalarplot!(vis2, subgridp, p3, clear = false, label = "ref sol")

            end

        end

        ##### Ionmonger solution
        scalarplot!(vis2, grid_intr, n_im , clear = false, marker ="", label = "", linestyle=:dash, color=:black)
        scalarplot!(vis2, grid_intr, p_im,  clear = false)
        scalarplot!(vis2, grid_intr, a_im,  clear = false, label = "IM (without reco)")

        scalarplot!(vis2, grid_intr, n_imIR , clear = false, marker ="", label = "", linestyle=:dot, color=:gray)
        scalarplot!(vis2, grid_intr, p_imIR,  clear = false)
        scalarplot!(vis2, grid_intr, a_imIR,  clear = false, label = "IM (with reco)", legend =:best, show = true)

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
    IV                   = zeros(0)
    biasValues           = zeros(0)

    # control.damp_initial = 0.5
    # control.damp_growth  = 1.21 # >= 1

    ## primary data for I-V scan protocol
    scanrate             = 0.1 * V/s
    number_tsteps        = 101
    endVoltage           = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend                 = endVoltage/scanrate
    tvalues              = range(0, stop = tend, length = number_tsteps)

    for istep = 2:number_tsteps

        t  = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep-1] # Time step size

        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: t = $(t),", " Δu  = ", Δu)
        end

        solve!(sol, inival, ctsys, control  = control, tstep = Δt)

        ## get I-V data
        current = get_current_val(ctsys, sol, inival, Δt)

        push!(IV, current)
        push!(biasValues, Δu)

        inival .= sol

    end # time loop

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

        if plotCTWithoutIntSpec
            if interfaceReco
                sol_ref  = readdlm("data/PSC-stationary-reference-sol-surface-reco.dat") # [coord sol_iphin sol_iphip sol_ipsi]
                IV_ref   = readdlm("data/PSC-stationary-reference-IV-surface-reco.dat")
            else
                sol_ref  = readdlm("data/PSC-stationary-reference-sol.dat") # [coord sol_iphin sol_iphip sol_ipsi]
                IV_ref   = readdlm("data/PSC-stationary-reference-IV.dat")
            end
        end

        vis3     = GridVisualizer(Plotter = PyPlot, layout=(2,1), size = (600,670), xlabel = "space [m]", ylabel = "potential [V]", fignumber=4)
        vis4     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "space [m]", ylabel = "density [m\$^{-3}\$]", fignumber=5)
        vis5     = GridVisualizer(Plotter = PyPlot, layout=(1,1), xlabel = "voltage [V]", ylabel = "current density [Am\$^{-2}\$]",fignumber=6)

        if interfaceSpecies

            phin_sol = views(sol, data.chargeCarrierList[iphin], subgridB, ctsys.fvmsys)
            phip_sol = views(sol, data.chargeCarrierList[iphip], subgridB, ctsys.fvmsys)
            psi_sol  = views(sol, data.index_psi, subgridB, ctsys.fvmsys)


            phinb_sol = views(sol, data.chargeCarrierList[iphinb], bgrid, ctsys.fvmsys)
            phipb_sol = views(sol, data.chargeCarrierList[iphipb], bgrid, ctsys.fvmsys)

            if ionicSpecies
                phia_sol = views(sol, data.chargeCarrierList[iphia], subgridB, ctsys.fvmsys)
            end
        else

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

            if ionicSpecies
                phia_soli = view(sol[iphia, :], subgridi)
            end

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

            if plotCTWithoutIntSpec
                scalarplot!(vis3[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best)
            end

            ##########################
            for i in eachindex(phin_sol)
                scalarplot!(vis3[2, 1], subgridB[i], phin_sol[i], clear = false, color=:green, linestyle=:solid, linewidth = 5)
                scalarplot!(vis3[2, 1], subgridB[i], phip_sol[i], clear = false, color=:red)

                if i == 3
                    scalarplot!(vis3[2, 1], subgridB[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                    scalarplot!(vis3[2, 1], subgridB[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$", color=:red)
                end
            end

            if ionicSpecies
                scalarplot!(vis3[2, 1], subgridB[regionIntrinsic], phia_sol[regionIntrinsic], clear = false, label = "\$ \\varphi_a \$", color=:gold)
            end

            # DA: current way out, when waiting for changes within ExtendableGrids and GridVisualize
            PyPlot.figure(4)
            PyPlot.plot(coord[icoordJ], phinb_sol, marker = "o", markersize = 12,  color =:darkgreen, label = "\$ \\bar{\\varphi}_n \$")
            PyPlot.plot(coord[icoordJ], phipb_sol, marker = "o", markersize = 12,  color =:darkred, label = "\$ \\bar{\\varphi}_p \$")
            println("value phin_b = ", phinb_sol)
            println("value phip_b = ", phipb_sol)

            if plotCTWithoutIntSpec
                scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
                scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best, show = true)
            end

        else
            scalarplot!(vis3[1, 1], subgridp, psi_solp,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis3[1, 1], subgridi, psi_soli,  clear = false, color=:blue, linewidth = 5)
            scalarplot!(vis3[1, 1], subgridn, psi_soln,  clear = false, color=:blue, linewidth = 5, label = "\$ \\psi \$")

            if plotCTWithoutIntSpec
                scalarplot!(vis3[1, 1], sol_ref[:, 1], sol_ref[:, 4], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                        legend =:best)
            end

            scalarplot!(vis3[2, 1], subgridp, phin_solp, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis3[2, 1], subgridp, phip_solp, clear = false, color=:red)
            scalarplot!(vis3[2, 1], subgridi, phin_soli, clear = false, color=:green, linestyle=:solid, linewidth = 5)
            scalarplot!(vis3[2, 1], subgridi, phip_soli, clear = false, color=:red)
            if ionicSpecies
                scalarplot!(vis3[2, 1], subgridi, phia_soli, clear = false, label = "\$ \\varphi_a \$", color=:gold)
            end
            scalarplot!(vis3[2, 1], subgridn, phin_soln, clear = false, label = "\$ \\varphi_n \$", color=:green)
            scalarplot!(vis3[2, 1], subgridn, phip_soln, clear = false, label = "\$ \\varphi_p \$", color=:red,legend =:best, show = true)

            if plotCTWithoutIntSpec
                scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 3], marker ="", clear = false, color=:black, linestyle=:dot, label = "")
                scalarplot!(vis3[2, 1], sol_ref[:, 1], sol_ref[:, 2], clear = false, color =:black, linestyle=:dot, label = "ref sol",
                            legend =:best, show = true)
            end

        end

        ## Ionmonger solution!

        psi_ETL_imIR  = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-psi-ETL-t-end.dat")     .-Vbi/2 .+ psi0_d .+ (1.2)/2)
        psi_intr_imIR = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-psi-intr-t-end.dat")    .-Vbi/2 .+ psi0_d .+ (1.2)/2)
        psi_HTL_imIR  = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-psi-HTL-t-end.dat")     .-Vbi/2 .+ psi0_d .+ (1.2)/2)

        a_imIR        = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-a-intr-t-end.dat"))
        p_imIR        = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-p-intr-t-end.dat"))
        n_imIR        = vec(readdlm("IM/with-surface-reco/0p1/IM-parameter-n-intr-t-end.dat"))
        ###################
        psi_ETL_im    = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-psi-ETL-t-end.dat")  .-Vbi/2 .+ psi0_d .+ (1.2)/2)
        psi_intr_im   = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-psi-intr-t-end.dat") .-Vbi/2 .+ psi0_d .+ (1.2)/2)
        psi_HTL_im    = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-psi-HTL-t-end.dat")  .-Vbi/2 .+ psi0_d .+ (1.2)/2)

        a_im          = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-a-intr-t-end.dat"))
        p_im          = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-p-intr-t-end.dat"))
        n_im          = vec(readdlm("IM/without-surface-reco/0p1/IM-parameter-n-intr-t-end.dat"))

        scalarplot!(vis3[1, 1], grid_ETL,  psi_ETL_im,  clear = false, marker ="", label ="", linestyle=:dash, color=:black, linewidth = 4)
        scalarplot!(vis3[1, 1], grid_intr, psi_intr_im, clear = false)
        scalarplot!(vis3[1, 1], grid_HTL,  psi_HTL_im,  clear = false, label = "IM (without reco)")

        scalarplot!(vis3[1, 1], grid_ETL,  psi_ETL_im,  clear = false, marker ="", label ="", linestyle=:dot, color=:gray, linewidth = 5)
        scalarplot!(vis3[1, 1], grid_intr, psi_intr_im, clear = false)
        scalarplot!(vis3[1, 1], grid_HTL,  psi_HTL_im,  clear = false, label = "IM (with reco)", legend =:best, show = true)

        ###############################################################################
        ##########                         Densities                         ##########
        ###############################################################################

        if interfaceSpecies

            if ionicSpecies
                scalarplot!(vis4, subgridB[regionIntrinsic], compute_densities(iphia, subgridB[regionIntrinsic][CellRegions][1], phia_sol[regionIntrinsic], psi_sol[regionIntrinsic]), clear = false, color=:gold, label ="\$ n_a \$", linewidth = 5, yscale=:log)
            end

            for i in eachindex(phin_sol)
                scalarplot!(vis4, subgridB[i], compute_densities(iphin, subgridB[i][CellRegions][1], phin_sol[i], psi_sol[i]), label ="", clear = false, color=:green, linewidth = 5, yscale=:log)
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

            if plotCTWithoutIntSpec
                n1 = compute_densities(iphin, regionDonor, sol_ref[1:icoord_n, 2], sol_ref[1:icoord_n, 4])
                p1 = compute_densities(iphip, regionDonor, sol_ref[1:icoord_n, 3], sol_ref[1:icoord_n, 4])

                n2 = compute_densities(iphin, regionIntrinsic, sol_ref[icoord_n:icoord_ni, 2], sol_ref[icoord_n:icoord_ni, 4])
                p2 = compute_densities(iphip, regionIntrinsic, sol_ref[icoord_n:icoord_ni, 3], sol_ref[icoord_n:icoord_ni, 4])

                n3 = compute_densities(iphin, regionAcceptor, sol_ref[icoord_ni:end, 2], sol_ref[icoord_ni:end, 4])
                p3 = compute_densities(iphip, regionAcceptor, sol_ref[icoord_ni:end, 3], sol_ref[icoord_ni:end, 4])

                scalarplot!(vis4, subgridB[1], n1, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
                scalarplot!(vis4, subgridB[1], p1, clear = false)
                scalarplot!(vis4, subgridB[2], n2, clear = false)
                scalarplot!(vis4, subgridB[2], p2, clear = false)
                scalarplot!(vis4, subgridB[3], n3, clear = false)
                scalarplot!(vis4, subgridB[3], p3, clear = false, label = "ref sol", legend =:best, show = true)
            end
        else

            scalarplot!(vis4, subgridp, compute_densities(iphin, regionAcceptor,  phin_solp, psi_solp), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis4, subgridp, compute_densities(iphip, regionAcceptor,  phip_solp, psi_solp), clear = false, color=:red)
            scalarplot!(vis4, subgridi, compute_densities(iphin, regionIntrinsic, phin_soli, psi_soli), clear = false, color=:green, linewidth = 5, yscale=:log)
            scalarplot!(vis4, subgridi, compute_densities(iphip, regionIntrinsic, phip_soli, psi_soli), clear = false, color=:red)
            scalarplot!(vis4, subgridn, compute_densities(iphin, regionDonor,     phin_soln, psi_soln), clear = false, color=:green, label ="\$ n_n \$", yscale=:log)
            scalarplot!(vis4, subgridn, compute_densities(iphip, regionDonor,     phip_soln, psi_soln), clear = false, label ="\$ n_p \$", color=:red)

            if ionicSpecies
                scalarplot!(vis4, subgridi, compute_densities(iphia, regionIntrinsic, phia_soli, psi_soli), clear = false, color=:gold, label ="\$ n_a \$", linewidth = 5, yscale=:log)
            end

            if plotCTWithoutIntSpec
                n1 = compute_densities(iphin, regionDonor, sol_ref[1:icoord_n, 2], sol_ref[1:icoord_n, 4])
                p1 = compute_densities(iphip, regionDonor, sol_ref[1:icoord_n, 3], sol_ref[1:icoord_n, 4])

                n2 = compute_densities(iphin, regionIntrinsic, sol_ref[icoord_n:icoord_ni, 2], sol_ref[icoord_n:icoord_ni, 4])
                p2 = compute_densities(iphip, regionIntrinsic, sol_ref[icoord_n:icoord_ni, 3], sol_ref[icoord_n:icoord_ni, 4])

                n3 = compute_densities(iphin, regionAcceptor, sol_ref[icoord_ni:end, 2], sol_ref[icoord_ni:end, 4])
                p3 = compute_densities(iphip, regionAcceptor, sol_ref[icoord_ni:end, 3], sol_ref[icoord_ni:end, 4])

                scalarplot!(vis4, subgridp, n1, clear = false, marker ="", label = "", linestyle=:dot, color=:black)
                scalarplot!(vis4, subgridp, p1, clear = false)
                scalarplot!(vis4, subgridi, n2, clear = false)
                scalarplot!(vis4, subgridi, p2, clear = false)
                scalarplot!(vis4, subgridn, n3, clear = false)
                scalarplot!(vis4, subgridn, p3, clear = false, label = "ref sol", legend =:best, show = true)
            end

        end

        ##### Ionmonger solution
        scalarplot!(vis4, grid_intr, n_im ,    clear = false, marker ="", label ="", linestyle=:dash, color=:black)
        scalarplot!(vis4, grid_intr, p_im,     clear = false)
        scalarplot!(vis4, grid_intr, a_im,     clear = false, label = "IM (without reco)")

        scalarplot!(vis4, grid_intr, n_imIR , clear = false, marker ="", label = "", linestyle=:dot, color=:gray)
        scalarplot!(vis4, grid_intr, p_imIR,  clear = false)
        scalarplot!(vis4, grid_intr, a_imIR,  clear = false, label = "IM (with reco)", legend =:best, show = true)

        ###############################################################################
        ##########                            IV                             ##########
        ###############################################################################


        IV_imIR        = readdlm("IM/with-surface-reco/0p1/IM-parameter-J.dat")
        IV_im         = readdlm("IM/without-surface-reco/0p1/IM-parameter-J.dat")

        if plotCTWithoutIntSpec
            scalarplot!(vis5, IV_ref[:, 1], IV_ref[:, 2], clear = false, label = "ref sol", color=:black, linestyle=:dot)
        end

        scalarplot!(vis5, biasValues,  abs.(IV),           clear = false, color=:green, linewidth = 5)
        scalarplot!(vis5, biasValues, abs.(IV_im[2:end]./(cm^2 .* 1.0e3)), clear = false, label = "IM (without reco)", color=:black, linestyle=:dot)
        scalarplot!(vis5, biasValues, abs.(IV_imIR[2:end]./(cm^2 .* 1.0e3)), clear = false, label = "IM (with reco)", color=:gray, linestyle=:dot, legend =:best, show = true)

    end

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, sol))/length(sol) # when using sparse storage, we get NaN values in solution
    return testval


end #  main



function test()
    testvalwithoutReco = -0.2988993117821689; testvalwithReco = -0.2988993117852727
    main(test = true, interfaceSpecies = true, leftInterface = true, interfaceReco=false) ≈ testvalwithoutReco && main(test = true, interfaceSpecies = true, leftInterface = true, interfaceReco=false) ≈ testvalwithReco

    # main(test = true, interfaceSpecies = false, leftInterface = true, interfaceReco=false) = -0.905305257878682
    # main(test = true, interfaceSpecies = false, leftInterface = true, interfaceReco=true)  = -0.892757917525781

end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
