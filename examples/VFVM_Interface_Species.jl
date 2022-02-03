#=

# 111: Interface Species
([source code](SOURCE_URL))

Nondimensionalized semiconductor device equations (with artificial doping)
in two setting:
a. Continuous quantities
b. Discontinuousquantities with additional Interfacequantities.

When varying reaction rate, the currents coincide.
=#

module VFVM_Interface_Species

using VoronoiFVM
using ExtendableGrids
using GridVisualize
using PyPlot
#using ChargeTransport


function main(;n=5, Plotter = PyPlot, tend = 20.0, unknown_storage=:sparse,
              reactionN1 = 1.e8, reactionP1 = 1.e8,     # "small" jumps
              reactionN2 = 5.e10, reactionP2 = 5.0e10)

    ################################################################################
    #### grid
    ################################################################################
    ## region numbers
    cm = 0.01
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

    coord_p_u               = collect(range(x0, h_pdoping/2, step=h_pdoping/(0.9*δ)))
    coord_p_g               = geomspace(h_pdoping/2, 
                                        h_pdoping, 
                                        h_pdoping/(1.2*δ), 
                                        h_pdoping/(1.2*δ), 
                                        tol=t)
    coord_i_g1              = geomspace(h_pdoping, 
                                        h_pdoping+h_intrinsic/k, 
                                        h_intrinsic/(7.1*δ), 
                                        h_intrinsic/(7.1*δ), 
                                        tol=t)
    coord_i_g2              = geomspace(h_pdoping+h_intrinsic/k, 
                                        h_pdoping+h_intrinsic,               
                                        h_intrinsic/(7.1*δ),    
                                        h_intrinsic/(7.1*δ), 
                                        tol=t)
    coord_n_g               = geomspace(h_pdoping+h_intrinsic,               
                                        h_pdoping+h_intrinsic+h_ndoping/2, 
                                        h_ndoping/(2.0*δ),   
                                        h_ndoping/(2.0*δ),      
                                        tol=t)
    coord_n_u               = collect(range(h_pdoping+h_intrinsic+h_ndoping/2, h_pdoping+h_intrinsic+h_ndoping, step=h_pdoping/(1.0*δ)))

    coord                   = glue(coord_p_u,coord_p_g,  tol=10*t) 
    icoord_p                = length(coord)
    coord                   = glue(coord,    coord_i_g1, tol=10*t)
    coord                   = glue(coord,    coord_i_g2, tol=10*t) 
    icoord_pi               = length(coord)
    coord                   = glue(coord,    coord_n_g,  tol=10*t)
    coord                   = glue(coord,    coord_n_u,  tol=10*t)
    grid                    = ExtendableGrids.simplexgrid(coord)
    numberOfNodes           = length(coord)

    ## cellmask! for defining the subregions and assigning region number (doping profiles do not intersect)
    cellmask!(grid, [0.0],                      [h_pdoping],                           regionAcceptor, tol = 1.0e-18)   # p-doped region   = 1
    cellmask!(grid, [h_pdoping],                [h_pdoping + h_intrinsic],             regionIntrinsic, tol = 1.0e-18)  # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic],  [h_pdoping + h_intrinsic + h_ndoping], regionDonor, tol = 1.0e-18)      # n-doped region   = 3

    ## bfacemask! for ``active'' boundary regions, i.e. internal interfaces. On the outer boundary regions, the 
    ## conditions will be formulated later
    bfacemask!(grid, [h_pdoping],               [h_pdoping],                           bregionJunction1)  # first  inner interface
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic],             bregionJunction2)  # second inner interface

    #gridplot(grid, Plotter = PyPlot, legend=:rt)

    ################################################################################
    #########  system
    ################################################################################

    sysC     = VoronoiFVM.System(grid, unknown_storage = unknown_storage)
    iphinC   = ContinuousQuantity(sysC, 1:numberOfRegions,    id = 1)
    iphipC   = ContinuousQuantity(sysC, 1:numberOfRegions,    id = 2)
    ipsiC    = ContinuousQuantity(sysC, 1:numberOfRegions,    id = 3)

    sysD     = VoronoiFVM.System(grid, unknown_storage = unknown_storage)
    iphinD   = DiscontinuousQuantity(sysD, 1:numberOfRegions, id = 1)
    iphipD   = DiscontinuousQuantity(sysD, 1:numberOfRegions, id = 2)
    iphin_b1 = InterfaceQuantity(sysD, bregionJunction1,      id = 3)
    iphip_b1 = InterfaceQuantity(sysD, bregionJunction1,      id = 4)
    ipsiD    = ContinuousQuantity(sysD, 1:numberOfRegions,    id = 7)

    NA  = [10.0, 0.0, 0.0];  ND = [0.0, 0.0, 10.0]

    function storageC!(f, u, node)

        etanC = - ( (u[iphinC] - u[ipsiC]) )
        etapC =   ( (u[iphipC] - u[ipsiC]) )

        f[iphinC] = - exp(etanC)
        f[iphipC] =   exp(etapC)
    
        f[ipsiC]  =  0.0

    end

    function storageD!(f, u, node)

        etanD = - ( (u[iphinD] - u[ipsiD]) )
        etapD =   ( (u[iphipD] - u[ipsiD]) )

        f[iphinD] = - exp(etanD)
        f[iphipD] =   exp(etapD)
    
        f[ipsiD]  =  0.0

    end
    #####################################
    function reactionC!(f, u, node)

        etanC     = - ( (u[iphinC] - u[ipsiC]) )
        etapC     =   ( (u[iphipC] - u[ipsiC]) )

        f[ipsiC]  = - (ND[node.region] - exp(etanC) + exp(etapC) - NA[node.region])
        ########################
        r0        = 1.0e-4
        recombC   = r0 * exp(etanC) * exp(etapC)

        f[iphinC] =  - recombC
        f[iphipC] =    recombC
    end

    function reactionD!(f, u, node)

        etanD    = - ( (u[iphinD] - u[ipsiD]) )
        etapD    =   ( (u[iphipD] - u[ipsiD]) )

        f[ipsiD] = - (ND[node.region] - exp(etanD) + exp(etapD) - NA[node.region])
        ########################
        r0        = 1.0e-4
        recombD   = r0 * exp(etanD) * exp(etapD)

        f[iphinD] =  - recombD
        f[iphipD] =    recombD
    end
    #####################################
    function fluxC!(f, u, node)

        f[ipsiC] =  - (u[ipsiC, 2] - u[ipsiC, 1])

        ########################
        bpC, bmC = fbernoulli_pm(-  (u[ipsiC, 2] - u[ipsiC, 1]) )

        etanC1 = - ( (u[iphinC, 1] - u[ipsiC, 1]) )
        etapC1 =   ( (u[iphipC, 1] - u[ipsiC, 1]) )
        etanC2 = - ( (u[iphinC, 2] - u[ipsiC, 2]) )
        etapC2 =   ( (u[iphipC, 2] - u[ipsiC, 2]) )

        f[iphinC]  =   (bmC * exp(etanC2) - bpC * exp(etanC1))
        f[iphipC]  = - (bpC * exp(etapC2) - bmC * exp(etapC1))
    end

    function fluxD!(f, u, node)

        f[ipsiD] =  - (u[ipsiD, 2] - u[ipsiD, 1])
 
        ########################
        bpD, bmD = fbernoulli_pm(-  (u[ipsiD, 2] - u[ipsiD, 1]) )

        etanD1 = - ( (u[iphinD, 1] - u[ipsiD, 1]) )
        etapD1 =   ( (u[iphipD, 1] - u[ipsiD, 1]) )

        etanD2 = - ( (u[iphinD, 2] - u[ipsiD, 2]) )
        etapD2 =   ( (u[iphipD, 2] - u[ipsiD, 2]) )

        f[iphinD]  =   (bmD * exp(etanD2) - bpD * exp(etanD1))
        f[iphipD]  = - (bpD * exp(etapD2) - bmD * exp(etapD1))
    end

    function breaction!(f, u, bnode)

        if bnode.region == bregionJunction1

            # left values
            nleft    = exp(- ( (u[iphinD, 1] - u[ipsiD]) ))
            pleft    = exp(  ( (u[iphipD, 1] - u[ipsiD]) ))

            # interface species
            n_interf = exp(- ( (u[iphin_b1]  - u[ipsiD]) )) 
            p_interf = exp(  ( (u[iphip_b1]  - u[ipsiD]) ))

            # right values
            nright   = exp(- ( (u[iphinD, 2] - u[ipsiD]) ))
            pright   = exp(  ( (u[iphipD, 2] - u[ipsiD]) ))
            ################

            # left and right reaction for n
            f[iphinD, 1] = reactionN1 * (nleft  - n_interf)
            f[iphinD, 2] = reactionN1 * (nright - n_interf) # minus from normal vector cancels out with sign from brackets

            
            # left and right reaction for p
            f[iphipD, 1] = reactionP1 * (pleft  - p_interf)
            f[iphipD, 2] = reactionP1 * (pright - p_interf) # minus from normal vector cancels out with sign from brackets

            @show f[iphipD, 1].value
            @show f[iphipD, 2].value

            # interface species reaction
            f[iphin_b1] =  - (f[iphinD, 1] + f[iphinD, 2]) # since its a reaction we need it with a minus
            f[iphip_b1] =  - (f[iphipD, 1] + f[iphipD, 2]) # since its a reaction we need it with a minus

        end

        if bnode.region == bregionJunction2

            # left values
            nleft    = exp(- ( (u[iphinD, 1] - u[ipsiD]) ))
            pleft    = exp(  ( (u[iphipD, 1] - u[ipsiD]) ))

            # right values
            nright   = exp(- ( (u[iphinD, 2] - u[ipsiD]) ))
            pright   = exp(  ( (u[iphipD, 2] - u[ipsiD]) ))
            ################ # same conclusions hold here true considering the signs

            f[iphinD, 1] =   reactionN2 *  (nleft - nright)
            f[iphinD, 2] = - reactionN2 *   (nleft - nright)

            # left and right reaction for p
            f[iphipD, 1] =   reactionP2 *  (pleft - pright)
            f[iphipD, 2] = - reactionP2 *   (pleft - pright)


        end
    end

    function bstorage!(f, u, bnode)

        f[ipsiD]  =  0.0

        if bnode.region == bregionJunction1

            etanb = - ( (u[iphin_b1] - u[ipsiD]) )
            etapb =   ( (u[iphip_b1] - u[ipsiD]) )
    
            f[iphin_b1] = - exp(etanb)
            f[iphip_b1] =   exp(etapb)
        
        end

    end

    #####################################
    physics!(sysC, VoronoiFVM.Physics(
        flux      = fluxC!,
        reaction  = reactionC!
    ))

    physics!(sysD, VoronoiFVM.Physics(
        flux      = fluxD!,
        reaction  = reactionD!,
        breaction = breaction!,
        bstorage  = bstorage!
    ))

    #####################################

    boundary_dirichlet!(sysC, iphinC, bregionAcceptor, 4.0)
    boundary_dirichlet!(sysC, iphipC, bregionAcceptor, 4.0)
    boundary_dirichlet!(sysC, ipsiC,  bregionAcceptor, 0.0)
    boundary_dirichlet!(sysC, iphinC, bregionDonor,    0.0)
    boundary_dirichlet!(sysC, iphipC, bregionDonor,    0.0)
    boundary_dirichlet!(sysC, ipsiC,  bregionDonor,    5.0)

    boundary_dirichlet!(sysD, iphinD, bregionAcceptor, 4.0)
    boundary_dirichlet!(sysD, iphipD, bregionAcceptor, 4.0)
    boundary_dirichlet!(sysD, ipsiD,  bregionAcceptor, 0.0)
    boundary_dirichlet!(sysD, iphinD, bregionDonor,    0.0)
    boundary_dirichlet!(sysD, iphipD, bregionDonor,    0.0)
    boundary_dirichlet!(sysD, ipsiD,  bregionDonor,    5.0)

    control                   = VoronoiFVM.NewtonControl()
    control.max_iterations    = 300
    control.verbose           = false
    control.tol_absolute      = 1.0e-10
    control.tol_relative      = 1.0e-10
    control.handle_exceptions = true
    control.tol_round         = 1.0e-10
    control.max_round         = 5
    control.damp_initial      = 0.05
    control.damp_growth       = 1.21 # >= 1
    control.max_round         = 5

    ################################################################################
    #########  Bias Loop
    ################################################################################

    ## Create a solution array
    inivalC = unknowns(sysC); inivalC .= 0.0
    solC    = unknowns(sysC)

    inivalD = unknowns(sysD); inivalD .= 0.0
    solD    = unknowns(sysD)


    biasval  = range(0, stop = 2.0, length = 10)
    Icspec   = zeros(0)
    Idspec   = zeros(0)

    for Δu in biasval

        boundary_dirichlet!(sysC, iphinC, bregionAcceptor, 4.0 + Δu)
        boundary_dirichlet!(sysC, iphipC, bregionAcceptor, 4.0 + Δu)
        boundary_dirichlet!(sysC, ipsiC,  bregionAcceptor, 0.0 + Δu)

        boundary_dirichlet!(sysD, iphinD, bregionAcceptor, 4.0 + Δu)
        boundary_dirichlet!(sysD, iphipD, bregionAcceptor, 4.0 + Δu)
        boundary_dirichlet!(sysD, ipsiD,  bregionAcceptor, 0.0 + Δu)

        println("Δu = ", Δu)

        solve!(solC, inivalC, sysC, control = control)
        inivalC .= solC
        solve!(solD, inivalD, sysD, control = control)
        inivalD .= solD

        ## get current
        factoryC = VoronoiFVM.TestFunctionFactory(sysC)
        tfC      = testfunction(factoryC, [1], [2])
        IC       = integrate(sysC, tfC, solC)

        factoryD = VoronoiFVM.TestFunctionFactory(sysD)
        tfD      = testfunction(factoryD, [1], [2])
        ID       = integrate(sysD, tfD, solD)

        valc = 0.0
        for ii = 1:length(IC)-1
            valc = valc + IC[ii]
        end

        vald = 0.0
        for ii = 1:length(ID)-1
            vald = vald + ID[ii]
        end

        push!(Icspec, valc)
        push!(Idspec, vald )

    end # bias loop
    

    ################################################################################
    #########  Plotting
    ################################################################################

    vis = GridVisualizer(Plotter = PyPlot, layout=(3,1))
    scalarplot!(vis[1, 1], grid, solC[iphinC, :], clear=false, label = "\$ \\varphi_n \$", color=:green)
    scalarplot!(vis[1, 1], grid, solC[iphipC, :], clear=false, label = "\$ \\varphi_p \$", color=:red, show=true)
    scalarplot!(vis[1, 1], grid, solC[ipsiC, :],  clear=false, label = "\$ \\psi \$",      color=:blue, show=true)
    Plotter.legend(fancybox = true, loc = "best", fontsize=11)
    Plotter.title("Solution (Continuous Setting)")

    subgrids = VoronoiFVM.subgrids(iphinD, sysD)
    phin_sol = VoronoiFVM.views(solD, iphinD, subgrids, sysD)
    phip_sol = VoronoiFVM.views(solD, iphipD, subgrids, sysD)
    psi_sol  = VoronoiFVM.views(solD, ipsiD, subgrids, sysD)

    for i = 1:length(phin_sol)
        scalarplot!(vis[2, 1], subgrids[i], phin_sol[i], clear = false, color=:green)
        scalarplot!(vis[2, 1], subgrids[i], phip_sol[i], clear = false, color=:red)
        scalarplot!(vis[2, 1], subgrids[i], psi_sol[i],  clear = false, color=:blue)
    end

    Plotter.title("Solution (Disontinuous Setting)")

    scalarplot!(vis[3, 1], biasval, Icspec, clear = false, label = "continuous", color=:green)
    scalarplot!(vis[3, 1], biasval, Idspec, clear = false, label = "discontinuous", color=:red)
    Plotter.legend(fancybox = true, loc = "best", fontsize=11)

    Plotter.legend(fancybox = true, loc = "best", fontsize=11)
    Plotter.title("IV Curves")

end # main

end # module
