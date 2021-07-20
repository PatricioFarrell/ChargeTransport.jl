
"""
$(SIGNATURES)

Method which determines with input parameters which inner interface model 
was chosen by user.

"""
function inner_interface_model(data::ChargeTransportData)

    countSurfaceReco = 0::Int64; countInterfaceCharge = 0::Int64

    # detect which interface model the user chooses by counting
    for ireg in 1:data.params.numberOfBoundaryRegions
        
        if     data.boundary_type[ireg] ==  interface_model_surface_recombination

            countSurfaceReco = countSurfaceReco + 1

        elseif data.boundary_type[ireg] == interface_model_ion_charge

            countInterfaceCharge = countInterfaceCharge + 1

        end

    end

     # build the system based on the input interface model
     if     countSurfaceReco > 0 # build surface_recombination based system

        return interface_model_surface_recombination

    elseif countInterfaceCharge > 0 # build ion interface charge system 

        # DA: currently for this case, since InterfaceQuantites is not well tested we stick with the no inferface case, i.e.
        #     the known case how we build the system. Since we did not change anything, the code in 
        #     Example107 works perfectly fine since the choice of species number is already set.
        return interface_model_none

    elseif countSurfaceReco + countInterfaceCharge == 0 # build system without further interface conditions

        return interface_model_none

    end

end

"""
$(SIGNATURES)

Method which determines with input parameters which inner interface model 
was chosen by user.

"""
function inner_interface_model(ctsys::ChargeTransportSystem)


    countSurfaceReco = 0::Int64; countInterfaceCharge = 0::Int64

    # detect which interface model the user chooses by counting
    for ireg in 1:ctsys.data.params.numberOfBoundaryRegions
        
        if     ctsys.data.boundary_type[ireg] ==  interface_model_surface_recombination

            countSurfaceReco = countSurfaceReco + 1

        elseif ctsys.data.boundary_type[ireg] == interface_model_ion_charge

            countInterfaceCharge = countInterfaceCharge + 1

        end

    end

     # build the system based on the input interface model
     if     countSurfaceReco > 0 # build surface_recombination based system

        return interface_model_surface_recombination

    elseif countInterfaceCharge > 0 # build ion interface charge system 

        # DA: currently for this case, since InterfaceQuantites is not well tested we stick with the no inferface case, i.e.
        #     the known case how we build the system. Since we did not change anything, the code in 
        #     Example107 works perfectly fine since the choice of species number is already set.
        return interface_model_none

    elseif countSurfaceReco + countInterfaceCharge == 0 # build system without further interface conditions

        return interface_model_none

    end

end

"""
$(SIGNATURES)

New system constructor which builds all necessary information needed based
on the input parameters concerning additional interface models.
The idea is that this new constructor will replace the current one.

"""
function ChargeTransportSystem2(grid, data ;unknown_storage)

    ctsys                 = ChargeTransportSystem()

    interface_model       = inner_interface_model(data)
    # here at this point a quantity based solving is built or not
    ctsys                 = build_system(grid, data, unknown_storage, interface_model)
    
    return ctsys

end

##########################################################
##########################################################
"""
$(SIGNATURES)

The core of the new system constructor. Here, the system for no additional interface model is build.

"""
function build_system(grid, data, unknown_storage, ::Type{interface_model_none})

    num_species  = data.params.numberOfCarriers + data.params.numberOfInterfaceCarriers + 1

    ctsys        = ChargeTransportSystem()

    ctsys.data   = data
    
    physics      = VoronoiFVM.Physics(data        = data,
                                      flux        = flux!,
                                      reaction    = reaction!,
                                      breaction   = breaction!,
                                      storage     = storage!,
                                      bstorage    = bstorage!
                                      )

    ctsys.fvmsys = VoronoiFVM.System(grid, physics, unknown_storage = unknown_storage)

    # for detection of number of species. In following versions we can simply delete num_species from physics initialization. 
    VoronoiFVM.increase_num_species!(ctsys.fvmsys, num_species) 

    return ctsys

end

"""
$(SIGNATURES)

The core of the new system constructor. Here, the system for additional
surface recombination at inner interfaces is build.

"""
function build_system(grid, data, unknown_storage, ::Type{interface_model_surface_recombination})

    ctsys        = ChargeTransportSystem()

    fvmsys = VoronoiFVM.System(grid, unknown_storage=unknown_storage)

    data.speciesQuantities[1] = DiscontinuousQuantity(fvmsys, 1:data.params.numberOfRegions; id = 1) # iphin
    data.speciesQuantities[2] = DiscontinuousQuantity(fvmsys, 1:data.params.numberOfRegions; id = 2) # iphip

    for icc in 3:data.params.numberOfCarriers
        data.speciesQuantities[icc] = ContinuousQuantity(fvmsys, data.enable_ion_vacancies; id = icc) # if ion vacancies are present
    end

    data.speciesQuantities[data.params.numberOfCarriers+1] = ContinuousQuantity(fvmsys, 1:data.params.numberOfRegions)    # last quantitiy is psi

    physics    = VoronoiFVM.Physics(data        = data,
                                    flux        = fluxDiscontqF!,
                                    reaction    = reactionDiscontqF!,
                                    breaction   = breactionDiscontqF!,
                                    storage     = storageDiscontqF!
                                    ### DA: insert additionally bstorage
                                    )

    # add the defined physics to system
    physics!(fvmsys, physics)

    # numberOfSpecies within system changes since we consider discontinuous quasi Fermi potentials
    data.params.numberOfSpeciesSystem = VoronoiFVM.num_species(fvmsys)

    ctsys.fvmsys = fvmsys
    ctsys.data   = data

    return ctsys

end

##########################################################
##########################################################

breactionDiscontqF!(f, u, bnode, data) = breactionDiscontqF!(f, u, bnode, data, data.boundary_type[bnode.region])


function breactionDiscontqF!(f, u, bnode, data, ::Type{interface_model_surface_recombination})

    ipsi = data.speciesQuantities[end]

    if data.calculation_type == inEquilibrium

        return emptyFunction()

    else
        params      = data.params
        paramsnodal = data.paramsnodal


        # Ei_a                = -4.05                 * eV   
        # Ei_i                = -4.60                 * eV   
        # Ei_d                = -5.00                 * eV   
    
        # EI                  = [Ei_a, Ei_i, Ei_d]
    
        # recombinationVelocity    = [1.0e1 1.0e7;
        #                             1.0e5 1.0e1]
        # recombinationTrapDensity = 0.5.* [trap_density!(1, 1, data, EI[1]) + trap_density!(1, 2, data, EI[2])   trap_density!(1, 2, data, EI[2]) + trap_density!(1, 3, data, EI[3]);
        # trap_density!(2, 1, data, EI[1]) + trap_density!(2, 2, data, EI[2])   trap_density!(2, 2, data, EI[2]) + trap_density!(2, 3, data, EI[3])]
                                

        # iphin = data.speciesQuantities[1]
        # iphip = data.speciesQuantities[2]
    
        ####################### idea 1 ###########################
        # for icc in data.speciesQuantities[1:2] # list of our charge carrier quantities

        #     etan1  = params.chargeNumbers[iphin.id] / params.UT * ( (u[iphin, 1] - u[ipsi]) + params.bandEdgeEnergy[iphin.id, bnode.cellregions[1]] / q ) # left
        #     etan2  = params.chargeNumbers[iphin.id] / params.UT * ( (u[iphin, 2] - u[ipsi]) + params.bandEdgeEnergy[iphin.id, bnode.cellregions[2]] / q ) # right

        #     etap1  = params.chargeNumbers[iphip.id] / params.UT * ( (u[iphip, 1] - u[ipsi]) + params.bandEdgeEnergy[iphip.id, bnode.cellregions[1]] / q ) # left
        #     etap2  = params.chargeNumbers[iphip.id] / params.UT * ( (u[iphip, 2] - u[ipsi]) + params.bandEdgeEnergy[iphip.id, bnode.cellregions[2]] / q ) # right

        #     n1     = (params.densityOfStates[iphin.id, bnode.cellregions[1]] + paramsnodal.densityOfStates[iphin.id, bnode.index])* data.F[iphin.id](etan1)
        #     n2     = (params.densityOfStates[iphin.id, bnode.cellregions[2]] + paramsnodal.densityOfStates[iphin.id, bnode.index])* data.F[iphin.id](etan2)

        #     p1     = (params.densityOfStates[iphip.id, bnode.cellregions[1]] + paramsnodal.densityOfStates[iphip.id, bnode.index])* data.F[iphip.id](etap1)
        #     p2     = (params.densityOfStates[iphip.id, bnode.cellregions[2]] + paramsnodal.densityOfStates[iphip.id, bnode.index])* data.F[iphip.id](etap2)
            
        #     average_n    = 0.5 * (n1 + n2)^(2/3)
        #     average_p    = 0.5 * (p1 + p2)^(2/3)
        #     average_phin = 0.5 * ( u[iphin, 1] - u[iphin, 2] )
        #     average_phip = 0.5 * ( u[iphip, 1] - u[iphip, 2] )

        #     exponentialTerm       = exp((q * average_phin - q  * average_phip ) / (kB * params.temperature))
        #     excessCarrierDensTerm = average_n * average_p * (1.0 - exponentialTerm)


        #     kernelSRH = 1.0 / (  1.0/recombinationVelocity[iphip.id, bnode.region-2] * (average_n + recombinationTrapDensity[iphin.id, bnode.region-2]) + 1.0/recombinationVelocity[iphin.id, bnode.region-2] * (average_p + recombinationTrapDensity[iphip.id, bnode.region-2]) )

        #     react2     = q * params.chargeNumbers[icc.id] *  kernelSRH *  excessCarrierDensTerm 

        #     f[icc, 1] =   react2
        #     f[icc, 2] = - react2
        #end

        ####################### idea 1 ###########################
        ##########################################################
        ####################### idea 2 ###########################

        # for icc in data.speciesQuantities[1:2] # list of our charge carrier quantities

        #     etan1  = params.chargeNumbers[iphin.id] / params.UT * ( (u[iphin, 1] - u[ipsi]) + params.bandEdgeEnergy[iphin.id, bnode.cellregions[1]] / q ) # left
        #     etan2  = params.chargeNumbers[iphin.id] / params.UT * ( (u[iphin, 2] - u[ipsi]) + params.bandEdgeEnergy[iphin.id, bnode.cellregions[2]] / q ) # right

        #     etap1  = params.chargeNumbers[iphip.id] / params.UT * ( (u[iphip, 1] - u[ipsi]) + params.bandEdgeEnergy[iphip.id, bnode.cellregions[1]] / q ) # left
        #     etap2  = params.chargeNumbers[iphip.id] / params.UT * ( (u[iphip, 2] - u[ipsi]) + params.bandEdgeEnergy[iphip.id, bnode.cellregions[2]] / q ) # right

        #     n1     = (params.densityOfStates[iphin.id, bnode.cellregions[1]] + paramsnodal.densityOfStates[iphin.id, bnode.index])* data.F[iphin.id](etan1)
        #     n2     = (params.densityOfStates[iphin.id, bnode.cellregions[2]] + paramsnodal.densityOfStates[iphin.id, bnode.index])* data.F[iphin.id](etan2)

        #     p1     = (params.densityOfStates[iphip.id, bnode.cellregions[1]] + paramsnodal.densityOfStates[iphip.id, bnode.index])* data.F[iphip.id](etap1)
        #     p2     = (params.densityOfStates[iphip.id, bnode.cellregions[2]] + paramsnodal.densityOfStates[iphip.id, bnode.index])* data.F[iphip.id](etap2)
            

        #     exponentialTerm1       = exp((q * u[iphin, 1] - q  * u[iphip, 1] ) / (kB * params.temperature))
        #     exponentialTerm2       = exp((q * u[iphin, 2] - q  * u[iphip, 2] ) / (kB * params.temperature))
        #     excessCarrierDensTerm1 = n1 * p1 * (1.0 - exponentialTerm1)
        #     excessCarrierDensTerm2 = n2 * p2 * (1.0 - exponentialTerm2)


        #     kernelSRH1 = 1.0 / (  1.0/recombinationVelocity[iphip.id, bnode.region-2] * (n1 + recombinationTrapDensity[iphin.id, bnode.region-2]) + 1.0/recombinationVelocity[iphin.id, bnode.region-2] * (p1 + recombinationTrapDensity[iphip.id, bnode.region-2]) )
        #     kernelSRH2 = 1.0 / (  1.0/recombinationVelocity[iphip.id, bnode.region-2] * (n2 + recombinationTrapDensity[iphin.id, bnode.region-2]) + 1.0/recombinationVelocity[iphin.id, bnode.region-2] * (p2 + recombinationTrapDensity[iphip.id, bnode.region-2]) )

        #     react1     = q * params.chargeNumbers[icc.id] *  kernelSRH1 *  excessCarrierDensTerm1 
        #     react2     = q * params.chargeNumbers[icc.id] *  kernelSRH2 *  excessCarrierDensTerm2 

        #     f[icc, 1] =   react1
        #     f[icc, 2] = - react2
        # end
        ####################### idea 2 ###########################

        for icc in data.speciesQuantities[1:2]

            eta1 = params.chargeNumbers[icc.id] / params.UT * ( (u[icc, 1] - u[ipsi]) + params.bandEdgeEnergy[icc.id, bnode.cellregions[1]] / q )
            eta2 = params.chargeNumbers[icc.id] / params.UT * ( (u[icc, 2] - u[ipsi]) + params.bandEdgeEnergy[icc.id, bnode.cellregions[2]] / q )

            dens1 = (params.densityOfStates[icc.id, bnode.cellregions[1] ] * data.F[icc.id](eta1))^(2/3)
            dens2 = (params.densityOfStates[icc.id, bnode.cellregions[2] ] * data.F[icc.id](eta2))^(2/3)


            # qF potentials 
            d         = [1.0e1 1.0e3;
                        1.0e7 1.0e1]
            react     = d[icc.id, bnode.region-2] * (u[icc, 1] - u[icc, 2])


            # densities 
            # d         = [1.0e-10 1.0e-10;
            #              1.0e-30 1.0e-15]
            # react     =  d[icc.id, bnode.region-2] * q * params.chargeNumbers[icc.id] * (dens1 - dens2)

            f[icc, 1] =   react
            f[icc, 2] = - react
        end

        ########    
       

        return f
    end

end
##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master reactionDiscontqF! function. This is the function which enters VoronoiFVM and hands over
reaction terms for concrete calculation type and bulk recombination model.

"""
reactionDiscontqF!(f, u, node, data) = reactionDiscontqF!(f, u, node, data, data.calculation_type)



function reactionDiscontqF!(f, u, node, data, ::Type{inEquilibrium})

    ipsi        = data.speciesQuantities[end]

    params      = data.params
    paramsnodal = data.paramsnodal

    for icc in data.speciesQuantities[1:end-1] # list of our quantities
        eta     = params.chargeNumbers[icc.id] / params.UT * ( (u[icc] - u[ipsi]) + params.bandEdgeEnergy[icc.id, node.region] / q )

        f[ipsi] = f[ipsi] - params.chargeNumbers[icc.id] * (params.doping[icc.id, node.region])  # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc.id] * params.densityOfStates[icc.id, node.region] * data.F[icc.id](eta) # add charge carrier
    
    end
    
    f[ipsi] = - data.λ1 * q * f[ipsi]

    for icc in data.speciesQuantities[1:end-1]
        f[icc] = u[icc] - 0.0
    end

    return f

end


recombination_kernel(data, ireg, iphin::VoronoiFVM.AbstractQuantity, iphip::VoronoiFVM.AbstractQuantity, n, p, ::Type{bulk_recombination_none}) = 0.0


function recombination_kernel(data, ireg, iphin::VoronoiFVM.AbstractQuantity, iphip::VoronoiFVM.AbstractQuantity,  n, p, ::Type{bulk_recombination_radiative})

    params = data.params

    return params.recombinationRadiative[ireg]

end


function recombination_kernel(data, ireg, iphin::VoronoiFVM.AbstractQuantity, iphip::VoronoiFVM.AbstractQuantity, n, p,::Type{bulk_recombination_trap_assisted})

    params = data.params

    kernelSRH = 1.0 / (  params.recombinationSRHLifetime[iphip.id, ireg] * (n + params.recombinationSRHTrapDensity[iphin.id, ireg]) + params.recombinationSRHLifetime[iphin.id, ireg] * (p + params.recombinationSRHTrapDensity[iphip.id, ireg]) )

    return kernelSRH

end


function recombination_kernel(data, ireg, iphin::VoronoiFVM.AbstractQuantity, iphip::VoronoiFVM.AbstractQuantity, n, p, ::Type{bulk_recombination_full})

    params = data.params

    # radiative recombination
    kernelRadiative = recombination_kernel(data, ireg, iphin, iphip, n, p, bulk_recombination_radiative)

    # SRH recombination
    kernelSRH       = recombination_kernel(data, ireg, iphin, iphip, n, p, bulk_recombination_trap_assisted)

    # Auger recombination
    kernelAuger     = (params.recombinationAuger[iphin.id, ireg] * n + params.recombinationAuger[iphip.id, ireg] * p)


    return kernelRadiative + kernelAuger + kernelSRH

end


# if using discontinuous Quantities, generation currently not tested!!!!
function reactionDiscontqF!(f, u, node, data, ::Type{outOfEquilibrium})

    ipsi        = data.speciesQuantities[end]

    params      = data.params
    paramsnodal = data.paramsnodal

    for icc in data.speciesQuantities[1:end-1] # list of our quantities

        eta     = params.chargeNumbers[icc.id] / params.UT * ( (u[icc] - u[ipsi]) + params.bandEdgeEnergy[icc.id, node.region] / q )

        f[ipsi] = f[ipsi] - params.chargeNumbers[icc.id] * (params.doping[icc.id, node.region])  # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc.id] * params.densityOfStates[icc.id, node.region] * data.F[icc.id](eta) # add charge carrier
    
    end
    
    f[ipsi] = - data.λ1 * q * f[ipsi]

    iphin = data.speciesQuantities[1]
    iphip = data.speciesQuantities[2]

    # rhs of continuity equations for electron and holes (bipolar reaction)
    etan                  = params.chargeNumbers[iphin.id] / params.UT * ( (u[iphin] - u[ipsi]) + params.bandEdgeEnergy[iphin.id, node.region] / q )
    etap                  = params.chargeNumbers[iphip.id] / params.UT * ( (u[iphip] - u[ipsi]) + params.bandEdgeEnergy[iphip.id, node.region] / q )
    n                     = (params.densityOfStates[iphin.id, node.region] + paramsnodal.densityOfStates[iphin.id, node.index])* data.F[iphin.id](etan)
    p                     = (params.densityOfStates[iphip.id, node.region] + paramsnodal.densityOfStates[iphip.id, node.index])* data.F[iphip.id](etap)
    exponentialTerm       = exp((q *u[iphin] - q  * u[iphip]) / (kB * params.temperature))
    excessCarrierDensTerm = n * p * (1.0 - exponentialTerm)

    for icc in data.speciesQuantities[1:2]

        # gives you the recombination kernel based on choice of user
        kernel = recombination_kernel(data, node.region, iphin, iphip, n, p, data.bulk_recombination_model)
            
        f[icc]          = q * params.chargeNumbers[icc.id] *  kernel *  excessCarrierDensTerm 
    end


    return f

end


##########################################################
##########################################################


fluxDiscontqF!(f, u, edge, data) = fluxDiscontqF!(f, u, edge, data, data.calculation_type)



function fluxDiscontqF!(f, u, edge, data, ::Type{inEquilibrium})

    ipsi    = data.speciesQuantities[end]

    params  = data.params

    dpsi    =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi] = - params.dielectricConstant[edge.region] * ε0 * dpsi

    return f
    
end


fluxDiscontqF!(f, u, edge, data, ::Type{outOfEquilibrium}) = fluxDiscontqF!(f, u, edge, data, data.flux_approximation)



"""
$(TYPEDSIGNATURES)

The classical Scharfetter-Gummel flux scheme. This also works for space-dependent band-edge energy, but
not for space-dependent effective DOS.

"""
function fluxDiscontqF!(f, u, edge, data, ::Type{ScharfetterGummel})

    ipsi         = data.speciesQuantities[end]

    params       = data.params
    paramsnodal  = data.paramsnodal
    
    dpsi         =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]      = - params.dielectricConstant[edge.region] * ε0 * dpsi

  
    for icc in data.speciesQuantities[1:end-1]

        j0       =  params.chargeNumbers[icc.id] * q *  params.UT * params.mobility[icc.id, edge.region] 

        etak     = params.chargeNumbers[icc.id] / params.UT * ( (u[icc, 1] - u[ipsi, 1]) + params.bandEdgeEnergy[icc.id, edge.region]/q )
        etal     = params.chargeNumbers[icc.id] / params.UT * ( (u[icc, 2] - u[ipsi, 2]) + params.bandEdgeEnergy[icc.id, edge.region]/q )

        bp, bm = fbernoulli_pm(params.chargeNumbers[icc.id] * dpsi / params.UT)

        f[icc] = - j0 * ( bm * params.densityOfStates[icc.id, edge.region]  * data.F[icc.id](etal) - bp * params.densityOfStates[icc.id, edge.region]  * data.F[icc.id](etak) ) 
        
    end
    
    return f

end


function fluxDiscontqF!(f, u, edge, data, ::Type{excessChemicalPotential})

    ipsi         = data.speciesQuantities[end]

    params       = data.params
    paramsnodal  = data.paramsnodal
    
    dpsi         =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]      = - params.dielectricConstant[edge.region] * ε0 * dpsi

  
    for icc in data.speciesQuantities[1:end-1]

        j0       =  params.chargeNumbers[icc.id] * q *  params.UT * params.mobility[icc.id, edge.region]

        etak     = params.chargeNumbers[icc.id] / params.UT * ( (u[icc, 1] - u[ipsi, 1]) + params.bandEdgeEnergy[icc.id, edge.region]/q )
        etal     = params.chargeNumbers[icc.id] / params.UT * ( (u[icc, 2] - u[ipsi, 2]) + params.bandEdgeEnergy[icc.id, edge.region]/q )

        Q        = params.chargeNumbers[icc.id] * ( dpsi/params.UT) + (etal - etak) - log(data.F[icc.id](etal)) + log(data.F[icc.id](etak) )

        bp, bm   = fbernoulli_pm(Q)

        f[icc] = - j0 * ( bm * params.densityOfStates[icc.id, edge.region]  * data.F[icc.id](etal) - bp * params.densityOfStates[icc.id, edge.region]  * data.F[icc.id](etak) ) 

    end

    return f

end
##########################################################
##########################################################
function breactionDiscontqF!(f, u, bnode, data, ::Type{ohmic_contact})

    params      = data.params
    paramsnodal = data.paramsnodal 

    ipsi        = data.speciesQuantities[end]

    # parameters
    α           = 1.0e-10                      # tiny penalty value
 
    innerRegion = bnode.cellregions[1] # tuple with [subRegionNumber 0]

    for icc in data.speciesQuantities[1:end-1]

        if typeof(icc) == VoronoiFVM.DiscontinuousQuantity{Int32}
            eta = params.chargeNumbers[icc.id] / params.UT * ( (u[icc.regionspec[innerRegion]] - u[ipsi]) + params.bBandEdgeEnergy[icc.id, bnode.region] / q )
        elseif typeof(icc) == VoronoiFVM.ContinuousQuantity{Int32}
            eta = params.chargeNumbers[icc.id] / params.UT * ( (u[icc] - u[ipsi]) + params.bBandEdgeEnergy[icc.id, bnode.region] / q )
        end

        f[ipsi] = f[ipsi] - params.chargeNumbers[icc.id] * ( params.bDoping[icc.id, bnode.region] ) # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc.id] * (params.bDensityOfStates[icc.id, bnode.region] + paramsnodal.densityOfStates[icc.id, bnode.index]) * data.F[icc.id](eta) # add charge carrier


    end
 
    f[ipsi] = - data.λ1 * 1 / α *  q * f[ipsi]

    
    return f

end

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Master storage! function. This is the function which enters VoronoiFVM and hands over
a storage term, if we consider transient problem.

"""
storageDiscontqF!(f, u, node, data) = storageDiscontqF!(f, u, node, data, data.model_type)

storageDiscontqF!(f, u, node, data, ::Type{model_stationary})  = emptyFunction()


"""
$(TYPEDSIGNATURES)

The storage term for time-dependent problems.
Currently, for the time-dependent current densities the implicit Euler scheme is used.
Hence, we have 

``f[n_\\alpha] =  z_\\alpha  q ∂_t n_\\alpha`` 

and for the electrostatic potential
``f[ψ] = 0``.

"""
function storageDiscontqF!(f, u, node, data, ::Type{model_transient})

    ipsi        = data.speciesQuantities[end]
    params      = data.params
    paramsnodal = data.paramsnodal

    for icc in data.speciesQuantities[1:end-1]
        eta     = params.chargeNumbers[icc.id] / params.UT * ( (u[icc] - u[ipsi]) + params.bandEdgeEnergy[icc.id, node.region] / q )
        f[icc]  = q * params.chargeNumbers[icc.id] * (params.densityOfStates[icc.id, node.region] + paramsnodal.densityOfStates[icc.id, node.index]) * data.F[icc.id](eta)

    end

    f[ipsi] = 0.0
    
    return f
end