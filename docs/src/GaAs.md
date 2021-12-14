Electronic Charge Carriers
================================

In both of the following examples we solve a system of partial differential equations for the description of the current flow in a bipolar three layer device (also called van Roosbroeck system)

```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( (p(\psi, \varphi_p) - N_A ) - (n(\psi, \varphi_n) - N_D) \Big),\\
	q \partial_t n(\psi, \varphi_n) -\nabla \cdot \mathbf{j}_n &= -qR(n,p), \\
	q \partial_t p(\psi, \varphi_p) + \nabla \cdot \mathbf{j}_p &= -qR(n,p)
\end{aligned}
```
for a given applied voltage.
The procedure will be the following
1. Initialize grid information
2. Initialize model information
3. Solve the problem in equilibrium
4. Solve the problem for an applied bias 

## Example 1: Solving the Stationary 1D Problem
We consider a three layer GaAs p-i-n device in one dimension. Here, we discuss [Example101_PIN](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Example101_PIN.jl).


### Step 1: Declare grid information
We have three layers and two external boundaries. We want to solve the system of equations on a uniform mesh with a local grid refinement.
Within this part we assign the indices of subregions and boundaries and declare from where to where the subregions are defined.

```julia
# region numbers
regionAcceptor          = 1          # p doped region
regionIntrinsic         = 2          # intrinsic region
regionDonor             = 3          # n doped region
regions                 = [regionAcceptor, regionIntrinsic, regionDonor]
numberOfRegions         = length(regions)

# boundary region numbers
bregionAcceptor         = 1
bregionDonor            = 2
bregions                = [bregionAcceptor, bregionDonor]
numberOfBoundaryRegions = length(bregions)

# grid
refinementfactor        = 2^(n-1)
h_pdoping               = 2 * μm
h_intrinsic             = 2 * μm
h_ndoping               = 2 * μm
coord                   = initialize_pin_grid(refinementfactor,
                                             h_pdoping,
                                             h_intrinsic,
                                             h_ndoping)

grid                    = simplexgrid(coord)

# ExtendableGrids.cellmask! for defining the subregions and assigning region number
cellmask!(grid, [0.0 * μm],[h_pdoping], regionAcceptor)
cellmask!(grid, [h_pdoping],[h_pdoping + h_intrinsic], regionIntrinsic)
cellmask!(grid, [h_pdoping + h_intrinsic],[h_pdoping + h_intrinsic + h_ndoping], regionDonor)
```

### Step 2: Declare model information
In this step, the indices of charge carriers need to be assigned. Automatically, the index numberOfCarriers 
$ + 1 $
is set for the electric potential. Additionally, we can choose here between several other model information e.g. such as the underlying statistics function or the recombination model. The possible choices are denoted above the respective variable.

```julia
# set indices of the quasi Fermi potentials
iphin                    = 1 # electron quasi Fermi potential
iphip                    = 2 # hole quasi Fermi potential
numberOfCarriers         = 2

# initialize Data instance and fill in data
data                     = Data(grid, numberOfCarriers)

# Following variable declares, if we want to solve stationary or transient problem
data.model_type          = model_stationary

# Following choices are possible for F: Boltzmann, FermiDiracOneHalfBednarczyk,
# FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
data.F                  .= Boltzmann

# Here, we need to specify which numbers are associated with electron and hole quasi Fermi
# potential. Further, the desired recombination processes can be chosen here. By default
# we use the stationary model for the SRH recombination.
data.bulk_recombination  = set_bulk_recombination(;iphin = iphin, iphip = iphip, 
                                                   bulk_recomb_Auger = true,
                                                   bulk_recomb_radiative = true,
                                                   bulk_recomb_SRH = true)

# Following choices are possible for boundary model: For contacts currently only ohmic_contact
# and schottky_contact are possible. For inner boundaries we have interface_model_none,
# interface_model_surface_recombination.
data.boundary_type[bregionAcceptor] = ohmic_contact                       
data.boundary_type[bregionDonor]    = ohmic_contact   
    
# Following choices are possible for the flux_discretization scheme: scharfetter_gummel,
# scharfetter_gummel_graded, excess_chemical_potential, excess_chemical_potential_graded,
# diffusion_enhanced, generalized_sg
data.flux_approximation             = excess_chemical_potential
```

Lastly, you are filling in your previously defined or externally read in parameter values.

```julia
# Params is a struct which contains all necessary physical parameters. If one wants to
# simulate space-dependent variable, one additionally needs to generate a ParamsNodal
# struct, see Example102.
params                                              = Params(grid, numberOfCarriers)

params.temperature                                  = T
params.UT                                           = (kB * params.temperature) / q
params.chargeNumbers[iphin]                         = -1
params.chargeNumbers[iphip]                         =  1

for ibreg in 1:numberOfBoundaryRegions   # boundary region data

    params.bDensityOfStates[iphin, ibreg]           = Nc
    params.bDensityOfStates[iphip, ibreg]           = Nv
    params.bBandEdgeEnergy[iphin, ibreg]            = Ec
    params.bBandEdgeEnergy[iphip, ibreg]            = Ev
end

for ireg in 1:numberOfRegions           # interior region data

    params.dielectricConstant[ireg]                 = εr

    # effective DOS, band-edge energy and mobilities
    params.densityOfStates[iphin, ireg]             = Nc
    params.densityOfStates[iphip, ireg]             = Nv
    params.bandEdgeEnergy[iphin, ireg]              = Ec
    params.bandEdgeEnergy[iphip, ireg]              = Ev
    params.mobility[iphin, ireg]                    = mun
    params.mobility[iphip, ireg]                    = mup

    # recombination parameters
    params.recombinationRadiative[ireg]             = Radiative
    params.recombinationSRHLifetime[iphin, ireg]    = SRH_LifeTime
    params.recombinationSRHLifetime[iphip, ireg]    = SRH_LifeTime
    params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
    params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity
    params.recombinationAuger[iphin, ireg]          = Auger
    params.recombinationAuger[iphip, ireg]          = Auger

end

# interior doping
params.doping[iphin, regionDonor]                   = Nd   
params.doping[iphin, regionIntrinsic]               = ni    
params.doping[iphip, regionIntrinsic]               = 0.0     
params.doping[iphip, regionAcceptor]                = Na

# boundary doping
params.bDoping[iphin, bregionDonor]                 = Nd     
params.bDoping[iphip, bregionAcceptor]              = Na 

# Region dependent params is now a substruct of data which is again a substruct of the
# system and will be parsed in next step.
data.params                                         = params

# in the last step, we initialize our system with previous data which is likewise dependent
# on the parameters. Important that this is in the end, otherwise our VoronoiFVMSys is
# not dependent on the data we initialized but rather on default data.
ctsys                                               = System(grid, data, 
                                                            unknown_storage=unknown_storage)
```

For this system we apply ohmic contacts. For the equilibrium calculations the applied voltage is zero.

```julia
set_ohmic_contact!(ctsys, bregionAcceptor, 0.0)
set_ohmic_contact!(ctsys, bregionDonor, 0.0)
```

### Step 3: Solve for equilibrium solution
With this call, the equilibrium solution will be calculated. Note that control refers to the Newton control
parameters and is a struct within VoronoiFVM.
```julia
solution         = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
initialGuess    .= solution 
```

### Step 4: Solve the stationary model with varying applied voltage 
Now, we increase the applied voltage, which is incorporated to the model by the quasi Fermi potential boundary conditions, and
solve the underlying problem for each new set of boundary conditions. Note that it is important to mark that we are now in outOfEqulibrium calculations.
```julia
ctsys.data.calculation_type  = outOfEquilibrium
maxBias                      = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
biasValues                   = range(0, stop = maxBias, length = 32)

for Δu in biasValues

    # set non equilibrium boundary conditions
    set_ohmic_contact!(ctsys, bregionAcceptor, Δu)

    solve!(solution, initialGuess, ctsys, control = control, tstep = Inf)

    initialGuess .= solution

end # bias loop
```

### Step 5: Postprocessing
By adding the following line to the previous loop
```julia
# get I-V data
current = get_current_val(ctsys, solution)
```
we have the possibility to calculate the total current.

Further, there are several different plotting routines which help to assess the quality of the numerical solution, see [ct_plotting.jl](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/src/ct_plotting.jl).

## Example 2: Adding a Nodal Dependent Doping

Now, instead of applying a region dependent doping it is possible to apply a nodal dependent one. (This is indeed also possible for other quantities, see the description of [ParamsNodal](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/ab0684293845859fb142ea69d786a88b597a8b67/src/ct_system.jl#L426).)
For this, go to previous Step 2, where you build your parameter set and adjust the doping initialization (code snippet from [Example102\_PIN\_nodal\_doping.jl](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Example102_PIN_nodal_doping.jl))

```julia
paramsnodal  = ParamsNodal(grid, numberOfCarriers)

# initialize the space dependent doping
NDoping           =   1.0e17  / cm^3
κ = 500.0
for icoord = 1:numberOfNodes
    paramsnodal.doping[icoord] = NDoping * 0.5 * ( 1.0  +  tanh( (0.1 - coord[icoord]/μm) *κ )  - ( 1.0 + tanh( (coord[icoord]/μm - 0.2) * κ )) )
end

# parse the substruct containg the nodal dependent parameters to the struct data
data.paramsnodal  = paramsnodal
```