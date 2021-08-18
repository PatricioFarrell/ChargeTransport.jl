Gallium Arsenide p-i-n Device Simulation
================================

In both of these examples we will solve the following system of partial differential equations for the description of the current flow in a bipolar GaAs p-i-n device

```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( (p(\psi, \varphi_p) - N_A ) - (n(\psi, \varphi_n) - N_D) \Big),\\
	\nabla \cdot \mathbf{j}_n &= qR(n,p), \\
	\nabla \cdot \mathbf{j}_p &= -qR(n,p).
\end{aligned}
```
The procedure will be the following
1. Initialize grid information
2. Initialize model information
3. Solve the problem in equilibrium
4. Solve the stationary problem for an applied bias 

## Example 1: Solving the Stationary 1D Problem
We consider a three layer p-i-n device in one dimension. We assume each layer to be 
$ 2 \mu\mathrm{m} $ 
long.


### Step 1: Declare your grid information
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

# cellmask! for defining the subregions and assigning region number
cellmask!(grid, [0.0 * μm],[h_pdoping], regionAcceptor)
cellmask!(grid, [h_pdoping],[h_pdoping + h_intrinsic], regionIntrinsic)
cellmask!(grid, [h_pdoping + h_intrinsic],[h_pdoping + h_intrinsic + h_ndoping], regionDonor)
```

### Step 2: Declare your model information
In this step, the indices of charge carriers need to be assigned. Automatically, the index numberOfCarriers 
$ + 1 $
is set for the electric potential. Additionally, we can choose here between several other model information e.g. such as the underlying statistics function or the recombination model. The possible choices are denoted above the respective variable.

```julia
# set indices of the quasi Fermi potentials
iphin              = 1 # electron quasi Fermi potential
iphip              = 2 # hole quasi Fermi potential
numberOfCarriers   = 2

# initialize ChargeTransportData instance and fill in data
data                                = ChargeTransportData(grid, numberOfCarriers)

#### declare here all necessary information concerning the model ###

# Following variable declares, if we want to solve stationary or transient problem
data.model_type                     = model_stationary

# Following choices are possible for F: Boltzmann, FermiDiracOneHalfBednarczyk,
# FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
data.F                             .= Boltzmann

# Here the user can specify, if they assume continuous or discontinuous charge carriers.
data.isContinuous[iphin]            = true
data.isContinuous[iphip]            = true

# Following choices are possible for recombination model: bulk_recomb_model_none,
# bulk_recomb_model_trap_assisted, bulk_recomb_radiative, bulk_recomb_full <: 
# bulk_recombination_model 
data.bulk_recombination             = set_bulk_recombination(iphin = iphin, iphip = iphip, 
                                             bulk_recombination_model = bulk_recombination)

# Following choices are possible for boundary model: For contacts currently only
# ohmic_contact and schottky_contact are possible. For inner boundaries we have
# interface_model_none, interface_model_surface_recombination, interface_model_ion_charge
# (distinguish between left and right).
data.boundary_type[bregionAcceptor] = ohmic_contact                       
data.boundary_type[bregionDonor]    = ohmic_contact   
    
# Following choices are possible for the flux_discretization scheme: ScharfetterGummel, 
# ScharfetterGummel_Graded, excessChemicalPotential, excessChemicalPotential_Graded, 
# diffusionEnhanced, generalized_SG
data.flux_approximation             = excessChemicalPotential
```

Lastly, you are filling in your previously defined or externally read in parameter values.

```julia
# Params is a struct which contains all necessary physical parameters. If one wants to
# simulate space-dependent variable, one additionally needs to generate a ParamsNodal
# struct, see Example102.
params                                              = ChargeTransportParams(grid,
                                                                            numberOfCarriers)

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
ctsys                                               = ChargeTransportSystem(grid, data, 
                                                            unknown_storage=unknown_storage)
```

For this system we apply ohmic contacts with a zero applied voltage in the first step

```julia
set_ohmic_contact!(ctsys, bregionAcceptor, 0.0)
set_ohmic_contact!(ctsys, bregionDonor, 0.0)
```

### Step 3: Solve for equilibrium solution

### Step 4: Solve the stationary model with varying applied voltage 

## Example 2: Adding a Nodal Dependent Doping

Now, instead of applying a region dependent doping it is possible to apply a nodal dependent one. For this, go to previous Step 2, where you build your parameter set and do

xxx
