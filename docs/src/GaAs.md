van Roosbroeck system
================================

In both of the following examples, we solve the van Roosbroeck equations, a system of partial differential equations which describe current flow in a bipolar multi layer device:

```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( (p(\psi, \varphi_p) - N_A ) - (n(\psi, \varphi_n) - N_D) \Big),\\
	q \partial_t n(\psi, \varphi_n) -\nabla \cdot \mathbf{j}_n &= -qR(n,p), \\
	q \partial_t p(\psi, \varphi_p) + \nabla \cdot \mathbf{j}_p &= -qR(n,p).
\end{aligned}
```
Ohmic contacts will be used as boundary conditions. We will proceed as follows

Step 1: Initialize grid

Step 2: Initialize physical model

Step 3: Solve the problem in equilibrium

Step 4: Solve the problem for an applied bias

## Example 1: Stationary 1D problem (region doping)
We consider a three-layer GaAs p-i-n device in one dimension. We will explain [the PIN example](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Ex101_PIN.jl) in
greater detail.


### Step 1: Initialize grid
We have three layers and two external boundaries. We would like to solve the van Roosbroeck system on a uniform mesh with local grid refinement. We declare subregions and external boundaries.

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

# ExtendableGrids.cellmask! assigns region numbers
cellmask!(grid, [0.0 * μm],[h_pdoping], regionAcceptor)
cellmask!(grid, [h_pdoping],[h_pdoping + h_intrinsic], regionIntrinsic)
cellmask!(grid, [h_pdoping + h_intrinsic],[h_pdoping + h_intrinsic + h_ndoping], regionDonor)
```

### Step 2: Initialize physical model
Next, we choose relevant physical models such as the underlying statistics function or the recombination model. Additional options are stated in the comments.
Furthermore, we define the charge carrier indices. The index for the electrostatic potential is set automatically to `numberOfCarriers + 1`.

```julia
# Set indices for the quasi Fermi potentials
iphin                  = 1    # electrons
iphip                  = 2    # holes
numberOfCarriers       = 2

# Initialize Data instance
data                   = Data(grid, numberOfCarriers)

# Solve the stationary problem instead of the transient one
data.modelType         = Stationary

# Choose statistical relation between density and qF potential
# options: Boltzmann, FermiDiracOneHalfBednarczyk,
#          FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
data.F                .= Boltzmann

# Enable/Disable recombination processes, the default is stationary SRH recombination.
data.bulkRecombination = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                 bulk_recomb_Auger = true,
                                                 bulk_recomb_radiative = true,
                                                 bulk_recomb_SRH = true)

# choose boundary models
# exterior boundaries: OhmicContact and SchottkyContact
# interior boundaries: InterfaceModelNone, InterfaceModelSurfaceReco.
data.boundaryType[bregionAcceptor] = OhmicContact
data.boundaryType[bregionDonor]    = OhmicContact

# choose flux discretization scheme: ScharfetterGummel ScharfetterGummelGraded,
# ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
data.fluxApproximation             = ExcessChemicalPotential
```

Next, we fill in pre-defined or externally read in parameter values.

```julia
# params contains all necessary physical parameters
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

# Initialize a ChargeTransport struct
data.params   = params
ctsys         = System(grid, data, unknown_storage=unknown_storage)
```

We use ohmic contacts and set the applied voltage to zero in equilibrium.

```julia
set_contact!(ctsys, bregionAcceptor, Δu = 0.0)
set_contact!(ctsys, bregionDonor,    Δu = 0.0)
```

### Step 3: Solve the problem in equilibrium
Solve the equilibrium. Note that `control` refers to the Newton control
parameters given in `VoronoiFVM`.
```julia
solution      = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
initialGuess .= solution
```

### Step 4: Solve the problem for an applied bias
Starting from the equilibrium solution, we increase the applied voltage. Note that it is important to set `OutOfEqulibrium`.
```julia
data.calculationType = OutOfEquilibrium
maxBias              = voltageAcceptor # bias at acceptor boundary
biasValues           = range(0, stop = maxBias, length = 32)

for Δu in biasValues
    set_contact!(ctsys, bregionAcceptor, Δu = Δu) # non equilibrium bc
    solve!(solution, initialGuess, ctsys, control = control, tstep = Inf)
    initialGuess .= solution
end
```

### Step 5: Postprocessing
By adding the following line to the previous loop
```julia
current = get_current_val(ctsys, solution)
```
we have the possibility to calculate the total current.

Moreover, there are several different plotting routines, see [ct_plotting.jl](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/src/ct_plotting.jl).

## Example 2: Stationary 1D problem (nodal doping)

Now, instead of using regionwise doping it is possible to apply a nodal doping. (This is indeed also possible for other physical parameters, see the description of [ParamsNodal](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/ab0684293845859fb142ea69d786a88b597a8b67/src/ct_system.jl#L426).)
For this, go to previous Step 2, where you build your parameter set and adjust the doping initialization (code snippet is from [this example](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Ex102_PIN_nodal_doping.jl))

```julia
paramsnodal = ParamsNodal(grid, numberOfCarriers)

# initialize the space dependent doping
NDoping = 1.0e17  / cm^3; κ = 500.0
for icoord = 1:numberOfNodes
    t1 = tanh( (0.1 - coord[icoord]/μm) *κ )
    t2 = 1.0 + tanh( (coord[icoord]/μm - 0.2) * κ )
    paramsnodal.doping[icoord] = NDoping * 0.5 * ( 1.0  +  t1  - t2 )
end

data.paramsnodal  = paramsnodal
```