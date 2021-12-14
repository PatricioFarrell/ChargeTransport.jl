Electronic and Ionic Charge Carriers
================================
In the following examples we will pay attention to specific types of drift-diffusion models for the description of the charge transport in perovskite solar cells, where we have electric and ionic charge carriers. Here, we assume to have three domains, denoted by 
$\mathbf{\Omega} = \mathbf{\Omega}_{\text{HTL}} \cup \mathbf{\Omega}_{\text{intr}} \cup \mathbf{\Omega}_{\text{ETL}}  $. 
The considered unknowns are the quasi Fermi potentials of electrons, holes and anion vacancies 
$\varphi_n, \varphi_p, \varphi_a$ 
and the electric potential 
$\psi$.
The underlying PDE reads, see [Abdel2021](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865)
```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( (p(\psi, \varphi_p) - N_A ) - (n(\psi, \varphi_n) - N_D) \Big),\\
	q \partial_t n(\psi, \varphi_n) - \nabla \cdot \mathbf{j}_n &= q\Bigl(G(\mathbf{x}) - R(n,p) \Bigr), \\
	q partial_t p(\psi, \varphi_p) + \nabla \cdot \mathbf{j}_p &= \Bigl(G(\mathbf{x}) - R(n,p) \Bigr),
\end{aligned}
``` 
for 
$\mathbf{x} \in \mathbf{\Omega}_{\text{HTL}} \cup  \mathbf{\Omega}_{\text{ETL}} $, $[0, t_F]$ and for $\mathbf{x} \in \mathbf{\Omega}_{\text{intr}} $, $[0, t_F]$ we have 
```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( p(\psi, \varphi_p)  - n(\psi, \varphi_n) + a(\psi, \varphi_a) - C_0 \Big),\\
q \partial_t n(\psi, \varphi_n)	- \nabla \cdot \mathbf{j}_n &= \Bigl(G(\mathbf{x}) - R(n,p) \Bigr), \\
	q \partial_t p(\psi, \varphi_p) + \nabla \cdot \mathbf{j}_p &= \Bigl(G(\mathbf{x}) - R(n,p) \Bigr),\\
	q \partial_t a(\psi, \varphi_a) + \nabla \cdot \mathbf{j}_a &= 0.
\end{aligned}
``` 

**General Information**
The extensions to the previous discussed simulation procedure in the previous example are the following

- another charge carrier, the anion vacancy, occurs (which is incorporated by the respective quasi Fermi potential)
- we allow jumps within the parameters entering the state equation
- the transient model is considered
- we allow a generation $G$ to be present

A quick survey on how to use our solver to adjust the input parameters such that these features can be simulated will be given in the following.

## Example 1: Solving the Stationary Problem with Graded Interfaces
We assume only electric charge carriers in this example. By default, we assume abrupt inner interfaces. If one wishes to simulate graded interfaces, where for example the effective density of states and the band-edge energy may vary, we refer to [Example105](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Example105_PSC_gradedFlux.jl) for more details.

First, we need to define two additional thin interface layers

```julia
# region numbers
regionDonor             = 1       # n doped region
regionJunction1         = 2
regionIntrinsic         = 3       # intrinsic region
regionJunction2         = 4
regionAcceptor          = 5       # p doped region
```
which need to be taken into account by the initialization of the grid.

Second, since we allow varying parameters within the thin interface layers, the flux discretization scheme needs to be adjusted and we need to construct a nodal dependent parameter struct

```julia
data.flux_approximation = scharfetter_gummel_graded

paramsnodal             = ParamsNodal(grid, numberOfCarriers)
```

Lastly, the respective parameters need to be graded. Currently, only a linear grading is implemented.

```julia
paramsnodal.bandEdgeEnergy[iphin, :]  = gradingParameter(paramsnodal.bandEdgeEnergy[iphin, :],
                                                        coord, regionTransportLayers, regionJunctions,
                                                        h, heightLayers, lengthLayers, EC)
```

## Example 2: A Linear I-V Measurement Scan Protocol
Here, the key parts of [Example106](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Example106_PSC_withIons_IVMeasurement.jl) are shortly summarized.

First, the charge carriers indices need to be extended since we assume here additional mobile anion vacancies
```julia
iphin                       = 2 # electron quasi Fermi potential
iphip                       = 1 # hole quasi Fermi potential
iphia                       = 3 # anion vacancy quasi Fermi potential
    
numberOfCarriers            = 3 # electrons, holes and anion vacancies
```
Another change is the choice of model_type since we consider here the dynamic problem. Further, we need to enable the ionic charge carriers on the specific active layers.
```julia
data.model_type             = model_transient
data.enable_ionic_carriers  = enable_ionic_carriers(ionic_carriers = [iphia], regions = [regionIntrinsic])
```

Arriving now at out of equilibrium calculations, we need to specify the scanrate and other information to set the time mesh. Currently, only linear scan protocols are predefined.

```julia
# primary data for I-V scan protocol
scanrate          = 1.0 * V/s
number_tsteps     = 31
endVoltage        = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary

# with fixed timestep sizes we can calculate the times
# a priori
tvalues           = set_time_mesh(scanrate, endVoltage, number_tsteps, type_protocol = linearScanProtocol)
```
Lastly, the time loop needs to be performed. Note that within the solve! method, we need to specify the time step.
```julia    
for istep = 2:number_tsteps
        
    t             = tvalues[istep]       # Actual time
    Δu            = t * scanrate         # Applied voltage 
    Δt            = t - tvalues[istep-1] # Time step size
        
    # Apply new voltage
    # set non equilibrium boundary conditions
    set_ohmic_contact!(ctsys, bregionAcceptor, Δu)

    # Solve time step problems with timestep Δt. initialGuess plays the role of the solution
    # from last timestep
    solve!(solution, initialGuess, ctsys, control  = control, tstep = Δt)


    initialGuess .= solution

end # time loop
```
## Example 3: Perovskite Solar Cell under I-V Scan Protocol and Illumination
Now, we can add to the previous calculation an illumination protocol. For this one needs to add

```julia
data.generation_model    = generation_uniform
```
and specify the uniform generation rate in each considered region, i.e.

```julia
for ireg in 1:numberOfRegions
    params.generationUniform[ireg]  = generationUniform[ireg]
end
```
where the input data is stored in generationUniform. Note that as Beer-Lambert generation is implemented, but yet not well-tested.
Further, we suggest to perform a time loop while increasing the generation rate and afterwards applying the scan protocol with a full generation due to numerical stability, see for this [Example107](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Example107_PSC_uniform_Generation.jl).

## Example 4: Solving a 2D Problem
Lastly, the code is capable of doing multi-dimensional calculations.

For a 2D mesh it is possible to use a structured grid via [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl), for this see [Example108](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Example108_PSC_2D_tensorGrid.jl).
But it is also possible to use the Julia wrapper [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl) to use Jonathan Richard Shewchuk's Triangle mesh generator, see [Example201 for the simulation on a rectangular grid](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Example201_PSC_2D_unstructuredGrid.jl) or [Example201 for a non-rectangular one](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Example201_2D_non_rectangularGrid.jl).

Lastly, with help of the [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl) wrapper, three dimensional meshes can be generated, see [Example202](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/examples/Example202_3D_grid.jl).