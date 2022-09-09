Interface models
================================

With help of `ChargeTransport.jl` you have the possibility to likewise simulate meaningful
interface models. For this, the quasi Fermi potentials do not necessarily need to be continuous.

In the following, we introduce briefly the main code snippets. We will illustrate this feature on the basis of the standard van Roosbroeck system.


Let $\mathbf{\Omega} \subseteq \mathbb{R}^d$,
$d \leq 3$, be an open, connected and bounded spatial domain with
$\mathbf{\Omega} = \left( \cup_{k } \overline{\mathbf{\Omega}}_{k}\right)^{\circ}
$,
$k \in \{\text{p}, \text{n}\}$, corresponding to the layers of a p-n device.
We denote the non-empty interface with codimension $1$ between both subdomains by
$\mathbf{\Sigma} = \partial \mathbf{\Omega}_{\text{p}} \cap \partial \mathbf{\Omega}_{\text{n}}$.
The model with $(\varphi_n, \varphi_p, \psi)$ as unknowns is given by

```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( C + z_p n_p(\psi, \varphi_p) + z_n n_n(\psi, \varphi_n) \Big),\\
	z_n q \partial_t n_n(\psi, \varphi_n) + \nabla \cdot \mathbf{j}_n &= - z_nqR(n_n,n_p), \\
	z_p q \partial_t n_p(\psi, \varphi_p) + \nabla \cdot \mathbf{j}_p &= - z_p qR(n_n,n_p).
\end{aligned}
```
In most applications following interface conditions for $\mathbf{x} \in \mathbf{\Sigma}$ are imposed

1. Continuity of
    1. quasi Fermi potentials
    2. electric potential
2. Conservation of
    1. electric current densities
    2. electric displacement flux

In the following we discuss the case, when 1.1., i.e. the continuity of quasi Fermi potentials, is not satisfied anymore.

Note that, it is important to identify within your grid that the interface is an actual internal boundary.

```julia
bregionJunction = 3

bfacemask!(grid, [h_ndoping], [h_ndoping + h_pdoping], bregionJunction)

```
Otherwise, you will run into issues.

Discontinuous quasi Fermi potentials
================================

To allow discontinuities you need to add the following lines.

```julia
data.isContinuous[iphin] = false
data.isContinuous[iphip] = false
```

**Do we likewise allow surface recombination for this case? If yes, how is it defined?
Current implementation:**

```julia
data.boundaryType[bregionJunction] = InterfaceModelDiscontqF
```

*To Do: The user does not necessarily need to type in InterfaceModelDiscontqF. This is an
information I can generate by my own internally.*

Possible (reasonable) choices here: InterfaceModelNone, InterfaceModelSurfaceReco


Further, we introduce the reactions rates for the discontinuity.

```julia
params.bReactDiscont[iphin, bregionJunction]     = 1.0e15 / s
params.bReactDiscont[iphip, bregionJunction]     = 1.0e15 / s
```

The larger these values, the more continuity can be observed.


Additional interface species
================================

Current code snippets needed

```julia
data.isContinuous[iphin] = false
data.isContinuous[iphip] = false
```

```julia
data.boundaryType[bregionJunction] = InterfaceModelDiscontqFInterfaceSpecies
```


**1. Is the interface species model, we introduce a specific one or the general model which shall
     always happen when we have discontinuous qF with additional interface species?
	 (need to know such that we have proper DataTypes)**


```julia
enable_interface_carriers!(data, bulkSpecies = [iphin, iphip],
                           interfaceSpecies = [iphinb, iphipb],
						   boundaryRegion = bregionJunction)
```
**2. does it make sense to allow one interfaceCarrier on both boundaries or should it be another carrier than?**


*1. Adjust wording 2. boundary either array or value 3. this should act as Jürgens enable_species!(). Then, carrier is added to a list!*
```julia
enable_interface_carrier!(data, bulkCarrier = iphin,
                         interfaceCarrier = iphinb,
						 boundaryRegion::Float64)
						 # or boundaryRegion::Array{Float64,1}
```


Visualization of effects
================================

**Some issues with VoronoiFVM.jl and ExtendableGrids.jl which need to be fixed.**
