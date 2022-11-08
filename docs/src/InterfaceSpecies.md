[Interface Species] (@id interfaceSpecies)
================================

With help of `ChargeTransport.jl` you have the possibility to likewise simulate present electric
interface charge carriers at the internal boundaries.

In the following, we introduce briefly the main code snippets along with the underlying continuous
charge transport model. We will illustrate this feature on the basis of the standard van Roosbroeck system.

Let $\mathbf{\Omega} \subseteq \mathbb{R}^d$,
$d \leq 3$, be an open, connected and bounded spatial domain divided into to appropriate subspaces $\mathbf{\Omega}_{-}, \mathbf{\Omega}_{+} $, corresponding to the layers of a p-n device.
We denote the non-empty interface with codimension $1$ between both subdomains by
$\mathbf{\Gamma} = \partial \mathbf{\Omega}_{-} \cap \partial \mathbf{\Omega}_{+}$.
The model with $(\psi, \varphi_n, \varphi_p)$ as unknowns is given by

```math
\begin{aligned}
	- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \Big( C + z_p n_p(\psi, \varphi_p) + z_n n_n(\psi, \varphi_n) \Big),\\
	z_n q \partial_t n_n(\psi, \varphi_n) + \nabla \cdot \mathbf{j}_n &= - z_nqR(n_n,n_p), \\
	z_p q \partial_t n_p(\psi, \varphi_p) + \nabla \cdot \mathbf{j}_p &= - z_p qR(n_n,n_p).
\end{aligned}
```
In most applications following interface conditions for $\mathbf{x} \in \mathbf{\Gamma}$ are imposed

1. Continuity of
    1. quasi Fermi potentials
    2. electric potential
2. Conservation of
    1. electric current densities
    2. electric displacement flux

Now, we discuss the case, the case of present interface charge carriers. Thus, all the previously mentioned conditions do not necessarily need to be true anymore.
For simplicity, we assume only a continuous electric potential.
The following figure illustrates the change in the set of unknowns.

![Interface Model code structure](images/interface-model-schematics.png)


Implementation
================================


Note that, it is important to identify within your grid that the interface is an actual internal boundary.

```julia
bregActive = 3

bfacemask!(grid, [h_ndoping], [h_ndoping + h_pdoping], bregActive)

```
Otherwise, you will run into issues.


Next, you need to enable the interface carrier. For this, we need information on the index of interface carrier, the corresponding bulk carrier and the active interface region

```julia
enable_interface_carrier!(data, bulkCarrier = iphin,
                          interfaceCarrier = iphinb,
                          bregions = [bregActive])
enable_interface_carrier!(data, bulkCarrier = iphip,
                          interfaceCarrier = iphipb,
                          bregions = [bregActive])
```

Further, we introduce the reaction coefficient (see paper)

```julia
params.bReactionCoefficient[iphin, bregionJunction]     = UT * mun * Nc/d
params.bReactionCoefficient[iphip, bregionJunction]     = UT * mup * Nv/d
```

Additional Information:
We can likewise use surface recombination. For this, adjust:

```julia
data.boundaryType[bregionJunction] = InterfaceRecombination # InterfaceNone
```

How to infer other possible surface effects and how to incorporate them into the code is explained [here](@ DeveloperSide).


Implementation (Interface Recombination)
================================
It would be nice to have a constructor like this


```julia
data.interfaceRecombination  = set_interface_recombination(;iphin = iphinb, iphip = iphipb,
                                                            bregions = [bregActive])
```
without having the user to pay so much attention concerning the boundaryType ...




Visualization of effects
================================

**Some issues with VoronoiFVM.jl and ExtendableGrids.jl which need to be fixed.**
