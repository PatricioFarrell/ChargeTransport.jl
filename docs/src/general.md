Introducing the Model
================================
This solver is based on the nonlinear partial differential equations solver [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl). As the name indicates, for the spatial discretization the finite volume method is used.

[VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) considers a continuous problem of \$ m \$ coupled partial differential equations on a domain \$ \\mathbf{\Omega} \\subset \mathbb{R}^d\$ which is
```math
\begin{aligned} 
\partial_t s(u) + \nabla \cdot  \mathbf{j} \Bigl(u, \nabla u\Bigr) + r(u) = f, 
\end{aligned}
```
where \$ s \$ denotes the storage term, $ \\mathbf{j}$ the flux, \$ r\$ the reaction and \$ f \$ the source term.
In the following we specify these general and rather abstract nonlinear, time-dependent transport-reaction equation.

## General Charge Transport Model
Let \$ n_\alpha \$ denote the density of charge carrier \$ \\alpha \$ and \$ \\psi \$ the electrostatic potential. For describing the charge transport of the carriers within a device we use drift-diffusion equations. We consider the Poisson equation which is coupled to continuity equations for each charge carrying species \$ \\alpha \$ . The charge transport model which our solver aims to simulate reads as follows
```math
\begin{aligned}
- \nabla \cdot (\varepsilon_s \nabla \psi) &= q \sum_{\alpha} (n_\alpha - C_\alpha),\\
z_\alpha q \partial_t n_\alpha +  \nabla\cdot \mathbf{j}_\alpha 
	&= 
	z_\alpha q	r_\alpha.
\end{aligned}
```
In analogy to semiconductor theory, we assume that for the problems, we consider that the charge carrier density
$n_{\alpha}$
can be linked to corresponding quasi Fermi potentials 
$\varphi_\alpha$
via a statistical relation 
```math
\begin{aligned}
n_\alpha = N_\alpha \mathcal{F}_\alpha \Bigl(\eta_\alpha(\psi, \varphi_\alpha) \Bigr), \quad \eta_\alpha = z_\alpha \frac{q (\varphi_\alpha - \psi) + E_\alpha}{k_B T},
\end{aligned}
```
where the physical quantities are clarified in the list of notations. Here, the electric current is defined as
```math
\begin{aligned}
    \mathbf{j}_\alpha 
	=
    - (z_\alpha)^2 q \mu_\alpha  
    n_\alpha  
    \nabla\varphi_\alpha
    ~
\end{aligned}
```
with the gradient of the quasi Fermi potentials as driving forces. With this definition of the electric current 
$\mathbf{j}_{\alpha}$
the set of unknowns is given by the quasi Fermi potentials
$\varphi_\alpha$ of each mobile species
$\alpha$
and the electrostatic potential 
$\psi$, i.e. we have
$\alpha + 1$
unknowns.

For consistency with literature, it is possible to rewrite the electric current mathematically equivalent in terms of densities. For this, we introduce the generalized Einstein relation
```math
\begin{aligned}
	D_\alpha =\mu_\alpha U_T g_\alpha \left( \eta_\alpha \right),
\end{aligned}
```
where $g$ is the nonlinear diffusion enhancement given by
```math
\begin{aligned}
g_\alpha(\eta_\alpha) = \frac{\mathcal{F_\alpha}(\eta_\alpha)}{\mathcal{F_\alpha}'(\eta_\alpha)}.
\end{aligned}
```
Mathematically, the diffusion enhancement can be seen as nonlinear, potential-dependent diffusion. With the help of the generalized Einstein relation  it is now possible to derive the electric currents in drift-diffusion form
```math
\begin{aligned}
	\mathbf{j}_\alpha  = -q \Bigl( z_\alpha D_\alpha  \nabla n_\alpha + z_\alpha^2 \mu_\alpha n_\alpha \nabla \psi \Bigr),
\end{aligned}
```
where the diffusion may be nonlinear.

## Concrete Example: Perovskite Solar Cells

A concrete example of a charge transport model is the one describing charge carriers in a perovskite cell device.

## List of Notations

| **symbol** | **physical quantity** |   |   |   |   | **symbol** | **physical quantity** |
| :---:         |     :---:      |          :---: |          :---: |          :---: |          :---: |          :---: |          :---: |
| \$ \\alpha \$   | mobile charge carrier     |      |      |      |      | \$ n_\\alpha \$    | charge carrier density of \$ \\alpha \$    |
| \$\\varepsilon_s\$     | dielectric permittivity       |      |      |      |      | \$ \\psi \$      | electrostatic potential      |
| \$ q \$     | elementary charge       |      |      |      |      | \$C_\\alpha\$      | doping or background charge w.r.t. \$ \\alpha \$      |
| \$ z_\\alpha \$     | charge number w.r.t. \$ \\alpha \$       |      |      |      |      | \$ \\mathbf{j}_\\alpha \$      | current density w.r.t. \$ \\alpha \$      |
| \$ r_\\alpha \$     | production/reaction rate w.r.t. \$ \\alpha \$       |      |      |      |      | \$ \\mathbf{j}_\\alpha \$      | current density w.r.t. \$ \\alpha \$      |
| \$ N_\\alpha \$     | effective density of states w.r.t. \$ \\alpha \$       |      |      |      |      | \$ \\mathcal{F}_\\alpha \$      | statistics function      |
| \$ \\varphi_\\alpha \$     | quasi Fermi potential w.r.t. \$ \\alpha \$       |      |      |      |      | \$ E_\\alpha \$      | band-edge energy w.r.t. \$ \\alpha \$      |
| \$ k_B \$     | Boltzmann constant       |      |      |      |      | \$ T \$      | temperature      |
| \$ \\mu_\\alpha \$     | mobility of carrier \$ \\alpha \$      |      |      |      |      | \$ D_\alpha \$      | diffusion coefficient      |
| \$ U_T \$     | thermal voltage      |      |      |      |      | \$ g_\alpha \$      | diffusion enhancement      |
|       |        |      |      |      |      |        |        |
|       |        |      |      |      |      |        |        |
| \$ \\text{V}_{\\alpha} \$     | ionic defect of mobile ionic species \$\\alpha \$      |      |      |      |      | \$ n, p, a \$      | density of electrons, holes and anion vacancies      |
| \$ N_{\\text{A}} \$     | acceptor doping      |      |      |      |      | \$ N_{\\text{D}} \$      | donor doping      |
| \$ C_0 \$     | background charge      |      |      |      |      | \$ G \$      | generation rate      |
| \$ R \$     | recombination rate      |      |      |      |      |       |       |

## Information on How to Use the Code
Following general information holds, when using this code:
- When introducing **continuous quasi Fermi potentials:** If there are $N$ numbers of species $\alpha$ with $\alpha = 1, ..., N$ it is assumed that the first $N-1$ species correspond to the charge carriers. The final species corresponds to the electrostatic potential $\psi$.
- When introducing (possible) **discontinuous quasi Fermi potentials**:


## History and literature
..........


