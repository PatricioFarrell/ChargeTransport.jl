
Mathematical Description of the Problem
================================
## Why are we interested in charge transport?
Understanding charge transport in a medium is one of the keys to understand the exact physical operations happening in this medium.
Where charge transport happens, there is also an electric current flowing. 

These charge transport mechanisms can be captured in mathematical models. The model itself sometimes cannot be solved analytically and numerical simulations help us to solve the underlying model.

Now, ongoing interdisciplinary research is made for improving existing electronic devices or providing new technologies -- sophisticated theoretical tools can help us to predict the behavior of these devices of interest.

Hence, on the one hand reliable and physically correct models and on the other hand physics preserving numerical techniques are powerful utensils for a potentially optimizations of device designs without actually building cost-worthy prototypes, just by understanding the charge transport.

As one might surmised, there is not only one correct way to go for the model and for the simulations tools. Hence, we introduce in the following the theoretical models we exclusively work with.

## Charge Transport Description via Drift-Diffusion Equations
For describing the charge transport of carriers within a device we use drift-diffusion equations. We consider the Poisson equation for the electric potential \$ \\psi \$
```math
\begin{aligned}
- \nabla \cdot (\varepsilon_s \nabla \psi) &= qC + q \sum_{\alpha} z_\alpha n_\alpha.
\end{aligned}
```
Here, 
$\varepsilon_s$
denotes the dielectric permittivity and $ q $ the elementary charge. The right-hand side of the Poisson equation is given by the space charge density which is the sum of charge carrier density
$z_\alpha n_\alpha$
multiplied by the respective charge number and some doping $ C $, if given.
This nonlinear equation is coupled to additional continuity equations for each charge carrier $ \alpha $, which describe the motion of free charge carriers in dependance of an electric field
```math
\begin{aligned}
z_\alpha q \partial_t n_\alpha +  \nabla\cdot \mathbf{j}_\alpha 
	&= 
	z_\alpha q	r_\alpha.
\end{aligned}
```
Here,
$\mathbf{j}_\alpha$
stands for the carrier current and $ r_\alpha $ for some production/reduction rate. In our applications we only allow specific recombination processes and generation as the aforementioned rate. Further comments on this can be found in 
[the information concerning the numerics](@ref system).

## Analogy to Semiconductor Theory
This package mainly focuses currently semiconductors as an application. Hence, we assume that the charge carrier density
$n_\alpha$
can be linked to corresponding quasi Fermi potentials $\varphi_\alpha$ via a statistical relation 
```math
\begin{aligned}
n_\alpha = N_\alpha \mathcal{F}_\alpha \Bigl(\eta_\alpha(\psi, \varphi_\alpha) \Bigr), \quad \eta_\alpha = z_\alpha \frac{q (\varphi_\alpha - \psi) + E_\alpha}{k_B T},
\end{aligned}
```
where the physical quantities are clarified in the list of notations. With this definition we can formulate the carrier current given by
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
with the gradient of the quasi Fermi potentials as driving forces. 

!!! note

    The unknowns here are always defined as the quasi Fermi potentials \$ \varphi_a \$ and the electric potential \$ \psi \$.




For a good overview on drift-diffusion models and semiconductor applications such as the underlying numerical simulation we recommend the following books and book chapter:

1. P. Farrell, D. H. Doan, M. Kantner, J. Fuhrmann, T. Koprucki, and N. Rotundo. “Drift-Diffusion Models”. In: Optoelectronic Device Modeling and Simulation: Fundamentals, Materials, Nanostructures, LEDs, and Amplifiers. CRC Press Taylor & Francis Group, 2017, pp. 733–771.
2. S. Selberherr. Analysis and Simulation of Semiconductor Devices. Springer-Verlag, 1984.
3. S. M. Sze and K. K. Ng. Physics of Semiconductor Devices. Wiley, 2006.


## List of Notations

| **symbol** | **physical quantity** |   |   |   |   | **symbol** | **physical quantity** |
| :---:         |     :---:      |          :---: |          :---: |          :---: |          :---: |          :---: |          :---: |
| \$ \\alpha \$   | mobile charge carrier     |      |      |      |      | \$ n_\\alpha \$    | charge carrier density of \$ \\alpha \$    |
| \$\\varepsilon_s\$     | dielectric permittivity       |      |      |      |      | \$ \\psi \$      | electrostatic potential      |
| \$ q \$     | elementary charge       |      |      |      |      | \$C\$      | doping/background charge      |
| \$ z_\\alpha \$     | charge number w.r.t. \$ \\alpha \$       |      |      |      |      | \$ r_\\alpha \$     | production/reaction rate w.r.t. \$ \\alpha \$       |      |      |      |      | \$ \\mathbf{j}_\\alpha \$      | current density w.r.t. \$ \\alpha \$      |
| \$ N_\\alpha \$     | effective density of states w.r.t. \$ \\alpha \$       |      |      |      |      | \$ \\mathcal{F}_\\alpha \$      | statistics function      |
| \$ \\varphi_\\alpha \$     | quasi Fermi potential w.r.t. \$ \\alpha \$       |      |      |      |      | \$ E_\\alpha \$      | band-edge energy w.r.t. \$ \\alpha \$      |
| \$ k_B \$     | Boltzmann constant       |      |      |      |      | \$ T \$      | temperature      |
| \$ \\mu_\\alpha \$     | mobility of carrier \$ \\alpha \$      |      |      |      |      |        |        |

