ChargeTransport.jl -- Simulating charge transport in semiconductors
================================

[![Build status](https://github.com/PatricioFarrell/ChargeTransport.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/PatricioFarrell/ChargeTransport.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://PatricioFarrell.github.io/ChargeTransport.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PatricioFarrell.github.io/ChargeTransport.jl/dev)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6275688.svg)](https://doi.org/10.5281/zenodo.6275688)


`ChargeTransport.jl` simulates charge transport in semiconductors. To this end, it discretizes
the semiconductor drift-diffusion equations via the Voronoi finite volume method as implemented in [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl).

### Special features

- heterostructures
- 1D, 2D and 3D simulations
- stationary and transient simulations
- IV curves and scan protocols
- an arbitrary amount of charge carriers may be added
- thermodynamically consistent, physics preserving numerical methods
- different charge carrier statistics per species (Boltzmann, Blakemore, Fermi-Dirac)

`ChargeTransport.jl` is a free software. For research purposes you may use it under the terms of the GNU Affero General Public License (AGPL). As a company you may contact any of the authors directly to obtain a commercial license. If you use this package in your work, please cite it according to [CITATION.bib](https://raw.githubusercontent.com/PatricioFarrell/ChargeTransport.jl/master/CITATION.bib)

### The following papers rely on `ChargeTransport.jl`

[1.] D. Abdel, P. Farrell and J. Fuhrmann. [Assessing the quality of the excess chemical potential flux scheme for degenerate semiconductor device simulation.](https://link.springer.com/article/10.1007/s11082-021-02803-4) Optical and Quantum Electronics 53 (163) (2021).

[2.] D. Abdel, P. Vágner, J. Fuhrmann and P. Farrell. [Modelling charge transport in perovskite solar cells: Potential-based and limiting ion depletion.](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865) Electrochimica Acta 390 (2021).

[3.] D. Abdel, C. Chainais-Hillairet, P. Farrell and M. Herda. [Numerical analysis of a finite volume scheme for charge transport in perovskite solar cells.](https://doi.org/10.48550/arXiv.2209.07934) preprint (2022).
