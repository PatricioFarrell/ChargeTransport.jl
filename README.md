ChargeTransport.jl -- A drift-diffusion solver 
================================

[![Build status](https://github.com/PatricioFarrell/ChargeTransport.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/PatricioFarrell/ChargeTransport.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://PatricioFarrell.github.io/ChargeTransport.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PatricioFarrell.github.io/ChargeTransport.jl/dev)


This package is a prototype for solving drift-diffusion equations for the simulation of charge transport in solids that centers physics preserving discretization schemes.

Concerning the spatial discretization we rely on a Voronoi finite volume method implemented within the solver [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl). Currently, we use for the time discretization an implicit Euler method.

The simulations in following papers are based on ChargeTransport.jl

1. D. Abdel, P. Farrell and J. Fuhrmann.[ Assessing the quality of the excess chemical potential flux scheme for degenerate semiconductor device simulation.](https://link.springer.com/article/10.1007/s11082-021-02803-4) Optical and Quantum Electronics 53 (163) (2021).

2. D. Abdel, P. Vágner, J. Fuhrmann and P. Farrell. [Modelling charge transport in perovskite solar cells: Potential-based and limiting ion depletion.](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865) Electrochimica Acta 390 (2021).
