

Got curious? -- Short introduction into the developer side (@id Developerside)
================================

## Test environment

in package environment, run the tests, see, if your changes have affected other well-tested examples in the example folder.


## Organisation of the code

### The charge carrier organization


##  Discontinuous quasi Fermi potentials


As shortly mentioned in [the discussion of interface species](@ref interfaceSpecies), you
have the possibility to likewise simulate meaningful interface models. For this, the quasi
Fermi potentials do not necessarily need to be continuous.

Currently, only present ionic charge carriers, trap carriers and interface carriers are tested.

But the code is capable to simulate much more.

To understand how to properly extend the software and to introduce other types of charge carriers
or further surface effects, we have to understand the basic organization of the code.



To allow discontinuities you need to add the following lines.

```julia
data.isContinuous[iphin] = false
data.isContinuous[iphip] = false
```
