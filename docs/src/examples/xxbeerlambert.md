```julia
using PyPlot
using DelimitedFiles

cm = 0.01
m = 1.0
s = 1.0


h_pdoping       = 3.00e-6 * cm + 1.0e-7 *cm
h_intrinsic     = 3.00e-5 * cm
h_ndoping       = 8.50e-6 * cm + 1.0e-7 *cm

h = h_pdoping + h_intrinsic +h_ndoping

grid1 = readdlm("Driftfusion-pedotpss-grid.dat")[220:510]
grid2 = readdlm("Driftfusion-pedotpss-grid.dat")[219:509]

Fph = (4.235553944562899e21)  / (m^2 * s)
α   = 2.078826021188540e7  / (m)

x1 = h_pdoping
x2 = h_pdoping + h_intrinsic

xx1 = x1-1.0e-9:1e-9:x2
xx2 = x1:1e-9:x2+1.0e-9


G = Array{Real,1}(undef, length(grid1) )
println(length(grid))
println(length(G))

G .= Fph * α * exp.( - α.*(grid1[:].-grid2[:]))

println(1.0e-6*G)

PyPlot.plot(grid1[:], 1.0e-6*G[:])
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

