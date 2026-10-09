module test_neutral_diffusion
using Test, ChargeTransport, ExtendableGrids, VoronoiFVM
function test()
    grid = simplexgrid(collect(range(0.0, 1.0; length = 51)))
    data = Data(grid, 3)
    data.params.chargeNumbers .= [-1, 1, 0]
    enable_iodine_species!(data; iodineSpecies = 3, regions = [1],
        referenceDensity = 1.0e18, diffusionCoefficient = 0.1)
    ct = ChargeTransport.System(grid, data)
    # Exercise the CT callbacks, solving only the neutral continuity equation.
    function flux(f, u, edge, data)
        tmp = zeros(eltype(u), 4)
        v = zeros(eltype(u), 4, 2)
        v[3, :] .= u[3, :]
        ChargeTransport.flux!(tmp, v, edge, data, OutOfEquilibrium)
        f[3] = tmp[3]
    end
    function storage(f, u, node, data)
        tmp = zeros(eltype(u), 4)
        v = zeros(eltype(u), 4)
        v[3] = u[3]
        ChargeTransport.storage!(tmp, v, node, data, OutOfEquilibrium)
        f[3] = tmp[3]
    end
    sys = VoronoiFVM.System(grid; data, flux, storage, species = [3])
    u = VoronoiFVM.unknowns(sys; inival = 0.0)
    x = grid[Coordinates][1, :]
    u[3, :] .= log.(1 .+ 0.5 .* cos.(pi .* x))
    density(f, u, node, data) = (f[3] = get_density!(u, node, data, 3))
    mass(v) = VoronoiFVM.integrate(sys, density, v)[3, 1]
    initial = mass(u)
    full = zeros(4, length(x))
    full[3, :] .= u[3, :]
    @test isapprox(integrated_density(ct; sol = full, icc = 3, ireg = 1), initial; rtol = 1e-14)
    @test abs(get_current_val(ct, full)) < 1e-12
    @test abs(get_current_val(ct, full, zeros(size(full)), 0.1)) < 1e-12
    control = ChargeTransport.SolverControl()
    control.verbose = false
    control.abstol = 1e-12
    control.reltol = 1e-12
    drift = 0.0
    for step in 1:100
        u = VoronoiFVM.solve(sys; inival = u, tstep = 0.2, control)
        drift = max(drift, abs(mass(u) / initial - 1))
    end
    spread = maximum(abs.(exp.(u[3, :]) .- 1))
    println("Neutral diffusion: relative mass drift = ", drift, "; final relative nonuniformity = ", spread)
    @test drift < 1e-10
    @test spread < 1e-7
    @test all(iszero, sys.boundary_factors[3, :])
    data.params.chargeNumbers[3] = 1
    @test_throws ArgumentError ChargeTransport.System(grid, data)
    return true
end
end
