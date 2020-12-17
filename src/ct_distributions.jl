"""
Distribution integrals
"""

using ForwardDiff

"""
The Boltzmann approximation exp(x)
"""
function Boltzmann(x::Real)
    exp(x)
end

"""
The Blakemore approximation 1/(exp(-x) + γ)
"""
function Blakemore(x::Real, γ::Real )
    1/(exp(-x) + γ)
end

"""
The Blakemore approximation 1/(exp(-x) + γ) with
γ = 0.27.
"""
function Blakemore(x::Real)
    Blakemore(x, 0.27)
end


"""
The Fermi-Dirac integral of order minus one.
"""
function FermiDiracMinusOne(x::Real)
    Blakemore(x, 1.0)
end


"""
The incomplete Fermi-Dirac integral of order 1/2, 
implemented according to Bednarczyk and Bednarczyk
"The Approximation of the Fermi-Dirac integral F_1/2()"
"""
function FermiDiracOneHalf(x::Real)

    a = x^4 + 33.6*x * (1.0 - 0.68*exp(-0.17*(x+1)^2)) + 50
    sqrt(pi) / ( 2 * (3/4*sqrt(pi) * a^(-3/8) + exp(-x) ) )

end

"""
Degenerate limit of incomplete Fermi-Dirac integral of order 1/2.
"""
function degenerateLimit(x)
        x < 0 ? NaN : 4/(3*sqrt(pi))*x^(3/2)
end

"""
Plot different distribution integrals
"""
function plotDistributions(;Plotter=nothing)

    Plotter.close()

    rcParams = Plotter.PyDict(Plotter.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"

    # println(PyPlot.matplotlib.rcParams)


    x = -5:0.01:10;

    Plotter.semilogy(x, ChargeTransportInSolids.FermiDiracOneHalf.(x), label="\$F_{1/2}\$");
    Plotter.semilogy(x, Boltzmann.(x), label="Boltzmann");
    Plotter.semilogy(x, ones(size(x))/0.27, "--", label="\$1/\\gamma=3.\\overline{703}\$", color=(0.6,0.6,0.6,1));
    Plotter.semilogy(x, Blakemore.(x), label="Blakemore (\$\\gamma=0.27\$)");   
    Plotter.semilogy(x, degenerateLimit.(x),label="degenerate limit");

    Plotter.xlabel("\$\\eta\$")
    Plotter.ylabel("\$\\mathcal{F}(\\eta)\$")
    Plotter.title("Distributions")
    Plotter.legend()

    Plotter.show();
end

"""
Plot diffusion enhancements
"""
function plotDiffusionEnhancements()

    PyPlot.close()

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"

    x = -5:0.01:10;

    f = ChargeTransportInSolids.FermiDiracOneHalf; df = x -> ForwardDiff.derivative(f,x)
    PyPlot.semilogy(x, f.(x)./df.(x), label="\$F_{1/2}\$");

    f = ChargeTransportInSolids.Boltzmann; df = x -> ForwardDiff.derivative(f,x)
    PyPlot.semilogy(x, f.(x)./df.(x), label="Boltzmann");

    f = ChargeTransportInSolids.Blakemore; df = x -> ForwardDiff.derivative(f,x)
    PyPlot.semilogy(x, f.(x)./df.(x), label="Blakemore (\$\\gamma=0.27\$)");

    f = ChargeTransportInSolids.degenerateLimit; df = x -> ForwardDiff.derivative(f,x)
    PyPlot.semilogy(x, f.(x)./df.(x), label="degenerate limit");

    PyPlot.xlabel("\$\\eta\$")
    PyPlot.ylabel("\$g(\\eta)\$")
    PyPlot.title("Diffusion Enhancements")
    PyPlot.legend()

    PyPlot.show();
end



export Boltzmann, Blakemore, FermiDiracMinusOne, FermiDiracOneHalf



