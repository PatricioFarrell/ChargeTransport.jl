"""
Distribution integrals.
"""

using ForwardDiff

"""
$(TYPEDSIGNATURES)


The Boltzmann approximation ``\\exp(x)``.
"""
function Boltzmann(x::Real)
    exp(x)
end

"""
$(TYPEDSIGNATURES)


The Blakemore approximation ``1/(\\exp(-x) + γ)`` with variable real scalar ``γ``, see 
[Blakemore1952, "The parameters of partially degenerate semiconductors"].

"""
function Blakemore(x::Real, γ::Real )
    1/(exp(-x) + γ)
end


"""
$(TYPEDSIGNATURES)


The Blakemore approximation ``1/(\\exp(-x) + γ)`` with
``γ = 0.27``, see 
[Blakemore1952, "The parameters of partially degenerate semiconductors"].
"""
function Blakemore(x::Real)
    Blakemore(x, 0.27)
end


"""
$(TYPEDSIGNATURES)


The Fermi-Dirac integral of order ``-1`` which reads 
``1/(\\exp(-x) + 1)``, see [Blakemore1982: "Approximations for Fermi-Dirac integrals, especially the function ``F_{1/2} (η)`` used to describe electron density in a semiconductor"].
"""
function FermiDiracMinusOne(x::Real)
    Blakemore(x, 1.0)
end


"""
$(TYPEDSIGNATURES)


The Fermi-Dirac integral of order ``0`` which reads 
``\\log(\\exp(x) + 1)``.
"""
function FermiDiracZero(x::Real)
    ex = exp(x)
    y  = 1+ex
    w  = y-1
    z  = w==0 ? ex : ex*log(y)/w
    return z
    #log( exp(x) + 1)
end


"""
$(TYPEDSIGNATURES)


The incomplete Fermi-Dirac integral of order 1/2, 
implemented according to [Bednarczyk1978, 
"The Approximation of the Fermi-Dirac integral ``F_{1/2}()``"].
"""
function FermiDiracOneHalfBednarczyk(x::Real)

    a = x^4 + 33.6*x * (1.0 - 0.68*exp(-0.17*(x+1)^2)) + 50
    sqrt(pi) / ( 2 * (3/4*sqrt(pi) * a^(-3/8) + exp(-x) ) )

end

"""
$(TYPEDSIGNATURES)


The incomplete Fermi-Dirac integral of order 1/2, 
implemented according to the software package TeSCA, see https://wias-berlin.de/software/index.jsp?lang=1&id=TeSCA.
"""
function FermiDiracOneHalfTeSCA(x::Real)
    if x < 1.6107
        ex = exp(x)
        y  = 1+ex
        w  = y-1
        z  = w==0 ? ex : ex*log(y)/w
#        z = log(1+ exp(x) )
        return ( 1 + 0.16 * z ) * z
    elseif 1.6107 <= x <= 344.7
        z = log( 1 + exp( x^(3/4)) )
        return 0.3258 - 0.0321 * z  + 0.7523 * z^2
    else
        z = x^(3/4)
        return 0.3258 - 0.0321 * z  + 0.7523 * z^2
    end

end


"""
$(TYPEDSIGNATURES)


Degenerate limit of incomplete Fermi-Dirac integral of order 1/2.
"""
function degenerateLimit(x)
        x < 0 ? NaN : 4/(3*sqrt(pi))*x^(3/2)
end

"""
$(TYPEDSIGNATURES)


Plot different distribution integrals.
"""
function plotDistributions(;Plotter=nothing)

    Plotter.close()

    rcParams = Plotter.PyDict(Plotter.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"

    # println(PyPlot.matplotlib.rcParams)


    x = -5:0.1:700;

    Plotter.semilogy(x, FermiDiracOneHalfBednarczyk.(x), label="\$F_{1/2}  \$ (Bednarczyk)");
    Plotter.semilogy(x, FermiDiracOneHalfTeSCA.(x), label="\$F_{1/2} \$ (TeSCA)");
    #Plotter.semilogy(x, Boltzmann.(x), label="Boltzmann");
    # Plotter.semilogy(x, ones(size(x))/0.27, "--", label="\$1/\\gamma=3.\\overline{703}\$", color=(0.6,0.6,0.6,1));
    # Plotter.semilogy(x, Blakemore.(x), label="Blakemore (\$\\gamma=0.27\$)");   
    # Plotter.semilogy(x, degenerateLimit.(x),label="degenerate limit");

    Plotter.xlabel("\$\\eta\$")
    Plotter.ylabel("\$\\mathcal{F}(\\eta)\$")
    Plotter.title("Distributions")
    Plotter.legend()

    Plotter.show();
end

"""
$(TYPEDSIGNATURES)


Plot diffusion enhancements.
"""
function plotDiffusionEnhancements()

    PyPlot.close()

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"

    x = -5:0.01:10;

    f = ChargeTransportInSolids.FermiDiracOneHalfBednarczyk; df = x -> ForwardDiff.derivative(f,x)
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



