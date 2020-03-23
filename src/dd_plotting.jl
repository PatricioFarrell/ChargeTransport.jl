"""
$(SIGNATURES)

Visualize doping and bDoping (x) to make sure they agree.
"""
function plotDoping(g::VoronoiFVM.AbstractGrid, data::DDFermiData)

    if length(g.coord[1]) != 1
        println("plotDoping is so far only implemented in 1D")
    end

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"

    colors = ["green", "red", "blue", "yellow"]
    styles = ["-", "--", "-.", ":"]

    # plot different doping values in interior
    for icc = 1:data.numberOfSpecies - 1 
        for i in 1:length(g.cellregions)

            # determine doping value in cell and number of cell nodes
            cellValue            = data.doping[g.cellregions[i],icc] 
            numberLocalCellNodes = length(g.cellnodes[:,i])

            # patch together cells
            PyPlot.semilogy(g.coord[g.cellnodes[:,i]], 
                        repeat(cellValue:cellValue,numberLocalCellNodes), 
                        color=colors[icc],
                        linewidth=3,
                        linestyle=styles[icc]);
        end

        # legend
        PyPlot.semilogy(NaN,NaN,color=colors[icc],linewidth=3,label="icc="*string(icc))
    end

    # plot different doping values on boundary
    for icc = 1: data.numberOfSpecies - 1 
        for i in 1:length(g.bfaceregions)

            # determine doping value in cell and number of cell nodes
            cellValue            = data.bDoping[g.bfaceregions[i],icc] 
            numberLocalCellNodes = length(g.bfacenodes[:,i])

            # patch together cells
            PyPlot.semilogy(g.coord[g.bfacenodes[:,i]], 
                        marker="x",
                        markersize=10,
                        repeat(cellValue:cellValue,numberLocalCellNodes), 
                        color=colors[icc]);
        end

    end

    PyPlot.xlabel("\$x\$")
    PyPlot.ylabel("\$N_{icc}\$")
    PyPlot.title("Doping")
    PyPlot.legend()

    PyPlot.show();
    PyPlot.figure()
end

"""
$(SIGNATURES)

Plot electroneutral potential.
"""
function plotElectroNeutralSolutionBoltzmann(grid, psi0)
        PyPlot.plot(grid.coord[:],psi0, label = "electroneutral potential (Boltzmann)", color="g", marker="o")
        PyPlot.xlabel("space [m]")
        PyPlot.ylabel("potential [V]")
        PyPlot.legend(loc="upper left")
        PyPlot.show()
        PyPlot.figure()
end


"""
$(SIGNATURES)

Plot electrostatic potential as well as the electron and hole quasi Fermi potentials.

"""
function plotSolution(sys, U0)
    dddata = VoronoiFVM.data(sys)

    PyPlot.clf()
    @views begin
        PyPlot.subplot(211)
        PyPlot.plot(sys.grid.coord[1,:], U0[3,:], label = "electrostatic potential", color="g", marker="o")
        PyPlot.plot(sys.grid.coord[1,:], U0[1,:], label = "quasi Fermi electron", color="b", marker="o", linestyle = "dashed")
        PyPlot.plot(sys.grid.coord[1,:], U0[2,:], label = "quasi Fermi hole", color="r", marker="o", linestyle = "dashdot")
        PyPlot.grid()
        PyPlot.xlabel("space [m]")
        PyPlot.ylabel("potential [V]")
        PyPlot.legend(loc="upper left")
        PyPlot.gcf()
    end

end

"""
$(SIGNATURES)

Plot the IV curve.

"""
function plotIV(biasValues,IV)
    PyPlot.subplot(212)
    PyPlot.plot(biasValues[1:length(IV)], IV) 
    PyPlot.grid()
    PyPlot.xlabel("bias [V]")
    PyPlot.ylabel("total current [A]")
    PyPlot.pause(1.0e-5)
end