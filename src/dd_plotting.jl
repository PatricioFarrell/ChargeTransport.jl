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

    # plot different doping values in interior
    for icc = 1:data.numberOfSpecies - 1 
        for i in 1:length(g.cellregions)

            # determine doping value in cell and number of cell nodes
            cellValue            = data.doping[g.cellregions[i],icc] 
            numberLocalCellNodes = length(g.cellnodes[:,i])

            # patch together cells
            PyPlot.plot(g.coord[g.cellnodes[:,i]], 
                        repeat(cellValue:cellValue,numberLocalCellNodes), 
                        color=colors[icc],
                        linewidth=3);
        end

        # legend
        PyPlot.plot(NaN,NaN,color=colors[icc],linewidth=3,label="icc="*string(icc))
    end

    # plot different doping values on boundary
    for icc = 1: data.numberOfSpecies - 1 
        for i in 1:length(g.bfaceregions)

            # determine doping value in cell and number of cell nodes
            cellValue            = data.bDoping[g.bfaceregions[i],icc] 
            numberLocalCellNodes = length(g.bfacenodes[:,i])

            # patch together cells
            PyPlot.plot(g.coord[g.bfacenodes[:,i]], 
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

Plot electrostatic potential, the electron and hole quasi Fermi potential as well as the IV curve.

"""
function plot_solution(sys, U0)
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