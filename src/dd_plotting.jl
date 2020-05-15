"""
$(SIGNATURES)

Plot densities of system.
"""
function plotDensities(grid, sys, U0, Δu )

    dddata        = data(sys)
    ipsi          = dddata.numberOfSpecies

    coord         = grid[Coordinates]
    bfaceregions  = grid[BFaceRegions]
    bfacenodes    = grid[BFaceNodes]
    cellregions   = grid[CellRegions]
    cellnodes     = grid[CellNodes]
    numberOfCoord = length(coord)

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"

    colors = ["green", "red", "blue", "yellow"]
    linestyles = ["-", ":", "--", "-."]

    #if length(coord[1]) != 1
    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    densities = Array{Real,2}(undef, dddata.numberOfSpecies-1, length(coord))

    for icc = 1:dddata.numberOfSpecies-1

        E = dddata.bBandEdgeEnergy[bfaceregions[1],icc] + dddata.bandEdgeEnergyNode[bfacenodes[1],icc]
        eta = dddata.chargeNumbers[icc] / dddata.UT * ( (dddata.contactVoltage[1]- U0[ipsi, 1]) + E / q )
        densities[icc, 1] = dddata.bDensityOfStates[1, icc] * dddata.F(eta)

        for i = 2:numberOfCoord-1
            # frage: wie am besten etaF() benutzen? -> etaF() hat als Eingeparameter etwas vom Typ VoronoiFVM.Node
            E   = dddata.bandEdgeEnergy[cellregions[i], icc] + dddata.bandEdgeEnergyNode[i, icc]
            eta = dddata.chargeNumbers[icc] / dddata.UT * ( (U0[icc, i]- U0[ipsi, i]) + E / q )
            densities[icc, i] = dddata.densityOfStates[cellregions[i], icc] * dddata.F(eta)

        end
        E = dddata.bBandEdgeEnergy[bfaceregions[2],icc] + dddata.bandEdgeEnergyNode[bfacenodes[2],icc]
        eta = dddata.chargeNumbers[icc] / dddata.UT * ( (dddata.contactVoltage[2]- U0[ipsi, numberOfCoord]) + E / q )
        densities[icc, numberOfCoord] = dddata.bDensityOfStates[2, icc] * dddata.F(eta)

    end


    PyPlot.clf()
    for icc = 1:dddata.numberOfSpecies-1
        PyPlot.plot(coord[1,:]./μm, densities[icc,:], label = " density (icc = $icc)", color = colors[icc], linewidth = 3, linestyle = "dashed")
    end
    PyPlot.grid()
    PyPlot.xlabel("space [\$\\mu m \$]")
    PyPlot.ylabel("density [\$\\frac{1}{m^3}\$]")
    PyPlot.legend(fancybox = true, loc = "best")
    PyPlot.title("bias \$\\delta U\$ = $Δu")
    PyPlot.pause(1.0e-0)

end

"""
$(SIGNATURES)

Plot energies of system (physical variant).
"""
function plotEnergies(grid, sys, U0, Δu)

    dddata        = data(sys)
    ipsi          = dddata.numberOfSpecies

    coord         = grid[Coordinates]
    bfaceregions  = grid[BFaceRegions]
    bfacenodes    = grid[BFaceNodes]
    cellregions   = grid[CellRegions]
    cellnodes     = grid[CellNodes]
    numberOfCoord = length(coord)

    #if length(coord[1]) != 1
    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    energies   = Array{Real,2}(undef, dddata.numberOfSpecies-1, length(coord))
    fermiLevel = Array{Real,2}(undef, dddata.numberOfSpecies-1, length(coord))

    colors = ["green", "red", "blue", "yellow"]
    linestyles = ["-", ":", "--", "-."]
# warum das - bei den Fermi level?
    for icc = 1:dddata.numberOfSpecies-1
        E = dddata.bBandEdgeEnergy[bfaceregions[1],icc] + dddata.bandEdgeEnergyNode[bfacenodes[1],icc]
        energies[icc, 1]   = E - q * U0[ipsi, 1]
        fermiLevel[icc, 1] = - q * U0[icc, 1]

        for i = 2:numberOfCoord-1
            # frage: wie am besten etaF() benutzen? -> etaF() hat als Eingeparameter etwas vom Typ VoronoiFVM.Node
            E   = dddata.bandEdgeEnergy[cellregions[i], icc] + dddata.bandEdgeEnergyNode[i, icc]
            energies[icc, i]   = E - q *U0[ipsi, i]
            fermiLevel[icc, i] = -q* U0[icc, i]

        end
        E = dddata.bBandEdgeEnergy[bfaceregions[2],icc] + dddata.bandEdgeEnergyNode[bfacenodes[2],icc]
        energies[icc, numberOfCoord]   = E - q * U0[ipsi, numberOfCoord]
        fermiLevel[icc, numberOfCoord] = - q * U0[icc, numberOfCoord]
    end


    PyPlot.clf()
    for icc = 1:dddata.numberOfSpecies-1
        PyPlot.plot(coord[1,:]./μm, energies[icc,:],
                    label = "\$E_i-\\psi\$ (icc = $icc)",
                    linewidth = 3,
                    color = colors[icc],
                    linestyle = linestyles[1])
        PyPlot.plot(coord[1,:]./μm, fermiLevel[icc,:],
                    label = "\$ q \\varphi_i\$ (icc = $icc)",
                    linewidth = 3,
                    color = colors[icc],
                    linestyle = linestyles[2])
    end

    PyPlot.grid()
    PyPlot.xlabel("space [\$\\mu m\$]")
    PyPlot.ylabel("Energy[\$eV\$]")
    PyPlot.legend(fancybox = true, loc = "best")
    PyPlot.title("bias \$\\delta U\$ = $Δu")
    PyPlot.pause(1.0e-0)

end



"""
$(SIGNATURES)
Plot band-edge energies.
"""

function plotEnergies(grid::ExtendableGrid, data::DDFermiData)
    coord       = grid[Coordinates]
    cellregions = grid[CellRegions]
    cellnodes   = grid[CellNodes]

    #if length(coord[1]) != 1
    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"

    colors = ["green", "red", "blue", "yellow"]
    styles = ["-", ":", "--", "-."]

    # plot different band-edge energies values in interior
    for icc = 1:data.numberOfSpecies - 1
        for i in 1:length(cellregions)
            # determine band-edge energy value in cell and number of cell nodes
            cellValue            = ( data.bandEdgeEnergy[cellregions[i],icc] + data.bandEdgeEnergyNode[i,icc] )/q
            numberLocalCellNodes = length(cellnodes[:,i])
            # patch together cells
            PyPlot.plot(coord[cellnodes[:,i]],
            repeat(cellValue:cellValue,numberLocalCellNodes),
            marker="x",
            color=colors[icc],
            linewidth=3,
            linestyle=styles[icc]);
        end

        # legend
        PyPlot.plot(NaN,NaN,color=colors[icc],linewidth=3,label="icc="*string(icc))
    end

    # plot different band-edge energy values on boundary
    bfaceregions = grid[BFaceRegions]
    bfacenodes = grid[BFaceNodes]
    for icc = 1: data.numberOfSpecies - 1
        for i in 1:length(bfaceregions)

            # determine band-edge energy value in cell and number of cell nodes
            cellValue            = (data.bBandEdgeEnergy[bfaceregions[i],icc] + data.bandEdgeEnergyNode[bfacenodes[i],icc])/q
            numberLocalCellNodes = length(bfacenodes[:,i])
            # patch together cells
            PyPlot.plot(coord[bfacenodes[:,i]],
            marker="x",
            markersize=10,
            repeat(cellValue:cellValue,numberLocalCellNodes),
            color=colors[icc]);
        end

    end

    PyPlot.xlabel("\$x\$")
    PyPlot.title("band-edge Energies")
    PyPlot.legend()
    PyPlot.show();
    PyPlot.figure()

end


"""
$(SIGNATURES)
Visualize doping and bDoping (x) to make sure they agree.
"""
function plotDoping(g::ExtendableGrid, data::DDFermiData)
    #todo_da: add following line
    coord  = g[Coordinates]
    #if length(g.coord[1]) != 1
    if length(coord[1]) != 1
        println("plotDoping is so far only implemented in 1D")
    end

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"

    colors = ["green", "red", "blue", "yellow"]
    styles = ["-",":", "--", "-."]

    # plot different doping values in interior
    #todo_da: add following line and delete dependency of cellregions from g
    cellregions = g[CellRegions]
    cellnodes = g[CellNodes]
    coord = g[Coordinates]
    for icc = 1:data.numberOfSpecies - 1
        for i in 1:length(cellregions)

            # determine doping value in cell and number of cell nodes
            #todo_da: deleted dependency on g
            cellValue            = data.doping[cellregions[i],icc]
            numberLocalCellNodes = length(cellnodes[:,i])

            # patch together cells
            PyPlot.semilogy(coord[cellnodes[:,i]],
            repeat(cellValue:cellValue,numberLocalCellNodes),
            color=colors[icc],
            linewidth=3,
            linestyle=styles[icc]);
        end

        # legend
        PyPlot.semilogy(NaN,NaN,color=colors[icc],linewidth=3,label="icc="*string(icc))
    end

    # plot different doping values on boundary
    #todo_da: following two lines added and delete dependency on g
    bfaceregions = g[BFaceRegions]
    bfacenodes = g[BFaceNodes]
    for icc = 1: data.numberOfSpecies - 1
        for i in 1:length(bfaceregions)

            # determine doping value in cell and number of cell nodes
            cellValue            = data.bDoping[bfaceregions[i],icc]
            numberLocalCellNodes = length(bfacenodes[:,i])

            # patch together cells
            PyPlot.semilogy(coord[bfacenodes[:,i]],
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
    savefig("pin-doping.eps")
    PyPlot.figure()
end

"""
$(SIGNATURES)
Plot electroneutral potential.
"""
function plotElectroNeutralSolutionBoltzmann(grid, psi0)
    #todo_da: add following line
    coord = grid[Coordinates]
    PyPlot.plot(coord[:],psi0, label = "electroneutral potential (Boltzmann)", color="g", marker="o")
    PyPlot.xlabel("space [m]")
    PyPlot.ylabel("potential [V]")
    PyPlot.legend(loc="upper left")
    PyPlot.show()
    PyPlot.figure()
end

"""
$(SIGNATURES)
Plot electrostatic potential as well as the electron and hole quasi-Fermi
potentials for fixed time and fixed boundary values.
    """

    function plotSolution(grid, sys, U0, Δu, time)
        dddata = VoronoiFVM.data(sys)
        coord  = grid[Coordinates]

        PyPlot.clf()
        @views begin
            PyPlot.plot(coord[1,:], U0[3,:], label = "electrostatic potential", color="g", marker="o")
            PyPlot.plot(coord[1,:], U0[1,:], label = "quasi-Fermi electron", color="b", marker="o", linestyle = "dashed")
            PyPlot.plot(coord[1,:], U0[2,:], label = "quasi-Fermi hole", color="r", marker="o", linestyle = "dashdot")
            PyPlot.grid()
            PyPlot.xlabel("space [m]")
            PyPlot.ylabel("potential [V]")
            PyPlot.legend(loc="upper left")
            PyPlot.title("applied bias = $Δu [V] and time = $time [s]")
            PyPlot.gcf()
        end

    end


    """
    $(SIGNATURES)
    Plot electrostatic potential as well as the electron and hole quasi-Fermi potentials in stationary case.
    """
    #todo_da:add dependency on grid
    function plotSolution(grid, sys, U0)
        dddata = VoronoiFVM.data(sys)
        coord = grid[Coordinates]

        PyPlot.clf()
        @views begin
            PyPlot.subplot(211)
            #todo_da changed sys.grid.coord into coord
            PyPlot.plot(coord[1,:], U0[3,:], label = "electrostatic potential", color="g", marker="o")
            PyPlot.plot(coord[1,:], U0[1,:], label = "quasi-Fermi electron", color="b", marker="o", linestyle = "dashed")
            PyPlot.plot(coord[1,:], U0[2,:], label = "quasi-Fermi hole", color="r", marker="o", linestyle = "dashdot")
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
        #PyPlot.subplot(212)
        PyPlot.plot(biasValues[1:length(IV)], IV)
        PyPlot.grid()
        PyPlot.xlabel("bias [V]")
        PyPlot.ylabel("total current [A]")
        PyPlot.pause(1.0e-5)
    end
