"""
$(SIGNATURES)

Plot densities of system.
"""
function plotDensities(grid, data, sol, bias)

    coord = grid[Coordinates]

    if length(coord[1]) != 1
        println("plotDensities is so far only implemented in 1D")
    end

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"
    colors = ["green", "red", "gold", "purple"]
    linestyles = ["-", ":", "--", "-."]
    densityNames  = ["n", "p", "a", "c"]

    densities = computeDensities(grid, data, sol)

    PyPlot.clf() 
    for icc = 1:data.numberOfSpecies-1
        PyPlot.semilogy(coord[1,:]./1, densities[icc,:], label = densityNames[icc], color = colors[icc], linewidth = 2, linestyle = "dashed")
    end
    PyPlot.grid()
    PyPlot.xlabel("space [\$ m \$]")
    PyPlot.ylabel("density [\$\\frac{1}{m^3}\$]")
    PyPlot.legend(fancybox = true, loc = "best")
    PyPlot.title("bias \$\\Delta u\$ = $bias")
    PyPlot.pause(0.00001)

end

"""
$(SIGNATURES)

Plot energies of system (physical variant).
"""
function plotEnergies(grid, data, sol, Δu)

    ipsi          = data.numberOfSpecies
    coord         = grid[Coordinates]

    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    colors = ["green", "red", "gold", "purple"]
    linestyles = ["-", ":", "--", "-."]
    labelBandEdgeEnergy = ["\$E_c-\\psi\$ ", "\$E_v-\\psi\$ ", "\$E_a-\\psi\$ ", "\$E_{cat}-\\psi\$ "]
    labelPotential = ["\$ - q \\varphi_n\$", "\$ - q \\varphi_p\$", "\$ - q \\varphi_a\$", "\$ - q \\varphi_c\$"]

    energies, fermiLevel = computeEnergies(grid, data, sol)


    for icc = 1:data.numberOfSpecies-1
        PyPlot.plot(coord[1,:]./1, energies[icc,:]./q,
                    label = labelBandEdgeEnergy[icc],
                    linewidth = 2,
                    color = colors[icc],
                    linestyle = linestyles[1])
        PyPlot.plot(coord[1,:]./1, fermiLevel[icc,:]./q,
                    label = labelPotential[icc],
                    linewidth = 2,
                    color = colors[icc],
                    linestyle = linestyles[2])
    end

    PyPlot.grid()
    PyPlot.xlabel("space [\$ m\$]")
    PyPlot.ylabel("energies [\$eV\$]")
    PyPlot.legend(fancybox = true, loc = "best")
    PyPlot.title("bias \$\\Delta u\$ = $Δu")
    PyPlot.pause(1.0e-5)

end



"""
$(SIGNATURES)
Plot band-edge energies.
"""

function plotEnergies(grid::ExtendableGrid, data::ChargeTransportData)
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

    colors = ["green", "red", "gold", "purple"]
    styles = ["-", ":", "--", "-."]
    densityNames  = ["\$ E_c\$", "\$ E_v \$", " \$ E_a \$", " \$ E_{cat}\$"]

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
    PyPlot.legend(fancybox = true, loc = "best")
    PyPlot.show();
    PyPlot.figure()

end


"""
$(SIGNATURES)
Visualize doping and bDoping (x) to make sure they agree.
"""
function plotDoping(g::ExtendableGrid, data::ChargeTransportData)
    coord  = g[Coordinates]
    if length(coord[1]) != 1
        println("plotDoping is so far only implemented in 1D")
    end

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"

    colors        = ["green", "red", "gold", "purple"]
    styles        = ["-",":", "--", "-."]
    densityNames  = ["n", "p", "a", "c"]
    

    # plot different doping values in interior
    cellregions = g[CellRegions]
    cellnodes   = g[CellNodes]
    coord       = g[Coordinates]
    for icc = 1:data.numberOfSpecies - 1
        for i in 1:length(cellregions)

            # determine doping value in cell and number of cell nodes
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
        PyPlot.semilogy(NaN,NaN,color=colors[icc],linewidth=3,label=densityNames[icc])
    end

    # plot different doping values on boundary
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
    PyPlot.legend(fancybox = true, loc = "best")

    PyPlot.show();
    PyPlot.figure()
end

"""
$(SIGNATURES)
Plot electroneutral potential.
"""
function plotElectroNeutralSolutionBoltzmann(grid, psi0)
    coord = grid[Coordinates]
    PyPlot.plot(coord[:],psi0, label = "electroneutral potential", color="g", marker="o")
    PyPlot.xlabel("space [m]")
    PyPlot.ylabel("potential [V]")
    PyPlot.legend(fancybox = true, loc = "best")
    PyPlot.show()
    PyPlot.figure()
end

"""
$(SIGNATURES)
Plot electrostatic potential as well as the electron and hole quasi-Fermi
potentials for fixed time and fixed boundary values.
"""

function plotSolution(coord, solution, Eref,  Δu)

    ipsi = size(solution)[1] # convention: psi is the last species

    colors        = ["green", "red", "gold", "purple"]
    linestyles    = ["--", "-.", "-", ":"]
    densityNames  = ["n", "p", "a", "c"]  
    PyPlot.clf() 

    PyPlot.plot(coord, (solution[ipsi,:] - Eref/q*ones(length(solution[ipsi,:]))), label = "\$\\psi\$", color="b")

    for icc in 1:ipsi-1
        PyPlot.plot(coord./1, solution[icc,:], label =  densityNames[icc], color= colors[icc], linestyle = linestyles[icc])
    end
            
    PyPlot.grid()
    PyPlot.xlabel("space [m]")
    PyPlot.ylabel("potential [V]")
    PyPlot.legend(fancybox = true, loc = "best")
    PyPlot.title("bias \$\\Delta u\$ = $Δu")
    PyPlot.gcf()

end


"""
$(SIGNATURES)
Plot electrostatic potential as well as the electron and hole quasi-Fermi potentials in stationary case.
"""
function plotSolution(coord, solution, Eref) # need to be dependent on Eref
    ipsi = size(solution)[1] # convention: psi is the last species
    
    colors        = ["green", "red", "yellow"]
    linestyles    = ["--", "-.", "-", ":"] 
    densityNames  = ["n", "p", "a", "c"]  
    PyPlot.clf()       
        
    PyPlot.plot(coord./1, solution[ipsi,:]-Eref/q*ones(length(solution[ipsi,:])), label = "\$\\psi\$", color="b")
                                                   
    for icc in 1:ipsi-1
    PyPlot.plot(coord./1, solution[icc,:], label = densityNames[icc], color= colors[icc], linestyle = linestyles[icc])
    end

    PyPlot.grid()
    PyPlot.xlabel("space [m]")
    PyPlot.ylabel("potential [V]")
    PyPlot.legend(fancybox = true, loc = "best")
    PyPlot.pause(1.0e-5)
end

"""
$(SIGNATURES)
Plot the IV curve.
"""
function plotIV(biasValues,IV, Δu)
    PyPlot.plot(biasValues[1:length(IV)], IV)
    PyPlot.grid()
    PyPlot.title("bias \$\\Delta u\$ = $Δu")
    PyPlot.xlabel("bias [V]")
    PyPlot.ylabel("total current [A]")
    PyPlot.pause(1.0e-5)
end
