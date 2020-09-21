"""
$(SIGNATURES)

Plot densities of system (with heterojunctions.)
Currently, only working for non-interfacial recombination.
"""
function plotDensities(grid, data, sol, bias)
    PyPlot.clf()

    if dim_space(grid) > 1
        println("ComputeDensities is so far only tested in 1D")
    end

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12
    rcParams["font.sans-serif"] = "Arial"
    colors = ["green", "red", "gold", "purple"]
    linestyles = ["-", ":", "--", "-."]
    densityNames  = ["n", "p", "a", "c"]

    ipsi      = data.numberOfSpecies

    cellnodes   = grid[CellNodes]
    cellregions = grid[CellRegions]
    coordinates = grid[Coordinates]
    for icc in 1:data.numberOfSpecies - 1

        # first cell
        u1    = sol[:, 1]
        u2    = sol[:, 2]
        ireg = cellregions[1]

        icc1 = computeDensities(u1, data, 1, 1, icc, ipsi, false) # breg = 1 since we are on the left boundary
        icc2 = computeDensities(u2, data, 2, ireg, icc, ipsi, true) 

        label_icc = densityNames[icc]
        
        PyPlot.semilogy([coordinates[1]./1, coordinates[2]./1], [icc1, icc2], label = label_icc, color = colors[icc], linewidth = 2) 

        for icell in 2:size(cellnodes,2) - 1
            in_region = true
            i1   = cellnodes[1,icell]
            i2   = cellnodes[2,icell]
            ireg = cellregions[icell]
            node = i1 

            u1    = sol[:, i1]
            u2    = sol[:, i2]
     
            icc1 = computeDensities(u1, data, i1, ireg, icc, ipsi, in_region)
            icc2 = computeDensities(u2, data, i2, ireg, icc, ipsi, in_region)
        
            PyPlot.semilogy([coordinates[i1]./1, coordinates[i2]./1], [icc1, icc2],  color = colors[icc], linewidth = 2) 
            
        end

        # last cell
        u1    = sol[:, end-1]
        u2    = sol[:, end]
        ireg = cellregions[end]
        node = cellnodes[2, end]

        icc1 = computeDensities(u1, data, node-1, ireg, icc, ipsi, true)
        icc2 = computeDensities(u2, data, node, 2, icc, ipsi, false) # breg = 2 since we are on the right boundary

        PyPlot.semilogy([coordinates[node-1]./1, coordinates[node]./1], [icc1, icc2], color = colors[icc], linewidth = 2) 
    end

    PyPlot.grid()
    PyPlot.xlabel("space [\$ m \$]")
    PyPlot.ylabel("density [\$\\frac{1}{m^3}\$]")
    PyPlot.legend(fancybox = true, loc = "best")
    #PyPlot.ylim((0.0, 1.0e24))
    PyPlot.title("bias \$\\Delta u\$ = $bias")
    PyPlot.pause(0.00001)

end


"""
$(SIGNATURES)

Plot energies of system (physical variant).
"""
function plotEnergies(grid, data, sol, Δu)
    PyPlot.clf()

    ipsi          = data.numberOfSpecies

    cellnodes   = grid[CellNodes]
    cellregions = grid[CellRegions]
    coord       = grid[Coordinates]

    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    colors = ["green", "red", "gold", "purple"]
    linestyles = ["-", ":", "--", "-."]
    labelBandEdgeEnergy = ["\$E_c-\\psi\$ ", "\$E_v-\\psi\$ ", "\$E_a-\\psi\$ ", "\$E_{cat}-\\psi\$ "]
    labelPotential = ["\$ - q \\varphi_n\$", "\$ - q \\varphi_p\$", "\$ - q \\varphi_a\$", "\$ - q \\varphi_c\$"]

    for icc in 1:data.numberOfSpecies - 1

        # first cell
        ireg = cellregions[1]

        E1          = data.bBandEdgeEnergy[1, icc] + data.bandEdgeEnergyNode[1, icc] # left boundary
        E2          = data.bandEdgeEnergy[1, icc] + data.bandEdgeEnergyNode[2, icc] 
        energy_icc1 = E1 - q * sol[ipsi, 1]
        energy_icc2 = E2 - q * sol[ipsi, 2]

        label_energy = labelBandEdgeEnergy[icc]
        PyPlot.plot([coord[1]./1, coord[2]./1], [energy_icc1, energy_icc2]./q, label = label_energy, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

        for icell in 2:size(cellnodes,2) - 1

            i1   = cellnodes[1,icell]
            i2   = cellnodes[2,icell]
            ireg = cellregions[icell]

            E1    = data.bandEdgeEnergy[ireg, icc] + data.bandEdgeEnergyNode[i1, icc]
            E2    = data.bandEdgeEnergy[ireg, icc] + data.bandEdgeEnergyNode[i2, icc]

            energy_icc1 = E1 - q * sol[ipsi, i1]
            energy_icc2 = E2 - q * sol[ipsi, i2]

            PyPlot.plot([coord[i1]./1, coord[i2]./1], [energy_icc1, energy_icc2]./q, linewidth = 2, color = colors[icc], linestyle = linestyles[1]) 
        end

        ireg = cellregions[end]
        node = cellnodes[2, end]

        E1          = data.bandEdgeEnergy[ireg, icc] + data.bandEdgeEnergyNode[node-1, icc] 
        E2          = data.bBandEdgeEnergy[2, icc] + data.bandEdgeEnergyNode[end, icc] # right boundary
        energy_icc1 = E1 - q * sol[ipsi, end-1]
        energy_icc2 = E2 - q * sol[ipsi, end]

        PyPlot.plot([coord[end-1]./1, coord[end]./1], [energy_icc1, energy_icc2]./q, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

        PyPlot.plot(coord[1,:]./1, - sol[icc,:], label = labelPotential[icc], linewidth = 2, color = colors[icc], linestyle = linestyles[2])
   
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
    densityNames  = ["\$\\varphi_n\$", "\$\\varphi_p\$", "\$\\varphi_a\$", "\$\\varphi_c\$"]  
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
    PyPlot.clf()
    ipsi = size(solution)[1] # convention: psi is the last species
    
    colors        = ["green", "red", "yellow"]
    linestyles    = ["--", "-.", "-", ":"] 
    densityNames  = ["\$\\varphi_n\$", "\$\\varphi_p\$", "\$\\varphi_a\$", "\$\\varphi_c\$"]  
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
    PyPlot.clf()
    PyPlot.plot(biasValues[1:length(IV)], IV)
    PyPlot.grid()
    PyPlot.title("bias \$\\Delta u\$ = $Δu")
    PyPlot.xlabel("bias [V]")
    PyPlot.ylabel("total current [A]")
    PyPlot.pause(1.0e-5)
end
