"""
$(SIGNATURES)

Plot densities of system (with heterojunctions.)
Currently, only working for non-interfacial recombination.
"""
function plotDensities(Plotter, grid, data, sol, bias)
    Plotter.clf()

    if dim_space(grid) > 1
        println("ComputeDensities is so far only tested in 1D")
    end

    rcParams                    = Plotter.PyDict(Plotter.matplotlib."rcParams")
    rcParams["font.size"]       = 12
    rcParams["font.sans-serif"] = "Arial"
    colors                      = ["green", "red", "gold", "purple"]
    linestyles                  = ["-", ":", "--", "-."]
    densityNames                = ["n", "p", "a", "c"]

    ipsi                        = data.numberOfCarriers + 1

    cellnodes                   = grid[CellNodes]
    cellregions                 = grid[CellRegions]
    coordinates                 = grid[Coordinates]
    for icc in 1:data.numberOfCarriers

        # first cell
        u1                      = sol[:, 1]
        u2                      = sol[:, 2]
        ireg                    = cellregions[1]

        icc1                    = computeDensities(u1, data, 1, 1, icc, ipsi, false) # breg = 1 since we are on the left boundary
        icc2                    = computeDensities(u2, data, 2, ireg, icc, ipsi, true) 

        label_icc               = densityNames[icc]
        
        Plotter.semilogy([coordinates[1]./1, coordinates[2]./1], [icc1, icc2], label = label_icc, color = colors[icc], linewidth = 2) 

        for icell in 2:size(cellnodes,2) - 1
            in_region = true
            i1        = cellnodes[1,icell]
            i2        = cellnodes[2,icell]
            ireg      = cellregions[icell]
            node      = i1 

            u1        = sol[:, i1]
            u2        = sol[:, i2]
     
            icc1      = computeDensities(u1, data, i1, ireg, icc, ipsi, in_region)
            icc2      = computeDensities(u2, data, i2, ireg, icc, ipsi, in_region)
        
            Plotter.semilogy([coordinates[i1]./1, coordinates[i2]./1], [icc1, icc2],  color = colors[icc], linewidth = 2)      
        end

        # last cell
        u1            = sol[:, end-1]
        u2            = sol[:, end]
        ireg          = cellregions[end]
        node          = cellnodes[2, end]

        icc1          = computeDensities(u1, data, node-1, ireg, icc, ipsi, true)
        icc2          = computeDensities(u2, data, node, 2, icc, ipsi, false) # breg = 2 since we are on the right boundary

        Plotter.semilogy([coordinates[node-1]./1, coordinates[node]./1], [icc1, icc2], color = colors[icc], linewidth = 2) 
    end

    Plotter.grid()
    Plotter.xlabel("space [\$m\$]")
    Plotter.ylabel("density [\$\\frac{1}{m^3}\$]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.title("bias \$\\Delta u\$ = $bias")
    Plotter.pause(0.00001)

end

function plotDensitiesIKZ(Plotter, grid, data, sol, bias)
    Plotter.clf()

    if dim_space(grid) > 1
        println("ComputeDensities is so far only tested in 1D")
    end

    rcParams                    = Plotter.PyDict(Plotter.matplotlib."rcParams")
    rcParams["font.size"]       = 12
    rcParams["font.sans-serif"] = "Arial"
    colors                      = ["green", "red", "gold", "purple"]
    linestyles                  = ["-", ":", "--", "-."]
    densityNames                = ["n",  "p", "density (defects of \$ \\mathrm{Ti}_\\mathrm{Sr})\$)", "a", "c"]

    ipsi                        = data.numberOfCarriers + 1

    cellnodes                   = grid[CellNodes]
    cellregions                 = grid[CellRegions]
    coordinates                 = grid[Coordinates]
    for icc in 1:data.numberOfCarriers

        # first cell
        u1                      = sol[:, 1]
        u2                      = sol[:, 2]
        ireg                    = cellregions[1]

        icc1                    = computeDensities(u1, data, 1, 1, icc, ipsi, false) # breg = 1 since we are on the left boundary
        icc2                    = computeDensities(u2, data, 2, ireg, icc, ipsi, true) 

        label_icc               = densityNames[icc]
        
        Plotter.semilogy([coordinates[1]./1, coordinates[2]./1], [icc1, icc2], label = label_icc, color = colors[icc], linewidth = 2) 

        for icell in 2:size(cellnodes,2) - 1
            in_region = true
            i1        = cellnodes[1,icell]
            i2        = cellnodes[2,icell]
            ireg      = cellregions[icell]
            node      = i1 

            u1        = sol[:, i1]
            u2        = sol[:, i2]
     
            icc1      = computeDensities(u1, data, i1, ireg, icc, ipsi, in_region)
            icc2      = computeDensities(u2, data, i2, ireg, icc, ipsi, in_region)
        
            Plotter.semilogy([coordinates[i1]./1, coordinates[i2]./1], [icc1, icc2],  color = colors[icc], linewidth = 2)      
        end

        # last cell
        u1            = sol[:, end-1]
        u2            = sol[:, end]
        ireg          = cellregions[end]
        node          = cellnodes[2, end]

        icc1          = computeDensities(u1, data, node-1, ireg, icc, ipsi, true)
        icc2          = computeDensities(u2, data, node, 2, icc, ipsi, false) # breg = 2 since we are on the right boundary

        Plotter.semilogy([coordinates[node-1]./1, coordinates[node]./1], [icc1, icc2], color = colors[icc], linewidth = 2) 
    end

    Plotter.grid()
    Plotter.xlabel("space [\$m\$]")
    Plotter.ylabel("density [\$\\frac{1}{m^3}\$]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.title("bias \$\\Delta u\$ = $bias")
    Plotter.tight_layout()
    Plotter.pause(0.00001)

end


"""
$(SIGNATURES)

Plot energies of system (physical variant).
"""
function plotEnergies(Plotter, grid, data, sol, Δu)
    Plotter.clf()

    ipsi                = data.numberOfCarriers + 1

    cellnodes           = grid[CellNodes]
    cellregions         = grid[CellRegions]
    coord               = grid[Coordinates]

    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    colors              = ["green", "red", "gold", "purple"]
    linestyles          = ["-", ":", "--", "-."]
    labelBandEdgeEnergy = ["\$E_c-\\psi\$ ", "\$E_v-\\psi\$ ", "\$E_a-\\psi\$ ", "\$E_{cat}-\\psi\$ "]
    labelPotential      = ["\$ - q \\varphi_n\$", "\$ - q \\varphi_p\$", "\$ - q \\varphi_a\$", "\$ - q \\varphi_c\$"]

    for icc in 1:2
        # first cell
        ireg         = cellregions[1]
        E1           = data.bBandEdgeEnergy[icc, 1] + data.bandEdgeEnergyNode[icc, 1] # left boundary
        E2           = data.bandEdgeEnergy[icc, 1]  + data.bandEdgeEnergyNode[icc, 2] 
        energy_icc1  = E1 - q * sol[ipsi, 1]
        energy_icc2  = E2 - q * sol[ipsi, 2]
        label_energy = labelBandEdgeEnergy[icc]

        Plotter.plot([coord[1]./1, coord[2]./1], [energy_icc1, energy_icc2]./q, label = label_energy, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

        for icell in 2:size(cellnodes,2) - 1
            i1          = cellnodes[1,icell]
            i2          = cellnodes[2,icell]
            ireg        = cellregions[icell]

            E1          = data.bandEdgeEnergy[icc, ireg] + data.bandEdgeEnergyNode[icc, i1]
            E2          = data.bandEdgeEnergy[icc, ireg] + data.bandEdgeEnergyNode[icc, i2]

            energy_icc1 = E1 - q * sol[ipsi, i1]
            energy_icc2 = E2 - q * sol[ipsi, i2]

            Plotter.plot([coord[i1]./1, coord[i2]./1], [energy_icc1, energy_icc2]./q, linewidth = 2, color = colors[icc], linestyle = linestyles[1]) 
        end

        ireg        = cellregions[end]
        node        = cellnodes[2, end]
        E1          = data.bandEdgeEnergy[icc, ireg] + data.bandEdgeEnergyNode[icc, node-1] 
        E2          = data.bBandEdgeEnergy[icc, 2] + data.bandEdgeEnergyNode[icc, end] # right boundary
        energy_icc1 = E1 - q * sol[ipsi, end-1]
        energy_icc2 = E2 - q * sol[ipsi, end]

        Plotter.plot([coord[end-1]./1, coord[end]./1], [energy_icc1, energy_icc2]./q, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

        Plotter.plot(coord[1,:]./1, - sol[icc,:], label = labelPotential[icc], linewidth = 2, color = colors[icc], linestyle = linestyles[2])
   
   end
   
   Plotter.grid()
   Plotter.xlabel("space [\$m\$]")
   Plotter.ylabel("energies [\$eV\$]")
   Plotter.legend(fancybox = true, loc = "best")
   Plotter.tight_layout()
   Plotter.title("bias \$\\Delta u\$ = $Δu")
   Plotter.pause(1.0e-5)

end


function plotEnergiesIKZ(Plotter, grid, data, sol, Δu)
    Plotter.clf()

    ipsi                = data.numberOfCarriers + 1

    cellnodes           = grid[CellNodes]
    cellregions         = grid[CellRegions]
    coord               = grid[Coordinates]

    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    colors              = ["green", "red", "gold", "purple"]
    linestyles          = ["-", ":", "--", "-."]
    labelBandEdgeEnergy = ["\$E_c-\\psi\$ ", "\$E_v-\\psi\$ ", "\$E_a-\\psi\$ ", "\$E_{cat}-\\psi\$ "]
    labelPotential      = ["\$ - q \\varphi_n\$", "\$ - q \\varphi_p\$", "\$ - q \\varphi_a\$", "\$ - q \\varphi_c\$"]

    for icc = 1:data.numberOfCarriers
        # first cell
        ireg         = cellregions[1]
        E1           = data.bBandEdgeEnergy[icc, 1] + data.bandEdgeEnergyNode[icc, 1] # left boundary
        E2           = data.bandEdgeEnergy[icc, 1]  + data.bandEdgeEnergyNode[icc, 2] 
        energy_icc1  = E1 - q * sol[ipsi, 1]
        energy_icc2  = E2 - q * sol[ipsi, 2]
        label_energy = labelBandEdgeEnergy[icc]

        Plotter.plot([coord[1]./1, coord[2]./1], [energy_icc1, energy_icc2]./q , label = label_energy, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

        for icell in 2:size(cellnodes,2) - 1
            i1          = cellnodes[1,icell]
            i2          = cellnodes[2,icell]
            ireg        = cellregions[icell]

            E1          = data.bandEdgeEnergy[icc, ireg] + data.bandEdgeEnergyNode[icc, i1]
            E2          = data.bandEdgeEnergy[icc, ireg] + data.bandEdgeEnergyNode[icc, i2]

            energy_icc1 = E1 - q * sol[ipsi, i1]
            energy_icc2 = E2 - q * sol[ipsi, i2]

            Plotter.plot([coord[i1]./1, coord[i2]./1], [energy_icc1, energy_icc2]./q, linewidth = 2, color = colors[icc], linestyle = linestyles[1]) 
        end

        ireg        = cellregions[end]
        node        = cellnodes[2, end]
        E1          = data.bandEdgeEnergy[icc, ireg] + data.bandEdgeEnergyNode[icc, node-1] 
        E2          = data.bBandEdgeEnergy[icc, 2] + data.bandEdgeEnergyNode[icc, end] # right boundary
        energy_icc1 = E1 - q * sol[ipsi, end-1]
        energy_icc2 = E2 - q * sol[ipsi, end]

        Plotter.plot([coord[end-1]./1, coord[end]./1], [energy_icc1, energy_icc2]./q, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

        Plotter.plot(coord[1,:]./1, - sol[icc,:], label = labelPotential[icc], linewidth = 2, color = colors[icc], linestyle = linestyles[2])
    end
   
   Plotter.grid()
   Plotter.xlabel("space [\$m\$]")
   Plotter.ylabel("energies [\$eV\$]")
   Plotter.legend(fancybox = true, loc = "best")
   Plotter.title("bias \$\\Delta u\$ = $Δu")
   Plotter.tight_layout()
   Plotter.pause(1.0e-5)

end



"""
$(SIGNATURES)
Plot band-edge energies.
"""

function plotEnergies(Plotter, grid::ExtendableGrid, data)
    coord       = grid[Coordinates]
    cellregions = grid[CellRegions]
    cellnodes   = grid[CellNodes]

    #if length(coord[1]) != 1
    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    rcParams                    = Plotter.PyDict(Plotter.matplotlib."rcParams")
    rcParams["font.size"]       = 12
    rcParams["font.sans-serif"] = "Arial"
    colors                      = ["green", "red", "gold", "purple"]
    styles                      = ["-", ":", "--", "-."]
    densityNames                = ["\$ E_c\$", "\$ E_v \$", " \$ E_a \$", " \$ E_{cat}\$"]

    # plot different band-edge energies values in interior
    for icc = 1:data.numberOfCarriers
        for i in 1:length(cellregions)
            # determine band-edge energy value in cell and number of cell nodes
            cellValue            = ( data.bandEdgeEnergy[icc, cellregions[i]] + data.bandEdgeEnergyNode[icc, i] )/q
            numberLocalCellNodes = length(cellnodes[:,i])
            # patch together cells
            Plotter.plot(coord[cellnodes[:,i]],
                        repeat(cellValue:cellValue,numberLocalCellNodes),
                        marker="x",
                        color=colors[icc],
                        linewidth=3,
                        linestyle=styles[icc]);
        end

        Plotter.plot(NaN, NaN, color=colors[icc], linewidth = 3, label = "icc="*string(icc)) # legend
    end

    # plot different band-edge energy values on boundary
    bfaceregions = grid[BFaceRegions]
    bfacenodes   = grid[BFaceNodes]

    for icc = 1: data.numberOfCarriers

        for i in 1:length(bfaceregions)
            # determine band-edge energy value in cell and number of cell nodes
            cellValue            = (data.bBandEdgeEnergy[icc, bfaceregions[i]] + data.bandEdgeEnergyNode[icc, bfacenodes[i]])/q
            numberLocalCellNodes = length(bfacenodes[:,i])

            # patch together cells
            Plotter.plot(coord[bfacenodes[:,i]],
                        repeat(cellValue:cellValue,numberLocalCellNodes),
                        marker="x",
                        markersize=10,
                        color=colors[icc]);
        end

    end
    Plotter.xlabel("\$x\$")
    Plotter.title("band-edge energies")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.show();
    Plotter.figure()

end


"""
$(SIGNATURES)
Visualize doping and bDoping (x) to make sure they agree.
"""
function plotDoping(Plotter, g::ExtendableGrid, data)

    coord       = g[Coordinates]
    cellregions = g[CellRegions]
    cellnodes   = g[CellNodes]
    coord       = g[Coordinates]

    if length(coord[1]) != 1
        println("plotDoping is so far only implemented in 1D")
    end

    rcParams                    = Plotter.PyDict(Plotter.matplotlib."rcParams")
    rcParams["font.size"]       = 12
    rcParams["font.sans-serif"] = "Arial"
    colors                      = ["green", "red", "gold", "purple"]
    styles                      = ["-",":", "--", "-."]
    densityNames                = ["n", "p", "a", "c"]
    

    # plot different doping values in interior

    for icc = 1:data.numberOfCarriers

        for i in 1:length(cellregions)
            # determine doping value in cell and number of cell nodes
            cellValue            = data.doping[icc, cellregions[i]]
            numberLocalCellNodes = length(cellnodes[:,i])

            # patch together cells
            Plotter.semilogy(coord[cellnodes[:,i]],
                            repeat(cellValue:cellValue,numberLocalCellNodes),
                            color=colors[icc],
                            linewidth=3,
                            linestyle=styles[icc]);
        end
        Plotter.semilogy(NaN, NaN, color = colors[icc], linewidth = 3, label = densityNames[icc]) # legend

    end

    # plot different doping values on boundary
    bfaceregions = g[BFaceRegions]
    bfacenodes   = g[BFaceNodes]

    for icc = 1: data.numberOfCarriers

        for i in 1:length(bfaceregions)
            # determine doping value in cell and number of cell nodes
            cellValue            = data.bDoping[icc, bfaceregions[i]]
            numberLocalCellNodes = length(bfacenodes[:,i])

            # patch together cells
            Plotter.semilogy(coord[bfacenodes[:,i]],
                            repeat(cellValue:cellValue,numberLocalCellNodes),
                            marker="x",
                            markersize=10,
                            color=colors[icc]);
        end

    end
    Plotter.xlabel("\$x\$")
    Plotter.ylabel("\$N_{icc}\$")
    Plotter.title("Doping")
    Plotter.legend(fancybox = true, loc = "best")
    #Plotter.show();
end

"""
$(SIGNATURES)
Plot electroneutral potential.
"""

function plotElectroNeutralSolutionBoltzmann(Plotter, grid, psi0)
    coord = grid[Coordinates]
    Plotter.plot(coord[:],psi0, label = "electroneutral potential", color="g", marker="o")
    Plotter.xlabel("space [m]")
    Plotter.ylabel("potential [V]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.show()
end

"""
$(SIGNATURES)
Plot electrostatic potential as well as the electron and hole quasi-Fermi
potentials for fixed time and fixed boundary values.
"""

function plotSolution(Plotter, coord, solution, Eref,  Δu)

    ipsi = size(solution)[1] # convention: psi is the last species

    colors        = ["green", "red", "gold", "purple"]
    linestyles    = ["--", "-.", "-", ":"]
    densityNames  = ["\$\\varphi_n\$", "\$\\varphi_p\$", "\$\\varphi_a\$", "\$\\varphi_c\$"]  

    Plotter.clf() 
    Plotter.plot(coord, (solution[ipsi,:] + Eref/q*ones(length(solution[ipsi,:]))), label = "\$\\psi\$", color="b")

    for icc in 1:ipsi-1
        Plotter.plot(coord./1, solution[icc,:], label =  densityNames[icc], color= colors[icc], linestyle = linestyles[icc])
    end
            
    Plotter.grid()
    Plotter.xlabel("space [m]")
    Plotter.ylabel("potential [V]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.title("bias \$\\Delta u\$ = $Δu")
    Plotter.gcf()

end

function plotSolutionIKZ(Plotter, coord, solution, Eref,  Δu)

    ipsi = size(solution)[1] # convention: psi is the last species

    colors        = ["green", "red", "gold", "purple"]
    linestyles    = ["--", "-.", "-", ":"]
    densityNames  = ["\$\\varphi_n\$", "\$\\varphi_p\$", "qF potential Ti-antisite defect", "\$\\varphi_a\$", "\$\\varphi_c\$"]  

    Plotter.clf() 
    Plotter.plot(coord, (solution[ipsi,:] + Eref/q*ones(length(solution[ipsi,:]))), label = "\$\\psi\$", color="b")

    for icc in 1:ipsi-1
        Plotter.plot(coord./1, solution[icc,:] , label =  densityNames[icc], color= colors[icc], linestyle = linestyles[icc])
    end
            
    Plotter.grid()
    Plotter.xlabel("space [m]")
    Plotter.ylabel("potential [V]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.title("bias \$\\Delta u\$ = $Δu")
    Plotter.gcf()

end

"""
$(SIGNATURES)
Plot electrostatic potential as well as the electron and hole quasi-Fermi potentials in stationary case.
"""

function plotSolution(Plotter, coord, solution, Eref) # need to be dependent on Eref
    Plotter.clf()
    ipsi = size(solution)[1] # convention: psi is the last species
    
    colors        = ["green", "red", "yellow"]
    linestyles    = ["--", "-.", "-", ":"] 
    densityNames  = ["\$\\varphi_n\$", "\$\\varphi_p\$", "\$\\varphi_a\$", "\$\\varphi_c\$"]  
    Plotter.clf()       
        
    Plotter.plot(coord./1, solution[ipsi,:]-Eref/q*ones(length(solution[ipsi,:])), label = "\$\\psi\$", color="b")
                                                   
    for icc in 1:ipsi-1
        Plotter.plot(coord./1, solution[icc,:], label = densityNames[icc], color= colors[icc], linestyle = linestyles[icc])
    end

    Plotter.grid()
    Plotter.xlabel("space [m]")
    Plotter.ylabel("potential [V]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.pause(1.0e-5)
end

"""
$(SIGNATURES)
Plot the IV curve.
"""
function plotIV(Plotter, biasValues,IV, Δu)
    Plotter.plot(biasValues[1:length(IV)], IV)
    Plotter.grid()
    Plotter.title("bias \$\\Delta u\$ = $Δu")
    Plotter.xlabel("bias [V]")
    Plotter.ylabel("total current [A]")
    Plotter.pause(1.0e-5)
end
