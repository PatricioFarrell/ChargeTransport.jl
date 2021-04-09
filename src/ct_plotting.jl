"""
$(TYPEDSIGNATURES)
Plotting routine, where the charge carrier densities are depicted
in dependence of space. The case of heterojunctions is tested, but yet
multidimensional calculations are missing.
Currently, for a matching legend, we need the following order:

index 1: electrons as charge carrier with the corresponding density ``n``,

index 2: holes as charge carrier with the corresponding density ``p``,

index 3: anion vacancies as charge carrier with the corresponding density ``a``,

index 4: cation vacancies as charge carrier with the corresponding density ``c``.

One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plotDensities(Plotter, grid, data, sol, title, ;plotGridpoints=false)
    Plotter.clf()

    if dim_space(grid) > 1
        error("ComputeDensities is so far only tested in 1D")
    end

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

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

        Plotter.semilogy([coordinates[1]./1, coordinates[2]./1], 1.0e-6 .*[icc1, icc2], marker = marker, label = label_icc, color = colors[icc], linewidth = 2) #multiplying by 1.0e-6 gives us the densities in cm^(-3)

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
        
            Plotter.semilogy([coordinates[i1]./1, coordinates[i2]./1], 1.0e-6 .*[icc1, icc2], marker = marker, color = colors[icc], linewidth = 2) #multiplying by 1.0e-6 gives us the densities in cm^(-3)     
        end

        # last cell
        u1            = sol[:, end-1]
        u2            = sol[:, end]
        ireg          = cellregions[end]
        node          = cellnodes[2, end]

        icc1          = computeDensities(u1, data, node-1, ireg, icc, ipsi, true)
        icc2          = computeDensities(u2, data, node, 2, icc, ipsi, false) # breg = 2 since we are on the right boundary

        Plotter.semilogy([coordinates[node-1]./1, coordinates[node]./1], 1.0e-6 .*[icc1, icc2], marker = marker, color = colors[icc], linewidth = 2) #multiplying by 1.0e-6 gives us the densities in cm^(-3)
    end

    Plotter.grid()
    Plotter.xlabel("space [\$m\$]")
    Plotter.ylabel("density [\$\\frac{1}{cm^3}\$]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.title(title)
    Plotter.tight_layout()
    Plotter.pause(0.001)

end

"""
$(TYPEDSIGNATURES)

With this method it is possible to plot the energies

``E_\\alpha - q \\psi \\quad \\text{w.r.t. space.}``

The case of heterojunctions is tested, but yet
multidimensional calculations are missing. 
Currently, for a matching legend, we need the following order:

index 1: electrons as charge carrier with the corresponding density ``n``,

index 2: holes as charge carrier with the corresponding density ``p``,

index 3: anion vacancies as charge carrier with the corresponding density ``a``,

index 4: cation vacancies as charge carrier with the corresponding density ``c``.

One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plotEnergies(Plotter, grid, data, sol, title, ;plotGridpoints=false)
    Plotter.clf()

    ipsi                = data.numberOfCarriers + 1

    cellnodes           = grid[CellNodes]
    cellregions         = grid[CellRegions]
    coord               = grid[Coordinates]

    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    colors              = ["green", "red", "gold", "purple"]
    linestyles          = ["-", ":", "--", "-."]
    labelBandEdgeEnergy = ["\$E_c-q\\psi\$ ", "\$E_v-q\\psi\$ ", "\$E_a-q\\psi\$ ", "\$E_{cat}-q\\psi\$ "]
    labelPotential      = ["\$ - q \\varphi_n\$", "\$ - q \\varphi_p\$", "\$ - q \\varphi_a\$", "\$ - q \\varphi_c\$"]

    for icc in 1:2
        # first cell
        ireg         = cellregions[1]
        E1           = data.bBandEdgeEnergy[icc, 1] + data.bandEdgeEnergyNode[icc, 1] # left boundary
        E2           = data.bandEdgeEnergy[icc, 1]  + data.bandEdgeEnergyNode[icc, 2] 
        energy_icc1  = E1 - q * sol[ipsi, 1]
        energy_icc2  = E2 - q * sol[ipsi, 2]
        label_energy = labelBandEdgeEnergy[icc]

        Plotter.plot([coord[1]./1, coord[2]./1], [energy_icc1, energy_icc2]./q, marker = marker, label = label_energy, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

        for icell in 2:size(cellnodes,2) - 1
            i1          = cellnodes[1,icell]
            i2          = cellnodes[2,icell]
            ireg        = cellregions[icell]

            E1          = data.bandEdgeEnergy[icc, ireg] + data.bandEdgeEnergyNode[icc, i1]
            E2          = data.bandEdgeEnergy[icc, ireg] + data.bandEdgeEnergyNode[icc, i2]

            energy_icc1 = E1 - q * sol[ipsi, i1]
            energy_icc2 = E2 - q * sol[ipsi, i2]

            Plotter.plot([coord[i1]./1, coord[i2]./1], [energy_icc1, energy_icc2]./q, marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[1]) 
        end

        ireg        = cellregions[end]
        node        = cellnodes[2, end]
        E1          = data.bandEdgeEnergy[icc, ireg] + data.bandEdgeEnergyNode[icc, node-1] 
        E2          = data.bBandEdgeEnergy[icc, 2] + data.bandEdgeEnergyNode[icc, end] # right boundary
        energy_icc1 = E1 - q * sol[ipsi, end-1]
        energy_icc2 = E2 - q * sol[ipsi, end]

        Plotter.plot([coord[end-1]./1, coord[end]./1], [energy_icc1, energy_icc2]./q, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

        Plotter.plot(coord[1,:]./1, - sol[icc,:], label = labelPotential[icc], marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[2])
   
   end
   
   Plotter.grid()
   Plotter.xlabel("space [\$m\$]")
   Plotter.ylabel("energies [\$eV\$]")
   Plotter.legend(fancybox = true, loc = "best")
   Plotter.title(title)
   Plotter.tight_layout()
   Plotter.pause(1.0e-5)

end

"""
$(SIGNATURES)
With this method it is possible to depict the band-edge energies ``E_\\alpha ``. 
This can be useful for debugging when dealing with heterojunctions.
"""
function plotEnergies(Plotter, grid::ExtendableGrid, data)
    coord       = grid[Coordinates]
    cellregions = grid[CellRegions]
    cellnodes   = grid[CellNodes]

    #if length(coord[1]) != 1
    if length(coord[1]) != 1
        error("plotEnergies is so far only implemented in 1D")
    end

    colors      = ["green", "red", "gold", "purple"]
    styles      = ["-", ":", "--", "-."]
    EnergyNames = ["\$ E_c\$", "\$ E_v \$", " \$ E_a \$", " \$ E_{cat}\$"]

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

        Plotter.plot(NaN, NaN, color=colors[icc], linewidth = 3, linestyle = styles[icc], label = EnergyNames[icc]) # legend
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
    Plotter.xlabel("space [\$m\$]")
    Plotter.ylabel("energies [\$eV\$]")
    Plotter.title("Band-edge energies \$ E_α\$")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.tight_layout()
    Plotter.show()

end


"""
$(TYPEDSIGNATURES)
Possibility to plot the considered doping. This is especially useful 
for making sure that the interior and the boundary doping agree.
"""
function plotDoping(Plotter, g::ExtendableGrid, data)

    coord       = g[Coordinates]
    cellregions = g[CellRegions]
    cellnodes   = g[CellNodes]
    coord       = g[Coordinates]

    if length(coord[1]) != 1
        error("plotDoping is so far only implemented in 1D")
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
            Plotter.plot(coord[cellnodes[:,i]],
            1.0e-6 .*repeat(cellValue:cellValue,numberLocalCellNodes),
                            color=colors[icc],
                            linewidth=3,
                            linestyle=styles[icc]); #multiplying by 1.0e-6 gives us the densities in cm^(-3) 
        end
        Plotter.plot(NaN, NaN, color = colors[icc], linewidth = 3, label = densityNames[icc]) # legend

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
            Plotter.plot(coord[bfacenodes[:,i]],
            1.0e-6 .*repeat(cellValue:cellValue,numberLocalCellNodes),
                            marker="x",
                            markersize=10,
                            color=colors[icc]);
        end

    end
    Plotter.yscale("symlog")
    Plotter.xlabel("space [\$m\$]")
    Plotter.ylabel("Doping [\$\\frac{1}{cm^3}\$]")
    Plotter.title("Doping values for charge carriers")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.tight_layout()
    #Plotter.show();
end

"""
$(TYPEDSIGNATURES)
Plotting routine for depicting the electroneutral potential.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plotElectroNeutralSolutionBoltzmann(Plotter, grid, psi0, ;plotGridpoints=false)

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    coord = grid[Coordinates]
    Plotter.plot(coord[:],psi0, label = "electroneutral potential \$ ψ_0 \$", color="b", marker= marker)
    Plotter.xlabel("space [m]")
    Plotter.ylabel("potential [V]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.tight_layout()
    Plotter.show()
end

"""
$(TYPEDSIGNATURES)
Method for plotting the solution vectors: the electrostatic potential ``\\psi``
as well as the charge carriers.
The case of heterojunctions is tested, but yet
multidimensional calculations are missing. 
Currently, for a matching legend, we need the following order:

index 1: electrons as charge carrier with the corresponding density ``n``,

index 2: holes as charge carrier with the corresponding density ``p``,

index 3: anion vacancies as charge carrier with the corresponding density ``a``,

index 4: cation vacancies as charge carrier with the corresponding density ``c``.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plotSolution(Plotter, coord, solution, Eref, title, ;plotGridpoints=false)

    if size(solution)[1] > 4
        ipsi = 4
    else
        ipsi = size(solution)[1] # convention: psi is the last species
    end

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end


    colors        = ["green", "red", "gold", "purple"]
    linestyles    = ["--", "-.", "-", ":"]
    densityNames  = ["\$\\varphi_n\$", "\$\\varphi_p\$", "\$\\varphi_a\$", "\$\\varphi_c\$"]  

    Plotter.clf() 
    Plotter.plot(coord, (solution[ipsi,:] + Eref/q*ones(length(solution[ipsi,:]))), marker = marker, label = "\$\\psi\$", color="b")

    for icc in 1:ipsi-1
        Plotter.plot(coord./1, solution[icc,:], label =  densityNames[icc], marker = marker, color= colors[icc], linestyle = linestyles[icc])
    end
            
    Plotter.grid()
    Plotter.xlabel("space [m]")
    Plotter.ylabel("potential [V]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.title(title)
    Plotter.tight_layout()
    Plotter.gcf()

end

function plotSolution(Plotter, grid, solution, Eref, agrid, t, Δu)


    # Create a visualizer. Works with Plots (fast once compiled) and PyPlot
    visualizer = p = GridVisualizer(Plotter = Plotter, layout = (1,1) )

    ipsi = size(solution)[1] # convention: psi is the last species

    colors        = ["green", "red", "gold", "purple"]
    linestyles    = ["--", "-.", "-", ":"]
    densityNames  = ["\$\\varphi_n\$", "\$\\varphi_p\$", "\$\\varphi_a\$", "\$\\varphi_c\$"]  
    
    Plotter.clf() 
    scalarplot!(p[1,1], grid, (solution[ipsi,:] + Eref/q*ones(length(solution[ipsi,:]))), label = "\$\\psi\$", color="b",  marker = "x", 
    title="time \$ t =\$ $t, bias \$\\Delta u\$ = $Δu", clear = true)

    for icc in 1:2
        scalarplot!(p[1,1], grid, solution[icc,:], label =  densityNames[icc], color= colors[icc], linestyle = linestyles[icc], clear = false)
    end

    for icc in 3:ipsi-1
        scalarplot!(p[1,1], agrid , view(solution[icc,:], agrid), label =  densityNames[icc], color= colors[icc], linestyle = linestyles[icc], clear = false)
    end

    reveal(p)

    Plotter.grid()
    Plotter.xlabel("space [m]")
    Plotter.ylabel("potential [V]")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.title("time \$ t =\$ $t, bias \$\\Delta u\$ = $Δu")
    Plotter.tight_layout()
    Plotter.gcf()
end

"""
$(TYPEDSIGNATURES)
Method for showing the total current in dependence of the applied voltage.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plotIV(Plotter, biasValues,IV, Δu, ;plotGridpoints=false)

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    Plotter.plot(biasValues[1:length(IV)], IV, marker = marker)
    Plotter.grid()
    Plotter.title("bias \$\\Delta u\$ = $Δu")
    Plotter.xlabel("bias [V]")
    Plotter.ylabel("total current [A]")
    Plotter.tight_layout()
    Plotter.pause(1.0e-5)
end
