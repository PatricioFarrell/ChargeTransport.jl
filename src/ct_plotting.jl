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
function plot_densities(Plotter, grid, data, sol, title, ;plotGridpoints=false)

    params      = data.params
    paramsnodal = data.paramsnodal

    Plotter.clf()

    if dim_space(grid) > 2
        error("plot_densities! is so far only tested in 1D and 2D")
    end

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    colors       = Array{String, 1}(undef, 4)
    linestyles   = Array{String, 1}(undef, 4)
    densityNames = Array{String, 1}(undef, 4)

    iphin     = data.bulk_recombination.iphin
    iphip     = data.bulk_recombination.iphip
    # DA: Caution!!! This can cause problems when plotting!!
    ionic_vac = 3:4#data.enable_ion_vacancies.ionic_vacancies

    colors[iphin]     = "green";            linestyles[iphin]     = "-";          densityNames[iphin]     = "n"
    colors[iphip]     = "red";              linestyles[iphip]     = ":";          densityNames[iphip]     = "p"
    colors[ionic_vac] = ["gold", "purple"]; linestyles[ionic_vac] = ["--", "-."]; densityNames[ionic_vac] = ["a", "c"]

    # colors                      = ["green", "red", "gold", "purple"]
    # linestyles                  = ["-", ":", "--", "-."]
    # densityNames                = ["n", "p", "a", "c"]

    ipsi                        = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1

    cellnodes                   = grid[CellNodes]
    bfacenodes                  = grid[BFaceNodes]
    cellregions                 = grid[CellRegions]
    bfaceregions                = grid[BFaceRegions]
    coordinates                 = grid[Coordinates]
    
    for icc in 1:params.numberOfCarriers

        # first cell
        u1                      = sol[:, 1]
        u2                      = sol[:, 2]
        ireg                    = cellregions[1]

        icc1                    = compute_densities!(u1, data, 1, 1, icc, ipsi, false) # breg = 1 since we are on the left boundary
        icc2                    = compute_densities!(u2, data, 2, ireg, icc, ipsi, true) 

        label_icc               = densityNames[icc]

        if dim_space(grid) == 1
            
            Plotter.semilogy([coordinates[1]./1, coordinates[2]./1], 1.0e-6 .*[icc1, icc2], marker = marker, label = label_icc, color = colors[icc], linewidth = 2) #multiplying by 1.0e-6 gives us the densities in cm^(-3)
    
            for icell in 2:size(cellnodes,2) - 1
                in_region = true
                i1        = cellnodes[1,icell]
                i2        = cellnodes[2,icell]
                ireg      = cellregions[icell]
                node      = i1 

                u1        = sol[:, i1]
                u2        = sol[:, i2]
     
                icc1      = compute_densities!(u1, data, i1, ireg, icc, ipsi, in_region)
                icc2      = compute_densities!(u2, data, i2, ireg, icc, ipsi, in_region)
        
                Plotter.semilogy([coordinates[i1]./1, coordinates[i2]./1], 1.0e-6 .*[icc1, icc2], marker = marker, color = colors[icc], linewidth = 2) #multiplying by 1.0e-6 gives us the densities in cm^(-3)     
            end

            # last cell
            u1            = sol[:, end-1]
            u2            = sol[:, end]
            ireg          = cellregions[end]
            node          = cellnodes[2, end]

            icc1          = compute_densities!(u1, data, node-1, ireg, icc, ipsi, true)
            icc2          = compute_densities!(u2, data, node, 2, icc, ipsi, false) # breg = 2 since we are on the right boundary

            Plotter.semilogy([coordinates[node-1]./1, coordinates[node]./1], 1.0e-6 .*[icc1, icc2], marker = marker, color = colors[icc], linewidth = 2) #multiplying by 1.0e-6 gives us the densities in cm^(-3)

        elseif dim_space(grid) == 2

            # inner regions
            for ibcell in 1:length(bfaceregions)

                in_region = false
                i1        = bfacenodes[1,ibcell]
                i2        = bfacenodes[2,ibcell]
                ibreg      = bfaceregions[ibcell]

                u1        = sol[:, i1]
                u2        = sol[:, i2]
     
                icc1      = compute_densities!(u1, data, i1, ibreg, icc, ipsi, in_region)
                icc2      = compute_densities!(u2, data, i2, ibreg, icc, ipsi, in_region)
    
                Plotter.surf([grid[XCoordinates][i1], grid[XCoordinates][i2]], [grid[YCoordinates][i1], grid[YCoordinates][i2]], [log(icc1), log(icc2)] );
            end
        
        end

    end

    Plotter.grid()
    Plotter.xlabel("space [\$m\$]")
    Plotter.ylabel("density [\$\\frac{1}{cm^3}\$]")
    Plotter.legend(fancybox = true, loc = "best", fontsize=11)
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
function plot_energies(Plotter, grid, data, sol, title, ;plotGridpoints=false)

    params      = data.params
    paramsnodal = data.paramsnodal

    Plotter.clf()

    ipsi                = params.numberOfCarriers  + params.numberOfInterfaceCarriers + 1

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

    colors              = Array{String, 1}(undef, 4)
    linestyles          = Array{String, 1}(undef, 4)
    labelBandEdgeEnergy = Array{String, 1}(undef, 4)
    labelPotential      = Array{String, 1}(undef, 4)

    iphin     = data.bulk_recombination.iphin
    iphip     = data.bulk_recombination.iphip
    # DA: Caution!!! This can cause problems when plotting!!
    ionic_vac = 3:4#data.enable_ion_vacancies.ionic_vacancies

    colors[iphin]     = "green";            linestyles[iphin]     = "-";          labelBandEdgeEnergy[iphin]     = "\$E_c-q\\psi\$"; labelPotential[iphin] = "\$ - q \\varphi_n\$"
    colors[iphip]     = "red";              linestyles[iphip]     = ":";          labelBandEdgeEnergy[iphip]     = "\$E_v-q\\psi\$"; labelPotential[iphip] = "\$ - q \\varphi_n\$"
    colors[ionic_vac] = ["gold", "purple"]; linestyles[ionic_vac] = ["--", "-."]; labelBandEdgeEnergy[ionic_vac] = ["\$E_a-q\\psi\$ ", "\$E_{cat}-q\\psi\$ "]; labelPotential[ionic_vac] = [ "\$ - q \\varphi_a\$", "\$ - q \\varphi_c\$"]

    for icc in [iphin, iphip]
        # first cell
        ireg         = cellregions[1]
        E1           = params.bBandEdgeEnergy[icc, 1] + paramsnodal.bandEdgeEnergy[icc, 1] # left boundary
        E2           = params.bandEdgeEnergy[icc, 1]  + paramsnodal.bandEdgeEnergy[icc, 2] 
        energy_icc1  = E1 - q * sol[ipsi, 1]
        energy_icc2  = E2 - q * sol[ipsi, 2]
        label_energy = labelBandEdgeEnergy[icc]

        Plotter.plot([coord[1]./1, coord[2]./1], [energy_icc1, energy_icc2]./q, marker = marker, label = label_energy, linewidth = 2, color = colors[icc], linestyle = linestyles[1])

        for icell in 2:size(cellnodes,2) - 1
            i1          = cellnodes[1,icell]
            i2          = cellnodes[2,icell]
            ireg        = cellregions[icell]

            E1          = params.bandEdgeEnergy[icc, ireg] + paramsnodal.bandEdgeEnergy[icc, i1]
            E2          = params.bandEdgeEnergy[icc, ireg] + paramsnodal.bandEdgeEnergy[icc, i2]

            energy_icc1 = E1 - q * sol[ipsi, i1]
            energy_icc2 = E2 - q * sol[ipsi, i2]

            Plotter.plot([coord[i1]./1, coord[i2]./1], [energy_icc1, energy_icc2]./q, marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[1]) 
        end

        ireg        = cellregions[end]
        node        = cellnodes[2, end]
        E1          = params.bandEdgeEnergy[icc, ireg] + paramsnodal.bandEdgeEnergy[icc, node-1] 
        E2          = params.bBandEdgeEnergy[icc, 2] + paramsnodal.bandEdgeEnergy[icc, end] # right boundary
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
function plot_energies(Plotter, grid::ExtendableGrid, data)

    params      = data.params
    paramsnodal = data.paramsnodal 

    coord       = grid[Coordinates]
    cellregions = grid[CellRegions]
    cellnodes   = grid[CellNodes]

    if length(coord[1]) != 1
        error("plotEnergies is so far only implemented in 1D")
    end

    colors      = Array{String, 1}(undef, 4)
    styles      = Array{String, 1}(undef, 4)
    EnergyNames = Array{String, 1}(undef, 4)

    iphin     = data.bulk_recombination.iphin
    iphip     = data.bulk_recombination.iphip
    # DA: Caution!!! This can cause problems when plotting!!
    ionic_vac = 3:4#data.enable_ion_vacancies.ionic_vacancies

    colors[iphin]     = "green";            styles[iphin]     = "-";          EnergyNames[iphin]     = "\$E_c\$";
    colors[iphip]     = "red";              styles[iphip]     = ":";          EnergyNames[iphip]     = "\$E_v\$"
    colors[ionic_vac] = ["gold", "purple"]; styles[ionic_vac] = ["--", "-."]; EnergyNames[ionic_vac] = ["\$E_a\$ ", "\$E_{cat}\$ "]

    # plot different band-edge energies values in interior
    for icc = 1:params.numberOfCarriers
        for i in 1:length(cellregions)
            # determine band-edge energy value in cell and number of cell nodes
            cellValue            = ( params.bandEdgeEnergy[icc, cellregions[i]] + paramsnodal.bandEdgeEnergy[icc, i] )/q
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

    for icc = 1: params.numberOfCarriers

        for i in 1:length(bfaceregions)
            # determine band-edge energy value in cell and number of cell nodes
            cellValue            = (params.bBandEdgeEnergy[icc, bfaceregions[i]] + paramsnodal.bandEdgeEnergy[icc, bfacenodes[i]])/q
            numberLocalCellNodes = length(bfacenodes[:,i])

            # patch together cells
            Plotter.plot(coord[bfacenodes[:,i]],
                        repeat(cellValue:cellValue,numberLocalCellNodes),
                        marker="x",
                        markersize=10,
                        color=colors[icc]);
        end

    end

    Plotter.grid()
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
function plot_doping(Plotter, g::ExtendableGrid, data::ChargeTransportData)

    params      = data.params
    coord       = g[Coordinates]
    cellregions = g[CellRegions]
    cellnodes   = g[CellNodes]
    coord       = g[Coordinates]

    if length(coord[1]) != 1
        error("plotDoping is so far only implemented in 1D")
    end

    
    colors       = Array{String, 1}(undef, 4)
    styles       = Array{String, 1}(undef, 4)
    densityNames = Array{String, 1}(undef, 4)

    iphin     = data.bulk_recombination.iphin
    iphip     = data.bulk_recombination.iphip
    # DA: Caution!!! This can cause problems when plotting!!
    ionic_vac = 3:4#data.enable_ion_vacancies.ionic_vacancies

    colors[iphin]     = "green";            styles[iphin]     = "-";          densityNames[iphin]     = "n";
    colors[iphip]     = "red";              styles[iphip]     = ":";          densityNames[iphip]     = "p"
    colors[ionic_vac] = ["gold", "purple"]; styles[ionic_vac] = ["--", "-."]; densityNames[ionic_vac] = ["a", "c"]
    

    # plot different doping values in interior
    for icc = 1:params.numberOfCarriers

        for i in 1:length(cellregions)
            # determine doping value in cell and number of cell nodes
            cellValue            = params.doping[icc, cellregions[i]]
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

    for icc = 1: params.numberOfCarriers

        for i in 1:length(bfaceregions)
            # determine doping value in cell and number of cell nodes
            cellValue            = params.bDoping[icc, bfaceregions[i]]
            numberLocalCellNodes = length(bfacenodes[:,i])

            # patch together cells
            Plotter.plot(coord[bfacenodes[:,i]],
            1.0e-6 .*repeat(cellValue:cellValue,numberLocalCellNodes),
                            marker="x",
                            markersize=10,
                            color=colors[icc]);
        end

    end

    Plotter.grid()
    Plotter.yscale("symlog")
    Plotter.xlabel("space [\$m\$]")
    Plotter.ylabel("Doping [\$\\frac{1}{cm^3}\$]")
    Plotter.title("Doping values for charge carriers")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.tight_layout()

end

"""
Plot doping for nodal dependent doping
"""
function plot_doping(Plotter, g::ExtendableGrid, paramsnodal::ChargeTransportParamsNodal)

    coord  = g[Coordinates]

    Plotter.plot(coord[:], 1.0e-6 .*paramsnodal.doping[:], color = "green", marker = "x")

    Plotter.grid()
    Plotter.yscale("symlog")
    Plotter.xlabel("space [\$m\$]")
    Plotter.ylabel("Doping [\$\\frac{1}{cm^3}\$]")
    Plotter.title("Doping values for charge carriers")
    Plotter.tight_layout()


end

"""
$(TYPEDSIGNATURES)
Plotting routine for depicting the electroneutral potential.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plot_electroNeutralSolutionBoltzmann(Plotter, grid, psi0, ;plotGridpoints=false)

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    coord = grid[Coordinates]
    
    Plotter.grid()
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
function plot_solution(Plotter, grid, data, solution, title, ;plotGridpoints=false)

    coord                 = grid[Coordinates]'
    ipsi                  = data.params.numberOfCarriers + data.params.numberOfInterfaceCarriers + 1
    
    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    colors       = Array{String, 1}(undef, 4)
    linestyles   = Array{String, 1}(undef, 4)
    densityNames = Array{String, 1}(undef, 4)

    iphin     = data.bulk_recombination.iphin
    iphip     = data.bulk_recombination.iphip
    # DA: Caution!!! This can cause problems when plotting!!
    ionic_vac = 3:4#data.enable_ion_vacancies.ionic_vacancies

    colors[iphin]     = "green";            linestyles[iphin]     = "--";       densityNames[iphin]     = "\$\\varphi_n\$"
    colors[iphip]     = "red";              linestyles[iphip]     = "-.";       densityNames[iphip]     = "\$\\varphi_p\$"
    colors[ionic_vac] = ["gold", "purple"]; linestyles[ionic_vac] = ["-", ":"]; densityNames[ionic_vac] = ["\$\\varphi_a\$", "\$\\varphi_c\$"]

    Plotter.clf() 
    Plotter.plot(coord, (solution[ipsi,:] + data.params.Eref/q*ones(length(solution[ipsi,:]))), marker = marker, label = "\$\\psi\$", color="b")

    for icc in 1:data.params.numberOfCarriers
        Plotter.plot(coord./1, solution[icc,:], label =  densityNames[icc], marker = marker, color= colors[icc], linestyle = linestyles[icc])
    end
            
    Plotter.grid()
    Plotter.xlabel("space [m]")
    Plotter.ylabel("potential [V]")
    Plotter.legend(fancybox = true, loc = "best", fontsize=11)
    Plotter.title(title)
    Plotter.tight_layout()
    Plotter.gcf()

end

function plot_solution(Plotter, grid, solution, Eref, agrid, t, Δu)


    # Create a visualizer. Works with Plots (fast once compiled) and PyPlot
    p = GridVisualizer(Plotter = Plotter, layout = (1,1) )

    ipsi = data.params.numberOfCarriers + data.params.numberOfInterfaceCarriers + 1

    colors       = Array{String, 1}(undef, 4)
    linestyles   = Array{String, 1}(undef, 4)
    densityNames = Array{String, 1}(undef, 4)

    iphin     = data.bulk_recombination.iphin
    iphip     = data.bulk_recombination.iphip
    # DA: Caution!!! This can cause problems when plotting!!
    ionic_vac = 3:4#data.enable_ion_vacancies.ionic_vacancies

    colors[iphin]     = "green";            linestyles[iphin]     = "--";       densityNames[iphin]     = "\$\\varphi_n\$"
    colors[iphip]     = "red";              linestyles[iphip]     = "-.";       densityNames[iphip]     = "\$\\varphi_p\$"
    colors[ionic_vac] = ["gold", "purple"]; linestyles[ionic_vac] = ["-", ":"]; densityNames[ionic_vac] = ["\$\\varphi_a\$", "\$\\varphi_c\$"]
    
    Plotter.clf() 
    scalarplot!(p[1,1], grid, (solution[ipsi,:] + Eref/q*ones(length(solution[ipsi,:]))), label = "\$\\psi\$", color="b",  marker = "x", 
    title="time \$ t =\$ $t, bias \$\\Delta u\$ = $Δu", clear = true)

    for icc in [iphin, iphip]
        scalarplot!(p[1,1], grid, solution[icc,:], label =  densityNames[icc], color= colors[icc], linestyle = linestyles[icc], clear = false)
    end

    for icc in 3:data.params.numberOfCarriers
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
function plot_IV(Plotter, biasValues,IV, Δu, ;plotGridpoints=false)

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
