
"""
$(TYPEDSIGNATURES)
Method which can be used to construct the arrays parsed to the plotting routines for labeling.
The description for electrons and holes are predefined. If one wishes to extend by labels for,
e.g. mobile ionic carriers or traps, this can be done within the main file.

"""
function set_plotting_labels(data)

    label_energy   = Array{String, 2}(undef, 2, data.params.numberOfCarriers) # band-edge energies and potential
    label_BEE      = Array{String, 1}(undef, data.params.numberOfCarriers)    # band-edge energie parameters
    label_density  = Array{String, 1}(undef, data.params.numberOfCarriers)
    label_solution = Array{String, 1}(undef, data.params.numberOfCarriers)

    # indices (∈ IN) of electron and hole quasi Fermi potentials specified by user
    iphin       = data.bulkRecombination.iphin # integer index of φ_n
    iphip       = data.bulkRecombination.iphip # integer index of φ_p

    ## for electrons
    label_energy[1, iphin] = "\$E_c-q\\psi\$"; label_energy[2, iphin] = "\$ - q \\varphi_n\$"; label_BEE[iphin] = "\$E_c\$"
    label_density[iphin]   = "n";              label_solution[iphin]  = "\$ \\varphi_n\$"

    ## for holes
    label_energy[1, iphip] = "\$E_v-q\\psi\$"; label_energy[2, iphip] = "\$ - q \\varphi_p\$"; label_BEE[iphip] = "\$E_v\$"
    label_density[iphip]   = "p";              label_solution[iphip]  = "\$ \\varphi_p\$"


    return label_solution, label_density, label_energy, label_BEE
end

"""
$(TYPEDSIGNATURES)
Plotting routine, where the charge carrier densities are depicted
in dependence of space. The case of heterojunctions is tested, but yet
multidimensional plottings are not included.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.

"""
function plot_densities(Plotter, ctsys, solution, title, label_density, ;plotGridpoints=false)

    Plotter.clf()

    grid            = ctsys.fvmsys.grid
    data            = ctsys.fvmsys.physics.data
    numberOfRegions = grid[NumCellRegions]

    if dim_space(grid) > 1
        error("plot_densities is so far only tested in 1D")
    end

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    params     = data.params
    colors     = ["green", "red", "gold", "purple", "orange"]

    for icc in 1:params.numberOfCarriers

        # grids = Array{ExtendableGrid, 1}(undef, numberOfRegions)
        # nicc  = Array{Array{Float64, 1}, 1}(undef, numberOfRegions)

        ## first region for label
        label_icc = label_density[icc]
        subg      = subgrid(grid, [1])
        ncc       = get_density(solution, 1, ctsys, icc)

        Plotter.semilogy(subg[Coordinates]', 1.0e-6 .*ncc, marker = marker, label = label_icc, color = colors[icc], linewidth = 2)

        ## additional regions
        for ireg in 2:numberOfRegions
            subg = subgrid(grid, [ireg])
            ncc  = get_density(solution, ireg, ctsys, icc)

            ## Note that this implies a 1D plot, for multidimensional plots, you may work with
            ## GridVisualize.jl or write your own code.
            Plotter.semilogy(subg[Coordinates]', 1.0e-6 .*ncc, marker = marker, color = colors[icc], linewidth = 2)
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
multidimensional plottings are not included.

One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.

"""
function plot_energies(Plotter, ctsys, sol, title, label_energy, ;plotGridpoints=false)

    Plotter.clf()

    grid        = ctsys.fvmsys.grid
    data        = ctsys.fvmsys.physics.data
    params      = data.params
    paramsnodal = data.paramsnodal
    ipsi        = data.index_psi
    cellnodes   = grid[CellNodes]
    cellregions = grid[CellRegions]
    coord       = grid[Coordinates]

    if length(coord[1]) != 1
        println("plotEnergies is so far only implemented in 1D")
    end

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    colors         = ["green", "red", "gold", "purple", "orange"]
    linestyles     = ["-", ":", "--", "-.", "-"]

    for icc in data.electricCarrierList
        # first cell
        ireg         = cellregions[1]
        E1           = params.bBandEdgeEnergy[icc, 1] + paramsnodal.bandEdgeEnergy[icc, 1] # left boundary
        E2           = params.bandEdgeEnergy[icc, 1]  + paramsnodal.bandEdgeEnergy[icc, 2]
        energy_icc1  = E1 - q * sol[ipsi, 1]
        energy_icc2  = E2 - q * sol[ipsi, 2]

        Plotter.plot([coord[1]./1, coord[2]./1], [energy_icc1, energy_icc2]./q, marker = marker, label = label_energy[1, icc], linewidth = 2, color = colors[icc], linestyle = linestyles[1])

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

        Plotter.plot(coord[1,:]./1, - sol[icc,:], label = label_energy[2, icc], marker = marker, linewidth = 2, color = colors[icc], linestyle = linestyles[2])

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
function plot_energies(Plotter, ctsys, label_BEE)

    grid        = ctsys.fvmsys.grid
    data        = ctsys.fvmsys.physics.data
    params      = data.params
    paramsnodal = data.paramsnodal

    coord       = grid[Coordinates]
    cellregions = grid[CellRegions]
    cellnodes   = grid[CellNodes]

    if length(coord[1]) != 1
        error("plot_energies is so far only implemented in 1D")
    end

    colors      = ["green", "red", "gold", "purple", "orange"]
    linestyles  = ["-", ":", "--", "-.", "-"]

    # plot different band-edge energies values in interior
    for icc = 1:params.numberOfCarriers
        for i in eachindex(cellregions)
            # determine band-edge energy value in cell and number of cell nodes
            cellValue            = ( params.bandEdgeEnergy[icc, cellregions[i]] + paramsnodal.bandEdgeEnergy[icc, i] )/q
            numberLocalCellNodes = length(cellnodes[:,i])
            # patch together cells
            Plotter.plot(coord[cellnodes[:,i]],
                        repeat(cellValue:cellValue,numberLocalCellNodes),
                        marker="x",
                        color=colors[icc],
                        linewidth=3,
                        linestyle=linestyles[icc]);
        end

        Plotter.plot(NaN, NaN, color=colors[icc], linewidth = 3, linestyle = linestyles[icc], label = label_BEE[icc]) # legend
    end

    # plot different band-edge energy values on boundary
    bfaceregions = grid[BFaceRegions]
    bfacenodes   = grid[BFaceNodes]

    for icc = 1: params.numberOfCarriers

        for i in eachindex(bfaceregions)
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
    Plotter.title("Band-edge energies \$ E_\\alpha \$")
    Plotter.legend(fancybox = true, loc = "best")
    Plotter.tight_layout()
    Plotter.show()

end


"""
$(TYPEDSIGNATURES)
Possibility to plot the considered doping. This is especially useful
for making sure that the interior and the boundary doping agree.

"""
function plot_doping(Plotter, ctsys, label_density)

    g           = ctsys.fvmsys.grid
    data        = ctsys.fvmsys.physics.data
    params      = data.params
    coord       = g[Coordinates]
    cellregions = g[CellRegions]
    cellnodes   = g[CellNodes]
    coord       = g[Coordinates]

    if length(coord[1]) != 1
        error("plot_doping is so far only implemented in 1D")
    end

    colors       = ["green", "red", "gold", "purple", "orange"]
    linestyles   = ["-", ":", "--", "-.", "-"]

    # plot different doping values in interior
    for icc = 1:params.numberOfCarriers

        for i in eachindex(cellregions)
            # determine doping value in cell and number of cell nodes
            cellValue            = params.doping[icc, cellregions[i]]
            numberLocalCellNodes = length(cellnodes[:,i])

            # patch together cells (multiplying by 1.0e-6 gives us the densities in cm^(-3))
            Plotter.plot(coord[cellnodes[:,i]],
            1.0e-6 .*repeat(cellValue:cellValue,numberLocalCellNodes),
                            color=colors[icc],
                            linewidth=3,
                            linestyle=linestyles[icc]);
        end
        # legend
        Plotter.plot(NaN, NaN, color = colors[icc], linewidth = 3, label = label_density[icc])

    end

    # plot different doping values on boundary
    bfaceregions = g[BFaceRegions]
    bfacenodes   = g[BFaceNodes]

    for icc = 1: params.numberOfCarriers

        for i in eachindex(bfaceregions)
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
Plot doping for nodal dependent doping.
"""
function plot_doping(Plotter, g::ExtendableGrid, paramsnodal::ParamsNodal)

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
multidimensional plottings are not included.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plot_solution(Plotter, ctsys, solution, title, label_solution, ;plotGridpoints=false)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    if dim_space(grid) > 1
        error("plot_solution is so far only tested in 1D")
    end

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    coord        = grid[Coordinates]'
    ipsi         = data.index_psi

    colors       = ["green", "red", "gold", "purple", "orange"]
    linestyles   = ["-", ":", "--", "-.", "-"]

    Plotter.clf()
    Plotter.plot(coord, solution[ipsi,:], marker = marker, label = "\$\\psi\$", color="b", linewidth= 3)
    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOn
        ipsiStandard = data.barrierLoweringInfo.ipsiStandard
        Plotter.plot(coord, solution[ipsiStandard,:], marker = marker, label = "\$\\psi\$ (Schottky contacts)", color="black", linestyle = ":", linewidth= 3)
    end

    # electrons and holes
    for icc ∈ data.electricCarrierList
        Plotter.plot(coord./1, solution[icc,:], label =  label_solution[icc], marker = marker, color= colors[icc], linestyle = linestyles[1], linewidth= 3)
    end

    for icc ∈ data.ionicCarrierList
        # DA: could be modified by using subgrids from ExtendableGrids.
        # this is needed to only plot present ionic charge carrier in respective defined regions
        regions    = grid[CellRegions]

        subregions = zeros(Int64, 0)
        for ix in 1:length(icc.regions)
            subreg = findall(x -> x == icc.regions[ix], regions)
            append!(subregions, subreg)
            push!(subregions, subregions[end]+1)
        end
        subgrid    = coord[subregions]

        icc        = icc.ionicCarrier

        Plotter.plot(subgrid./1, solution[icc, subregions], label =  label_solution[icc], marker = marker, color= colors[icc], linestyle = linestyles[1], linewidth= 3)
    end

    for icc ∈ data.trapCarrierList
        # this is needed to only plot present ionic charge carrier in respective defined regions
        regions    = grid[CellRegions]

        subregions = zeros(Int64, 0)
        for ix in 1:length(icc.regions)
            subreg = findall(x -> x == icc.regions[ix], regions)
            append!(subregions, subreg)
            push!(subregions, subregions[end]+1)
        end
        subgrid    = coord[subregions]

        icc = icc.trapCarrier

        Plotter.plot(subgrid./1, solution[icc, subregions], label =  label_solution[icc], marker = marker, color= colors[icc], linestyle = linestyles[1], linewidth= 3)
    end

    Plotter.grid()
    Plotter.xlabel("space [m]")
    Plotter.ylabel("potential [V]")
    Plotter.legend(fancybox = true, loc = "best", fontsize=11)
    Plotter.title(title)
    Plotter.tight_layout()
    Plotter.gcf()

end

function plot_solution(Plotter, grid, solution, agrid, t, Δu, label_solution)

    # Create a visualizer. Works with Plots (fast once compiled) and PyPlot
    p = GridVisualizer(Plotter = Plotter, layout = (1,1) )

    ipsi         = data.index_psi

    colors       = ["green", "red", "gold", "purple", "orange"]
    linestyles   = ["-", ":", "--", "-.", "-"]

    Plotter.clf()
    scalarplot!(p[1,1], grid, solution[ipsi,:], label = "\$\\psi\$", color="b", marker = "x", title="time \$ t =\$ $t, bias \$\\Delta u\$ = $Δu", clear = true)

    for icc in [iphin, iphip]
        scalarplot!(p[1,1], grid, solution[icc,:], label = label_solution[icc], color= colors[icc], linestyle = linestyles[icc], clear = false)
    end

    for icc in 3:data.params.numberOfCarriers
        scalarplot!(p[1,1], agrid , view(solution[icc,:], agrid), label = label_solution[icc], color= colors[icc], linestyle = linestyles[icc], clear = false)
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
Method for showing the total current.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plot_IV(Plotter, biasValues,IV, title, ;plotGridpoints=false)

    if plotGridpoints == true
        marker = "o"
    else
        marker = ""
    end

    Plotter.plot(biasValues[1:length(IV)], IV, marker = marker)
    Plotter.grid()
    Plotter.title(title)
    Plotter.xlabel("bias [V]")
    Plotter.ylabel("total current [A]")
    Plotter.tight_layout()
    Plotter.pause(1.0e-5)
end
