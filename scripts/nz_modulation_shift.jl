using DataFrames, DrWatson

include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

if ~(@isdefined df)
    df = collect_results(datadir("sims", "gcorrect"))
    data_num = size(df)[1]
    depths = zeros(data_num)
    shifts = zeros(data_num)
    
    alpha_qm = -0.00124
    
    for ii in range(1, data_num)
    
        H_eigen = df[ii, "solution"]
        zz = df[ii, "zz"]
        depths[ii] = df[ii, "depth"]
    
        center_indx = find_center_index(H_eigen, zz, 1)
        ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
        ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])
    
        sinkz_nz0 = get_sine_sq_exp(ψ_nz0, zz)
        sinkz_nz1 = get_sine_sq_exp(ψ_nz1, zz)
        
        shifts[ii] = -alpha_qm*depths[ii]*(sinkz_nz1 - sinkz_nz0)
    end
    sortedindx = sortperm(depths)
    depths = depths[sortedindx]
    shifts = shifts[sortedindx]
    shifts = shifts[depths .> 10]
    depths = depths[depths .> 10]

else
    print("df exist, skip loading \n")
end


fig = plot(
    depths, 
    shifts/ustrip(fclock) ,
    yticks=[0, 3, 6]*1e-17, 
    ylim=(-0.5e-17, 6e-17), 
    dpi=300,
    lw=2,
    marker=:circle,
    ylabel="Differential shift", 
    xlabel="Lattice depth "*L"\mathrm{(E_r)}",
    xticks=[0, 30, 50, 100, 200, 300],
    # xticklabels=["", "", "","","","","",""],
    # size=(400, 300),
    label="WS",
    left_margin = 10Plots.px,
    right_margin=-10Plots.px
)
plot!(depths, -sqrt.(depths)*alpha_qm/ustrip(fclock), lw=2, label=L"\propto\sqrt{u}", legend=:bottomright,)

fig2 = plot(depths, shifts/ustrip(fclock) + sqrt.(depths)*alpha_qm/ustrip(fclock), 
    ylim=(0e-19, 5e-19), 
    lw=2, 
    # ylabel="Differential shift", 
    yticks=[-2, 0, 2, 4]*1e-19,
    xticks=[0, 30, 50, 100, 200, 300],
    xlabel="Lattice depth "*L"\mathrm{(E_r)}",
    label="Difference",
    right_margin = 5Plots.px,
    left_margin=-10Plots.px,
    bottom_margin=10Plots.px,
    legend=true
    ) 
figtot = plot(fig, fig2, layout=(1, 2), size=(700, 300), link=:x,)
# Plots.pdf(figtot, plotsdir("nz_modulation_lock_simulation.pdf"))
