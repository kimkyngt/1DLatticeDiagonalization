using DataFrames, DrWatson

include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

df = collect_results(datadir("sims"))

data_num = size(df)[1]
depths = zeros(data_num)
shifts = zeros(data_num)

alpha_qm = -0.00125
beta = -0.446e-6

for ii in range(1, data_num)

    H_eigen = df[ii, "solution"]
    zz = df[ii, "zz"]
    depths[ii] = df[ii, "depth"]

    center_indx = find_center_index(H_eigen, zz, 1)
    ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
    ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])

    sinkz_nz0 = get_sine_sq_exp(ψ_nz0, zz)
    sinkz_nz1 = get_sine_sq_exp(ψ_nz1, zz)
    
    shifts[ii] = -alpha_qm*depths[ii]*(sinkz_nz1 - sinkz_nz0) + beta*depths[ii]*(2*sqrt(depths[ii]) -3)
end
I = sortperm(depths)
depths = depths[I]
shifts = shifts[I]

fig = plot(
    depths, 
    shifts/ustrip(fclock) ,
    yticks=[0, 1, 2, 3, 4, 5, 6]*1e-17, 
    ylim=(-0.5e-17, 5e-17), 
    dpi=300,
    lw=2,
    marker=:circle,
    legend=false,
    ylabel="Differential shift", 
    xlabel="Lattice depth (Er)",
    xticks=0:50:300,
    # size=(400, 300),
)
Plots.pdf(fig, plotsdir("nz_modulation_lock_simulation.pdf"))
