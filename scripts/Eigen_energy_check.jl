using DrWatson, DataFrames
include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

# # Extract ground state energy
# if ~(@isdefined df)
#     df = collect_results(datadir("sims", "gcorrect"))
# else
#     print("df exist, skip loading \n")
# end
# E_ground = zeros(size(df)[1], 3) # get
# depths = zeros(size(df)[1])
# for ii in range(1, size(df)[1])
#     E, depth = get_axial_eigen_energy(df, ii)
#     E_ground[ii, :] = E[1:3]
#     depths[ii] = depth
# end
# wsave(datadir("eigen_energy.jld2"), Dict("eigen_energy"=>E_ground, "depth"=>depths))

@time data = wload(datadir("ground_energy.jld2"))
E_ground = data["E"]
depths = data["depth"]

fig1 = plot(depths, depths + E_ground, st=:scatter, lc=1, label="WS state, "*L"n_z=0", msw=0, markersize=4, size=(400, 300))
plot!(depths[depths .> 10], depths[depths .> 10] + data["E1"][depths .> 10], st=:scatter, lc=0, label="WS state, "*L"n_z=1", markersize=4, msw=0)
xx = range(1, 300)
ν_z = 2*xx.^(1/2)
ν_rec = 1
plot!(xx, (0+1/2)*ν_z .- 0.5*ν_rec*(0^2 + 0 + 1/2), lw=1, lc=:black, label="Harmonic model")
plot!(xx, (1+1/2)*ν_z .- 0.5*ν_rec*(1^2 + 1 + 1/2), lw=1, lc=:black,label=false)
plot!(xlabel="Lattice depth "*L"(\mathrm{E_r})", ylabel="Eigen energy "*L"(\mathrm{E_r})", legend=:topleft,xscale=:linear, yscale=:linear)

# Difference
ν_z = 2*depths.^(1/2)
ν_rec = 1

fig2 = plot(depths, depths + E_ground - ((0+1/2)*ν_z .- 0.5*ν_rec*(0^2 + 0 + 1/2)), st=:scatter, lc=1, label="WS state, "*L"n_z=0", msw=0, markersize=4, size=(400, 300))
plot!(depths[depths .> 10], depths[depths .> 10] + data["E1"][depths .> 10] - ((1+1/2)*ν_z[depths .> 10]  .- 0.5*ν_rec*(1^2 + 1 + 1/2)), st=:scatter, lc=0, label="WS state, "*L"n_z=1", markersize=4, msw=0)

plot!(xlabel="Lattice depth "*L"(\mathrm{E_r})", ylabel="Eigen energy "*L"(\mathrm{E_r})", legend=:topleft,xscale=:linear, yscale=:linear, legend_location=:bottomright)

plot(fig1, fig2, layout=(2, 1), size=(400, 500))

# Plots.pdf(plotsdir("Eigen_energy_vs_depths.pdf"))