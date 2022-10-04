using DrWatson, Interpolations, QuadGK, SpecialFunctions, Plots, LaTeXStrings
include(srcdir("tools.jl"))
rabi_data = wload(datadir("rabi_carrier.jld2"))

depth = rabi_data["depth"]
rabi = rabi_data["rabi_carrier"]
p = sortperm(depth)
depth = depth[p]
rabi = rabi[:, p]

# Make interpolation function from the numerical data

rabi_0 = LinearInterpolation(depth, rabi[1, :], extrapolation_bc=0)
rabi_1 = LinearInterpolation(depth[depth .> 7.2], rabi[2, :][depth .> 7.2], extrapolation_bc=0)
rabi_2 = LinearInterpolation(depth[depth .> 16], rabi[3, :][depth .> 16], extrapolation_bc=0)

#  Check the interpolation function
function check_rabi_freqs(;kwargs...)
    depths = range(5, 30, length=10000)
    fig = plot(depths, rabi_0.(depths), label=L"n_z=0")
    plot!(depths[rabi_1.(depths) .∉ 0], rabi_1.(depths)[rabi_1.(depths) .∉ 0], label=L"n_z=1")
    # plot!(depths, rabi_2.(depths), label=L"n_z=2")
    plot!(
        xlabel="Lattice depth (Eᵣ)",
        ylabel="Fractional rabi frequency",
        # minorticks=10,
        grid=:true, 
        size=(400, 300),
        legend=:bottomright
        ;kwargs...)
    return fig
end

# plotly()
fig = check_rabi_freqs(
    # xscale=:log10, 
    yscale=:linear, 
# title="Testing interpolation"
)
# scatter!(depth[depth .> 7.2], rabi[2, :][depth .> 7.2])
# scatter!(depth[depth .> 16], rabi[3, :][depth .> 16])


# Plots.pdf(
#     fig, 
#     # plotsdir("Rabi_interpolation.pdf")
#     )

# Plots.pdf(
#     check_rabi_freqs(;xlims=(10, 16), ylims=(0.05, 0.3), title="Testing interpolation", legend=:bottomright), 
#     plotsdir("Rabi_interpolation.pdf")
#     )
fig
