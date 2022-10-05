using DrWatson, Interpolations, Plots, LaTeXStrings

include(srcdir("tools.jl"))
rabi_data = wload(datadir("rabi_frequencies.jld2"))

depth = rabi_data["depth"]
rabi = rabi_data["rabi_freqs"]
p = sortperm(depth)
depth = depth[p]
fig = plot(depth, rabi[1, p].*rabi[2, p], label=L"\Omega_0 \Omega_1")
plot!(depth, rabi[1, p].*rabi[3, p], label=L"\Omega_0 \Omega_2")
plot!(depth, rabi[1, p].*rabi[4, p], label=L"\Omega_0 \Omega_3")
plot!(depth, rabi[1, p], label=L"\Omega_0", ls=:dash)
plot!(depth, rabi[2, p], label=L"\Omega_1", ls=:dash)
plot!(depth, rabi[3, p], label=L"\Omega_2", ls=:dash)
plot!(depth, rabi[4, p], label=L"\Omega_3", ls=:dash)
plot!(
    yscale=:log10, 
    xscale=:log10, 
    ylims=(1e-2, 1.2),
    xlabel="Lattice depth (Eᵣ)",
    minorticks=10
    )

Plots.pdf(fig, plotsdir("line_pulling_WS", "Rabi_frequency_products.pdf"))

Δ = 868
T_rabi = 1.2
Ω_0 = 1/T_rabi * 0.5
fig2 = plot(depth, Ω_0^2*rabi[1, p].*rabi[2, p]/Δ/ustrip(fclock), label=L"\Omega_0 \Omega_1/\Delta")
plot!(depth, Ω_0^2*rabi[1, p].*rabi[3, p]/Δ/ustrip(fclock), label=L"\Omega_0 \Omega_2/\Delta")
plot!(depth, Ω_0^2*rabi[1, p].*rabi[4, p]/Δ/ustrip(fclock), label=L"\Omega_0 \Omega_3 /\Delta")
plot!(
    title="T_rabi = "*string(T_rabi),
    yscale=:log10, 
    # xscale=:log10, 
    ylims=(1e-20, 3e-18),
    xlims=(2.5, 12),
    ylabel="(Maximal) Line pulling from Ω₀Ωᵢ/Δ",
    xlabel="Lattice depth (Eᵣ)",
    yminorticks=10
    )
fig2
Plots.pdf(fig2, plotsdir("line_pulling_WS", "Line_pulling_estimation_"*savename(Dict([("T_rabi", T_rabi)]))*".pdf"))
Plots.png(fig2, plotsdir("line_pulling_WS", "Line_pulling_estimation_"*savename(Dict([("T_rabi", T_rabi)]))*".png"))
