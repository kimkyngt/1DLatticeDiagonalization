using DrWatson, Interpolations, QuadGK, LsqFit, Plots, LaTeXStrings
include(srcdir("tools.jl"))
rabi_data = wload(datadir("rabi_carrier.jld2"))

depth = rabi_data["depth"]
rabi = rabi_data["rabi_carrier"]
p = sortperm(depth)
depth = depth[p]
rabi = abs.(rabi[:, p])

# Make interpolation function from the numerical data

rabi_0 = LinearInterpolation(depth, rabi[1, :], extrapolation_bc=0)
rabi_1 = LinearInterpolation(depth[depth .> 7.2], rabi[2, :][depth .> 7.2], extrapolation_bc=0)
rabi_2 = LinearInterpolation(depth[depth .> 16], rabi[3, :][depth .> 16], extrapolation_bc=0)

#  Check the interpolation function
function check_rabi_freqs(;kwargs...)
    depths = range(1, 300, length=10000)
    fig = plot(depths, rabi_0.(depths), label=L"n_z=0")
    plot!(depths, rabi_1.(depths), label=L"n_z=1")
    plot!(depths, rabi_2.(depths), label=L"n_z=2")
    plot!(
        xlabel="Lattice depth (Eᵣ)",
        ylabel="Fractional rabi frequency",
        minorticks=5,
        ;kwargs...)
    return fig
end

# Parameters
U_0 = 15
T_r = get_Tr(U_0)*1.2
w_0 = 260e-6 # cavity waist
rmax = w_0 # in 

rr = range(-rmax, rmax, length=1000)

fig_rabi = plot(rr*1e6, get_rho(rr, U_0, T_r*0.6, 1)/get_rho(0, U_0, T_r*0.6, 1), label=L"\rho_{n_z=1}(r)", color=2, fill=0, alpha=0.3, lw=0, 
title=string(U_0)*"Er"
)
plot!(rr*1e6, get_rho(rr, U_0, T_r, 0)/get_rho(0, U_0, T_r, 0), label=L"\rho_{n_z=0}(r)", color=1, fill=0, alpha=0.3, lw=0)
omega_0 = rabi_0.(U_0*exp.(-2*rr.^2/w_0^2))
omega_1 = rabi_1.(U_0*exp.(-2*rr.^2/w_0^2))
plot!(rr[omega_0 .∉ 0]*1e6, omega_0[omega_0 .∉ 0],label=L"Ω_{n_z = 0}(r)", color=1)
plot!(rr[omega_1 .∉ 0]*1e6, omega_1[omega_1 .∉ 0],label=L"Ω_{n_z = 1}(r)", color=2, ylims=(-0.01, 1.01), xlabel="r (μm)")

fig_potential = plot(rr*1e6, get_Unz(rr, U_0, w_0, 0), label=L"U_{n_z=0}(r)", color=1, bottom_margin=20Plots.px, left_margin=20Plots.px)
plot!(rr*1e6, get_Unz(rr, U_0, w_0, 1), label=L"U_{n_z=1}(r)", color=2)
hline!([0], color=:grey, ls=:dash, label="")
plot!(rr*1e6, -U_0*exp.(-2*rr.^2/w_0^2), label=L"U(r)", legend=:bottomright, ylabel="Energy (Eᵣ)", color=:grey,)

# fig_tot = plot(fig_rabi, fig_potential, layout=(1 ,2), size=(800, 300) , xlabel="r (μm)")

fig_tot = plot(fig_rabi, size=(400, 300))

sim_params = Dict(
    "depth" => U_0,
    "temp" => ustrip(T_r),
)
Plots.pdf(fig_tot, plotsdir("Rabi_inhomogeneity", savename(sim_params, "pdf")
))

fig_tot
