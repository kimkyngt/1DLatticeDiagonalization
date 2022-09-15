using DrWatson, Interpolations, QuadGK, LsqFit, Plots, LaTeXStrings
include(srcdir("tools.jl"))
rabi_data = wload(datadir("rabi_carrier.jld2"))

depth = rabi_data["depth"]
rabi = rabi_data["rabi_carrier"]
p = sortperm(depth)
depth = depth[p]
rabi = rabi[:, p]

# Make interpolation function from the numerical data
rabi_0 = LinearInterpolation(depth, rabi[1, :], extrapolation_bc=Flat())
rabi_1 = LinearInterpolation(depth[depth .> 7.2], rabi[2, :][depth .> 7.2], extrapolation_bc=Flat())
rabi_2 = LinearInterpolation(depth[depth .> 16], rabi[3, :][depth .> 16], extrapolation_bc=Flat())

#  Check the interpolation function
function check_rabi_freqs(;kwargs...)
    depths = range(1, 300, length=1000)
    fig = plot(depths, rabi_0.(depths), label=L"n_z=0")
    plot!(depths, rabi_1.(depths), label=L"n_z=1")
    plot!(depths, rabi_2.(depths), label=L"n_z=2")
    plot!(
        xlabel="Lattice depth (Eᵣ)",
        ylabel="Fractional rabi frequency"
        ;kwargs...)
    return fig
end


check_rabi_freqs(;xscale=:log10, yscale=:log10, title="Testing interpolation", legend=:bottomright)


# Plots.pdf(
#     check_rabi_freqs(;xscale=:log10, yscale=:log10, title="Testing interpolation", legend=:bottomright), 
#     plotsdir("Rabi_interpolation.pdf")
#     )

# Plots.pdf(
#     check_rabi_freqs(;xlims=(10, 16), ylims=(0.05, 0.3), title="Testing interpolation", legend=:bottomright), 
#     plotsdir("Rabi_interpolation.pdf")
#     )

function get_rho(r, U_0, T_r, nz)
    """r is number in meter, U_0 is number in Er, T_r is u"K" """
    # ω̄ = 2*π*2.44u"Hz"*(2*√(U_0) - 1/2) # for the ground state
    A = get_Unz(r, U_0, w_0, nz)*uconvert(NoUnits, Er/k_B/T_r)
    return exp.(-A)
end


# Parameters
U_0 = 13
T_r = get_Tr(U_0)
w_0 = 260e-6 # cavity waist
rmax = w_0 # in 

rr = range(-rmax, rmax, length=1000)

fig_rabi = plot(rr*1e6, get_rho(rr, U_0, T_r*0.6, 1)/get_rho(0, U_0, T_r*0.6, 1), label=L"\rho_{n_z=1}(r)", color=2, fill=0, alpha=0.3, lw=0, title=string(U_0)*"Er")
plot!(rr*1e6, get_rho(rr, U_0, T_r, 0)/get_rho(0, U_0, T_r, 0), label=L"\rho_{n_z=0}(r)", color=1, fill=0, alpha=0.3, lw=0)
plot!(rr*1e6, rabi_0.(U_0*exp.(-2*rr.^2/w_0^2)),label=L"Ω_{n_z = 0}", color=1)
plot!(rr*1e6, rabi_1.(U_0*exp.(-2*rr.^2/w_0^2)),label=L"Ω_{n_z = 1}", color=2, ylims=(-0.01, 1.01))

fig_potential = plot(rr*1e6, get_Unz(rr, U_0, w_0, 0), label=L"U_{n_z=0}(r)", color=1)
plot!(rr*1e6, get_Unz(rr, U_0, w_0, 1), label=L"U_{n_z=1}(r)", color=2)
hline!([0], color=:grey, ls=:dash, label="")
plot!(rr*1e6, -U_0*exp.(-2*rr.^2/w_0^2), label=L"U(r)", legend=:bottomright, ylabel="Lattice depth (Er)", color=:grey,)

fig_tot = plot(fig_rabi, fig_potential, layout=(2 ,1), size=(400, 500) , 
xlabel="x (μm)")

sim_params = Dict(
    "depth" => U_0,
    "temp" => ustrip(T_r),
)
Plots.pdf(fig_tot, plotsdir("Rabi_inhomogeneity", savename(sim_params, "pdf")
))

fig_tot
