using DrWatson, Interpolations, QuadGK
include(srcdir("tools.jl"))
include(srcdir("initialize.jl"))
function _get_rho(r, U_0, T_r; nz=0)
    """r is number in meter, U_0 is number in Er, T_r is u"K" """
    # ω̄ = 2*π*2.44u"Hz"*(2*√(U_0) - 1/2) # for the ground state
    A = get_Unz_sag(r, U_0, w_0, nz) * uconvert(NoUnits, Er / k_B / T_r)
    rho = exp(-A)
    return rho
end
function get_Tr(U_0)
    31e-9*(U_0 - 2.2)^(0.58)*1u"K"
end
# Parameters
U_0 = 551.7*(0.54 - 0.0017 - 0.52)
T_r_nz0 = get_Tr(U_0)
# T_r_nz0 = 17u"nK"
T_r_nz1 = T_r_nz0 / 1.4

w_0 = 260e-6 # cavity waist
rmax = w_0
rmin = -rmax
νr = uconvert(u"Hz", sqrt(U_0 * Er / m87Sr / pi^2 / (w_0 * 1u"μm")^2))

rr = range(rmin, rmax, length=400)
fig = plot(rr * 1e6, -U_0 * exp.(-2 * rr .^ 2 / w_0^2), label="-U(x)/Eᵣ",  size=(400, 300), ls=:dash, lc=:black, lw=1)
plot!(rr * 1e6, -U_0 * exp.(-2 * rr .^ 2 / w_0^2) .+ ustrip(m87Sr * g_0 * sin(θtilt) / Er) * rr, label="",  size=(600, 400), lc=:black, lw=2)
# plot!(rr * 1e6, get_Unz_sag(rr, U_0, w_0, 0), label="nz=0", lc=1)
# plot!(rr * 1e6, get_Unz_sag(rr, U_0, w_0, 1), label="nz=1", lc=2)
# plot!(rr*1e6, get_Unz_sag(rr, U_0, w_0, 2), label="nz=2")
# plot!(rr*1e6, get_Unz_sag(rr, U_0, w_0, 3), label="nz=3")


rho_nz0 = [_get_rho(r, U_0, T_r_nz0) for r in rr] / _get_rho(0, U_0, T_r_nz0)
rho_nz1 = [_get_rho(r, U_0, T_r_nz1, nz=1) for r in rr] / _get_rho(0, U_0, T_r_nz1, nz=1)
plot!(
    rr * 1e6, U_0 /1.2 * rho_nz0 .- U_0,
    xlabel="x" * " (μm)", ylabel="U(x)/Eᵣ", label=L"\rho^{n_z = 0}", ls=:solid, fill=-U_0, alpha=0.5, lw=2, color=1,
    legend=:best)

# hline!([0], color=:grey, ls=:dash, label=:none)
plot!(xlims=(rmin * 1e6, rmax * 1e6), legend=false, ylims=(-11, 0), size=(350, 200), grid=false)
# annotate!(0, -5, string(round(u"nK",T_r_nz0)), grid=false)

# saving
d = (U_0=U_0, T_r=ustrip(u"nK", T_r_nz0), rmax=rmax)
# Plots.pdf(fig, plotsdir("Unz", savename("Unz", d, "pdf")))
fig


