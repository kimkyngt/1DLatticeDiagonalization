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
    if U_0 > 15
        return 42 * √(U_0) * u"nK"
    else
        return (-45.2 + 14.1 * U_0) * u"nK"
    end
end
# Parameters
U_0 = 3
T_r_nz0 = get_Tr(U_0 * 1.2)*0 
T_r_nz1 = T_r_nz0 / 1.4

w_0 = 260e-6 # cavity waist
rmax = w_0
rmin = -rmax
νr = uconvert(u"Hz", sqrt(U_0 * Er / m87Sr / pi^2 / (w_0 * 1u"μm")^2))

rr = range(rmin, rmax, length=400)
fig = plot(rr * 1e6, -U_0 * exp.(-2 * rr .^ 2 / w_0^2), label="-U(x)/Eᵣ", minorticks=5, size=(400, 300), ls=:dash, lc=:black, lw=1)
plot!(rr * 1e6, -U_0 * exp.(-2 * rr .^ 2 / w_0^2) .+ ustrip(m87Sr * g_0 * sin(θtilt) / Er) * rr, label="-U(x)/Eᵣ + mgsin(θ)x", minorticks=5, size=(600, 400), lc=:black, lw=1)
plot!(rr * 1e6, get_Unz_sag(rr, U_0, w_0, 0), label="nz=0", lc=1)
plot!(rr * 1e6, get_Unz_sag(rr, U_0, w_0, 1), label="nz=1", lc=2)
# plot!(rr*1e6, get_Unz_sag(rr, U_0, w_0, 2), label="nz=2")
# plot!(rr*1e6, get_Unz_sag(rr, U_0, w_0, 3), label="nz=3")


rho_nz0 = [_get_rho(r, U_0, T_r_nz0) for r in rr] / _get_rho(0, U_0, T_r_nz0)
rho_nz1 = [_get_rho(r, U_0, T_r_nz1, nz=1) for r in rr] / _get_rho(0, U_0, T_r_nz1, nz=1)
plot!(
    rr * 1e6, U_0 / 3 * rho_nz0 .+ get_Unz_sag(0, U_0, w_0, 0),
    xlabel=L"\rho" * " (μm)", ylabel=L"U_{n_z}(\rho)~(E_{rec})", label=L"\rho^{n_z = 0}", ls=:solid, fill=get_Unz_sag(0, U_0, w_0, 0), alpha=0.5, lw=1, color=3,
    legend=:best)
# plot!(
#     rr*1e6, U_0/3*rho_nz0 .+ get_Unz_sag(0, U_0, w_0, 0), 
#     xlabel=L"\rho"*" (μm)", ylabel=L"U_{n_z}(\rho)", label=L"\rho^{n_z = 0}, no g", ls=:solid, fill=get_Unz_sag(0, U_0, w_0, 0), alpha=0.5, lw=0,
#     legend=:best)

# plot!(rr*1e6, U_0/3*rho_nz1 .+ get_Unz_sag(0, U_0, w_0, 1), 
#         xlabel=L"\rho"*" (μm)", ylabel=L"U_{n_z}(\rho)", label=L"\rho^{n_z = 1}", ls=:solid, fill=get_Unz_sag(0, U_0, w_0, 1), alpha=0.5, lw=0,
#         )

# plot!(rr*1e6, U_0/3*(rho_nz0 - rho_nz1), 
#         xlabel=L"\rho"*" (μm)", ylabel=L"U_{n_z}(\rho)", label=L"\Delta \rho", ls=:solid, fill=0, alpha=0.5, lw=1,
#         legend=:best)

plot!(xlabel=L"\rho" * " (μm)", ylabel=L"U_{n_z}(\rho)", label=L"\Delta \rho", ls=:solid, fill=0, alpha=0.5, lw=1, legend=:bottomright)


hline!([0], color=:grey, ls=:dash, label=:none)
plot!(title=string(round(U_0)) * "Er", xlims=(rmin * 1e6, rmax * 1e6))

# saving
d = (U_0=U_0, T_r=ustrip(u"nK", T_r_nz0), rmax=rmax)
Plots.pdf(fig, plotsdir("Unz", savename("Unz", d, "pdf")))
fig


