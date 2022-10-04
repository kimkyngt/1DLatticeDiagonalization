using DrWatson, Interpolations, QuadGK, LsqFit, Plots, LaTeXStrings
include(srcdir("tools.jl"))

energy_data = wload(datadir("ground_energy.jld2"))
rabi_data = wload(datadir("rabi_carrier.jld2"))

depth = rabi_data["depth"]
rabi = rabi_data["rabi_carrier"]
p = sortperm(depth)
depth = depth[p]
rabi = rabi[:, p]


rabi_0 = LinearInterpolation(depth, rabi[1, :], extrapolation_bc=Flat())
rabi_1 = LinearInterpolation(depth[depth .> 7.2], rabi[2, :][depth .> 7.2], extrapolation_bc=Flat())
rabi_2 = LinearInterpolation(depth[depth .> 16], rabi[3, :][depth .> 16], extrapolation_bc=Flat())


function integrand(U_0, T_r, r, t)
    get_rho(r, U_0, T_r)*sin(rabi_0(U_0*exp(-2*r^2/w_0^2))*t/2)^2
end

function get_pe_avg(U_0, T_r, t)
    pe, err = quadgk(r -> r*integrand(U_0, T_r, r, t), 0, rmax)
    return pe
end

function get_normalziation(U_0, T_r)
    N, err = quadgk(r -> r*get_rho(r, U_0, T_r), 0, rmax)
    return N
end

# Parameters
U_0 = 5
T_r = 40u"nK"#80u"nK"
w_0 = 260e-6 # cavity waist
rmax = w_0*2 # in 

# fig1 density
rr = range(0, rmax, length=100)
fig_density = plot(rr*1e6, [get_rho(r, U_0, T_r) for r in rr] / get_rho(0, U_0, T_r), xlabel="r (μm)", label="ρ")
plot!(rr*1e6, rabi_0.(U_0*exp.(-2*rr.^2/w_0^2)), color=2, ylims=(0, 1.1) , label="Ω/Ω₀")
plot!(rr*1e6, exp.(-2*rr.^2/w_0^2), label="-U(r)/U₀", legend=:right)
annotate!(150, 1.05, string(round(U_0, digits=2))*"Er, "*string(round(ustrip(T_r), digits=2))*"nK", 9)

# fig2 rabi flopping
tt = collect(range(0, 40π, length=100))
pebar = [get_pe_avg(U_0, T_r, t) for t in tt] / get_normalziation(U_0, T_r)
fig_rabi = plot(tt/π, pebar, st=:scatter, label="Sim")
model(t, p) = p[1] .+ p[2]*exp.(-p[3] * t) .* sin.(p[4]*t/2).^2
p0 = [0.0, 1.0, 0.2, 0.5]
fit = curve_fit(model, tt, pebar, p0)
plot!(tt/π, model(tt, fit.param), label="Fit", legend=true, ylabel="Excitation fraction", xlabel="Time (Ω₀t/π)")
annotate!(2, 1.05, L"$Ω/Ω_{fit} = $"*string(round(rabi_0(U_0)/fit.param[4], digits=2)), 9)

# println(fit.param)
fig = plot(fig_density, fig_rabi, size= (600, 300), bottom_margin = 10Plots.px, ylim=(0, 1.1))

cond=(U_0 = U_0, Tr = ustrip(u"nK", T_r))
# Plots.pdf(fig, plotsdir("rabi_flopping_sim", savename(cond, "pdf")))
# Check normalization
# integral, err = quadgk(x -> get_rho(x, U_0, T_r), 0, 1e-2, )
fig
