using DrWatson
using Plots, LaTeXStrings
include(srcdir("lineshape.jl"))

T_rabi = 0.38
omega0 = π/T_rabi
param1 = Dict([
    ("Omega0" , omega0),
    ("Omega1" , omega0*0.1),
    ("Omega2" , omega0*0.01),
    ("detuning_max", omega0*8),
    ("detuning_num", 100),
    ("T_rabi", T_rabi),
    ("Delta" , 2pi*868),
    ("initial_condition", "pure"),
    ("N_lattice", 2*10+1),
    ("reltol", 1e-9),
])

param2 = copy(param1)
param2["initial_condition"] = "plus"

result1 = get_spectrum(param1)
result2 = get_spectrum(param2)
@tagsave(datadir("lineshape", "lineshape_comparison", savename(param1, "jld2")), result1)
@tagsave(datadir("lineshape", "lineshape_comparison", savename(param2, "jld2")), result2)

fig1 = plot(result1["detunings"]/(2pi), result1["exc_frac"], lw=2, label=param1["initial_condition"], color=1)
vline!([sum(result1["detunings"]/(2pi) .* result1["exc_frac"])/sum(result1["exc_frac"])], ls=:dash, color=1, label=param1["initial_condition"]*"COM", lw=1)

plot!(result2["detunings"]/(2pi), result2["exc_frac"], lw=2, label=param2["initial_condition"], xlabel="Hz", title="T_rabi="*string(T_rabi), color=2)
vline!([sum(result2["detunings"]/(2pi) .* result2["exc_frac"])/sum(result2["exc_frac"])], ls=:dash, color=2, label=param2["initial_condition"]*"COM", lw=1)

fig2 = plot(result1["detunings"]/(2pi), result1["exc_frac"] - result2["exc_frac"], lw=2, label="Δ", xlabel="Hz", title="initial_condition: pure vs plus")
fig = plot(fig1, fig2, size=(800, 300), margins=20Plots.px)
Plots.pdf(fig, plotsdir("lineshape", "lineshape_comparison", savename(param2)))
fig