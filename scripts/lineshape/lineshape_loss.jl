using DrWatson
using Plots, LaTeXStrings
include(srcdir("lineshape.jl"))

T_rabi = 0.02
omega0 = π/T_rabi
param1 = Dict([
    ("Omega0" , omega0),
    ("Omega1" , omega0*0.1),
    ("Omega2" , omega0*0.01),
    ("detuning_max", omega0*2),
    ("detuning_num", 50),
    ("T_rabi", T_rabi),
    ("Delta" , 2pi*868),
    ("initial_condition", "pure"),
    ("N_lattice", 2*10+1),
    ("reltol", 1e-9),
    ("gamma", 0.2),
    ("initial_spin", "up"),
])

param2 = copy(param1)
param2["gamma"] = 0.0

result1 = get_spectrum_loss(param1)
result2 = get_spectrum_loss(param2)

## Plotting 
fig1 = plot(result1["detunings"]/(2pi), result1["exc_frac"], 
    lw=2, 
    label="γ = "*string(param1["gamma"]), 
    color=1
    )
plot!(
    result2["detunings"]/(2pi), result2["exc_frac"], 
    lw=2, label="γ = "*string(param2["gamma"]), 
    xlabel="Hz", title="T_rabi="*string(T_rabi), color=2
    )

fig2 = plot(result1["detunings"]/(2pi), result1["exc_frac"] - result2["exc_frac"], lw=2, label="Δ", xlabel="Hz", title="initial_condition: loss vs no loss")

# Discriminator and line pulling
modulation_depth = 0.8/T_rabi
frac_to_detuning1, diff_frac1 = get_discriminator(result1, modulation_depth)
frac_to_detuning2, diff_frac2= get_discriminator(result2, modulation_depth)
fig3 = plot(diff_frac1, frac_to_detuning1(diff_frac1)/2/pi, lw=2)
plot!(diff_frac2, frac_to_detuning2(diff_frac2)/2/pi,  
    lw=2, legend=false, 
    xlabel=L"\delta \rho^{ee}", 
    ylabel="Frequency (Hz)", 
    title="Line pulling (1E-18): "*string(round(frac_to_detuning1(0)/ustrip(fclock)*1e18/2/pi, digits=1))*", "*string(round(frac_to_detuning2(0)/ustrip(fclock)*1e18/2/pi, digits=1)),
    marker=:circle,
    markerstrokewith=0,
    )

fig = plot(fig1, fig2, fig3, layout=(2, 2), size=(800, 600), margins=20Plots.px)
# # Plots.pdf(fig, plotsdir("lineshape", "lineshape_comparison", savename(param2)))
# # @tagsave(datadir("lineshape", "lineshape_comparison", savename(param1, "jld2")), result1)
# # @tagsave(datadir("lineshape", "lineshape_comparison", savename(param2, "jld2")), result2)


# fig