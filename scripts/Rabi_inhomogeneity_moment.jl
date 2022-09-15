using DrWatson, Interpolations, QuadGK, LsqFit, Plots, LaTeXStrings
include(srcdir("tools.jl"))
rabi_data = wload(datadir("rabi_carrier.jld2"))

depth = rabi_data["depth"]
rabi = rabi_data["rabi_carrier"]
p = sortperm(depth)
depth = depth[p]
rabi = rabi[:, p]

# Make interpolation function from the numerical data
rabi_0 = linear_interpolation(depth, rabi[1, :], extrapolation_bc=Flat())
rabi_1 = linear_interpolation(depth[depth .> 7.2], rabi[2, :][depth .> 7.2], extrapolation_bc=Flat())
rabi_2 = linear_interpolation(depth[depth .> 16], rabi[3, :][depth .> 16], extrapolation_bc=Flat())

# rabi_1 = extrapolate(interpolate(depth[depth .> 7.2], rabi[2, :][depth .> 7.2], Gridded(Linear())), extrapolation_bc=Flat())

#  Check the interpolation function
function check_rabi_freqs(;kwargs...)
    depths = range(1, 300, length=10000)
    fig = line(depths, rabi_0.(depths), label=L"n_z=0")
    plot!(depths, rabi_1.(depths), label=L"n_z=1")
    plot!(depths, rabi_2.(depths), label=L"n_z=2")
    plot!(
        xlabel="Lattice depth (Eᵣ)",
        ylabel="Fractional rabi frequency"
        ;kwargs...
        )
    # scatter!(depth, rabi[2, :], )
    return fig
end
check_rabi_freqs(;xscale=:log10, yscale=:log10, title="Testing interpolation", legend=:bottomright)

# function get_rho(r, U_0, T_r, nz)
#     """r is number in meter, U_0 is number in Er, T_r is u"K" """
#     # ω̄ = 2*π*2.44u"Hz"*(2*√(U_0) - 1/2) # for the ground state
#     A = get_Unz(r, U_0, w_0, nz)*uconvert(NoUnits, Er/k_B/T_r)
#     return exp.(-A)
# end


# # Parameters
# U_0 = 15
# T_r = get_Tr(U_0)
# w_0 = 260e-6 # cavity waist
# rmax = w_0 # in 

# y = rabi_1.(U_0*exp.(-2*rr.^2/w_0^2))
# DyDx = (y - circshift(y, 1))/(rr[2] - rr[1])

# rr = range(-rmax, rmax, length=1000)
# plot(rr*1e6, y,label=L"Ω_{n_z = 1}(r)", grid=true)
# plot!(rr*1e6, DyDx/maximum(DyDx),label=L"Ω'_{n_z = 1}(r)",)
