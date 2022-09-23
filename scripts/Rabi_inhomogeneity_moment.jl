using DrWatson, Plots, LaTeXStrings, LsqFit, Polynomials, Interpolations
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

# Fit Rabi frequency to a series of polynomial to extract Rabi frequency for the 2nd order radial sideband.
U_0 = 15
w_0 = 260e-6
rmax = w_0 # in 
rmaxfit = w_0/3
rr = range(-rmax, rmax, length=1000)
ys = @. rabi_1(-get_U(U_0, rr, w_0))
rr_tofit = rr[ys .∉ 0]
ys_tofit = ys[ys .∉ 0]

ys_tofit = ys_tofit[abs.(rr_tofit) .< rmaxfit]
rr_tofit = rr_tofit[abs.(rr_tofit) .< rmaxfit]

n = 10
function model(x, p)
    return sum([p[i]x.^(2*(i-1)) for i in range(1, n)])
end
p0 = [float(1.0) for ii in range(1, n)]
fitresult = curve_fit(model, rr_tofit, ys_tofit, p0)
print(fitresult.param)
plot(rr[ys .∉ 0]*1e6, ys[ys .∉ 0])
plot!(rr_tofit*1e6, model(rr_tofit, fitresult.param ))