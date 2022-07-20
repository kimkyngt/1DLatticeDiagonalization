using DrWatson, Interpolations, LsqFit

include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

#load and sort
rabi_calc = wload(datadir("rabi_frequencies.jld2"))
depths = rabi_calc["depth"]
rabi_freqs = rabi_calc["rabi_freqs"]
p = sortperm(depths)
depths = depths[p]
rabi_freqs = rabi_freqs[:, p]

# cutoff for nz = 1
depths_nz0 = depths
rabi_freqs_nz0 = rabi_freqs[1:4, :]
depths_nz1 = depths[depths .> 7.9]
rabi_freqs_nz1 = rabi_freqs[5:8, depths .> 7.9]

# check loaded data
fig1 = plot()
for ii=1:4
    plot!(depths_nz0, abs2.(rabi_freqs_nz0[ii, :]), 
    xscale=:log10,  
    yscale=:log10,
    ylim=(1e-4, 2)
    )
end
fig2 = plot()
for ii=1:4
    plot!(depths_nz1, abs2.(rabi_freqs_nz1[ii, :]),
    xscale=:log10, 
    yscale=:log10,
    ylim=(1e-4, 2),
    )
end
fig = plot(fig1, fig2, layout=(2, 1))

xs = depths_nz1
A = abs2.(rabi_freqs_nz1[1, :])
inter_func = LinearInterpolation(xs, A, extrapolation_bc=Line())
function func(volt_Tr, p)
    """eff_depth(Er) = OFFSET(Er) + SLOPE(Er/V)*voltage(V) - kB*Tr (Er)"""
    # voltage = volt_Tr[1]
    # Tr = volt_Tr[2]
    eff_depth = p[1] .+ p[2]*volt_Tr[1,:] .- ustrip(kB/Er)*volt_Tr[2,:]
    return inter_func(eff_depth)
end

# Load experimental data
data1 = JSON.parse(JSON.parsefile(datadir("exp_pro", "rabi_freqs_20220625.json")))
data2 = JSON.parse(JSON.parsefile(datadir("exp_pro", "rabi_freqs_20220626.json")))
data = merge(data1, data2)
open(datadir("rabi_frequencies_data.json"),"w") do f
    JSON.print(f,data) 
end
# Fitting test
voltage = (0.54 .+ data["nz1 ws0"]["depth mV"])
Tr = [ustrip(get_Tr(data["nz1 ws0"]["depth Erec"][ii], "nz1")) for ii in range(1, length(voltage))]
volt_Tr = [[voltage[ii], Tr[ii]] for ii in range(1, length(voltage))]

Ω0 = maximum(data["nz0 ws0"]["rabi frequency"])
y = abs2.(data["nz1 ws0"]["rabi frequency"]/Ω0)

p0 = [0, 0]
fit = curve_fit(func, volt_Tr, y, p0)