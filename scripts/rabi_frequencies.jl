using DrWatson, DataFrames

include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))

df = collect_results(datadir("sims"))


depths = zeros(size(df)[1])
rabi_freqs = zeros(Complex, 8, size(df)[1])

for ii in range(1, size(df)[1])
    depths[ii] = df[ii, "depth"]
    rabi_freqs[:, ii] = get_rabi_frequency(df, ii)
end

fig = plot(dpi=300)
for ii in range(1,8)
    plot!(depths, abs2.(rabi_freqs[ii, :]) , 
        seriestype=:line, 
        yscale=:log10, 
        xscale=:log10,
        xlims=(8, 30),
        ylims=(1e-4, 1),    
        xticks=(5:2:30, 5:2:30)
    )
end
fig


### nz=1 only
fig = plot(dpi=300)
for ii in range(5,8)
    plot!(depths, abs2.(rabi_freqs[ii, :]) , 
        seriestype=:line, 
        yscale=:log10, 
        xscale=:log10,
        xlims=(8, 30),
        ylims=(1e-4, 1),    
        xticks=([9, 10, 20, 30], [9, 10, 20, 30])
    )
end
fig
