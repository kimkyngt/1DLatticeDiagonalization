using DrWatson, DataFrames

include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))

if ~(@isdefined df)
    df = collect_results(datadir("sims"))
else
    print("df exist, skip loading \n")
end


depths = zeros(size(df)[1])
rabi_freqs = zeros(Complex, 8, size(df)[1])

for ii in range(1, size(df)[1])
    depths[ii] = df[ii, "depth"]
    rabi_freqs[:, ii] = get_rabi_frequency(df, ii)
end

wsave(datadir("sims", "rabi_frequencies.jld2"), Dict([("depth", depths), ("rabi_freqs", rabi_freqs)]))

fig1 = plot()
for ii in range(1, 4)
    plot!(depths, abs2.(rabi_freqs[ii, :]) , 
        seriestype=:line, 
        yscale=:log10, 
        xscale=:log10,
        xlims=(5, 60),
        ylims=(1e-4, 2),    
        xticks=(3:10:60, 3:10:60),
        color="black",
        title="n_z=0",
    )
end
fig2 = plot()
for ii in range(5,8)
    plot!(depths, abs2.(rabi_freqs[ii, :]) , 
        seriestype=:line, 
        yscale=:log10, 
        xscale=:log10,
        xlims=(5, 60),
        ylims=(1e-4, 2),    
        xticks=([8, 9, 10, 13,20, 30], [8, 9, 10, 13, 20, 30]),
        color="black", 
        title="n_z = 1"
    )
end
plot(fig1, fig2, layout=(2, 1))

# nz1
function fit_to_data()
    fig1 = plot()

    depthoffset = 0.95
    depthscale = 1.00
    Tr = get_Tr.(depths, "nz1")
    depths_eff = depthoffset .+ depthscale*depths + 1*(BoltzmannConstant*Tr/Er)
    for ii in range(5,8)
        plot!(depths_eff, abs2.(rabi_freqs[ii, :]) , 
            seriestype=:line, 
            yscale=:log10, 
            xscale=:log10,
            xlims=(7, 60),
            ylims=(1e-4, 2),    
            xticks=([8, 9, 10, 20, 30], [8, 9, 10, 20, 30]),
            color="black", 
            title="n_z = 1"
        )
    end
    overlay_data!("nz1")

    Tr = get_Tr.(depths, "nz0")
    depths_eff = depthoffset .+ depthscale*depths + 1*(BoltzmannConstant*Tr/Er)
    fig2 = plot()
    for ii in range(1, 4)
        plot!(depths_eff, abs2.(rabi_freqs[ii, :]) , 
            seriestype=:line, 
            yscale=:log10, 
            xscale=:log10,
            xlims=(5, 60),
            ylims=(1e-4, 2),    
            xticks=(3:5:60, 3:5:60),
            color="black",
            title="n_z=0",
        )
    end
    overlay_data!("nz0")
    fig = plot(fig1, fig2, layout=(2, 1), legend=false, size=(600, 600))
    return fig
end
fig = fit_to_data()
