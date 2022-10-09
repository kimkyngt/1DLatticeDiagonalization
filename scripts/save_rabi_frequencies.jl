using DrWatson, DataFrames

include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

if ~(@isdefined df)
    df = collect_results(datadir("sims", "gcorrect"))
else
    print("df exist, skip loading \n")
end


depths = zeros(size(df)[1])
rabi_freqs = zeros(Complex, 8, size(df)[1])

for ii in range(1, size(df)[1])
    depths[ii] = df[ii, "depth"]
    rabi_freqs[:, ii] = get_rabi_frequency(df, ii)
end

# Hand fittings

# # nz1
# function fit_to_data()
#     fig1 = plot()

#     depthoffset = 0.95
#     depthscale = 1.00
#     Tr = get_Tr.(depths, "nz1")
#     depths_eff = depthoffset .+ depthscale*depths + 1*(BoltzmannConstant*Tr/Er)
#     for ii in range(5,8)
#         plot!(depths_eff, abs2.(rabi_freqs[ii, :]) , 
#             seriestype=:line, 
#             yscale=:log10, 
#             xscale=:log10,
#             xlims=(7, 60),
#             ylims=(1e-4, 2),    
#             xticks=([8, 9, 10, 20, 30], [8, 9, 10, 20, 30]),
#             color="black", 
#             title="n_z = 1"
#         )
#     end
#     overlay_data!("nz1")

#     Tr = get_Tr.(depths, "nz0")
#     depths_eff = depthoffset .+ depthscale*depths + 1*(BoltzmannConstant*Tr/Er)
#     fig2 = plot()
#     for ii in range(1, 4)
#         plot!(depths_eff, abs2.(rabi_freqs[ii, :]) , 
#             seriestype=:line, 
#             yscale=:log10, 
#             xscale=:log10,
#             xlims=(5, 60),
#             ylims=(1e-4, 2),    
#             xticks=(3:5:60, 3:5:60),
#             color="black",
#             title="n_z=0",
#         )
#     end
#     overlay_data!("nz0")
#     fig = plot(fig1, fig2, layout=(2, 1), legend=false, size=(600, 600))
#     return fig
# end
# fig = fit_to_data()


# Save data
tosave = Dict([("depth", depths), ("rabi_freqs", rabi_freqs)])
wsave(datadir("rabi_frequencies.jld2"), tosave)
open(datadir("rabi_frequencies.json"),"w") do f
    JSON.print(f,tosave) 
end


# # carrier only 
rabi_carrier = zeros(Complex, 3, size(df)[1])

for ii in range(1, size(df)[1])
    depths[ii] = df[ii, "depth"]
    rabi_carrier[:, ii] = get_rabi_carrier_frequency(df, ii)
end

tosave = Dict([("depth", depths), ("rabi_carrier", rabi_carrier)])
wsave(datadir("rabi_carrier.jld2"), tosave)
open(datadir("rabi_carrier.json"),"w") do f
    JSON.print(f,tosave) 
end