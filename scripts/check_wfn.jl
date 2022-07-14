using DataFrames, DrWatson

include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))

df = collect_results(datadir("sims", "gcorrect"))

draw_wfn(df, findfirst(df[:, "depth"] .== 8.5), 0)

plot_eigen_spectrum(df,findfirst(df[:, "depth"] .== 8.5) )