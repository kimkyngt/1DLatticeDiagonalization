using DataFrames, DrWatson

include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))

df = collect_results(datadir("sims"))

draw_wfn(df, findfirst(df[:, "depth"] .== 7.5), 0)

