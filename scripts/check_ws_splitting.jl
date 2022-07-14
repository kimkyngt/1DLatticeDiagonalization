using DrWatson
include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))

plotly()
using DelimitedFiles
excs = readdlm(datadir("exp_raw", "WS_3_8_excs.txt"))
freq = readdlm(datadir("exp_raw", "WS_3_8_freq.txt"))
plot(freq, excs)
