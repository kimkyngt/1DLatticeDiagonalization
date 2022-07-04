using JSON, DrWatson, Plots


fig = plot(dpi=300)
data = JSON.parse(JSON.parsefile(datadir("exp_pro", "rabi_freqs_20220625.json")))
Ω0 = maximum(data["nz0 ws0"]["rabi frequency"])

for key in keys(data)
    print(key)
    plot!(data[key]["depth Erec"],abs2.(data[key]["rabi frequency"]/Ω0), marker=:circle, xscale=:log10, yscale=:log10, label=key, lw=0)
end

data = JSON.parse(JSON.parsefile(datadir("exp_pro", "rabi_freqs_20220626.json")))
for key in keys(data)
    print(key)
    plot!(data[key]["depth Erec"],abs2.(data[key]["rabi frequency"]/Ω0), marker=:circle, xscale=:log10, yscale=:log10, label=key, lw=0)
end

fig