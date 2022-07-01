using JSON, DrWatson

data = JSON.parse(JSON.parsefile(datadir("exp_pro", "rabi_freqs_20220625.json")))

plot(data["nz1 ws1"])