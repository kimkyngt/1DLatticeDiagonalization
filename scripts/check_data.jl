using DataFrames

df = collect_results(datadir("sims"))
print(df.depth)