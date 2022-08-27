using DrWatson, DataFrames, Trapz, Plots

include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

if ~(@isdefined df)
    df = collect_results(datadir("sims", "gcorrect"))
else
    print("df exist, skip loading \n")
end
sort!(df, [:depth])
# Preallocations
Ndata = size(df)[1]
depths = zeros(Ndata)
powers = [0, 1, 2, 4, 6]
moments_nz0 = zeros(Ndata, length(powers))
moments_nz1 = zeros(Ndata, length(powers))

# Compute matrix elements
@time for ii in range(1, Ndata)
    H_eigen = df[ii, "solution"]
    center_indx = find_center_index(H_eigen, df[ii, "zz"], 1)
    ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
    ψ_nz0 = ψ_nz0/sqrt(trapz(df[ii, "zz"], abs2.(ψ_nz0)))
    ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])
    ψ_nz1 = ψ_nz1/sqrt(trapz(df[ii, "zz"], abs2.(ψ_nz1)))

    for jj in range(1, length(powers)) 
        moments_nz0[ii, jj] = trapz(df[ii, "zz"], get_periodic_zpower(df[ii, "zz"], powers[jj]) .* abs2.(ψ_nz0)) 
        moments_nz1[ii, jj] = trapz(df[ii, "zz"], get_periodic_zpower(df[ii, "zz"], powers[jj]) .* abs2.(ψ_nz1)) 
    end
    depths[ii] = df[ii, "depth"]
end
indx = sortperm(depths)

fig_nz1 = plot(depths[depths .> 8], abs.(moments_nz1[depths .> 8, :]), 
xticks= ([ 7, 15, 30, 100], string.([7, 15, 30, 100])),
# st=:scatter, 
xscale=:log10, 
xlims=(2, 200),
# ylims=(-0.1, 1.1),
# yscale=:log10,
legend=:best,
label=[L"(kz)^0" L"(kz)^1" L"(kz)^2" L"(kz)^4" L"(kz)^6"], ylabel="Expectation value", 
title=L"n_z=1",
grid=:true,
xlabel="Lattice depth (Eᵣ)")

fig_nz0 = plot(depths[depths .> 3], abs.(moments_nz0[depths .> 3, :]), xticks= ([ 3, 6, 10, 30, 100], string.([ 3, 6, 10, 30, 100])),
xlims=(2, 200),
# ylims=(-0.1, 1.1),
# st=:scatter, 
xscale=:log10, 
# yscale=:log10,
legend=:best,
grid=:true,
title=L"n_z=0",
label=[L"(kz)^0" L"(kz)^1" L"(kz)^2" L"(kz)^4" L"(kz)^6"], ylabel="Expectation value", 
xlabel="Lattice depth (Eᵣ)")

fig = plot(fig_nz0, fig_nz1, layout=(2, 1), size=(500, 500), margins=10Plots.px,)
Plots.pdf(fig, plotsdir("moments", "moments_axial.pdf"))
Plots.png(fig, plotsdir("moments", "moments_axial.png"))
fig