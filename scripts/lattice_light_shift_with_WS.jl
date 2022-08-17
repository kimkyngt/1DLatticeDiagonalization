using DataFrames, DrWatson

include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

if ~(@isdefined df)
    df = collect_results(datadir("sims", "gcorrect"))
else
    print("df exist, skip loading \n")
end

# Preallocations
Ndata = size(df)[1]
depths = zeros(Ndata)
coskz_sq = zeros(Ndata)
sinkz_sq = zeros(Ndata)

# Compute matrix elements
@time for ii in range(1, Ndata)
    H_eigen = df[ii, "solution"]
    zz = df[ii, "zz"]

    center_indx = find_center_index(H_eigen, zz, 1)
    ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])

    depths[ii] = df[ii, "depth"]
    sinkz_sq[ii] = get_sine_sq_exp(ψ_nz0, zz)
    coskz_sq[ii] = get_cos_sq_exp(ψ_nz0, zz)
end
indx = sortperm(depths)


# Light shift coefficients
alpha_qm = -0.00124
dalpha_E1dnu = 1.855e-11
delta_E1 = -00e6
alpha_E1 = dalpha_E1dnu*delta_E1

# Plot harmonic ligt shift curves
dnu_harmonic = (alpha_E1 - alpha_qm)*0.5*sqrt.(depths) .- alpha_E1*depths
dnu_WS = -alpha_E1*coskz_sq.*depths - alpha_qm*sinkz_sq.*depths

plt = plot(depths[indx], dnu_harmonic[indx]/ustrip(fclock), st=:line, label="Harmonic")
plot!(depths[indx], dnu_WS[indx]/ustrip(fclock), st=:line,ls=:dash,  label="WS")
plot!(
    ylim=(0, 1.5e-17), xlim=(0, 50), 
    yticks=collect(0:0.5:1.5)*1e-17,
    xlabel="Lattice depth (Eᵣ)",
    ylabel=L"\Delta \nu_{LS}",
    legend=:true
    )

plt2 = plot(depths[indx], (dnu_harmonic[indx] - dnu_WS[indx])/ustrip(fclock), label="Harmonic - WS", ylims=(-1e-19, 1e-19), legend=:true, yticks=[-1, 0, 1]*1e-19, xlabel="Lattice depth (Eᵣ)", xscale=:log10)

fig = plot(plt, plt2, size=(700, 300), margin=10Plots.px)
Plots.pdf(plt, plotsdir("Lightshift_harmonic_WS_comparison.pdf"))

