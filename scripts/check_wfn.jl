using DataFrames, DrWatson

include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

if ~(@isdefined df)
    df = collect_results(datadir("sims", "gcorrect"))
else
    print("df exist, skip loading \n")
end


# plot_eigen_spectrum(df,findfirst(df[:, "depth"] .== 30))

function draw_wfn_script(df, data_indx, siteindx;kwargs...)
    H_eigen = df[data_indx, "solution"]
    zz = df[data_indx, "zz"]
    depth = df[data_indx, "depth"]
    mg = uconvert(NoUnits, m87Sr*g_n*λ813/2/Er/π)
    center_indx = find_center_index(H_eigen, zz, siteindx)
    ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
    E0 = real(H_eigen.values[center_indx[1]])
    ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])
    E1 = real(H_eigen.values[center_indx[2]])

    fig = plot(zz/π, -depth*cos.(zz).^2 .+ mg*zz, label="lattice", color=:black, lw=1, alpha=0.7)
    
    wfn_scaling = 17
    plot!(
        zz/π, wfn_scaling*real.(ψ_nz0) .+ E0,
        fill= + E0, 
        lw=0,
        # xlabel="Lattice site", 
        label="nz = 0", 
        color=palette(:tab20c)[1], 
        alpha=1,
        )
    plot!(zz/π, wfn_scaling*real.(ψ_nz1) .+ E1,
        fill = E1,
        lw=0,
        alpha=1,
        # title=df[data_indx, "depth"], 
        xlim=[-3.5, 3.5], 
        # xlabel="Lattice site",
        label="nz = 1", 
        color=palette(:tab20c)[5], 
    )
    plot!(;kwargs...)
    
    return fig
end

draw_wfn_script(df, findfirst(df[:, "depth"] .== 15), 0, legend=:false, 
# ylims=(-18, 2), 
size=((3+3/8)*96*2/3, (3+3/8)*96*3/4), grid=:false, axis=:off, ticks=:false,ylabel="Energy", font="helvetica")
# Plots.pdf(plotsdir("wfn_in_lattice.pdf"))