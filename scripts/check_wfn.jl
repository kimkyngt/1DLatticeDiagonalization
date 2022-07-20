using DataFrames, DrWatson

include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

if ~(@isdefined df)
    df = collect_results(datadir("sims", "gcorrect"))
else
    print("df exist, skip loading \n")
end


draw_wfn(df, findfirst(df[:, "depth"] .== 50), 0, size=(400, 250), grid=:true, legend=:true, yticks=false)
# plot_eigen_spectrum(df,findfirst(df[:, "depth"] .== 30))

function draw_wfn_script(df, data_indx, siteindx;kwargs...)
    H_eigen = df[data_indx, "solution"]
    zz = df[data_indx, "zz"]
    center_indx = find_center_index(H_eigen, zz, siteindx)
    ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
    ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])

    fig = plot(zz/π, 0.15*cos.(kclock/k813*zz), label="clock laser", color=palette(:tab10)[4], lw=1, alpha=0.5)
    
    plot!(
        zz/π, real.(ψ_nz0),
        fill=0, 
        xlim=[-10, 10], 
        xlabel="Lattice site", 
        label="nz = 1", 
        color=palette(:tab10)[1], 
        alpha=0.7,
        )
    plot!(zz/π, real.(ψ_nz1),
        fill=0,
        alpha=0.7,
        title=df[data_indx, "depth"], 
        xlim=[-6, 6], 
        xlabel="Lattice site",
        label="nz = 0", 
        color=palette(:tab10)[2], 
        ;kwargs...
    )
    
    return fig
end

draw_wfn_script(df, findfirst(df[:, "depth"] .== 50), 0, size=(400, 250), grid=:true, legend=:true, yticks=false)
