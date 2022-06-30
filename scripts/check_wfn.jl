using DataFrames

include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))

df = collect_results(datadir("sims"))

function draw_wfn(df, data_indx)
    H_eigen = df[data_indx, "solution"]
    zz = df[data_indx, "zz"]
    center_indx = find_center_index(H_eigen, zz)
    ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
    ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])
    fig = plot(zz/π, abs2.(ψ_nz0), dpi=300, 
        title=df[data_indx, "depth"], 
        xlim=[-5, 5]
    )
    plot!(zz/π, abs2.(ψ_nz1))
    return fig
end
draw_wfn(df, findfirst(df[:, "depth"] .== 11))

