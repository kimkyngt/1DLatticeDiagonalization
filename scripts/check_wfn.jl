using DataFrames, DrWatson

include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

if ~(@isdefined df)
    df = collect_results(datadir("sims", "gcorrect"))
else
    print("df exist, skip loading \n")
end

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

    max_lattice = 0.5
    fig = plot(zz[abs.(zz/π) .< max_lattice]/π, -depth*cos.(zz[abs.(zz/π) .< max_lattice]).^2 .+ 0*mg*zz[abs.(zz/π) .< max_lattice], label="lattice", color=:black, lw=2, alpha=0.7)
    # plot!(zz[abs.(zz/π) .< max_lattice]/π, -depth*cos.(zz[abs.(zz/π) .< max_lattice]).^2 .+ 0*mg*zz[abs.(zz/π) .< max_lattice], label="lattice", color=:black, lw=2, alpha=0.7)
    
    wfn_scaling = 70
    plot!(
        zz[abs.(zz/π) .< max_lattice]/π, wfn_scaling*abs2.(ψ_nz0[abs.(zz/π) .< max_lattice]) .+ E0,
        fill= + E0, 
        lw=1,
        # xlabel="Lattice site", 
        label="nz = 0", 
        color=1, 
        alpha=1,
        )
    plot!(zz[abs.(zz/π) .< max_lattice]/π, wfn_scaling*abs2.(ψ_nz1[abs.(zz/π) .< max_lattice]) .+ E1,
        fill = E1,
        lw=1,
        alpha=1,
        # title=df[data_indx, "depth"], 
        # xlim=[-5.5, 5.5], 
        # xlabel="Lattice site",
        label="nz = 1", 
        color=palette(:tab20c)[5], 
    )
    plot!(;kwargs...)
    
    return fig
end

depth = 40
fig = draw_wfn_script(df, findfirst(df[:, "depth"] .== depth), 0, legend=:false, 
# ylims=(-18, 2), 
size=((3+3/8)*96*1.5, (3+3/8)*96*1), grid=:false, axis=:off, ticks=:false,
# ylabel="Energy", 
font="helvetica",
margins=20Plots.px
)
# Plots.pdf(plotsdir("wavefunctions", string(depth)*" Er.pdf"))

# fig = plot_eigen_spectrum(df,findfirst(df[:, "depth"] .== 10), xlims=[-3.5, 3.5], legend=false, grid=false, size=(400, 400), ylim=(-12, 2), axis=false, xticks=false, yticks=false)
# Plots.pdf(plotsdir("schematics", "lattice_schematics.pdf"))
fig