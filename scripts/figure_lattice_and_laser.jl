using DrWatson, DataFrames
include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

# if ~(@isdefined df)
#     df = collect_results(datadir("sims", "gcorrect"))
# else
#     print("df exist, skip loading \n")
# end


# function find_wfn(df, data_indx, siteindx;kwargs...)
#     H_eigen = df[data_indx, "solution"]
#     zz = df[data_indx, "zz"]
#     center_indx = find_center_index(H_eigen, zz, siteindx)
#     ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
#     ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])
#     return ψ_nz0, ψ_nz1, zz/π
# end

# depth = 13
# ψ_nz0, ψ_nz1, lattice_sites = find_wfn(df, findfirst(df[:, "depth"] .== depth), 0)

zz = range(-1.5, 1.5, length=300)
A = 0.2
B = 0.1
U = -0.3*(cos.(π*zz).^2)
mgz = B*zz

# lattice potential
plot(zz, mgz + U, color=:black, alpha=0.7)
offset = 0.5
plot!(zz, mgz + U .+ offset, color=:black,alpha=0.7)
elevelalpha = 1
# Energy levels
plot!([-1.35, -0.72], [-0.2, -0.2], lw=1.5, color=palette(:tab20c)[5], alpha=elevelalpha)
plot!([-1.23, -0.82], [-0.3, -0.3], lw=1.5, color=palette(:tab20c)[1], alpha=elevelalpha)
plot!([-1.35, -0.72], [-0.2, -0.2] .+ offset, lw=1.5, color=palette(:tab20c)[5], alpha=elevelalpha)
plot!([-1.23, -0.82], [-0.3, -0.3] .+ offset, lw=1.5, color=palette(:tab20c)[1], alpha=elevelalpha)

plot!([-1.35, -0.72] .+ 1, [-0.2, -0.2] .+ B, lw=1.5, color=palette(:tab20c)[5], alpha=elevelalpha)
plot!([-1.23, -0.82] .+ 1, [-0.3, -0.3] .+ B, lw=1.5, color=palette(:tab20c)[1], alpha=elevelalpha)
plot!([-1.35, -0.72] .+ 1, [-0.2, -0.2] .+ B .+ offset, lw=1.5, color=palette(:tab20c)[5], alpha=elevelalpha)
plot!([-1.23, -0.82] .+ 1, [-0.3, -0.3] .+ B .+ offset, lw=1.5, color=palette(:tab20c)[1], alpha=elevelalpha)

plot!([-1.35, -0.72] .+ 2, [-0.2, -0.2] .+ 2*B, lw=1.5, color=palette(:tab20c)[5], alpha=elevelalpha)
plot!([-1.23, -0.82] .+ 2, [-0.3, -0.3] .+ 2*B, lw=1.5, color=palette(:tab20c)[1], alpha=elevelalpha)
plot!([-1.35, -0.72] .+ 2, [-0.2, -0.2] .+ 2*B .+ offset, lw=1.5, color=palette(:tab20c)[5], alpha=elevelalpha)
plot!([-1.23, -0.82] .+ 2, [-0.3, -0.3] .+ 2*B .+ offset, lw=1.5, color=palette(:tab20c)[1], alpha=elevelalpha)

# annotation
annotate!(-1.3, -0.32, (L"n_z = 0", 8, :right, palette(:tab20c)[1]))
annotate!(-1.4, -0.22, (L"n_z = 1", 8, :right, palette(:tab20c)[5]))

# annotate!(-1.45, -0.2, (L"\left|g, n_z = 1\right\rangle", 8, :right, palette(:tab20c)[5]))
# annotate!(-1.45, -0.3+offset, (L"\left|e, n_z = 0\right\rangle", 8, :right, palette(:tab20c)[1]))
# annotate!(-1.45, -0.2+offset, (L"\left|e, n_z = 1\right\rangle", 8, :right, palette(:tab20c)[5]))
annotate!(-1.5, -0.1+offset, ("³P₀", 8, :right, :black))
annotate!(-1.5, -0.1, ("¹S₀", 8, :right, :black))

# arrows
plot!([-0.15, -0.15], [-0.3, -0.2+offset] .+ B, arrow=true, color=:blue, lw=1)
annotate!(-0.22, 0.0, Plots.text("BSB", 6, "helvetica", color=:blue, rotation=90))

plot!([-0., -0.0], [-0.3, -0.3+offset].+ B, arrow=true, color=:black, lw=1)
annotate!(0.08, 0.0.+ B, Plots.text("Carrier", 6, "helvetica", color=:black, rotation=90))

plot!([-0, 1], [-0.3, -0.3+offset+B].+ B, arrow=true, color=:black, lw=1)
annotate!(0.7, 0.05.+ B, Plots.text("WS+1", 6, "helveitica",color=:black, rotation=56))
annotate!(-1.3, 0.05.+ B, Plots.text(L"\lambda_{\mathrm{lat}}/2", 6,color=:black, rotation=0))

# fine tune
plot!(axis=false, legend=false, grid=false, ticks=false, xlims=(-2.5, 2), size = ((3+3/8)*96, (3+3/8)*96*2.3/4), background_color = RGBA(1,1,1,0))

Plots.pdf(plotsdir("lattice_diagram.pdf"))
plot!()