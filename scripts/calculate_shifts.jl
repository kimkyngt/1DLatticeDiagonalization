H_eigen = soln["solution"]
zz = soln["zz"]

center_indx = find_center_index(H_eigen)
ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])
fig_wfn = plot(zz, ψ_nz0, label=L"n_z=0", dpi=300)
plot!(zz, ψ_nz1,
    xlabel=L"kz",
    label=L"n_z=1"
)

sinkz_nz0 = get_sine_sq_exp(ψ_nz0)
sinkz_nz1 = get_sine_sq_exp(ψ_nz1)
shift = 0.0012*depth*(sinkz_nz1 - sinkz_nz0)/ustrip(fclock)


# fig_spectrum = plot_eigen_spectrum(H_eigen)
# Plots.pdf(fig_spectrum, plotsdir("spectrums.pdf"))
