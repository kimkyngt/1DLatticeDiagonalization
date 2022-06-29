using DataFrames

include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))

df = collect_results(datadir("sims"))

data_num = size(df)[1]
depths = zeros(data_num)
shifts = zeros(data_num)

for ii in range(1, data_num)

    H_eigen = df[ii, "solution"]
    zz = df[ii, "zz"]
    depths[ii] = df[ii, "depth"]

    center_indx = find_center_index(H_eigen, zz)
    ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
    ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])

    sinkz_nz0 = get_sine_sq_exp(ψ_nz0, zz)
    sinkz_nz1 = get_sine_sq_exp(ψ_nz1, zz)
    
    shifts[ii] = 0.0012*depths[ii]*(sinkz_nz1 - sinkz_nz0)/ustrip(fclock)
end


# fig_spectrum = plot_eigen_spectrum(H_eigen)
# Plots.pdf(fig_spectrum, plotsdir("spectrums.pdf"))
