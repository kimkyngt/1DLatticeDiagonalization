using DrWatson, JSON, Plots, DataFrames
include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

# Experimental data
json_string = read(datadir("exp_raw", "axial_scan_data.json"), String)
data = JSON.parse(json_string)
detuning = data["nz0"]["detuning"]
rho_ee = data["nz0"]["exc_frac"]
p = sortperm(detuning)
plot(detuning[p]*1e-3, rho_ee[p])


# Calculate Rabi frequency ratios 
if ~(@isdefined df)
    df = collect_results(datadir("sims", "gcorrect"))
else
    print("df exist, skip loading \n")
end

function get_rabi_freq_BSB(df, data_indx)
    """Get list of rabi frequency of BSB"""
    H_eigen = df[data_indx, "solution"]
    zz = df[data_indx, "zz"]
    dz = zz[2] - zz[1]
    center_indx = find_center_index(H_eigen, zz, 0)
    rabi_freqs = zeros(Complex, 4)

    ψ_nz0 = H_eigen.vectors[:, center_indx[1]]
    ψ_nz1 = H_eigen.vectors[:,center_indx[2]]
    expikr_ψ = exp.(im*zz*(kclock/k813)) .* ψ_nz0
    for ii in range(1, 4)
        ψ_nz1_ws = circshift(ψ_nz1, (ii-1)*round(π/dz))
        rabi_freqs[ii+4] = ψ_nz1_ws' * expikr_ψ
    end

    return rabi_freqs
end

get_rabi_freq_BSB(df, findfirst(df[:, "depth"] .== 10))