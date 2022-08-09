using DrWatson, JSON, Plots, DataFrames
include(srcdir("initialize.jl"))
include(srcdir("tools.jl"))

# Experimental data
json_string = read(datadir("exp_raw", "axial_scan_data.json"), String)
data = JSON.parse(json_string)
detuning = data["nz0"]["detuning"]
rho_ee = data["nz0"]["exc_frac"]
p = sortperm(detuning)


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
    rabi_freqs = zeros(Complex, 6)

    ψ_nz0 = H_eigen.vectors[:, center_indx[1]]
    ψ_nz1 = H_eigen.vectors[:,center_indx[2]]
    expikr_ψ = exp.(im*zz*(kclock/k813)) .* ψ_nz0
    for ii in range(1, length(rabi_freqs))
        ψ_nz1_ws = circshift(ψ_nz1, (ii-1)*round(π/dz))
        rabi_freqs[ii] = ψ_nz1_ws' * expikr_ψ
    end

    # Get energy gap
    Egap = H_eigen.values[center_indx[2]] - H_eigen.values[center_indx[1]]
    return rabi_freqs, abs(Egap)
end

Omega, Egap = get_rabi_freq_BSB(df, findfirst(df[:, "depth"] .== 20))

# plot(detuning[p]*1e-3, rho_ee[p])
plot!([Egap*ustrip(Er/h)*1e-3 + 867.69*(ii-1)*1e-3 for ii in range(1, length(Omega))], [abs(Omega[ii]) for ii in range(1, length(Omega))], st=:stem)
