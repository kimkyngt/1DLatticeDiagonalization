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
    rabi_freqs = zeros(Complex, 15)

    ψ_nz0 = H_eigen.vectors[:, center_indx[1]]
    ψ_nz1 = H_eigen.vectors[:,center_indx[2]]
    expikr_ψ = exp.(im*zz*(kclock/k813)) .* ψ_nz0
    for ii in range(1, length(rabi_freqs))
        ψ_nz1_ws = circshift(ψ_nz1, (ii- 8)*round(π/dz))
        rabi_freqs[ii] = ψ_nz1_ws' * expikr_ψ
    end

    # Get energy gap
    Egap = H_eigen.values[center_indx[2]] - H_eigen.values[center_indx[1]]
    return rabi_freqs, abs(Egap)
end

function get_axial_sideband_shape(detuning, Tr, U_0)
    """
    detuning: number, Hz
    Tr: number, K
    U_0: number, Er
    offset: number, Hz
    """
    function heaviside(t)
        0.5 * (sign(t) + 1)
    end
    # Equation (A4) from Blatt2009PRA
    νrec = ustrip(Er/h)
    νz = 2*νrec*sqrt(U_0)
    γ̃ = νz - νrec*1 # from nz = 0 to nz = 1
    α = γ̃/νrec * νz * ustrip(h/k_B) / Tr
    σ = (1 - detuning/γ̃) * exp(-α*(1-detuning/γ̃)) * heaviside(γ̃-detuning)
    return σ
end
# paremeters
Tr = 145e-9
U_0 = 10.4
Omega, Egap = get_rabi_freq_BSB(df, findfirst(df[:, "depth"] .== U_0))

delta = range(-30e3, 30e3, length=6001)
νrec = ustrip(Er/h)
νz = 2*νrec*sqrt(U_0)
γ̃ = νz - νrec*1

yfit = zeros(length(delta), 15)
for ii in range(1, length(Omega))
    offset = -867.69*(ii-8)
    yfit[:, ii] = abs.(Omega[ii])*[get_axial_sideband_shape(d - offset, Tr, U_0) for d in delta]
end
y = sum(yfit, dims=2)

# plot and save plot
fig = plot(detuning[p]*1e-3, rho_ee[p], label="data")
plot!(delta/1e3, (y/maximum(y))*0.16, 
    label="fit", lw=1.5
)
plot!(delta/1e3, yfit*(1/maximum(y))*0.16, 
    label="", lw=1, color=:grey
)
plot!(
    xlabel="Clock detuning (kHz)", 
    ylabel="Excitation fraction",
    xlims=(10, 23), ylims=(-0.02, 0.2), 
    margin=4Plots.px
)

savefig(fig, plotsdir("axial_sideband", "10Er_fit.pdf"))
fig

# save data
datadict = Dict([
    ("data", Dict(
        "detuning"=>detuning[p],
        "frac"=>rho_ee[p],
    )),
    "theory"=>Dict(
        "detuning"=>delta,
        "frac"=>yfit*(1/maximum(y))*0.16
    )
]
)
open(datadir("axial_scan_fit.json"),"w") do f
    JSON.print(f, datadict)
end