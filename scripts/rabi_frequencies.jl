using DrWatson, DataFrames

include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))

df = collect_results(datadir("sims"))

function draw_wfn(df, data_indx)
    H_eigen = df[data_indx, "solution"]
    zz = df[data_indx, "zz"]
    center_indx = find_center_index(H_eigen, zz)
    ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
    ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])
    fig = plot(zz/π, real.(ψ_nz0), dpi=300, 
        title=df[data_indx, "depth"], 
        xlim=[-5, 5]
    )
    plot!(zz/π, real.(ψ_nz1))
    return fig
end

draw_wfn(df, findfirst(df[:, "depth"] .== 10))

function get_rabi_frequency(df, data_indx)
    H_eigen = df[data_indx, "solution"]
    zz = df[data_indx, "zz"]
    center_indx = find_center_index(H_eigen, zz)

    rabi_freqs = zeros(Complex, 8)
    ψ_nz0 = H_eigen.vectors[:, center_indx[1]]
    for ii in range(1, 4)
        ψ_nz0_ws = H_eigen.vectors[:, center_indx[1]+ii-1]
        rabi_freqs[ii] = ψ_nz0' *(exp.(im*zz).*ψ_nz0_ws)
    end
    
    ψ_nz1 = H_eigen.vectors[:, center_indx[2]]
    for ii in range(1, 4)
        ψ_nz1_ws = H_eigen.vectors[:, center_indx[2]+ii-1]
        rabi_freqs[ii+4] = ψ_nz1' *(exp.(im*zz).*ψ_nz1_ws)
    end

    return rabi_freqs
end

depths = zeros(size(df)[1])
rabi_freqs = zeros(Complex, 8, size(df)[1])

for ii in range(1, size(df)[1])
    depths[ii] = df[ii, "depth"]
    rabi_freqs[:, ii] = get_rabi_frequency(df, ii)
end

fig = plot(dpi=300)
for ii in range(1, 8)
    plot!(depths, abs2.(rabi_freqs[ii, :]) , 
        seriestype=:line, 
        yscale=:log10, 
        xscale=:log10,
        xlims=(8, 30),
        ylims=(1e-4, 1),    
        xticks=(5:5:30, 5:5:30)
    )
end
fig


plot_eigen_spectrum(df, findfirst(df[:, "depth"] .== 9))