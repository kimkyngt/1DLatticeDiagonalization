# Tools for diagonalization
using Plots, JSON, LaTeXStrings
function get_mean_position(ψ, zz)
    """get quantum averaged position"""
	return real(ψ'*(zz.*ψ))
end

function get_U(z, ρ, U0)
    """Give potential in a units of Er"""
    (-ustrip(U0)*(cos(z))^2*exp(-2*ρ^2/w0^2) .+ 867/868.23*ustrip(m87Sr*g_n/k813)*z)/ustrip(Er)
end

function get_U(U_0, ρ, w_0)
    """get_U(U_0, ρ, w_0)"""
    -U_0*exp.(-2*ρ.^2/w_0^2)
end

function get_Unz(r, U0, w_0, nz)
    """Give potential in a units of Er"""
    return -U0*(exp.(-r.^2/w_0.^2) .- sqrt(1/U0)*(nz+1/2)).^2 .+ (nz+1/2)^2  .- 0.5*(nz^2+nz+1/2) 
end

function get_Unz_sag(r, U0, w_0, nz, θ=θtilt)
    """Give potential with gravitational sagging in a units of Er"""
    return -U0*(exp.(-r.^2/w_0.^2) .- sqrt(1/U0)*(nz+1/2)).^2 .+ (nz+1/2)^2  .- 0.5*(nz^2+nz+1/2) .+ ustrip(m87Sr*g_0*sin(θ)/Er)*r
end

function find_center_index(soln, zz, siteindx)
    """Find x=0 eigen state's index"""
	zexpt = [get_mean_position(soln.vectors[:,ii], zz) for ii in range(1, size(soln.vectors)[1])]
	centered_indices = findall(abs.(zexpt) .< pi/2)
	return centered_indices
end

function draw_wfn(df, data_indx, siteindx;kwargs...)
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
        xlim=[-5, 5], 
        xlabel="Lattice site",
        label="nz = 0", 
        color=palette(:tab10)[2], 
        ;kwargs...
    )
    
    return fig
end


function plot_eigen_spectrum(df, ii)
    """Plot eigenvalue specturm to check the calculation"""
    F = df[ii, "solution"]
    zz = df[ii, "zz"]
    maxbandnum = 3
    numsites = df[ii, "numsites"]
    U_0 = df[ii, "depth"]
	fig =plot(
        zz/π, [get_U(z, 0u"m", U_0*Er) for z in zz], 
        label=L"U(z, r=0)",
        color=:black,
        alpha=0.7,
        size=(600, 300)
    )
    plot!(
        [get_mean_position(real.(F.vectors[:,ii]), zz) for ii in range(1, numsites*maxbandnum)]/π,
        real.(F.values[1:numsites*maxbandnum]),
		marker=:hline, markerstrokewidth=3, 
        markersize=6,
        st = :scatter,
        color=1,
        label=L"E_n(\langle z \rangle)"
    )
    plot!(
        legend=:bottomright,
        xlabel="Lattice site, "*L"2z / \lambda_L",
        ylabel="Energy (Eᵣ)",
        grid=:true
        )
    return fig
end

function get_sine_sq_exp(ψ, zz)
    return ψ'*((sin.(zz).^2).*ψ)
end

function get_cos_sq_exp(ψ, zz)
    return ψ'*((cos.(zz).^2).*ψ)
end

function get_sine_sq_exp_density(ψ, zz)
    return abs2.(ψ)'*(sin.(zz).^2)
end


function get_rabi_frequency(df, data_indx)
    """Get list of rabi frequency. nz=0 WS0, 1, 2, 3 and nz=1 WS0, 1, 2, 3"""
    H_eigen = df[data_indx, "solution"]
    zz = df[data_indx, "zz"]
    dz = zz[2] - zz[1]
    center_indx = find_center_index(H_eigen, zz, 0)

    rabi_freqs = zeros(Complex, 8)

    ψ_nz0 = H_eigen.vectors[:, center_indx[1]]
    expikr_ψ = exp.(im*zz*(kclock/k813)) .* ψ_nz0
    for ii in range(1, 4)
        # shift centered wavefunction to get different order ws states.
        ψ_nz0_ws = circshift(ψ_nz0, (ii-1)*round(π/dz))
        rabi_freqs[ii] = ψ_nz0_ws' * expikr_ψ
    end
    
    ψ_nz1 = H_eigen.vectors[:,center_indx[2]]
    expikr_ψ = exp.(im*zz*(kclock/k813)) .* ψ_nz1
    for ii in range(1, 4)
        ψ_nz1_ws = circshift(ψ_nz1, (ii-1)*round(π/dz))
        rabi_freqs[ii+4] = ψ_nz1_ws' * expikr_ψ
    end

    return rabi_freqs
end


function get_Tr(U_0)
    if U_0>15
        return 42*√(U_0)*u"nK"
    else
        return (-45.2+14.1*U_0)*u"nK"
    end
end

function get_rabi_carrier_frequency(df, data_indx)
    """Get list of carrier rabi frequencies up to nz=2"""
    H_eigen = df[data_indx, "solution"]
    zz = df[data_indx, "zz"]
    center_indx = find_center_index(H_eigen, zz, 0)
    function calc_Omega(ψ)
        ψ' *(exp.(im*zz*(kclock/k813)) .*ψ) /(ψ'*ψ)
    end

    rabi_freqs = [calc_Omega(H_eigen.vectors[:, center_indx[i]]) for i=[1, 2, 3]]

    return rabi_freqs
end

function get_rho(r::AbstractArray, U_0, T_r, nz=0)
    """r is number in meter, U_0 is number in Er, T_r is u"K" """
    # ω̄ = 2*π*2.44u"Hz"*(2*√(U_0) - 1/2) # for the ground state
    U = get_Unz(r, U_0, w_0, nz)
    U[U .> 0] .= Inf
    A = U*uconvert(NoUnits, Er/k_B/T_r)
    return exp.(-A)
end

function get_rho(r::Number, U_0, T_r, nz=0)
    """r is number in meter, U_0 is number in Er, T_r is u"K" """
    # ω̄ = 2*π*2.44u"Hz"*(2*√(U_0) - 1/2) # for the ground state
    A = get_Unz(r, U_0, w_0, nz)*uconvert(NoUnits, Er/k_B/T_r)
    if A < 0
        rho = exp.(-A)
    else
        rho = 1e-5
    end
    return rho
end

function overlay_data!(filter_string="nz1")
    data = JSON.parse(JSON.parsefile(datadir("exp_pro", "rabi_freqs_20220625.json")))
    Ω0 = maximum(data["nz0 ws0"]["rabi frequency"])

    for key in keys(data)
        if key[1:3] == filter_string
            plot!(data[key]["depth Erec"],abs2.(data[key]["rabi frequency"]/Ω0), marker=:circle, xscale=:log10, yscale=:log10, label=key, lw=0)
        end
    end

    data = JSON.parse(JSON.parsefile(datadir("exp_pro", "rabi_freqs_20220626.json")))
    for key in keys(data)
        if key[1:3] == filter_string
            plot!(data[key]["depth Erec"],abs2.(data[key]["rabi frequency"]/Ω0), marker=:circle, xscale=:log10, yscale=:log10, label=key, lw=0)
        end
    end

end

function get_Tr(u, state)
    """get measured, fitted radial temperature, returns in K"""
    if state == "nz0"
        Tr = 3.18e-8*u^0.65*1u"K"
    elseif state == "nz1"
        Tr = 1.53e-8*u^0.65*1u"K"
    end
    
    return Tr
end

function get_axial_eigen_energy(df, data_indx; num_state=5)
    """Axial eigen energy returns in Er"""
    H_eigen = df[data_indx, "solution"]
    zz = df[data_indx, "zz"]
    depth = df[data_indx, "depth"]
    center_indx = find_center_index(H_eigen, zz, 0)

    axial_energy = zeros(num_state)
    for ii in range(1, num_state)
        axial_energy[ii] = real(H_eigen.values[center_indx[ii]])
    end
    return axial_energy, depth
end

function get_periodic_zpower(zz, power)
    (zz .- pi*floor.(zz/pi .+ 1/2) ).^power
end

print("tools.jl imported \n")