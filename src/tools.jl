# Tools for diagonalization
using Plots
function get_mean_position(ψ, zz)
    """get quantum averaged position"""
	return real(ψ'*(zz.*ψ))
end

function get_rabi_freq(ψ1, ψ2)
    """Calculate Rabi frequency with normalized wavefunctions"""
	return ψ1'*(exp.(im*zz*kclock/k813).*ψ2)
end

function get_U(z, ρ, U0)
    """Give potential in a units of Er"""
    (-ustrip(U0)*(cos(z))^2*exp(-2*ρ^2/w0^2) + ustrip(m87Sr*g/k813)*z)/ustrip(Er)
end

function find_center_index(soln, zz)
    """Find x=0 eigen state's index"""
	zexpt = [get_mean_position(soln.vectors[:,ii], zz) for ii in range(1, size(soln.vectors)[1])]
	centered_indices = findall(abs.(zexpt) .< 0.1)
	return centered_indices
end

function plot_eigen_spectrum(F)
    maxbandnum = 3
	fig = plot([get_mean_position(real.(F.vectors[:,ii]), zz) for ii in range(1, numsites*maxbandnum)],
    real.(F.values[1:numsites*maxbandnum]),
		marker=:hline,markerstrokewidth=3, 
        markersize=6,lw=0, label=L"\langle z \rangle"*" vs Eigen energy")
	
	plot!(
        real.([z.data[ii, ii] for ii in 2:(num_z-1)]), 
        real.([U.data[ii, ii] for ii in 2:(num_z-1)]), 
        xlabel = L"k_{lat}z", ylabel = L"U (E_r)", 
        label="Potential to solve"
        )
	
	plot!(
        title="Depth = "*string(round(U_0/Er, digits=2))*" Er", 
	legend=:best
    )
    return fig
end

function get_sine_sq_exp(ψ, zz)
    return ψ'*((sin.(zz).^2).*ψ)
end
