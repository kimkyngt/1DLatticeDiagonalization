using DrWatson, QuantumOptics, OrdinaryDiffEq, ProgressMeter

function get_spectrum(params::Dict)
    """
    Get a rabi lineshape. 
    See notebooks/WS_ladder_clock_test.jl
    """
	@unpack Omega0, Omega1, Omega2, detuning_max, detuning_num, T_rabi, Delta, initial_condition, N_lattice, reltol = params
	ϕ = 813*pi/698
	N_cutoff = N_lattice - 1
	b_fock = FockBasis(N_cutoff)
	b_spin = SpinBasis(1//2)
	b = b_fock ⊗ b_spin

	# Fundamental operators
	n = number(b_fock)
	sz = sigmaz(b_spin)
	sx = sigmax(b_spin)
	sy = sigmay(b_spin)

	function get_laser_Hamiltonian(Ω_k, k)
		H = 0*n⊗sz
		for m in range(0, N_cutoff-k)
			H += (
				(fockstate(b_fock, m)⊗dagger(fockstate(b_fock, m+k))
				+ fockstate(b_fock, m+k)⊗dagger(fockstate(b_fock, m))) 
				⊗ ( real(Ω_k*exp(im*m*ϕ))*sx - imag(Ω_k*exp(im*m*ϕ))*sy )
			)
		end
		# Filter double counting case
		if k == 0
			H = 0.5*H
		end
		return 0.5*H
	end
	
	function get_H_tot(delta)
		# Constructing the Hamiltonian
		Hatom = Delta*n⊗one(b_spin) - 0.5*delta*one(b_fock)⊗(one(b_spin) + sz)
		Hlaser = get_laser_Hamiltonian(Omega0, 0) + get_laser_Hamiltonian(Omega1, 1) + get_laser_Hamiltonian(Omega2, 2)
		H_tot = Hatom+Hlaser
		return H_tot
	end


	function get_exc_fac(detuning)
		# Solve dynamics and get final excitation fraction

        if initial_condition == "pure"
            # println("Intial condition: pure")
            ψ_0 = normalize!(
                Ket(b_fock, vcat([0 for ii in 1:N_cutoff/2], 1, [0 for ii in 1:N_cutoff/2],)) ⊗ spindown(b_spin)
            )
        elseif initial_condition =="plus"
            # println("Intial condition: plus")
            ψ_0 = normalize!(
                Ket(b_fock, vcat([0 for ii in 1:N_cutoff/2], 1, 1, [0 for ii in 2:N_cutoff/2],)) ⊗ spindown(b_spin)
            )
        end
		tspan = range(0, T_rabi, length=2)
		tout, ψ_t = timeevolution.schroedinger(tspan, ψ_0, get_H_tot(detuning), alg = DP5(), reltol=reltol, maxiters=1e9)

		exc_frac = expect(one(b_fock)⊗(sz+one(b_spin))/2, ψ_t[end])
		return exc_frac
	end


	rho_ee = zeros(detuning_num)
	detunings = range(-detuning_max, detuning_max, detuning_num)
	@showprogress 1 "Computing spectrum..." for ii in 1:detuning_num
		rho_ee[ii] = get_exc_fac(detunings[ii])
	end
	result = Dict([
		("parameters", params),
		("detunings", detunings),
		("exc_frac", rho_ee)
	])
	return result
end

# function compute_maximum_linepulling(result::Dict, lockpoint)
# 	@unpack detunings, exc_fac = result
# 	minimum(abs.(exc_fac[detunings < 0] .- lockpoint))
# end

