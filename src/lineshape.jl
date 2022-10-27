using DrWatson, QuantumOptics, OrdinaryDiffEq, ProgressMeter, Interpolations

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

function get_discriminator(result::Dict, modulation_depth::Number)
	# modulation depth in Hz
	@unpack detunings, exc_frac, parameters = result
	@unpack T_rabi = parameters
	dx = detunings[2] - detunings[1]
	Nshift = Int(round(2*pi*modulation_depth/dx))
	diff_frac = circshift(exc_frac, Nshift) - circshift(exc_frac, -Nshift)
	diff_frac = diff_frac[abs.(detunings) .< 2*modulation_depth]
	detunings = detunings[abs.(detunings) .< 2*modulation_depth]
	p = sortperm(detunings)
	frac_to_detuning = LinearInterpolation(diff_frac[p], detunings[p], extrapolation_bc=0)
	# _detunings = range(minimum(detunings)*0.9, maximum(detunings)*0.9, length=100)

	return frac_to_detuning, diff_frac
end

function get_spectrum_loss(params::Dict)
    """
    Get a rabi lineshape with a decay within the system.
    See notebooks/WS_ladder_clock_loss_test.jl
    """
	@unpack Omega0, Omega1, Omega2, detuning_max, detuning_num, T_rabi, Delta, initial_spin, N_lattice, reltol , gamma = params
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

	function get_J()
		# Jump operator for the loss
		J = fockstate(b_fock, 0)⊗dagger(fockstate(b_fock, 1))
		for jj in 1:(N_cutoff-1)
			J += fockstate(b_fock, jj)⊗dagger(fockstate(b_fock, jj+1))
		end
		J = J⊗one(b_spin)
	end

	function get_exc_fac(detuning)
		# Solve dynamics and get final excitation fraction

        if initial_spin == "down"
            # println("Intial condition: pure")
            ψ_0 = normalize!(
                Ket(b_fock, vcat([0 for ii in 1:N_cutoff/2], 1, [0 for ii in 1:N_cutoff/2],)) ⊗ spindown(b_spin)
            )
        elseif initial_spin =="up"
            # println("Intial condition: plus")
            ψ_0 = normalize!(
                Ket(b_fock, vcat([0 for ii in 1:N_cutoff/2], 1, [0 for ii in 1:N_cutoff/2],)) ⊗ spinup(b_spin)
            )
        end
		tspan = range(0, T_rabi, length=2)
		tout, ψ_t = timeevolution.master(tspan, ψ_0, get_H_tot(detuning), [sqrt(gamma)*get_J()], alg = DP5(), reltol=reltol, maxiters=1e9)

		exc_frac = expect(one(b_fock)⊗(sz+one(b_spin))/2, ψ_t[end])
		return abs(exc_frac)
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


