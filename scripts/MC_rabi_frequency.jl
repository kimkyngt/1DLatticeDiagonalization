using DrWatson, QuantumOptics
include(srcdir("initialize.jl"))

numsites = 22
zmax = numsites*π/2 # units in k_lat z
num_z = 2^8
U_0 = 50*Er
# basis generation
bz = PositionBasis(-zmax, zmax, num_z); zz = samplepoints(bz)
bpz = MomentumBasis(bz); ppz = samplepoints(bpz)

# transformation operators
Tzpz = transform(bz, bpz)
Tpzz = transform(bpz, bz)

# Basic operators
z  = position(bz)
Pz = momentum(bpz)

# Hamiltonian construction
## Kinetic part
Hkin = Pz^2
Hkin_FFT = dense(LazyProduct(Tzpz, Hkin, Tpzz))

function get_U(z, ρ, U0)
    """Give potential in a units of Er"""
    (-ustrip(U0)*(cos(z))^2*exp(-2*ρ^2/w0^2) + ustrip(m87Sr*g/k813)*z)/ustrip(Er)
end


function get_Unz_rho(ρ)
	# ρ = 100e-6u"m"
	U = z*0
		# Hard boundary condition
	U.data[1, 1] = 1000000
	U.data[end, end] = 1000000
	for ii = 2:(num_z-1)
		U.data[ii, ii] = get_U(z.data[ii, ii], ρ, U_0)
	end

	H_tot = dense(LazySum(Hkin_FFT, U))
	F = eigen(H_tot.data)
	return F.values
end


ρlist = range(0u"m", w0/√2 *2.5, length = 30)
Unz = [get_Unz_rho(ρ) for ρ in ρlist]
plot()
for ii in 1:5
    y = [ real(Unz[jj][ii]) for jj in 1:length(ρlist) ]
    plot!(ρlist*(√2 /w0), y, lw=2)
end
plot!(title="Beloy Fig.2", xlims=(0, 2.5), ylims=(-50, 10), legend=false, minorticks=5)

