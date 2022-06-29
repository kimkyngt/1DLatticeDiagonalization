using DrWatson
print(projectname())
include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))

	
# basis generation
function diagonalize_lattice(sim_params::Dict)
    @unpack numsites, numz, depth, radius = sim_params
    zmax = numsites*π/2
    bz = PositionBasis(-zmax, zmax, numz); zz = samplepoints(bz)
    bpz = MomentumBasis(bz); ppz = samplepoints(bpz)

    # Basic operators
    z  = position(bz)
    Pz = momentum(bpz)

    # transformation operators
    Tzpz = transform(bz, bpz)
    Tpzz = transform(bpz, bz)

    # Hamiltonian construction
        # Kinetic part
    Hkin = Pz^2
    Hkin_FFT = LazyProduct(Tzpz, Hkin, Tpzz)
    Hkin_FFT = dense(Hkin_FFT)

    # lattice potential
    U = z*0
    U_0 = depth*Er
        # Hard boundary condition
    U.data[1, 1] = 1000000
    U.data[end, end] = 1000000
    for ii = 2:(numz-1)
        U.data[ii, ii] = get_U(z.data[ii, ii], radius*1u"m", U_0)
    end

    H_tot = dense(LazySum(Hkin_FFT, U))
    @time H_eigen = eigen(H_tot.data)

    output_d = copy(sim_params)
    output_d["solution"] = H_eigen
    output_d["zz"] = zz

    return output_d
end

# Computation parameters
sim_params = Dict{String, Any}(
    "numsites" => 50, 
    "numz" => 2^10,     
    "depth" => 9,
    "radius" => 0,
)

soln = diagonalize_lattice(sim_params)
H_eigen = soln["solution"]
zz = soln["zz"]

center_indx = find_center_index(H_eigen)
ψ_nz0 = real.(H_eigen.vectors[:, center_indx[1]])
ψ_nz1 = real.(H_eigen.vectors[:, center_indx[2]])
fig_wfn = plot(zz, ψ_nz0, label=L"n_z=0", dpi=300)
plot!(zz, ψ_nz1,
    xlabel=L"kz",
    label=L"n_z=1"
)

sinkz_nz0 = get_sine_sq_exp(ψ_nz0)
sinkz_nz1 = get_sine_sq_exp(ψ_nz1)
shift = 0.0012*depth*(sinkz_nz1 - sinkz_nz0)/ustrip(fclock)


# fig_spectrum = plot_eigen_spectrum(H_eigen)
# Plots.pdf(fig_spectrum, plotsdir("spectrums.pdf"))

@tagsave(datadir("sims", savename(sim_params, "jld2")), soln)
firstsim = readdir(datadir("sims"))[1]
wload(datadir("sims", firstsim))
