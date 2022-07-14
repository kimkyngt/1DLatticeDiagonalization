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

# Depth scan

depths =  [12.5:0.5:30;31:5:300]
print(length(depths))
for depth in depths
    sim_params = Dict{String, Any}(
        "numsites" => 41, 
        "numz" => 2000,     
        "depth" => depth,
        "radius" => 0,
    )
    soln = diagonalize_lattice(sim_params)
    @tagsave(datadir("sims", "gcorrect", savename(sim_params, "jld2")), soln)
    print(string(depth)*" Er done...")
end
print("Finished")
# firstsim = readdir(datadir("sims"))[1]
# wload(datadir("sims", firstsim))
