using DrWatson, QuantumOptics, LinearAlgebra
include(srcdir("tools.jl"))
include(srcdir("initialize.jl"))

# basis generation
function diagonalize_Unz_Cartesian(sim_params::Dict)
    @unpack xmax, numx, depth, nz = sim_params
    bx = PositionBasis(-xmax, xmax, numx) # x in units of 1/k813
    xx = samplepoints(bx)
    bpx = MomentumBasis(bx)

    # Basic operators
    x  = position(bx)
    Px = momentum(bpx)

    # transformation operators
    Txpx = transform(bx, bpx)
    Tpxx = transform(bpx, bx)

    # Hamiltonian construction
        # Kinetic part
    Hkin = Px^2
    Hkin_FFT = LazyProduct(Txpx, Hkin, Tpxx)
    Hkin_FFT = dense(Hkin_FFT)

    # lattice potential
    U = x*0
    U_0 = depth

    #     # Hard boundary condition
    # U.data[1, 1] = 1000000
    # U.data[end, end] = 1000000
    for ii = 1:numx
        U.data[ii, ii] = get_Unz(x.data[ii, ii], U_0, w0*k813, nz) # x in units of w_0
    end

    H_tot = dense(LazySum(Hkin_FFT, U))
    @time H_eigen = eigen(H_tot.data)

    output_d = copy(sim_params)
    output_d["solution"] = H_eigen
    output_d["xx"] = xx
    output_d["potential"] = U.data.nzval
    return output_d
end

# Depth scan
# depths =  [12.5:0.5:30;31:5:300]
depth = 12
# for depth in depths
    sim_params = Dict{String, Any}(
        "xmax" => k813*w0, # in units of 1/k813
        "numx" => 2000,     
        "depth" => depth,
        "nz" => 0,
    )
    soln = diagonalize_Unz_Cartesian(sim_params)
    # @tagsave(datadir("sims", "gcorrect", savename(sim_params, "jld2")), soln)
    print(string(depth)*" Er done...")
# end
print("Finished")


# Plot results
plot(soln["xx"][2:end-1]/k813, real.(soln["potential"][2:end-1]),xunit=u"μm", label="U(x)", ylabel="Eᵣ")

indx = 250
E = real.(soln["solution"].values[indx])
y = real.(soln["solution"].vectors[:, indx])
plot!(soln["xx"]/k813, abs(E)/3*y/maximum(y) .+ E, 
xlim=(-100, 100),
label="wavefunction",
)
hline!([real(soln["solution"].values[1])], label="U(0)")
hline!([real(soln["solution"].values[1])+1], label="U(0) + 1Er")
