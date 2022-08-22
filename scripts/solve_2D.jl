using DrWatson, QuantumOptics, Plots, LinearAlgebra

# Parameters
w0 = 260e-6
rmax = 3
zmax = 3

nr = 32
br = PositionBasis(-rmax, rmax, nr)

nz = 32
bz = PositionBasis(-zmax, zmax, nz)

bpr = MomentumBasis(br)
bpz = MomentumBasis(bz)

r = position(br)⊗one(bz)
z = one(br)⊗position(bz)
Pr = momentum(bpr)⊗one(bpz)
Pz = one(bpr)⊗momentum(bpz)

compbx = br⊗bz
compbp = bpr⊗bpz

Txp = transform(compbx, compbp)
Tpx = transform(compbp, compbx)

Hkin = Pr^2/2 + Pz^2/2 # kinetic energy in momentum space
@time Hkin_FFT = LazyProduct(Txp, Hkin, Tpx)# lazy tensor for the split-step method
Hkin_FFT = dense(Hkin_FFT)
U = 0.1*r^2 + 0.01*z^2 + 100*r^2*z^2 + 10*z^4

H = Hkin_FFT + U

@time E = eigen(collect(H.data))
plt = [heatmap(transpose(reshape(abs2.(E.vectors[:, ii]), (nr, nz))), colorbar=false, axis=false)
 for ii in range(1, 16)]
plot(plt..., layout=(4, 4), size=(1000, 1000))
