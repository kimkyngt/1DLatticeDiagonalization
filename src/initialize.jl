using PhysicalConstants, PhysicalConstants.CODATA2018, Unitful, UnitfulRecipes, LaTeXStrings, Plots, LinearAlgebra, QuantumOptics, OhMyREPL

λ813 = 813.4280e-9u"m"
k813 = 2*π/λ813
α813 = 286.0 * 4π*ε_0*a_0^3
fclock = 429_228_004_229_873.0u"Hz"
λclock = uconvert(u"m", c_0/fclock)
kclock = 2π/λclock
Er = uconvert(u"J", ħ^2 * k813^2 / 2 / m87Sr)
w0 = 260e-6u"m"

θtilt = acos(867.69/867.75) # tilt angle of the lattice

print("initialize.jl imported \n")