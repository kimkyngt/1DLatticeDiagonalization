using PhysicalConstants, PhysicalConstants.CODATA2018, Unitful, UnitfulRecipes, LaTeXStrings, Plots, LinearAlgebra, QuantumOptics, OhMyREPL

gr(dpi=300)

# c = SpeedOfLightInVacuum
# h = PlanckConstant
# kB = BoltzmannConstant
# e = ElementaryCharge 
# ħ = h/(2*π) 
# m_u = AtomicMassConstant
# m87Sr = 86.9088774970*m_u # Wikipedia
# me = ElectronMass
# a0 = BohrRadius
# μB = BohrMagneton
# Hartree = ħ^2/me/a0^2
# ϵ0 = VacuumElectricPermittivity
# μ0 = VacuumMagneticPermeability
# atm = StandardAtmosphere
# g = StandardAccelerationOfGravitation


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