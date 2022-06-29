using PhysicalConstants, PhysicalConstants.CODATA2018, Unitful, UnitfulRecipes, LaTeXStrings, Plots, LinearAlgebra, QuantumOptics	

c = SpeedOfLightInVacuum
h = PlanckConstant
kB = BoltzmannConstant
e = ElementaryCharge 
ħ = h/(2*π) 
m_u = AtomicMassConstant
m87Sr = 86.9088774970*m_u # Wikipedia
me = ElectronMass
a0 = BohrRadius
μB = BohrMagneton
Hartree = ħ^2/me/a0^2
ϵ0 = VacuumElectricPermittivity
μ0 = VacuumMagneticPermeability
atm = StandardAtmosphere
g = StandardAccelerationOfGravitation

λ813 = 813.4280e-9u"m"
k813 = 2*π/λ813
α813 = 286.0 * 4π*ϵ0*a0^3
klat = 2π/λ813
fclock = 429_228_004_229_873.0u"Hz"
λclock = uconvert(u"m", c/fclock)
kclock = 2π/λclock
Er = uconvert(u"J", ħ^2 * klat^2 / 2 / m87Sr)
w0 = 260e-6u"m"
