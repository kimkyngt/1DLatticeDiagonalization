
using DrWatson
print(projectname()*" activated \n")
include(srcdir("intialize.jl"))
include(srcdir("tools.jl"))


function get_reduction_factor_sag(V0)
    """V0 in units of Er = no units"""
    Tr = 22.52e-9*V0^0.62 * 1u"K"
    ωR = uconvert(u"m*s^-1",√(4*V0*Er/m87Sr))/w0
    A  = uconvert(u"m^-2", 1/2*m87Sr*ωR^2/kB/Tr)
    Ap = 2/w0^2 
    B = uconvert(u"m^-1", m87Sr*9.796u"m/s^2"*sin(θtilt)/(kB*Tr))
    ζ = (1/(1+Ap/A)) * exp(-B^2*Ap/(4*(A+A')*A))
   return ζ
end

function get_reduction_factor(V0)
    """V0 in units of Er = no units"""
    Tr = 22.52e-9*V0^0.62 * 1u"K"
    ωR = uconvert(u"m*s^-1",√(4*V0*Er/m87Sr))/w0

    A  = uconvert(u"m^-2", 1/2*m87Sr*ωR^2/kB/Tr)
    Ap = 2/w0^2 
    ζ = (1/(1+Ap/A))
   return ζ
end

function get_center_reduction(V0)
    """V0 in units of Er = no units"""
    ζ = exp( -uconvert(NoUnits, (m87Sr*9.796u"m/s^2"*sin(θtilt)*w0)^2/8/Er^2) /(V0- √(V0))^2)
   return ζ
end

V0s  = range(4, 100, length=1000)
plot(V0s, [get_reduction_factor_sag(V0) for V0 in V0s], xscale=:log10, label = L"\zeta"*" with sagging", lw=2)
plot!(V0s, [get_reduction_factor(V0) for V0 in V0s], xscale=:log10, label = L"\zeta", lw=2, xlabel="Lattice depth (Er)", legend=:bottomright)

plot!(V0s, [get_center_reduction(V0) for V0 in V0s], xscale=:log10, label = "Center reduction", lw=2, xlabel="Lattice depth (Er)", legend=:bottomright)

