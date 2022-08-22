using DrWatson, Interpolations, QuadGK, LsqFit, Plots, LaTeXStrings

function get_rho(r, U_0, T_r)
    """r is number in meter, U_0 is number in Er, T_r is u"K" """
    # ω̄ = 2*π*2.44u"Hz"*(2*√(U_0) - 1/2) # for the ground state
    A = -(U_0*exp(-2*(r^2/w_0^2)) - sqrt(U_0)*exp(-r^2/w_0^2) + 1/4)*uconvert(NoUnits, Er/k_B/T_r)
    rho = exp(-A)
    return rho
end
function get_Tr(U_0)
    if U_0>15
        return 42*√(U_0)*u"nK"
    else
        return (-45.2+14.1*U_0)*u"nK"
    end
end
# Parameters
U_0 = 4
T_r = get_Tr(U_0)*1.0
νr = uconvert(u"Hz", sqrt(U_0*Er/m87Sr/pi^2/(w_0*1u"μm")^2))
w_0 = 260e-6 # cavity waist
rmax = 300e-6 # in 

rr = range(0, rmax, length=100)
plot(rr*1e6, -U_0*exp.(-2*rr.^2/w_0^2), label="-U(r)/U₀")
plot!(rr*1e6, -(U_0*exp.(-2*(rr.^2/w_0^2)) .- 2*(0+1/2)sqrt(U_0)*exp.(-rr.^2/w_0^2) .+ (0^2 + 0 + 1/2)*1/2), label="nz=0")
plot!(rr*1e6, -(U_0*exp.(-2*(rr.^2/w_0^2)) .- 2*(1+1/2)sqrt(U_0)*exp.(-rr.^2/w_0^2) .+ (1^2 + 1 + 1/2)*1/2), label="nz=1")
plot!(rr*1e6, -(U_0*exp.(-2*(rr.^2/w_0^2)) .- 2*(2+1/2)sqrt(U_0)*exp.(-rr.^2/w_0^2) .+ (2^2 + 2 + 1/2)*1/2), label="nz=2")
# plot!(rr*1e6, -(U_0*exp.(-2*(rr.^2/w_0^2)) .- 2*(3+1/2)sqrt(U_0)*exp.(-rr.^2/w_0^2) .+ (3^2 + 3 + 1/2)*1/2), label="nz=3")
plot!(rr*1e6, U_0/2*[get_rho(r, U_0, T_r) for r in rr] / get_rho(0, U_0, T_r), xlabel="r (μm)", ylabel="Er", label="ρ", ls=:solid, fill=0, legend=:best)
hline!([0], color=:grey, ls=:dash, label=:none)
annotate!(150, 1.05, string(round(U_0))*"Er", 9, xlims=(0, rmax*1e6))


