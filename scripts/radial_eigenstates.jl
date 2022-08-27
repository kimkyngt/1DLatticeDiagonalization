using DrWatson
include(srcdir("tools.jl"))

# Parameters
U_0 = 5
T_r_nz0 = get_Tr(U_0*1.2)
T_r_nz1 = T_r_nz0 /1.4

w_0 = 260e-6 # cavity waist
rmax = w_0
rmin = -rmax
νr = uconvert(u"Hz", sqrt(U_0*Er/m87Sr/pi^2/(w_0*1u"μm")^2))

rr = range(rmin, rmax, length=400)
fig = plot(rr*1e6, -U_0*exp.(-2*rr.^2/w_0^2), label="-U(r)/Eᵣ", minorticks=5, size=(600, 400))
plot!(rr*1e6, get_Unz(rr, U_0, w_0, 0), label="nz=0")
plot!(rr*1e6, get_Unz(rr, U_0, w_0, 1), label="nz=1")
plot!(rr*1e6, get_Unz(rr, U_0, w_0, 2), label="nz=2")
plot!(rr*1e6, get_Unz(rr, U_0, w_0, 3), label="nz=3")