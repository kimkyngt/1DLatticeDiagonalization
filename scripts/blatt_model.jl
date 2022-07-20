using DrWatson
include(srcdir("initialize.jl"))

function get_qnx(nx, Tr=2e-6u"K", nu=450u"s^-1")
    zx = exp(-uconvert(NoUnits, h*nu/kB/Tr))
    return (1-zx)*zx^nx
end

function get_pe(Ω, t, Tr, nu=450u"s^-1")
    zx = exp(-uconvert(NoUnits, h*nu/kB/Tr))
    ϕ = 2*π*t*Ω*exp(-ηx^2/2)*exp(-ηz^2/2)
    return 1/2 + (1-zx)/2*(zx*cos(ϕ*(1-ηx^2)) - cos(ϕ))/(1 + zx^2 - 2*zx*cos(ϕ*ηx^2))
end


# Reproducing Figure 4 of the paper. 
        ## Parameters
ts = range(0, 0.08, length=10000)*1u"s"
ηz = 0.24
ηx = 10e-3*3.2

    ##Plot
plot(ts*1e3, [get_pe(59u"s^-1", t, 1e-6*1u"K") for t in ts],
    label="59 Hz, 1 uK"
)
plot!(ts*1e3, [get_pe(76u"s^-1", t, 3e-6*1u"K") for t in ts],
    title="FIG. 4", 
    label="76 Hz, 3 uK",
    legend=:true,
)
