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
U_0 = 30
T_r = get_Tr(U_0*1.1)
w_0 = 260e-6 # cavity waist
rmax = 300e-6 # in 

normalization = quadgk(r -> r*get_rho(r, U_0, T_r), 0, rmax)[1]
powers = [0, 2, 4, 6, 8, 10, ]
@time moments = [
    quadgk(r -> r^powers[ii]*r*get_rho(r, U_0, T_r), 0, rmax)[1] for ii in range(1, length(powers))
    ] / normalization
expansions = moments .* [(2/w_0)^powers[k] / factorial(k) for k in range(1, length(powers))]


fig = scatter(powers/2,  expansions, legend=false, xticks=powers/2, xlabel=L"k", ylabel=L"\left\langle \frac{1}{k!}\left(\frac{2r^2}{w_0^2}\right)^k \right\rangle"*", Thermal averaged",
title=string(U_0)*" Er, "*string(round(u"nK", T_r)),
yscale=:log10,
# yticks=[1, 0.1, 0.01, 0.001, 0.0001],
yminorticks=10,
grid=true,
)

d=(U_0=U_0, T_r=T_r, rmax=rmax)
Plots.pdf(fig, plotsdir("moments", savename("Moments", d, "pdf")))

fig