import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import *

# Simulation parameters
atoms_per_site = 13
Nsim = 500
chi = -1e-4 
f0 = 0.5
readout_noise_amplitude = 1
photon_noise_factor = np.sqrt(1/2.8) #sqrt(atom/photon)

# Average cloud shape
zz = np.arange(0, 512)
zz = zz - 200

Natom_peak = atoms_per_site*6.04/(813e-3 / 2) # peak atom number for a pixel at the center of the cloud
sigma_cloud = 0.5e-3/6.04e-6 
density_cloud = Natom_peak*np.exp(-(zz**2)/(2*sigma_cloud**2)) # Averaged cloud density. atoms per pixel. 


def linearfit(x, a, b):
    return a*x + b

def liear_odr(B, x):
    return B[0]*x + B[1]

def div0( a, b ):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide( a, b )
        c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c

def readout_noise(amplitude):
    return np.random.randn(512)*amplitude

def photon_noise(N, factor):
    """factor: sqrt(electron or atom/photon)"""
    return np.array(np.random.randn(512)*np.sqrt(N)*factor)*1.0

# sampling
def run_sim(Nsim, chi, f0, density_cloud):
    f0_ptbypt = np.zeros(Nsim)
    f0err_ptbypt = np.zeros(Nsim)
    Ne_mean = np.zeros(len(density_cloud))
    Ng_mean = np.zeros(len(density_cloud))
    for ii in range(Nsim):

        # load atoms
        Ntot_true = np.array([np.random.poisson(density_cloud[ii]) for ii in range(len(density_cloud))])

        # excite the atoms
        fraction = f0 + chi*Ntot_true
        Ne = np.array([np.random.binomial(Ntot_true[ii], fraction[ii]) for ii in range(len(Ntot_true))])
        Ng = Ntot_true - Ne

        # add detection noise
        Ne = Ne + photon_noise(Ne, photon_noise_factor)
        Ne = Ne + readout_noise(readout_noise_amplitude)
        Ng = Ng + photon_noise(Ng, photon_noise_factor)
        Ng = Ng + readout_noise(readout_noise_amplitude)
        Nbg = readout_noise(readout_noise_amplitude)
        Ntot = Ng + Ne - 2*Nbg
        frac_measured = div0(Ne - Nbg, Ntot)

        # Density shift correction
        # popt, pcov = curve_fit(linearfit, Ntot[Ntot>0], frac_measured[Ntot>0], sigma=1/np.sqrt(Ntot[Ntot>0]), p0=[0, 0.5])
        # perr = np.sqrt(np.diag(pcov))
        # f0_ptbypt[ii] = popt[1]
        # f0err_ptbypt[ii] = perr[1]

        goodindx = (Ntot>4*6.04/(813e-3 / 2)) 
        linear = Model(liear_odr)
        mydata = RealData(
             Ntot[goodindx], frac_measured[goodindx], 
            sx=1/np.sqrt(Ntot[goodindx]), # SEM from shot noise
            sy=np.sqrt(frac_measured[goodindx]*(1-frac_measured[goodindx]))/np.sqrt(Ntot[goodindx]), # SEM from QPN
            )
        myodr = ODR(mydata, linear, beta0=[0, 0.5])
        myoutput = myodr.run()
        f0_ptbypt[ii] = myoutput.beta[1]
        f0err_ptbypt[ii] = myoutput.sd_beta[1]

        Ne_mean = Ne_mean + Ne - Nbg
        Ng_mean = Ng_mean + Ng - Nbg

    Ne_mean = Ne_mean/Nsim
    Ng_mean = Ng_mean/Nsim
    # Processing with averaged density profile.
    Nmean = Ne_mean + Ng_mean
    fmean = div0(Ne_mean, Nmean)
    popt_mean = curve_fit(linearfit, Nmean[Nmean>0], fmean[Nmean>0], sigma=1/np.sqrt(Nmean[Nmean>0]))[0]
    f0_mean = popt_mean[1]

    return f0_ptbypt, f0err_ptbypt, f0_mean, Nmean, fmean, Ntot, frac_measured, fraction

f0_ptbypt, f0err_ptbypt, f0_mean, Nmean, fmean, Ntot, frac_measured, fraction= run_sim(Nsim, chi, f0, density_cloud)

# For testing single shot
goodindx = (Ntot>0)*(frac_measured>=0)
popt, pcov = curve_fit(linearfit, Ntot[goodindx], frac_measured[goodindx], sigma=1/np.sqrt(Ntot[goodindx]))
perr = np.sqrt(np.diag(pcov))

# Plot
fig, axes = plt.subplots(2, 2, figsize=(8, 6))

axes[0, 0].plot(Ntot/(6.04/(813e-3 / 2)), label="$N(x)$")
axes[0, 0].plot(density_cloud/(6.04/(813e-3 / 2)), label="$\\bar{N}(x)$")
axes[0, 0].legend()
axes[0, 0].set(xlabel="pixel", ylabel="Atoms per site")

axes[0, 1].plot(frac_measured, label="measured fraction")
axes[0, 1].plot(fraction, label="fraction")
axes[0, 1].legend()
axes[0, 1].set(xlabel="pixel", ylabel=r"$\rho^{ee}$", ylim=[0, 1])

axes[1, 0].plot(Ntot[goodindx], frac_measured[goodindx], '.', label="")
axes[1, 0].plot(Ntot, linearfit(Ntot, *popt), color="C1", label="Linear fit")
axes[1, 0].fill_between(Ntot, y1=linearfit(Ntot, *(popt+perr)), y2=linearfit(Ntot, *(popt-perr)), alpha=0.7, color="C1")
axes[1, 0].axhline(f0, color="red", label="true f0")
axes[1, 0].legend()
axes[1, 0].set(ylabel=r"$\rho^{ee}$", xlabel="Atoms per pixel", ylim=[0, 1])

axes[1, 1].hist(f0_ptbypt, color="black", histtype="step", label="pt-by-pt")
axes[1, 1].axvline(f0, color="red", label="true")
axes[1, 1].axvline(np.average(f0_ptbypt), color="black", lw=2, ls="--", label="pt-by-pt $\mu$")
axes[1, 1].axvline(f0_mean, color="blue", ls="--", label="mean")
axes[1, 1].legend()
axes[1, 1].set(xlabel="Excitation fraction")
fig.tight_layout()
plt.show()