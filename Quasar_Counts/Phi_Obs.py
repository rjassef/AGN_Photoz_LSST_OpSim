import numpy as np
import astropy.units as u
from astropy.constants import L_sun
from scipy.interpolate import interp1d

def get_Lfrac_lam(Lfrac, Lstar_10, qlf):
    """
    This function returns L_lam(L)/L_lam(Lstar). This function is only valid for UV/optical wavelengths, were we assume the conversion factors are just proportional to the B-band conversion.

    Parameters
    ----------

    Lfrac: numpy array
        Values of L/Lstar for which to calculate Lfrac_lam = L_lam/L_lam(Lstar)

    Lstar_10: float
        Value of Lstar in units of 10^10 Lsun.

    qlf: QLF object
        QLF being used.

    """
    D = np.tile(qlf.c_B*Lstar_10**qlf.k_B, [len(Lfrac),1])
    Lfrac_2D = np.tile(Lfrac, [len(qlf.c_B),1]).T
    return np.sum(D,axis=1)/np.sum(D*Lfrac_2D**(qlf.k_B-1),axis=1)

def jacobian(Lfrac, Lstar_10, qlf):
    """
    This function returns the jacobian dlog L / dlog L_lam.

    Parameters
    ----------

    Lfrac: numpy array
        Values of L/Lstar for which to calculate Lfrac_lam = L_lam/L_lam(Lstar)

    Lstar_10: float
        Value of Lstar in units of 10^10 Lsun.

    qlf: QLF object
        QLF being used.

    """
    D = np.tile(qlf.c_B*Lstar_10**qlf.k_B, [len(Lfrac),1])
    Lfrac_2D = np.tile(Lfrac, [len(qlf.c_B),1]).T
    return np.sum(-D*Lfrac_2D**qlf.k_B,axis=1) / np.sum(D*(qlf.k_B -1)*Lfrac_2D**qlf.k_B,axis=1)
    #return np.sum(D*(1.+qlf.k_B)*Lfrac_2D**qlf.k_B, axis=1)/np.sum(D*Lfrac_2D**qlf.k_B, axis=1)



def get_phi_lam_obs(z, qlf, lLfrac_lam_obs_min, lLfrac_lam_obs_max, lam_eff_filter):

    """
    This is the main function of this module. For a given redshift, it returns the observed qlf at a given observed wavelength equal to the effective wavelength of the filter used / (1+z).

    Parameters
    ----------

    z: float
        Redshift.

    qlf: QLF object.
        QLF model being used.

    lLfrac_lam_obs_min: float
        Lower limit of log L_lam_obs / L_lam(Lstar) to produce the observed QLF.

    lLfrac_lam_obs_max: float
        Upper limit of log L_lam_obs / L_lam(Lstar) to produce the observed QLF.

    lam_eff_filter: float
        Effective wavelength of the filter being used.

    """

    #Start by getting the value of Lstar in units of 10^10 Lsun, which will be useful later on.
    Lstar = 10.**(qlf.log_Lstar(z))*qlf.Lstar_units
    Lstar_10 = (Lstar/(1e10*L_sun)).to(1.).value

    #Set the grid in bolometric L/Lstar.
    lLfrac_min = -3.0
    lLfrac_max =  3.0#10.0
    dlLfrac    =  0.01
    lLfrac     = np.arange(lLfrac_min,lLfrac_max,dlLfrac)
    Lfrac      = 10.**lLfrac

    #Get the bolometric QLF evaluated in the grid of Lfrac.
    phi_bol = qlf.phi_bol_Lfrac(Lfrac, z)

    #Transform the bolometric QLF to the intrinsic luminosity QLF in the band. We assume that the bolometric correction in all bands of interest is proportional to the one in the B-band, as is done in the Hopkins07 provided code.
    phi_lam = phi_bol*jacobian(Lfrac, Lstar_10, qlf)
    Lfrac_lam    = get_Lfrac_lam(Lfrac, Lstar_10, qlf)
    lLfrac_lam   = np.log10(Lfrac_lam)
    #dlLfrac_lam  = dlLfrac/jacobian(Lfrac, Lstar_10, qlf)

    #Since there is a natural dispersion to the bolometric corrections, we convolve phi_lam with the uncertainty function to take it into account.
    phi_lam_2D        = np.tile(phi_lam, (len(phi_lam), 1))
    sigma             = qlf.get_sigma(Lfrac, Lstar_10)
    lLfrac_lam_sig    = lLfrac_lam
    sigma_2D          = np.tile(sigma, (len(sigma), 1))
    lLfrac_lam_2D     = np.tile(lLfrac_lam, (len(lLfrac_lam), 1))
    lLfrac_lam_sig_2D = np.tile(lLfrac_lam_sig, (len(lLfrac_lam), 1)).T

    p = (2.*np.pi)**-0.5 * sigma_2D**-1 * np.exp( -0.5*( (lLfrac_lam_sig_2D - lLfrac_lam_2D)/sigma_2D)**2)

    phi_lam_sig = np.sum(phi_lam_2D*p * dlLfrac, axis=1)

    #The next step is to convolve with the obscuration function. The issue here is that the observed luminosity in the band is a function of the intrinsic luminosity and the obscuration.
    lNH_min = 20.
    lNH_max = 26.
    dlNH    = 0.01
    lNH     = np.arange(lNH_min, lNH_max, dlNH)

    #This is the output grid.
    nlLfrac_lam_obs    = 100
    dlLfrac_lam_obs    = (lLfrac_lam_obs_max-lLfrac_lam_obs_min)/nlLfrac_lam_obs
    if dlLfrac_lam_obs > 0.1:
        dlLfrac_lam_obs    = 0.1
    lLfrac_lam_obs     = np.arange(lLfrac_lam_obs_min, lLfrac_lam_obs_max + 0.1*dlLfrac_lam_obs, dlLfrac_lam_obs)

    #Now, get the obscuration function in the observed band.
    ltheta_fact = 0.4*qlf.dgr(z).to(u.cm**2).value*1e22 * qlf.xi(lam_eff_filter/(1.+z))
    ltheta = 10.**(lNH-22) * ltheta_fact
    ltheta_2D = np.tile(ltheta, [len(lLfrac_lam_obs), 1])

    #For each NH, we will need to evaluate the unreddened QLF at different luminosities. So let's go ahead and make the array of the lLfracs_lam_sig at which we will need to evaluate.
    lLfrac_lam_sig_eval_2D = np.tile(lLfrac_lam_obs, [len(lNH), 1]).T + ltheta_2D

    #This function is the probability of a certain NH for a given bolometric luminosity. We will need to recalculate Lfrac to make sure it shows the correct bolometric luminosity at the evaluation point.
    f_NH = qlf.fNH(lNH, Lfrac, Lstar_10, z).T
    #print(np.sum(f_NH*dlNH, axis=1))

    #Now, for every NH, we will need to interpolate/extrapolate the QLF to evaluate it in the needed places.
    phi_lam_sig_f_NH_eval_2D = np.zeros((len(lLfrac_lam_obs),len(lNH)))
    for i in range(len(lNH)):
        #phi_lam_sig_f_NH_interp = interp1d(lLfrac_lam_sig, phi_lam_sig*f_NH[:,i], bounds_error=False, fill_value = 0.)
        log_phi_lam_sig_f_NH_interp = interp1d(lLfrac_lam_sig[5:-5], np.log10(phi_lam_sig[5:-5].value*f_NH[5:-5,i]+1e-32), kind='linear', fill_value = 'extrapolate')
        #phi_lam_sig_f_NH_eval_2D[:,i] = phi_lam_sig_f_NH_interp(lLfrac_lam_sig_eval_2D[:,i])
        phi_lam_sig_f_NH_eval_2D[:,i] = 10.**(log_phi_lam_sig_f_NH_interp(lLfrac_lam_sig_eval_2D[:,i]))

    phi_lam_obs = np.sum(phi_lam_sig_f_NH_eval_2D * dlNH, axis=1) *  phi_lam_sig.unit

    # import matplotlib.pyplot as plt
    # plt.plot(lLfrac, qlf.phi_bol_Lfrac(10.**lLfrac, z))
    # plt.plot(lLfrac_lam, phi_lam)
    # plt.plot(lLfrac_lam_sig, phi_lam_sig)
    # plt.plot(lLfrac_lam_obs, phi_lam_obs)
    # #plt.xscale('log')
    # plt.yscale('log')
    # plt.show()

    return phi_lam_obs, dlLfrac_lam_obs*u.dex
