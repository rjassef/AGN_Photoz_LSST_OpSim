import numpy as np
from scipy.integrate import quad, dblquad
from scipy import interpolate
import astropy.units as u
from astropy.constants import L_sun, c
from astropy.table import Table
#from scipy.interpolate import interp1d

import os
import re
root_path = re.search("(.*/AGN_Photoz_LSST_OpSim)/*",os.getcwd()).group(1)

import sys
sys.path.append(root_path+"/Quasar_Counts/")
from Phi_Obs import get_phi_lam_obs

sys.path.append(root_path+"/../Shen_QLF_tools/pubtools/")
from utilities import return_qlf_in_band

#LSST Filter wavelengths.
lam_eff = {'LSSTu': 3751.36*u.AA,
           'LSSTg': 4741.64*u.AA,
           'LSSTr': 6173.23*u.AA,
           'LSSTi': 7501.62*u.AA,
           'LSSTz': 8679.19*u.AA,
           'LSSTy': 9711.53*u.AA,
          }

#Integrand. For a given redshift, and functions that determines the minimum and maximum Lfrac=L/L* observable at a given redshift, this function returns the differential number of quasars per unit redshift per sterradian.
def dN_dz(z, m_bright, m_faint_use_func, LSSTfilter, qlf, cosmo):

    #If the faint magnitude is brighter than the bright mangitude limit, it means that no object should be counted because L(m_bright)<L(Mi_lim)
    m_faint = m_faint_use_func(z)
    if m_faint<m_bright:
        return 0

    #Get the luminosity limits.
    DL = cosmo.luminosity_distance(z)
    nu_obs = (c/lam_eff[LSSTfilter])
    fnu_0  = 3631*u.Jy
    f_L_fact = np.log10((4.*np.pi*DL**2 * nu_obs * fnu_0).to(u.erg/u.s).value)
    lL_lam_obs_min = f_L_fact - 0.4*m_faint
    lL_lam_obs_max = f_L_fact - 0.4*m_bright

    #Set the grid in which we will want the QLF.
    nlL_lam_obs    = 100
    dlL_lam_obs    = (lL_lam_obs_max-lL_lam_obs_min)/nlL_lam_obs
    if dlL_lam_obs > 0.1:
        dlL_lam_obs    = 0.1
    lL_lam_obs     = np.arange(lL_lam_obs_min, lL_lam_obs_max + 0.1*dlL_lam_obs, dlL_lam_obs)

    #We estimate numerically the observed QLF starting from the bolometric QLF. Integrate it assuming a constant value per bin.
    nu_rest = (nu_obs*(1.0+z)).to(u.Hz).value
    lL_lam_obs_S20, lphi_lam_obs_S20 = return_qlf_in_band(z, nu_rest, model=qlf.model)
    lphi_lam_obs = interpolate.interp1d(lL_lam_obs_S20, lphi_lam_obs_S20)
    dN_dVc = np.sum(10.**(lphi_lam_obs(lL_lam_obs))*dlL_lam_obs)*u.Mpc**-3

    #Calculate the differential comoving volume term.
    dVc_dz = cosmo.differential_comoving_volume(z)

    #Return the integrand of eqn. (1) in the comments below.
    return (dN_dVc*dVc_dz).to(u.sr**-1).value

###

"""
This is the main function of this script. For a given LSST band magnitude upper and lower limits m_bright, m_faint and a redshift range, it returns the number of quasars in that magnitude range in that redshift range.

The number of quasars is then:

N_QSO = area * integral_{zmin}^{zmax} dN/dz * dz                            (1)

For simplicity, we split the term dN/dz = dN/dVc * dVc/dz, where

dVc/dz = differential comoving volume element                               (2)

and

dN/dVc = integral_{Lfrac_min(z)}^{Lfrac_max(z)} dn/dLfrac dLfrac            (3)

which is now carried out by Shen et al. (2020) pubtools package.

We carry out the outer integration using scipy.integrate.quad. The first arument given to quad is the function we want to integrate, followed by the integration limits. Additionally, we can give a list of arguments using the args keyword that are passed to the function as well. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html for further details.

Parameter
---------
zmin : float
    Lower redshift boundary.

zmax : float
    Upper redshft boundary.

m_bright : float or -np.inf
    Bright mangitude limit.

m_faint : float
    Faint magntiude limit.

LSSTfilter : str
    LSST filter for the magnitude limits. Can be LSSTu, LSSTg, LSSTr, LSSTi, LSSTz or LSSTy

qlf : QLF object
    QLF object. Needs to have the method dndLfrac(Lfrac,z).

area : float, optional
    Sky area covered in sterradians. Default: 4pi.

mstar_data: astropy.table object, optional
    Astropy Table with the magnitudes of an Lstar quasar as a function of redshift, in the same format as that produced by the mstar.vandenberk.py script. If not provided, the mstar_z.vandenberk.dat table will be read.

Mi_lim: float, optional
    Faint absolute i-band magnitude to which to integrate the QLF.

cosmo: astropy.cosmology object, optional
    Astropy cosmology object. If not provided, astropy.cosmology.Plnack13 is used.

"""
def Nqso(zmin, zmax, m_bright, m_faint, LSSTfilter, qlf, area=4.*np.pi*u.sr, mstar_data=None, Mi_lim=np.inf, cosmo=None):

    #Set the cosmology to use if not provided.
    if cosmo is None:
        from astropy.cosmology import Planck13 as cosmo

    #If not provided as input, start by reading the table of observed magnitudes as a function of redshift, created by the script mstar.vandenberk.py. See mstar_z.vandenberk.png for illustration.
    if mstar_data is None:
        mstar_data = Table.read("../mstar_z.vandenberk.dat", format='ascii')

    #Only use the redshift range on which the template overlaps with the filter curves. NOTE: Outside this redshift range, the function will trigger an exception and return None.
    cond = (~np.isnan(mstar_data[LSSTfilter]))
    z_mstar = mstar_data['Redshift'][cond]
    mstar   = mstar_data[LSSTfilter][cond]

    #We want to estimate the number of quasars down to a faint limit of m_faint or Mi_lim, whichever is brightest.
    Mstar_i          = mstar_data['M_i'][cond]

    #Functions to calculate m_faint_use = min ( m_faint, m_i(Mi_lim) ) at an arbitrary redshift.
    Mstar_i_interp = interpolate.interp1d(z_mstar, Mstar_i)
    mstar_interp   = interpolate.interp1d(z_mstar, mstar)
    m_faint_use_func = lambda z: np.where(m_faint>mstar_interp(z)+(Mi_lim-Mstar_i_interp(z)), mstar_interp(z)+(Mi_lim-Mstar_i_interp(z)), m_faint)

    #Carry out the calculation.
    return area.to(u.sr).value * quad(dN_dz, zmin, zmax, args=(m_bright, m_faint_use_func, LSSTfilter, qlf, cosmo))[0] #, epsabs=1e-1, epsrel=1e-3)[0]
