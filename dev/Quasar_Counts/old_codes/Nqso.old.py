import numpy as np
from scipy.integrate import quad, dblquad
from scipy import interpolate
import astropy.units as u
from astropy.constants import L_sun
from astropy.table import Table

#LSST Filter wavelengths. Will be useful for the obscuation function we will need to set up.
lam_eff = {'LSSTu': 3751.36*u.AA,
           'LSSTg': 4741.64*u.AA,
           'LSSTr': 6173.23*u.AA,
           'LSSTi': 7501.62*u.AA,
           'LSSTz': 8679.19*u.AA,
           'LSSTy': 9711.53*u.AA,
          }

#Integrand. See comments.
def dp_dlLfrac_dlogNH(lLfrac, log_NH, gamma1, gamma2, z, qlf, lLx_44_func):
    if lLfrac > 10.0:
        return 0
    Lfrac = 10.**(lLfrac)
    lLx_44 = lLx_44_func(Lfrac)
    return qlf.fNH(log_NH, lLx_44, z)*(Lfrac**gamma1 + Lfrac**gamma2)**-1

#Integrand. For a given redshift, and functions that determines the minimum and maximum Lfrac=L/L* observable at a given redshift, this function returns the differential number of quasars per unit redshift per sterradian.
def dN_dz(z,lLfrac_min_func_0,lLfrac_max_func_0,LSSTfilter,qlf,cosmo):

    #Calculate phi_star_prime, make sure it is in the right units, and strip the units.
    phi_star = 10.**(qlf.log_phi_star(z)) * qlf.phi_star_units
    phi_star = phi_star.to(u.dex**-1 * u.Mpc**-3).value

    #These are the limits in log L/L* for NH=0.
    lLfrac_min = lLfrac_min_func_0(z)
    lLfrac_max = lLfrac_max_func_0(z)
    if lLfrac_min>lLfrac_max:
        return 0

    #The key here is that we need to transform those limits in log L/L* to any absorption. This factor is the one that transforms the luminosities.
    obscuration_fact = 0.4*qlf.dgr.to(u.cm**2).value*1e22 * qlf.xi(lam_eff[LSSTfilter]/(1.+z))

    lLfrac_min_func = lambda log_NH: lLfrac_min + 10.**(log_NH-22)*obscuration_fact
    lLfrac_max_func = lambda log_NH: lLfrac_max + 10.**(log_NH-22)*obscuration_fact

    c1x, k1x, c2x, k2x = 10.83, 0.28, 6.08, -0.020
    Lstar_10 = 10.**(qlf.log_Lstar(z)-10.)
    conv = 44.-10.-np.log10(L_sun.to(u.erg/u.s).value)
    lLx_44_func = lambda Lfrac: conv-np.log10( c1x*(Lfrac*Lstar_10)**(k1x-1) + c2x*(Lfrac*Lstar_10)**(k2x-1) )

    #dN_dVc = dblquad(dp_dLfrac_dlogNH, 20., 25., Lfrac_min_func, Lfrac_max_func, args=(alpha,beta,z,qlf))[0]
    dN_dVc = dblquad(dp_dlLfrac_dlogNH, 20., 25., lLfrac_min_func, lLfrac_max_func, args=(alpha,beta,z,qlf, lLx_44_func))[0]
    dN_dVc *= phi_star_prime
    #print(dN_dVc)
    #input()

    #Calculate the differential comoving volume term.
    dVc_dz = cosmo.differential_comoving_volume(z).to(u.Mpc**3/u.sr).value

    #Return the integrand of eqn. (1) in the comments below.
    return dN_dVc*dVc_dz

###

"""
This is the main function of this script. For a given LSST band magnitude upper and lower limits m_bright, m_faint and a redshift range, it returns the number of quasars in that magnitude range in that redshift range.

The qlf method needs to be able to return dn/dLfrac, i.e., the differential number density of quasars per unit comoving volume at a certain redshift and luminosity fraction Lfrac=L/L*(z).

The number of quasars is then:

N_QSO = area * integral_{zmin}^{zmax} dN/dz * dz                            (1)

For simplicity, we split the term dN/dz = dN/dVc * dVc/dz, where

dVc/dz = differential comoving volume element                               (2)

and

dN/dVc = integral_{Lfrac_min(z)}^{Lfrac_max(z)} dn/dLfrac dLfrac            (3)

Now, we could call qlf.dndLfrac, as done in number_of_quasars_fixed_mlim.py, but that it is very inneficient. Note that:

dn/dLfrac = phi_star_prime(z) / (Lfrac**-alpha(z) + Lfrac**-beta(z))        (4)

so we can make the integration a lot more efficient if we do not calculate each parameter that depends on z only each time we call dn/dLfrac for a fixed redshift. We define then

dp/dLfrac = dn/dLfrac / phi_star_prime = (Lfrac**-alpha + Lfrac**-beta)**-1 (5)

So that

dN/dVc = phi_star_prime * integral_{Lfrac_min(z)}^{Lfrac_max(z)} dp/dLfrac dLfrac                                                                      (6)

We carry out all the integrations using scipy.integrate.quad. The first arument given to quad is the function we want to integrate, followed by the integration limits. Additionally, we can give a list of arguments using the args keyword that are passed to the function as well. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html for further details.

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

    #In order to calculate dN/dVc as defined above, we need to create functions that returns Lfrac_min and Lfrac_max for an arbitrary redshift given the magnitude limits m_bright and m_faint. We do this by using the scipy.interpolate.interp1d module. In short, given a list of x and y=f(x) points, this module creates a function that returns f(xnew) for any value of xnew (within the boundaries of x) by interpolating between the original points. The original list of x and y data points are read from the file created by the mstar.vandenber.py script.

    #Set the cosmology to use if not provided.
    if cosmo is None:
        from astropy.cosmology import Planck13 as cosmo

    #If not provided as input, start by reading the table of observed magnitudes as a function of redshift, created by the script mstar.vandenberk.py. See mstar_z.vandenberk.png for illustration.
    if mstar_data is None:
        mstar_data = Table.read("mstar_z.vandenberk.dat", format='ascii')

    #Only use the redshift range on which the template overlaps with the filter curves. NOTE: Outside this redshift range, the function will trigger an exception and return None.
    cond = (~np.isnan(mstar_data[LSSTfilter]))
    z_mstar = mstar_data['Redshift'][cond]
    mstar   = mstar_data[LSSTfilter][cond]

    #We want to estimate the number of quasars down to a faint limit of m_faint or Mi_lim, whichever is brightest.
    Mstar_i          = mstar_data['M_i'][cond]


    lLfrac_min_mfaint = -0.4*(m_faint-mstar)
    lLfrac_min_Mi     = -0.4*(Mi_lim-Mstar_i)
    lLfrac_min_table  = np.where(lLfrac_min_Mi>lLfrac_min_mfaint, lLfrac_min_Mi, lLfrac_min_mfaint)

    #obscuration_func = lambda log_NH, z: 10.**(log_NH+np.log10(qlf.dgr.to(u.cm**2).value)) * qlf.xi(lam_eff[LSSTfilter]/(1.+z))

    #Create the interpolating function that returns the value of Lfrac_min=Lmin/L* for an apparent magnitude limit m_faint at an arbitrary redshift z.
    #Lfrac_min_func_0 = interpolate.interp1d(z_mstar, Lfrac_min_table)
    lLfrac_min_func_0 = interpolate.interp1d(z_mstar, lLfrac_min_table)

    #Now create the interpolating function that returns the value of Lfrac_max=Lmax/L* for an apparent magnitude limit m_bright at an arbitrary redshift z. A special case here is that m_bright can be specified as -np.inf to indicate no upper limit to the integration. In that case, we simply create a lambda function that for any value of redshift returns L/L* = infinite (np.inf).
    if np.isinf(m_bright):
        #Lfrac_max_func = lambda z : np.inf
        lLfrac_max_func = lambda z : np.inf
    else:
        #Lfrac_max_func_0 = interpolate.interp1d(z_mstar, 10.**(-0.4*(m_bright-mstar)))
        lLfrac_max_func_0 = interpolate.interp1d(z_mstar, -0.4*(m_bright-mstar))

    #Exit if the redshift range requested exceeds the redshift range allowed by the template.
    if zmin < np.min(z_mstar) or zmax > np.max(z_mstar):
        print("Redshift range outside of template boundaries.")
        print("Template only allows {0}<z<{1}".format(np.min(z_mstar), np.max(z_mstar)))
        return np.nan

    #Carry out the calculation.
    #return area.to(u.sr).value * quad(dN_dz, zmin, zmax, args=(Lfrac_min_func_0, Lfrac_max_func_0, obscuration_func, qlf, cosmo), epsabs=1e-3, epsrel=1e-3)[0]
    return area.to(u.sr).value * quad(dN_dz, zmin, zmax, args=(lLfrac_min_func_0, lLfrac_max_func_0, LSSTfilter, qlf, cosmo), epsabs=1e-3, epsrel=1e-3)[0]
