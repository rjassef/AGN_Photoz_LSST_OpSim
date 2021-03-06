#!/usr/bin/env python

import numpy as np
from scipy.integrate import quad
from scipy import interpolate
from astropy.cosmology import Planck13 as cosmo
import astropy.units as u
from astropy.table import Table

import Shen20

#Integrand. For a given Lfrac=L/L* and redshift z, the function returns the differential number of quasars per unit comoving volume per unit Lfrac. This is just basically a shorthand to strip the units from the qlf.dndLfrac method, as scipy.integrate.quad does not accept units.
def dn_dLfrac(Lfrac,z,qlf):
    return qlf.dndLfrac(Lfrac,z).to(u.Mpc**-3).value

#Integrand. For a given redshift, and functions that determines the minimum and maximum Lfrac=L/L* observable at a given redshift, this function returns the differential number of quasars per unit redshift per sterradian.
def dN_dz(z,Lfrac_min_func,Lfrac_max_func,qlf):
    dVc_dz = cosmo.differential_comoving_volume(z).to(u.Mpc**3/u.sr).value
    dN_dVc = quad(dn_dLfrac,Lfrac_min_func(z),Lfrac_max_func(z),args=(z,qlf))[0]
    return dN_dVc*dVc_dz

###

"""
This is the main function of this script. For a given LSST band magnitude upper and lower limits m_bright, m_faint and a redshift range, it returns the number of quasars in that magnitude range in that redshift range.

The qlf method needs to be able to return dn/dLfrac, i.e., the differential number density of quasars per unit comoving volume at a certain redshift and luminosity fraction Lfrac=L/L*(z).

The number of quasars is then:

N_QSO = area * integral_{zmin}^{zmax} dN/dz * dz

For simplicity, we split the term dN/dz = dN/dVc * dVc/dz, where

dVc/dz = differential comoving volume element

and

dN/dVc = integral_{Lfrac_min(z)}^{Lfrac_max(z)} dn/dLfrac dLfrac

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

"""
def Nqso(zmin,zmax,m_bright,m_faint,LSSTfilter,qlf,area=4.*np.pi):

    #In order to calculate dN/dVc as defined above, we need to create functions that returns Lfrac_min and Lfrac_max for an arbitrary redshift given the magnitude limits m_bright and m_faint. We do this by using the scipy.interpolate.interp1d module. In short, given a list of x and y=f(x) points, this module creates a function that returns f(xnew) for any value of xnew (within the boundaries of x) by interpolating between the original points. The original list of x and y data points are read from the file created by the mstar.vandenber.py script.

    #Star by reading the table of observed magnitudes as a function of redshift, created by the script mstar.vandenberk.py. See mstar_z.vandenberk.png for illustration
    mstar_data = Table.read("mstar_z.vandenberk.dat", format='ascii')

    #Only use the redshift range on which the template overlaps with the filter curves. NOTE: Outside this redshift range, the function will trigger an exception and return None.
    cond = (~np.isnan(mstar_data[LSSTfilter]))
    z_mstar = mstar_data['Redshift'][cond]
    mstar   = mstar_data[LSSTfilter][cond]

    #Create the interpolating function that returns the value of Lfrac_min=Lmin/L* for an apparent magnitude limit m_faint at an arbitrary redshift z.
    Lfrac_min_func = interpolate.interp1d(z_mstar, 10.**(-0.4*(m_faint-mstar)))

    #Now create the interpolating function that returns the value of Lfrac_max=Lmax/L* for an apparent magnitude limit m_bright at an arbitrary redshift z. A special case here is that m_bright can be specified as -np.inf to indicate no upper limit to the integration. In that case, we simply create a lambda function that for any value of redshift returns L/L* = infinite (np.inf).
    if np.isinf(m_bright):
        Lfrac_max_func = lambda z : np.inf
    else:
        Lfrac_max_func = interpolate.interp1d(z_mstar, 10.**(-0.4*(m_bright-mstar)))

    #Exit if the redshift range requested exceeds the redshift range allowed by the template.
    if zmin < np.min(z_mstar) or zmax > np.max(z_mstar):
        print("Redshift range outside of template boundaries.")
        print("Template only allows {0}<z<{1}".format(np.min(z_mstar), np.max(z_mstar)))
        return None

    #Carry out the calculation.
    return area * quad(dN_dz, zmin, zmax, args=(Lfrac_min_func, Lfrac_max_func, qlf))[0]

###

#Create the QLF object.
qlf = Shen20.QLF(model="B")

#Find the number of quasars in between 18 and 19 mag in the g-band between z=0.5 and z=0.6
print("Number of quasars with 18<g<19 and 0.5<redshift<0.6:")
print("{0:,.0f}".format(Nqso(0.50, 0.60, 18.0, 19.0, 'LSSTg', qlf)))
print()

#Find the number of quasars brighter than i=25 in the redshift range 1 to 2 over the entire sky.
print("Number of quasars with i<25 and 1<redshift<2:")
print("{0:,.0f}".format(Nqso(1.0, 2.0, -np.inf, 25.0, 'LSSTi', qlf)))
print()
