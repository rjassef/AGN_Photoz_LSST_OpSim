#!/usr/bin/env python

import numpy as np
from scipy.integrate import quad
from scipy import interpolate
from astropy.cosmology import Planck13 as cosmo
import astropy.units as u
from astropy.table import Table

import Shen20

#Integrand. For a given Lfrac=L/L* and redshift z, the function returns the differential number of quasars per unit comoving volume per unit Lfrac.
def dN_dVc_dLfrac(Lfrac,z,qlf):
    return qlf.dndLfrac(Lfrac,z).to(u.Mpc**-3).value

#Integrand. For a given redshift, and a function that determines the minimum Lfrac=L/L* observable, this function returns the differential number of quasars per unit redshift per sterradian.
def dN_dz(z,Lfrac_min,Lfrac_max,qlf):
    dVc_dz = cosmo.differential_comoving_volume(z).to(u.Mpc**3/u.sr).value
    dN_dVc = quad(dN_dVc_dLfrac,Lfrac_min(z),Lfrac_max(z),args=(z,qlf))[0]
    return dN_dVc*dVc_dz

###

"""
This is the main function of this script. For a given LSST band magnitude upper and lower limits m_bright, m_faint and a redshift range, it returns the number of quasars in that magnitude range in that redshift range.

The qlf method needs to be able to return dn/dLfrac, i.e., the differential number density of quasars per unit comoving volume at a certain redshift and luminosity fraction Lfrac=L/L*(z).

The number of quasars is then:

N_QSO = area * integral_{zmin}^{zmax} dN/dVc * dVc/dz * dz

and

dN/dVc = integral_{Lfrac_min}^{Lfrac_max} dn/dLfrac dLfrac


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

    #Read the table of observed magnitudes as a function of redshift, created by the script mstar.vandenberk.py.
    mstar_data = Table.read("mstar_z.vandenberk.dat", format='ascii')

    #Only use redshift ranges on which the template overlaps with the filter curves.
    cond = mstar_data[LSSTfilter]<99.
    z_mstar = mstar_data['Redshift'][cond]
    mstar   = mstar_data[LSSTfilter][cond]

    #Create the interpolating function that returns the value of Lfrac=L/L* for an apparent magnitude limit m_lim at redshift z.
    Lfrac_min = interpolate.interp1d(z_mstar, 10.**(-0.4*(m_faint-mstar)))
    if np.isinf(m_bright):
        Lfrac_max = lambda z : np.inf
    else:
        Lfrac_max = interpolate.interp1d(z_mstar, 10.**(-0.4*(m_bright-mstar)))

    #Exit if the redshift range requested exceeds the redshift range allowed by the template.
    if zmin < np.min(z_mstar) or zmax > np.max(z_mstar):
        print("Redshift range outside of template boundaries.")
        print("Template only allows {0}<z<{1}".format(np.min(z_mstar), np.max(z_mstar)))
        return None

    #Carry out the calculation.
    return area * quad(dN_dz, zmin, zmax, args=(Lfrac_min, Lfrac_max, qlf))[0]

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
