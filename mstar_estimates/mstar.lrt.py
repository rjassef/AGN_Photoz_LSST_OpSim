#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.cosmology import Planck13 as cosmo
from astropy.table import Table

import Shen20
import lrt
import SED_Model

#Start by getting L1450 for an Lstar quasar.
qlf = Shen20.QLF(model="B")

#Read the photometry zero point.
jyzero = np.loadtxt("bandmag.dat",usecols=[2])
filters = np.genfromtxt("bandmag.dat",usecols=[0],dtype='U')

#Get m* for a range in z and all LSST bands.
zmin = 0.01
zmax = 6.0
dlogz = 0.01
zs = 10.**(np.arange(np.log10(1+zmin), np.log10(1+zmax), dlogz)) - 1
mstar = np.zeros((len(zs),len(filters)))
Lnu1450_unnorm = None
for k,z in enumerate(zs):

    #Create the qso object.
    qso = SED_Model.lrt_model()
    qso.zspec = z
    qso.comp = np.zeros(4)
    qso.comp[0] = 1.0
    qso.ebv = 0
    qso.igm = 1.0

    #Get the unnormalized observed magnitudes.
    qso.get_model_fluxes()
    mag_unnorm = -2.5*np.log10(qso.jymod/jyzero)

    #Get the qso Lnu_1450
    if Lnu1450_unnorm is None:
        Lnu1450_unnorm = lrt.get_lnu(qso.comp,qso.zspec,0.145)*u.erg/u.s/u.Hz

    #Get L1450 for an L* quasar.
    Lbol = 10.**(qlf.log_Lstar(z))*qlf.Lstar_units
    nu_1450 = c/(1450.*u.AA)
    Lnu_1450 = (qlf.L1450(Lbol)/nu_1450).to(u.erg/u.s/u.Hz)

    #Normalize the magnitudes
    norm = -2.5*np.log10(Lnu_1450/Lnu1450_unnorm)
    mstar[k] = mag_unnorm + norm
    cond = lrt.cal1.lbar[:len(filters)]/(1.+qso.zspec) < 0.09
    mstar[k, cond] = 99.


output = np.zeros((len(zs),len(filters)+1))
output[:,0] = zs
output[:,1:]  = mstar
np.savetxt("mstar_z.lrt.dat",output,fmt='%15.5f')

import matplotlib.pyplot as plt
for j,filter in enumerate(filters):
    cond = mstar[:,j]<99.
    plt.plot(zs[cond],mstar[cond,j],label=filter)
plt.legend()
plt.xscale('log')
plt.ylim([13.,25.])
plt.xlabel('Redshift')
plt.ylabel('Observed magnitude of L* quasar (AB)')
plt.title('LRT AGN template, Shen et al. QLF')
plt.savefig('mstar_z.lrt.png')
