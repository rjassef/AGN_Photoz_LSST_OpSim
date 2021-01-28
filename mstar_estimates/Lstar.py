#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.cosmology import Planck13 as cosmo
from astropy.table import Table

import os
os.environ['PYSYN_CDBS'] = "."
import warnings
warnings.simplefilter("ignore")
import pysynphot as S

import Shen20


#Start by getting L1450 for an Lstar quasar.
qlf = Shen20.QLF(model="B")
z = 2.5
Lbol = 10.**(qlf.log_Lstar(z))*qlf.Lstar_units
nu_1450 = c/(1450.*u.AA)
Lnu_1450 = (qlf.L1450(Lbol)/nu_1450).to(u.erg/u.s/u.Hz)

#Now, load the vanden Brek composite.
#Read the spectrum.
qso_spec = Table.read("vandenberk_composite.txt",format='ascii.cds')
qso = S.ArraySpectrum(wave=qso_spec['Wave'], waveunits='angstrom', flux=qso_spec['FluxD']*1e-17, fluxunits='flam')

#Create the 1450 filter. Since the filter will cover the rest-frame 1450 angstrom flux density, we need it to be redshifted with the spectrum.
lam_rest_b1450 = np.arange(1350., 1550. , 1.)*u.AA
R_1450 = np.zeros(len(lam_rest_b1450))
R_1450[(lam_rest_b1450>=1400.*u.AA) & (lam_rest_b1450<=1500.*u.AA)] = 1

#Load the LSST filter curves
filters = ['u', 'g', 'r', 'i', 'z']
filtercurve = dict()
for filter in filters:
    data = np.loadtxt("LSST_LSST.{}.dat".format(filter), skiprows=1)
    filtercurve[filter] = S.ArrayBandpass(data[:,0], data[:,1], name="{}band".format(filter))


#Get m* for a range in z and all LSST bands.
zmin = 0.01
zmax = 6.0
dlogz = 0.01
zs = 10.**(np.arange(np.log10(1+zmin), np.log10(1+zmax), dlogz)) - 1
mstar = np.zeros((len(zs),len(filters)))
for k,z in enumerate(zs):

    #Redshift the qso template to redshift z
    qso_z = qso.redshift(z)

    #Setup the redshifted 1450A band.
    UV_band = S.ArrayBandpass(lam_rest_b1450.to(u.AA).value*(1.+z), R_1450, name='UV_band')

    #Estimate the observed flux at 1450A for an L* galaxy and renormalize the observed spectrum to match that.
    DL = cosmo.luminosity_distance(z)
    fnu_1450 = (Lnu_1450 * (1.+z) / (4.*np.pi*DL**2)).to(u.erg/u.s/u.cm**2/u.Hz)
    qso_z_renorm = qso_z.renorm(fnu_1450.value, 'fnu', UV_band)
    obs_1450 = S.Observation(qso_z_renorm, UV_band, binset=qso_z_renorm.wave)
    #print("QLF  fnu_1450 = {0:7.3e}\n".format(fnu_1450))
    #print("Spec fnu_1450 = {0:7.3e}\n".format(obs_1450.effstim('fnu')))

    #Finally, get the observed magnitudes.
    for j, filter in enumerate(filters):
        try:
            obs = S.Observation(qso_z_renorm, filtercurve[filter], binset=qso_z_renorm.wave)
            mstar[k,j] = obs.effstim('abmag')
        except (S.exceptions.PartialOverlap,S.exceptions.DisjointError):
            mstar[k,j] = 99.

output = np.zeros((len(zs),len(filters)+1))
output[:,0] = zs
output[:,1:]  = mstar
np.savetxt("mstar_z.dat",output,fmt='%15.5f')

import matplotlib.pyplot as plt
for j,filter in enumerate(filters):
    cond = mstar[:,j]<99.
    plt.plot(zs[cond],mstar[cond,j],label='lsst'+filter)
plt.legend()
plt.xscale('log')
plt.xlabel('Redshift')
plt.ylabel('Observed magnitude of L* quasar (AB)')
plt.title('vanden Berk et al. composite, Shen et al. QLF')
plt.savefig('mstar_z.png')
