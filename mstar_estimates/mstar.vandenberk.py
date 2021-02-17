#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.cosmology import Planck13 as cosmo
from astropy.table import Table

#This convoluted way of importing pysynphot is so that it does not generate a warning and take a long time to load. This happens because some tasks of pysynphot expect a number of files installed in the path 'PYSYN_CDBS' to work correctly. None of those tasks are being used in this script.
import os
os.environ['PYSYN_CDBS'] = "."
import warnings
warnings.simplefilter("ignore")
import pysynphot as S

#Module with the implementation of the Shen et al. (2020) QLF model.
import Shen20


"""
This script creates the table mstar_z.vandenberk.dat, which holds the apparent magnitude of an L* quasar (assuming the Shen et al. 2020 QLF) as a function of redshift for all the LSST bands, as well as the absolute magnitude at 1450A, M_1450.

All output magnitudes are in the AB system.

"""

#Create the QLF object.
qlf = Shen20.QLF(model="B")

#Now, load the vanden Berk composite. The file with the spectrum is the one provided with their journal article.
qso_spec = Table.read("vandenberk_composite.txt",format='ascii.cds')
qso = S.ArraySpectrum(wave=qso_spec['Wave'], waveunits='angstrom', flux=qso_spec['FluxD']*1e-17, fluxunits='flam')

#Create the UV band. Following the description in Table 1 of Shen et al. (2020), the filter is a box-shaped filter with a width of 100A centered at 1450A. Note that since we want the filter to cover the rest-frame 1450 angstrom flux density, we will need to redshift it with the spectrum.
lam_rest_b1450 = np.arange(1350., 1550. , 1.)*u.AA
R_1450 = np.zeros(len(lam_rest_b1450))
R_1450[(lam_rest_b1450>=1400.*u.AA) & (lam_rest_b1450<=1500.*u.AA)] = 1

#Load the LSST filter curves as pysynphot ArrayBandpass objects. These filtercurves were downloaded from SVO: http://svo2.cab.inta-csic.es/svo/theory//fps3/index.php?&mode=browse&gname=LSST&gname2=LSST
filters = ['u', 'g', 'r', 'i', 'z', 'y']
filtercurve = dict()
for filter in filters:
    data = np.loadtxt("LSST_LSST.{}.dat".format(filter), skiprows=1)
    filtercurve[filter] = S.ArrayBandpass(data[:,0], data[:,1], name="{}band".format(filter))

#Redshift grid.
zmin = 0.01
zmax = 6.0
#dlogz = 0.01
#zs = 10.**(np.arange(np.log10(1+zmin), np.log10(1+zmax), dlogz)) - 1
dz = 0.01
zs = np.arange(zmin, zmax, dz)

#For each redshift bin, calculate the apparent magnitude of an L* quasar in a each LSST band.
mstar = np.zeros((len(zs),len(filters)))
M_1450 = np.zeros(len(zs))
for k,z in enumerate(zs):

    #Redshift the qso template to redshift z.
    qso_z = qso.redshift(z)

    #Setup the redshifted UV band as an ArrayBandpass object.
    UV_band = S.ArrayBandpass(lam_rest_b1450.to(u.AA).value*(1.+z), R_1450, name='UV_band')

    #Get L1450 for an L* quasar.
    Lbol = 10.**(qlf.log_Lstar(z))*qlf.Lstar_units

    #Tranform the monochromatic luminosity into a luminosity density at 1450A.
    nu_1450 = c/(1450.*u.AA)
    Lnu_1450 = (qlf.L1450(Lbol)/nu_1450).to(u.erg/u.s/u.Hz)

    #Get the AB absolute magnitude in the UV band.
    M_1450[k] = -2.5*np.log10(Lnu_1450/(4.*np.pi*(10.*u.pc)**2)/(3631*u.Jy))

    #Transform the 1450A luminosity density of the L* quasar to the observed flux density.
    DL = cosmo.luminosity_distance(z)
    fnu_1450 = (Lnu_1450 * (1.+z) / (4.*np.pi*DL**2)).to(u.erg/u.s/u.cm**2/u.Hz)

    #Finally, renormalize the redshifted qso spectrum to have a flux_density of fnu_1450 in the UV band.
    qso_z_renorm = qso_z.renorm(fnu_1450.value, 'fnu', UV_band)
    obs_1450 = S.Observation(qso_z_renorm, UV_band, binset=qso_z_renorm.wave)

    #Finally, convolve the redshifted quasar template with the LSST filter curves to obtain the observed magnitudes of the L* quasar.
    for j, filter in enumerate(filters):
        try:
            obs = S.Observation(qso_z_renorm, filtercurve[filter], binset=qso_z_renorm.wave)
            mstar[k,j] = obs.effstim('abmag')
        except (S.exceptions.PartialOverlap,S.exceptions.DisjointError):
            #If the spectrum does not overlap or only partially overlaps with the filter curve, set the magnitude to 99. to indicate we could not estimate it.
            mstar[k,j] = 99.

#Print the apparent magnitudes of an L* quasar as a function of redshift.
output = np.zeros((len(zs),len(filters)+2))
output[:,0] = zs
output[:,1:len(filters)+1]  = mstar
output[:,-1] = M_1450
table_file = open("mstar_z.vandenberk.dat","w")
table_file.write("#Redshift")
for filter in filters:
    table_file.write("\tLSST{}".format(filter))
table_file.write("\tM_1450\n")
np.savetxt(table_file,output,fmt='%15.5f')
table_file.close()

#Make a plot of the magnitude of an L* quasar as a function of redshift for each of the LSST bands.
import matplotlib.pyplot as plt
for j,filter in enumerate(filters):
    cond = mstar[:,j]<99.
    plt.plot(zs[cond],mstar[cond,j],label='lsst'+filter)
plt.legend()
#plt.xscale('log')
plt.ylim([13.,25.])
plt.xlabel('Redshift')
plt.ylabel('Observed magnitude of L* quasar (AB)')
plt.title('vanden Berk et al. composite, Shen et al. QLF')
plt.savefig('mstar_z.vandenberk.png')
