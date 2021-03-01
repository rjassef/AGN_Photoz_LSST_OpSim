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

#Module with the implementation of the QLF model.
import Hopkins07


"""
This script creates the table mstar_z.vandenberk.dat, which holds the apparent magnitude of an L* quasar (assuming the Hopkins et al. 2007 QLF) as a function of redshift for all the LSST bands, as well as the absolute magnitude at B-bands, M_B.

All output magnitudes are in the AB system.

"""

#Create the QLF object.
qlf = Hopkins07.QLF(model="Full")

#Now, load the vanden Berk composite. The file with the spectrum is the one provided with their journal article.
qso_spec = Table.read("../../SEDs/vandenberk_composite.txt",format='ascii.cds')
qso = S.ArraySpectrum(wave=qso_spec['Wave'], waveunits='angstrom', flux=qso_spec['FluxD']*1e-17, fluxunits='flam')

#Read the B-band filter curve.
B_curve    = np.loadtxt("../../Filter_Curves/B_bessell.filter", skiprows=1)
lam_rest_B = B_curve[:,0] * u.AA
R_B        = B_curve[:,1]

#Load the LSST filter curves as pysynphot ArrayBandpass objects. These filtercurves were downloaded from SVO: http://svo2.cab.inta-csic.es/svo/theory//fps3/index.php?&mode=browse&gname=LSST&gname2=LSST
filters = ['u', 'g', 'r', 'i', 'z', 'y']
filtercurve = dict()
for filter in filters:
    data = np.loadtxt("../../Filter_Curves/LSST_LSST.{}.dat".format(filter), skiprows=1)
    filtercurve[filter] = S.ArrayBandpass(data[:,0], data[:,1], name="{}band".format(filter))

#Calculate the m_B - m_i color at z=0 so that we can easily transform M_B to M_i later in the code. Note that M_i for an L* quasar is a needed quantity to replicate Table 10.2 of https://www.lsst.org/sites/default/files/docs/sciencebook/SB_10.pdf. See code Table10_2.py for further details.
B_band  = S.ArrayBandpass(lam_rest_B.to(u.AA).value, R_B, name='B_band')
obs_B   = S.Observation(qso, B_band, binset=qso.wave)
obs_i   = S.Observation(qso, filtercurve['i'], binset=qso.wave)
#We define color_B_i = M_B - M_i
color_B_i = obs_B.effstim('abmag') - obs_i.effstim('abmag')

#Redshift grid.
zmin = 0.01
#zmax = 6.0
zmax = 7.0
#dlogz = 0.01
#zs = 10.**(np.arange(np.log10(1+zmin), np.log10(1+zmax), dlogz)) - 1
dz = 0.01
zs = np.arange(zmin, zmax, dz)

#For each redshift bin, calculate the apparent magnitude of an L* quasar in a each LSST band.
mstar = np.zeros((len(zs),len(filters)))
M_B = np.zeros(len(zs))
M_i = np.zeros(len(zs))
for k,z in enumerate(zs):

    #Redshift the qso template to redshift z. Note that qso is an ArraySpectrum pysynphot object, so we use the pysynphot.ArraySpectrum.redshift() method to create a redshifted version of the spectrum.
    qso_z = qso.redshift(z)

    #Setup the redshifted B band as a pysynphot.ArrayBandpass object.
    B_band = S.ArrayBandpass(lam_rest_B.to(u.AA).value*(1.+z), R_B, name='B_band')

    #Get the bolometric luminosity of an L* quasar.
    Lbol = 10.**(qlf.log_Lstar(z))*qlf.Lstar_units

    #Get the B monochromatic luminosity associated to the L* quasar given its bolometric luminosity. See Hopkins07.L_B() method documentation for further details.
    L_B = qlf.L_B(Lbol)

    #Tranform the monochromatic luminosity into a luminosity density at 4380A, the effective wavelength of the B-band used (see Table A2 of Bessell et al. 1998, A&A, 333, 231), which is defined as nu * Lnu (nu=4380A).
    nu_B  = c/(4380.*u.AA)
    Lnu_B = L_B/nu_B

    #Get the AB absolute magnitude in the B band.
    #The flux density, fnu, observed at a luminosity distance DL for a source of luminosity density Lnu is given by:
    #
    # Lnu = (4pi DL^2)/(1+z) * fnu
    #
    # This corresponds to eqn. (6) of Hogg et al. (2002, arXiv:astro-ph/0210394). See that reference for further discussion.
    #
    # So at 10pc, we can assume that the redshift is 0, and we have that the observed flux density would be:
    #
    # Fnu = Lnu / [4pi (10pc)^2]
    Fnu_at_10pc = Lnu_B / (4.*np.pi*(10.*u.pc)**2)
    #Now, the AB magnitude is just m = -2.5 log10 ( fnu / 3631 Jy) by definition. So the AB absolute magnitude in B is just:
    M_B[k] = -2.5*np.log10( Fnu_at_10pc / (3631*u.Jy) )

    #And we can easily get the absolute i band magnitude of the L* quasar by using the color term color_B_i = M_B - M_i estimated earlier.
    M_i[k] = M_B[k] - color_B_i

    #The task now is to renormalize our redshifted spectrum so that it has a BA monochromatic luminosity equal to that of an L* quasar.
    #Start by transforming the B luminosity density of the L* quasar to the observed flux density, using the equation listed above.
    DL = cosmo.luminosity_distance(z)
    fnu_B = (Lnu_B * (1.+z) / (4.*np.pi*DL**2) )

    #Transform fnu_B to be in units of erg / s /cm^2 / Hz and then strip the unit markers. This is important because the renorm method of the redshifted spectrum object expects this units for the flux density provided, but it will not accept astropy units, only floats.
    fnu_B = fnu_B.to(u.erg/u.s/u.cm**2/u.Hz)
    #Strip the units.
    fnu_B = fnu_B.value

    #Finally, renormalize the redshifted qso spectrum to have a flux_density of fnu_B in the B band. For this, we used the renorm method of the pysynphot.ArraySpectrum class.
    qso_z_renorm = qso_z.renorm(fnu_B, 'fnu', B_band)

    #Finally, convolve the redshifted quasar template with the LSST filter curves to obtain the observed magnitudes of the L* quasar.
    for j, filter in enumerate(filters):
        try:
            obs = S.Observation(qso_z_renorm, filtercurve[filter], binset=qso_z_renorm.wave)
            mstar[k,j] = obs.effstim('abmag')
        except (S.exceptions.PartialOverlap,S.exceptions.DisjointError):
            #If the spectrum does not overlap or only partially overlaps with the filter curve, set the magnitude to 99. to indicate we could not estimate it.
            mstar[k,j] = np.nan

#Print the apparent magnitudes of an L* quasar as a function of redshift.
output = np.zeros((len(zs),len(filters)+3))
output[:,0] = zs
output[:,1:len(filters)+1]  = mstar
output[:,-2] = M_B
output[:,-1] = M_i
table_file = open("mstar_z.vandenberk.dat","w")
table_file.write("#Redshift")
for filter in filters:
    table_file.write("\tLSST{}".format(filter))
table_file.write("\tM_B")
table_file.write("\tM_i\n")
np.savetxt(table_file,output,fmt='%15.5f')
table_file.close()

#Make a plot of the magnitude of an L* quasar as a function of redshift for each of the LSST bands.
import matplotlib.pyplot as plt
for j,filter in enumerate(filters):
    cond = (~np.isnan(mstar[:,j]))
    plt.plot(zs[cond],mstar[cond,j],label='lsst'+filter)
plt.legend()
#plt.xscale('log')
plt.ylim([13.,25.])
plt.xlabel('Redshift')
plt.ylabel('Observed magnitude of L* quasar (AB)')
plt.title('vanden Berk et al. composite, Hopkins et al. QLF')
plt.savefig('mstar_z.vandenberk.png')
