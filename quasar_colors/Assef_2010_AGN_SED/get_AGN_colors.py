#!/usr/bin/env python

import numpy as np
from astropy.cosmology import Planck13 as cosmo
from astropy.table import Table

import SED_Model
import lrt

#Setup the cosmology file.
cosmo_file = open("cosmo.dat","w")
cosmo_file.write("{}\n".format(cosmo.Om0)) #Omega_Matter
cosmo_file.write("{}\n".format(1.0-cosmo.Om0)) #Omega_Lambda
cosmo_file.write("{}\n".format(0.0)) #Omega_k
cosmo_file.write("{}\n".format(cosmo.H0.value)) #H0
cosmo_file.close()

#Read the photometry zero point.
jyzero = np.loadtxt("bandmag.dat",usecols=[2])
filters = np.genfromtxt("bandmag.dat",usecols=[0],dtype='U')

#Set the redshift grid.
zmin = 0.01
zmax = 6.0
dz   = 0.01
zs = np.arange(zmin,zmax,dz)

#Calculate the quasar template colors.
qso_table = np.zeros((len(zs),len(filters)))
qso_table[:,0] = zs
for k,z in enumerate(zs):

    #Create the qso object.
    qso = SED_Model.lrt_model()
    qso.zspec = z
    qso.comp = np.zeros(4)
    qso.comp[0] = 1.0
    qso.ebv = 0.0
    qso.igm = 1.0

    #Get the unnormalized observed magnitudes.
    qso.get_model_fluxes()

    #Transform them into magnitudes.
    mag = -2.5*np.log10(qso.jymod/jyzero)

    #Save the colors.
    for j in range(len(filters[:-1])):
        #Check first if the band actually overlaps with the SED at those wavelengths.
        if lrt.cal1.lbar[j]/(1.+qso.zspec) < 0.09:
            qso_table[k,j+1] = np.nan
        else:
            qso_table[k,j+1] = mag[j]-mag[j+1]


#Create a table and save it as ascii.
cato = open("Assef_2010_AGN_SED_LSST_colors.dat","w")
cato.write("# Redshift")
for j in range(len(filters[:-1])):
    cato.write(" {0}-{1}".format(filters[j][-1],filters[j+1][-1]))
cato.write("\n")
np.savetxt(cato,qso_table,fmt='%10.3f')
cato.close()
