#!/usr/bin/env python

import numpy as np
from astropy.table import Table
import os
import re

def Temple_colors(color,z,imag,ebv=0):

    #We only have two reddening values.
    if ebv==0:
        ebvname = ""
    elif ebv==0.1:
        ebvname = "_ebv01"
    else:
        print("Invalid reddening value.")
        return None

    #Get the path to this module. Necessary so that we can read the model files.
    module_directory = os.path.dirname(os.path.abspath(__file__))

    #Folder where we find the model files.
    model_root = module_directory+"/Matthew_Temple_Models/redshift_colours_imag"

    #First, try to see if the specific selection of magnitude and reddening exist already.
    try:
        fname = model_root+"{0:.0f}".format(imag*10)
        fname += ebvname+".dat"
        tab = Table.read(fname, format='ascii')
        zs = tab['Redshift']
        colors = tab[color]
    except FileNotFoundError:
        print("Invalid i-band magnitude")
        return None


    #Now, interpolate to get the right redshift.
    return np.interp(z,zs,colors)

def Assef10_colors(color,z):

    #Get the path to this module. Necessary so that we can read the model files.
    module_directory = os.path.dirname(os.path.abspath(__file__))

    #Folder where we find the model files.
    model_file = module_directory+"/Assef_2010_AGN_SED/Assef_2010_AGN_SED_LSST_colors.dat"

    tab = Table.read(model_file, format='ascii')
    zs = tab['Redshift']
    colors = tab[color]

    #Now, interpolate to get the right redshift.
    return np.interp(z,zs,colors)
