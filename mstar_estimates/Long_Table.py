#/usr/bin/env python

import numpy as np
from astropy.table import Table
import astropy.units as u
import multiprocessing as mp
from functools import partial

import Shen20
from Nqso import Nqso

def get_Nqso(z, mi, qlf, mstar_data, Mi_lim, mi_index_use):
    N = np.zeros((len(mi_index_use),len(z)-1))
    for j,i in enumerate(mi_index_use):
        for k in range(len(z[:-1])):
            N[j,k] = Nqso(z[k], z[k+1], mi[i], mi[i+1], 'LSSTi', qlf, mstar_data=mstar_data, Mi_lim=Mi_lim)
    return N


mstar_data = Table.read("mstar_z.vandenberk.dat", format='ascii')

#Create the QLF object.
qlf = Shen20.QLF(model="B")


z  = np.arange(0.1, 7.0, 0.1)
mi = np.arange(15.7, 26.3, 0.1)
#z  = np.arange(5.0, 6.0, 0.2)
#mi = np.arange(20.0, 22.0, 0.2)

Ncpu = mp.cpu_count()-1
mi_index_use = np.arange(len(mi)-1)
mi_index_use_split = np.array_split(mi_index_use, Ncpu)

Mi_lim = -20
Pool = mp.Pool(Ncpu)
func = partial(get_Nqso, z, mi, qlf, mstar_data, Mi_lim)

Output = Pool.map(func, mi_index_use_split)
Output = np.vstack(Output)
Pool.close()

cato = open("Long_Table.txt","w")
for i in range(len(mi)):
    cato.write("{0:5.1f} ".format(mi[i]))
cato.write("\n")
for k in range(len(z)):
    cato.write("{0:5.1f} ".format(z[k]))
cato.write("\n")

for i in range(len(mi[:-1])):
    for k in range(len(z[:-1])):
        cato.write("{0:12.2f}".format(Output[i,k]))
    cato.write("\n")

cato.close()
