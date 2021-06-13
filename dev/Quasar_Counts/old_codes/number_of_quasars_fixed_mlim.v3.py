#!/usr/bin/env python

import numpy as np

from Table_Nqso import Nqso as Table_Nqso

print("Using Table interpolation: ")

#Find the number of quasars in between 18 and 19 mag in the i-band between z=0.5 and z=0.6
print("Number of quasars with 18<i<19 and 0.5<redshift<0.6:")
print("{0:,.0f}".format(Table_Nqso(0.50, 0.60, 18.0, 19.0, 'LSSTi')))
print()

#Find the number of quasars brighter than i=25 in the redshift range 1 to 2 over the entire sky.
print("Number of quasars with i<25 and 1<redshift<2:")
print("{0:,.0f}".format(Table_Nqso(1.0, 2.0, 16.0, 25.0, 'LSSTi')))
print()

print("Number of quasars with i<25.15 and 1.7<redshift<2.3:")
print("{0:,.0f}".format(Table_Nqso(1.7, 2.3, 16.0, 25.15, 'LSSTi')))
print()

###

import Shen20
from Nqso import Nqso

print("Using full integrator: ")
qlf = Shen20.QLF(model="B")

#Find the number of quasars in between 18 and 19 mag in the i-band between z=0.5 and z=0.6
print("Number of quasars with 18<i<19 and 0.5<redshift<0.6:")
print("{0:,.0f}".format(      Nqso(0.50, 0.60, 18.0, 19.0, 'LSSTi', qlf)))
print()
#Find the number of quasars brighter than i=25 in the redshift range 1 to 2 over the entire sky.
print("Number of quasars with i<25 and 1<redshift<2:")
print("{0:,.0f}".format(      Nqso(1.0, 2.0, 16.0, 25.0, 'LSSTi', qlf)))
print()

print("Number of quasars with i<25.15 and 1.7<redshift<2.3:")
print("{0:,.0f}".format(      Nqso(1.7, 2.3, 16.0, 25.15, 'LSSTi', qlf)))
print()
