#!/usr/bin/env python

import numpy as np

import Shen20
from Nqso import Nqso

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
