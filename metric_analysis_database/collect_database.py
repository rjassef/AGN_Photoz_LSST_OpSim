#!/usr/bin/env python

import numpy as np
from astropy.table import Table

import os
import re
root_path = re.search("(.*/AGN_Photoz_LSST_OpSim)/*",os.getcwd()).group(1)

header = dict()

header["/uband_depth/Table_ug_band_depth.txt"] = [
    "run",
    "u_depth_3sig",
    "u_depth_3sig_rank",
    "g_depth_5sig",
    "g_depth_5sig_rank"
]

for z in ["1.0","2.0","3.0"]:
    header["/color_excess/Table_color_excess_z{}.txt".format(z)] = [
        "run",
        "cex_ug_{}".format(z),
        "cex_ug_{}_rank".format(z),
        "cex_gr_{}".format(z),
        "cex_gr_{}_rank".format(z),
        "cex_ri_{}".format(z),
        "cex_ri_{}_rank".format(z),
        "cex_iz_{}".format(z),
        "cex_iz_{}_rank".format(z),
        "cex_zy_{}".format(z),
        "cex_zy_{}_rank".format(z)
    ]

header["/mean_time_between_obs/Table_mean_nights_between_bands.txt"] = [
    "run",
    "nights_ug",
    "nights_ug_rank",
    "nights_gr",
    "nights_gr_rank",
    "nights_ri",
    "nights_ri_rank",
    "nights_iz",
    "nights_iz_rank",
    "nights_zy",
    "nights_zy_rank"
]

#Use the first medians table to create the master output table. This table should contain all the runs that other tables contain.
tables = list(header.keys())
master_table = Table.read(root_path+tables[0], format='ascii')
for j, col in enumerate(header[tables[0]]):
    master_table.rename_column('col{}'.format(j+1), col)

#Now, add the rest. Note that there is no certainty that all the runs are in each table, nor that they are in the same order.
for table in tables[1:]:
    mds_table = Table.read(root_path+table, format='ascii')
    for j, col in enumerate(header[table]):
        mds_table.rename_column('col{}'.format(j+1), col)

    for col in header[table][1:]:
        if col in master_table.colnames:
            print("Error: repeated column ",col)
            exit()
        master_table[col] = np.zeros(len(master_table['run']))

    for k, run in enumerate(mds_table['run']):
        i = np.argwhere(master_table['run']==run)[0][0]
        for j,col in enumerate(header[table][1:]):
            master_table[col][i] = mds_table[col][k]

master_table.write("metrics_database.fits", format='fits')
