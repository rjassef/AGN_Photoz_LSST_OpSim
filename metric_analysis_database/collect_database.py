#!/usr/bin/env python

import numpy as np
from astropy.table import Table

import os
import re
root_path = re.search("(.*/AGN_Photoz_LSST_OpSim)/*",os.getcwd()).group(1)

header = dict()

#First, add the the u-band depths analysis.
header["/uband_depth/Table_ug_band_depth.txt"] = [
    "run",
    "u_depth_3sig",
    "u_depth_3sig_rank",
    "u_depth_3sig_area",
    "g_depth_5sig",
    "g_depth_5sig_rank",
    "g_depth_5sig_area"
]

#Next, add the DDFs runs.
DDFs = ['AllDDFs', 'ECDFS', 'COSMOS', 'XMM-LSS', 'ELAISS1', 'EDFS']
for DDF in DDFs:
    for filter in ["u","g"]:
        sig = "5"
        if filter=="u":
            sig = "3"
        header["/DDFs_uband_depth/Table_DDF_{0}_band_depth.{1}.txt".format(filter, DDF)] = [
            "run",
            "{0}_{1}_depth_{2}sig".format(DDF,filter,sig),
            "{0}_{1}_depth_{2}sig_rank".format(DDF,filter,sig),
            "{0}_{1}_depth_{2}sig_area".format(DDF,filter,sig)
        ]


#Next, add the color excess runs.
filters = ["u","g","r","i","z","y"]
for z in ["1.0","2.0","3.0"]:
    tname = "/color_excess/Table_color_excess_z{}.txt".format(z)
    header[tname] = ["run"]
    for k in range(len(filters[:-1])):
        header[tname].append("cex_{0:s}{1:s}_{2}".format(filters[k], filters[k+1], z))
        header[tname].append("cex_{0:s}{1:s}_{2}_rank".format(filters[k], filters[k+1], z))
        header[tname].append("cex_{0:s}{1:s}_{2}_area".format(filters[k], filters[k+1], z))

#Add the mean nights runs.
tname = "/mean_time_between_obs/Table_mean_nights_between_bands.txt"
header[tname] = ["run"]
for k in range(len(filters[:-1])):
    header[tname].append("nights_{0}{1}".format(filters[k],filters[k+1]))
    header[tname].append("nights_{0}{1}_rank".format(filters[k],filters[k+1]))
    header[tname].append("nights_{0}{1}_area".format(filters[k],filters[k+1]))

#Finally add the number of QSOs detected.
tname = "/Quasar_Counts/Metric/Table_Nqso_i.txt"
header[tname] = [
    "run",
    "Nqso_i",
    "Nqso_i_rank",
    "Nqso_i_area"
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

master_table.write("metrics_database.fits", format='fits', overwrite=True)
