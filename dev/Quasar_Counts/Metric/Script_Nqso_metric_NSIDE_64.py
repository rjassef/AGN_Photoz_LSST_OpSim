import numpy as np
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
from lsst.sims.maf.metrics import BaseMetric

import os
import sys
sys.path.append("/home/idies/workspace/Storage/rjassef/persistent/LSST_OpSim/Scripts_NBs/")
from opsimUtils import *

sys.path.append("/home/idies/workspace/Storage/rjassef/persistent/AGN_Photoz_LSST_OpSim/")
from script_utils import find_completed_runs

from Nqso_metric import NqsoMetric

#We will use the same slicer for both bands. We use nside=64 to use less storage. Also, we add
#useCache=False again to deal with another warning.
#
# According to the DustMap script description: https://sims-maf.lsst.io/_modules/lsst/sims/maf/maps/dustMap.html using nside 64 instead of 128 should not be a problem.
NSIDE = 64
slicer = slicers.HealpixSlicer(nside=NSIDE, useCache=False)

#Also, all the metrics will use the same constraint of not looking at DDFs.
constraint = 'note not like "DD%"'

#Set up the metrics for succesive filter pairs.
Nqso_EM5 = list()
#filters = ['u', 'g', 'r', 'i', 'z', 'y']
filters = ['i']
for filter in filters:
    metric = NqsoMetric(filter)
    Nqso_EM5.append(metricBundles.MetricBundle(metric, slicer, constraint))

#Set up the bundleDict
bundleDict = dict()
for k,filter in enumerate(filters):
    bundleDict['Nqso_EM5_{}'.format(filter)] = Nqso_EM5[k]

#Set up folders.
your_username = "rjassef"
folder_mafoutput = "Nqso_EM5_{0:d}".format(NSIDE)
#outDir = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}'.format(your_username,folder_mafoutput)
outDir = '/home/idies/workspace/Temporary/{0}/scratch/MAFOutput/{1}'.format(your_username,folder_mafoutput)
if not os.path.exists(os.path.abspath(outDir)):
    os.mkdir(os.path.abspath(outDir))

#resultDbPath  = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}'.format(your_username, folder_mafoutput)
resultDbPath  = '/home/idies/workspace/Temporary/{0}/scratch/MAFOutput/{1}'.format(your_username, folder_mafoutput)
#metricDataPath = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}/MetricData/'.format(your_username, folder_mafoutput)
metricDataPath = '/home/idies/workspace/Temporary/{0}/scratch/MAFOutput/{1}/MetricData/'.format(your_username, folder_mafoutput)

#Find the list of completed runs.
n_metrics = len(Nqso_EM5)
completed_runs = find_completed_runs(n_metrics, resultDbPath, metricDataPath)

#Run for FBS 1.5
#FBS_version = "1.5"
for FBS_version in ["1.5", "1.6", "1.7"]:

    dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)

    opSimDbs, resultDbs = connect_dbs(dbDir, outDir)

    dbRuns = show_opsims(dbDir)
    for run in dbRuns:
        if run in completed_runs:
            continue
        for k in range(len(filters)):
            Nqso_EM5[k].setRunName(run)
        metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
                        opSimDbs[run], metricDataPath, resultDbs[run])
        metricGroup.runAll()

#Repeat for FBS 1.6
# FBS_version = "1.6"
# dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)

# opSimDbs, resultDbs = connect_dbs(dbDir, outDir)

# dbRuns = show_opsims(dbDir)
# for run in dbRuns:
#     if run in completed_runs:
#         continue
#     for k in range(len(filters)):
#         Nqso_EM5[k].setRunName(run)
#     metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
#                     opSimDbs[run], metricDataPath, resultDbs[run])
#     metricGroup.runAll()

# #Repeat for FBS 1.7
# FBS_version = "1.7"
# dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)

# opSimDbs, resultDbs = connect_dbs(dbDir, outDir)

# dbRuns = show_opsims(dbDir)
# for run in dbRuns:
#     if run in completed_runs:
#         continue
#     for k in range(len(filters)):
#         Nqso_EM5[k].setRunName(run)
#     metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
#                     opSimDbs[run], metricDataPath, resultDbs[run])
#     metricGroup.runAll()
