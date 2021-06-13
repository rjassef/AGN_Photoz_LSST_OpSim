import numpy as np
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
from lsst.sims.maf.metrics import BaseMetric

root_path = "/home/idies/workspace/Storage/rjassef/persistent/AGN_Photoz_LSST_OpSim/"

import os
import sys
sys.path.append(root_path+"../LSST_OpSim/Scripts_NBs/")
from opsimUtils import *

sys.path.append(root_path+"uband_depth")
from ExgalM5_with_cuts_AGN import ExgalM5_with_cuts_AGN

sys.path.append(root_path)
from script_utils import find_completed_runs


#We will use the same slicer and constraint for each metric. 
NSIDE=64
slicer = slicers.HealpixSlicer(nside=NSIDE)
constraint = 'note not like "DD%"' #remove DDFs

#Setup the metric bundles.
EM5 = list()
filters = ['u', 'g', 'r', 'i', 'z', 'y']
for filter in filters:
    metric = ExgalM5_with_cuts_AGN(lsstFilter=filter, metricName='EM5_{}'.format(filter))
    constraint_use = 'filter = "{}" and '.format(filter) + constraint
    EM5.append(metricBundles.MetricBundle(metric, slicer, constraint_use))

#Setup the bundle dictionary
bundleDict = dict()
for k,filter in enumerate(filters):
    bundleDict['EM5_{}'.format(filter)] = EM5[k]
    

#Setup the output folders.
your_username = "rjassef"
folder_mafoutput = "EM5_{0:d}_v2".format(NSIDE)
outDir = '/home/idies/workspace/Temporary/{0}/scratch/MAFOutput/{1}'.format(your_username, folder_mafoutput)
if not os.path.exists(os.path.abspath(outDir)):
    os.mkdir(os.path.abspath(outDir))

resultDbPath  = '/home/idies/workspace/Temporary/{0}/scratch/MAFOutput/{1}'.format(
    your_username, folder_mafoutput)
metricDataPath = '/home/idies/workspace/Temporary/{0}/scratch/MAFOutput/{1}/MetricData/'.format(
    your_username, folder_mafoutput)

#Find the list of completed runs.
n_metrics = len(EM5)
completed_runs = find_completed_runs(n_metrics, resultDbPath, metricDataPath)

#Run all OpSims in FBS 1.5, 1.6 and 1.7
FBS_versions = ["1.5","1.6","1.7"]
for FBS_version in FBS_versions:
    dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)
    opSimDbs, resultDbs = connect_dbs(dbDir, outDir)
    dbRuns = show_opsims(dbDir)
    for run in dbRuns:
        if run in completed_runs:
            continue
        for EM5_filt in EM5:
            EM5_filt.setRunName(run)
        metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
                        opSimDbs[run], metricDataPath, resultDbs[run])
        metricGroup.runAll()