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
sys.path.append("/home/idies/workspace/Storage/rjassef/persistent/AGN_Photoz_LSST_OpSim/")
from script_utils import find_completed_runs

sys.path.append("/home/idies/workspace/Storage/rjassef/persistent/LSST_OpSim/Scripts_NBs/")
from opsimUtils import *

#We are running a single metric here for every run. 
NSIDE=64
slicer = slicers.HealpixSlicer(nside=NSIDE)
constraint = 'visitExposureTime > 11'
metric = metrics.CountMetric(col='observationStartMJD', metricName='nvisitsLong')

#Setup the bundle. 
mb = metricBundles.MetricBundle(metric, slicer, constraint)
bundleDict = {'VisitCounts':mb}

your_username = "rjassef"
folder_mafoutput = "WFDfootprint_{0:d}".format(NSIDE)
outDir = '/home/idies/workspace/Temporary/{0}/scratch/MAFOutput/{1}'.format(your_username,folder_mafoutput)
if not os.path.exists(os.path.abspath(outDir)):
    os.mkdir(os.path.abspath(outDir))

resultDbPath = '/home/idies/workspace/Temporary/{0}/scratch/MAFOutput/{1}/'.format(
    your_username, folder_mafoutput)
metricDataPath = '/home/idies/workspace/Temporary/{0}/scratch/MAFOutput/{1}/MetricData/'.format(
    your_username, folder_mafoutput)

#Find the list of completed runs.
n_metrics = 1
completed_runs = find_completed_runs(n_metrics, resultDbPath, metricDataPath)

#Run for all the OpSim runs in FBS 1.5, 1.6 and 1.7.
FBS_versions = ["1.5","1.6","1.7"]
for FBS_version in FBS_versions:
    dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)
    opSimDbs, resultDbs = connect_dbs(dbDir, outDir)
    dbRuns = show_opsims(dbDir)
    for run in dbRuns:
        if run in completed_runs:
            continue
        mb.setRunName(run)
        metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
                        opSimDbs[run], metricDataPath, resultDbs[run])
        metricGroup.runAll()
