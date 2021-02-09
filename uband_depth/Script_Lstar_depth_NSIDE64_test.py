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

from ExgalM5_with_cuts_AGN import ExgalM5_with_cuts_AGN

#We will use the same slicer for both bands. We use nside=64 to use less storage. Also, we add 
#useCache=False again to deal with another warning.
#
# According to the DustMap script description: https://sims-maf.lsst.io/_modules/lsst/sims/maf/maps/dustMap.html using nside 64 instead of 128 should not be a problem.
NSIDE = 64
slicer_ug = slicers.HealpixSlicer(nside=NSIDE, useCache=False)

#Set up the MAF for u-band, 5sigma depth.
metric_u = ExgalM5_with_cuts_AGN(lsstFilter='u', metricName='ExgalM5_with_cuts_AGN_u')
constraint_u = 'filter = "u"'
constraint_u += ' and note not like "DD%"' # added so the sky plot won't saturate (remove DDFs)
EM5u = metricBundles.MetricBundle(metric_u, slicer_ug, constraint_u)

#Set up the MAF for g-band, 5sigma depth.
metric_g = ExgalM5_with_cuts_AGN(lsstFilter='g', metricName='ExgalM5_with_cuts_AGN_g')
constraint_g = 'filter = "g"'
constraint_g += ' and note not like "DD%"' # added so the sky plot won't saturate (remove DDFs)
EM5g = metricBundles.MetricBundle(metric_g, slicer_ug, constraint_g)

bundleDict = {'EM5u': EM5u, 'EM5g': EM5g}

#Setup folders 
your_username = "rjassef"
folder_mafoutput = "EM5_depths_{0:d}".format(NSIDE)

outDir = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}'.format(your_username,folder_mafoutput)
if not os.path.exists(os.path.abspath(outDir)):
    os.mkdir(os.path.abspath(outDir))
    
metricDataPath = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}/MetricData/'.format(
    your_username, folder_mafoutput)
resultDbPath  = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}'.format(
    your_username, folder_mafoutput)

#Find the list of completed runs.
n_metrics = len(list(bundleDict.keys()))
print(n_metrics)
completed_runs = find_completed_runs(n_metrics, resultDbPath, metricDataPath)

#Run for FBS 1.5 
FBS_version = "1.5"
dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)

opSimDbs, resultDbs = connect_dbs(dbDir, outDir)

dbRuns = show_opsims(dbDir)
for run in dbRuns:
    if run in completed_runs:
        continue
    EM5u.setRunName(run)
    EM5g.setRunName(run)
    metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
                    opSimDbs[run], metricDataPath, resultDbs[run])
    metricGroup.runAll()
    
#Repeat for FBS 1.6
FBS_version = "1.6"
dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)

opSimDbs, resultDbs = connect_dbs(dbDir, outDir)

dbRuns = show_opsims(dbDir)
for run in dbRuns:
    if run in completed_runs:
        continue
    EM5u.setRunName(run)
    EM5g.setRunName(run)
    metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
                    opSimDbs[run], metricDataPath, resultDbs[run])
    metricGroup.runAll()
