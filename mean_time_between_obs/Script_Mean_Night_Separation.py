#This is the same as in the Jupyter note book, but implemented as a script for running as a Job. 

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

#Setup the metric.
from MeanNightFilterSeparationMetric import MeanNightFilterSeparationMetric

#We will use the same slicer and constraint for each metric. Since there are no warnings about cache
#and minimum nside, we will use 64 to go faster.
slicer = slicers.HealpixSlicer(nside=64)
constraint = 'note not like "DD%"' # added so the sky plot won't saturate (remove DDFs)

Mean_Night_list = list()
filters = ['u', 'g', 'r', 'i', 'z', 'y']
for k in range(len(filters[:-1])):
    metric_use = MeanNightFilterSeparationMetric(filters[k], filters[k+1])
    constraint_use = constraint + ' and (filter = "{0:s}" or filter = "{1:s}")'.format(
        filters[k], filters[k+1])
    Mean_Night_list.append(metricBundles.MetricBundle(metric_use, slicer, constraint_use))
    
bundleDict = dict()
for k, Mean_Night in enumerate(Mean_Night_list):
    dict_name = "Mean_Night_{0:s}{1:s}".format(filters[k], filters[k+1])
    bundleDict[dict_name] = Mean_Night
    

#Setup folder names.
your_username = "rjassef"
folder_mafoutput = "Mean_Night_Filter_v2"
resultDbPath  = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}'.format(
    your_username, folder_mafoutput)
metricDataPath = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}/MetricData/'.format(
    your_username, folder_mafoutput)

#Find the list of completed runs.
n_metrics = len(Mean_Night_list)
completed_runs = find_completed_runs(n_metrics, resultDbPath, metricDataPath)

#Run for FBS v1.5
FBS_version = "1.5"
dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)
outDir = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}'.format(your_username,folder_mafoutput)

if not os.path.exists(os.path.abspath(outDir)):
    os.mkdir(os.path.abspath(outDir))
   
opSimDbs, resultDbs = connect_dbs(dbDir, outDir)

dbRuns = show_opsims(dbDir)
for run in dbRuns:
    if run in completed_runs:
        continue
    for Mean_Night in Mean_Night_list:
        Mean_Night.setRunName(run)
    metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
                    opSimDbs[run], metricDataPath, resultDbs[run])
    metricGroup.runAll()
    
#Now run for FBS v1.6
FBS_version = "1.6"
dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)

if not os.path.exists(os.path.abspath(outDir)):
    os.mkdir(os.path.abspath(outDir))
    
opSimDbs, resultDbs = connect_dbs(dbDir, outDir)

dbRuns = show_opsims(dbDir)
for run in dbRuns:
    for Mean_Night in Mean_Night_list:
        Mean_Night.setRunName(run)
    metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
                    opSimDbs[run], metricDataPath, resultDbs[run])
    metricGroup.runAll()
    
