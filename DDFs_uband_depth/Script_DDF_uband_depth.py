import numpy as np
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
from lsst.sims.maf.metrics import BaseMetric

import re
import os
import sys
sys.path.append("/home/idies/workspace/Storage/rjassef/persistent/LSST_OpSim/Scripts_NBs/")
from opsimUtils import *

sys.path.append("/home/idies/workspace/Storage/rjassef/persistent/AGN_Photoz_LSST_OpSim/uband_depth/")
from ExgalM5_with_cuts_AGN import ExgalM5_with_cuts_AGN

from DDF_utils import my_ddfInfo

####

#Function to construct the metric bundles for a given OpSim, DDF and filter.
def get_mb(filters, DDF_names, opsdb, slicer, run, metricDataPath):     
    
    EM5 = list()
    bundleDict = dict()
    
    for filter in filters:
        for DDF_name in DDF_names:
            
            metricName = 'ExgalM5_with_cuts_AGN_{0}_{1}'.format(DDF_name, filter)
            metric = ExgalM5_with_cuts_AGN(lsstFilter=filter,
                                   metricName=metricName)
           
            constraint = 'filter = "{}"'.format(filter)
            constraintName = '{}'.format(filter)
        
            if DDF_name == 'AllDDFs':
                constraint += ' and proposalId > 1'
                constraintName += "_and_proposalId_gt_1"
            else:
                propids = my_ddfInfo(opsdb, DDF_name)['proposalId']
                if propids is None:
                    continue
                constraint += ' and ('
                constraintName += "_and_"
                for i, propid in enumerate(propids):
                    if i>0:
                        constraint += ' or '
                        constraintName += "_or_"
                    constraint += 'proposalId = {}'.format(propid)
                    constraintName += "proposalId_{}".format(propid)
                constraint += ')'
                
            #Check if this metric has already been run.
            runName = re.sub("\.","_",run)
            fname = "{0}/{1}_{2}_{3}_HEAL.npz".format(metricDataPath, 
                                                      runName, metricName, 
                                                      constraintName)
            if os.path.exists(fname):
                continue
            
            EM5.append(metricBundles.MetricBundle(metric, slicer, constraint))
            bundleDict['EM5_{0}_{1}'.format(DDF_name, filter)] = EM5[-1]
    
    return EM5, bundleDict

####

#Function that runs all the metrics for all the OpSim, keeping track of the specific DDFs. 
def run_calcs(FBS_version, slicer, metricDataPath, filters, DDF_names):
    
    #OpSims folder. 
    dbDir = '/home/idies/workspace/lsst_cadence/FBS_{}/'.format(FBS_version)
    opSimDbs, resultDbs = connect_dbs(dbDir, outDir)
   
    #For each run, get all calculations.
    dbRuns = show_opsims(dbDir)
    for run in dbRuns:     

        print("Processing run ",run)
        
        #If OpSim run has no DDFs, then skip it. 
        opsim_propInfo = opSimDbs[run].fetchPropInfo()[1]
        #if len(opSimDbs[run].fetchPropInfo()[1]['DD']) == 0:
        if 'DD' not in opsim_propInfo.keys() or len(opsim_propInfo['DD'])==0:
            print("No DDFs in this run. Skipping to the next one.")
            continue
        
        #Set up the metric bundles for this run.
        EM5, bundleDict = get_mb(filters, DDF_names, opSimDbs[run], 
                                 slicer, run, metricDataPath)
        if len(EM5)==0:
            continue
            
        for EM5x in EM5:
            EM5x.setRunName(run)
        metricGroup = metricBundles.MetricBundleGroup(bundleDict,\
                        opSimDbs[run], metricDataPath, resultDbs[run])
        metricGroup.runAll()

####


#Use the same slicer for all. 
NSIDE = 64
my_slicer = slicers.HealpixSlicer(nside=NSIDE, useCache=False)

#Set the DDFs and filters to use.
DDF_names = ['AllDDFs', 'ECDFS', 'COSMOS', 'XMM-LSS', 'ELAISS1', 'EDFS']
filters = ['u','g']

#Set up the folders.
your_username = "rjassef"
folder_mafoutput = "DDFs_EM5_depths_{0:d}".format(NSIDE)

metricDataPath = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}/MetricData/'.format(
    your_username, folder_mafoutput)

resultDbPath  = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}'.format(
    your_username, folder_mafoutput)

outDir = '/home/idies/workspace/Storage/{0}/persistent/MAFOutput/{1}'.format(your_username,folder_mafoutput)

if not os.path.exists(os.path.abspath(outDir)):
    os.mkdir(os.path.abspath(outDir))
    
#Run for each FBS
run_calcs("1.5", my_slicer, metricDataPath, filters, DDF_names)
run_calcs("1.6", my_slicer, metricDataPath, filters, DDF_names)
run_calcs("1.7", my_slicer, metricDataPath, filters, DDF_names)