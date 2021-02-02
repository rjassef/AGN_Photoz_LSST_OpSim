import sys 
sys.path.append("/home/idies/workspace/Storage/rjassef/persistent/LSST_OpSim")
from Scripts_NBs.opsimUtils import *

def find_completed_runs(n_metrics, resultDbPath, metricDataPath, bundleDicts=None):
    
    # get a dictionary of resultDb from the results directory.
    resultDbs = getResultsDbs(resultDbPath)
    
    #Load the metrics that have been completed if that dictionary has not been provided as input.
    if bundleDicts is None:
        bundleDicts = dict()
        for runName in resultDbs:
            bundleDicts[runName] = bundleDictFromDisk(resultDbs[runName], runName, metricDataPath)
    
    #Find the ones that have been completed so we do not run them again.
    completed_runs = list()
    for key in list(bundleDicts.keys()):
        if len(list(bundleDicts[key].keys()))==n_metrics:
            completed_runs.append(key)
    return completed_runs