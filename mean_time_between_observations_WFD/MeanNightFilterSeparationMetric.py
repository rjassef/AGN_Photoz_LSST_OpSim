#Modified from the metric written by Phil Marshall for lensing quasars, available as a contributed metric.

from lsst.sims.maf.metrics import BaseMetric
import numpy as np

__all__ = ['MeanNightSeparationMetric']

class MeanNightFilterSeparationMetric(BaseMetric):

    def __init__(self, filter1, filter2, nightCol='night', 
                 FilterCol='filter', **kwargs):
        """
        nightCol = the name of the column defining the visit night
        FilterCol = the name of the column with the filters.
        """
        self.nightCol = nightCol
        self.FilterCol = FilterCol
        self.filter1 = filter1
        self.filter2 = filter2
        
        metricName="MeanNightFilterSeparationMetric_"+filter1+"_"+filter2
        
        super(MeanNightFilterSeparationMetric, self).__init__(col=[self.nightCol, self.FilterCol],\
                                                        metricName=metricName, **kwargs)

    #For each element in array a, find the closest element in array b.
    def find_min(self, a, b):
        aa = np.tile(a,(len(b),1))
        bb = np.tile(b,(len(a),1)).T
        kmin = np.argmin(np.abs(aa-bb), axis=0)
        return b[kmin]
        
    def run(self, dataSlice, slicePoint=None):
          
        # Extract the nights for filter 1 and filter 2:
        condition_1 = (dataSlice[self.FilterCol] == self.filter1)
        condition_2 = (dataSlice[self.FilterCol] == self.filter2)
            
        # Find unique nights:
        nights_1 = np.atleast_1d(np.unique(dataSlice[self.nightCol][condition_1]))
        nights_2 = np.atleast_1d(np.unique(dataSlice[self.nightCol][condition_2]))

        #At least one night of observations is needed, so exit with a badval if no 
        #observations were conducted in a given band..
        if len(nights_1)==0 or len(nights_2)==0:
            return self.badval
        
        #For each night_1, find the closest night_2
        nights_2_closest_to_1 = self.find_min(nights_1, nights_2)
     
        #Find the difference in nights in absolute terms.
        dnights = np.abs(nights_1 - nights_2_closest_to_1)

        #Find the average number of nights.
        campaignMean = np.average(dnights)

        return campaignMean # in days

# ======================================================================