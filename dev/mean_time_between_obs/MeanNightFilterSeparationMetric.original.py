#Modified from the metric written by Phil Marshall for lensing quasars, available as a contributed metric.

from lsst.sims.maf.metrics import BaseMetric
import numpy as np

__all__ = ['MeanNightSeparationMetric']

class MeanNightFilterSeparationMetric(BaseMetric):

    def __init__(self, filter1, filter2, seasonCol='season', nightCol='night', 
                 FilterCol='filter', **kwargs):
        """
        seasonCol = the name of the column defining the season number
        nightCol = the name of the column defining the visit night
        FilterCol = the name of the column with the filters.
        """
        self.seasonCol = seasonCol
        self.nightCol = nightCol
        self.FilterCol = FilterCol
        self.filter1 = filter1
        self.filter2 = filter2
        
        metricName="MeanNightFilterSeparationMetric_"+filter1+"_"+filter2
        
        super(MeanNightFilterSeparationMetric, self).__init__(col=[self.seasonCol, \
                                                             self.nightCol, self.FilterCol],\
                                                        metricName=metricName, **kwargs)

    #For each element in array a, find the closest element in array b.
    def find_min(self, a, b):
        aa = np.tile(a,(len(b),1))
        bb = np.tile(b,(len(a),1)).T
        kmin = np.argmin(np.abs(aa-bb), axis=0)
        return b[kmin]
        
    def run(self, dataSlice, slicePoint=None):

        # Small loop over the seasons:
        uniqSeasons = np.unique(dataSlice[self.seasonCol])
        seasonMeans = np.array([])
        for k in uniqSeasons:
            
            # Extract this season's nights for filter 1 and filter 2:
            condition = (dataSlice[self.seasonCol] == k)
            condition_1 = (condition) & (dataSlice[self.FilterCol] == self.filter1)
            condition_2 = (condition) & (dataSlice[self.FilterCol] == self.filter2)
            
            # Find unique nights, and sort:
            nights_1 = np.atleast_1d(np.sort(np.unique(dataSlice[self.nightCol][condition_1])))
            nights_2 = np.atleast_1d(np.sort(np.unique(dataSlice[self.nightCol][condition_2])))
       
         
            # Check the number of observing nights this season:
            # print "season, nights: ",k,nights
            if len(nights_1) == 0 or len(nights_2) == 0:
                # print "    zero nights, continuing..."
                continue
            else:
                #For each night_1, find the closest night_2
                nights_2_closest_to_1 = self.find_min(nights_1, nights_2)

                #Find the difference in nights in absolute terms.
                dnights = np.abs(nights_1 - nights_2_closest_to_1)

                #Find the average number of nights.
                thisSeasonMean = np.average(dnights)

            seasonMeans = np.append(seasonMeans,thisSeasonMean)
            
        # Take average over seasons, defensively:
        # print "taking average over the following seasons' mean separations:",seasonMeans
        if len(seasonMeans) > 0:
            campaignMean = np.average(seasonMeans)
        else:
            campaignMean = 0.0
        # print "result =",campaignMean
        return campaignMean # in days

# ======================================================================