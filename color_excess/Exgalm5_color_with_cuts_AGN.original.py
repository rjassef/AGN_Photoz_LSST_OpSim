import numpy as np
import lsst.sims.maf.metrics as metrics
from lsst.sims.maf.metrics import BaseMetric

class Exgalm5_color_with_cuts_AGN(BaseMetric):
    
    def __init__(self, lsstFilters, m5Col='fiveSigmaDepth', \
                 extinction_cut=1.0, \
                 metricName='Exgalm5_color_with_cuts_AGN', func=np.mean, **kwargs):
       
        #Fail initialization if there are less or more than 2 filters.
        if len(lsstFilters)!=2:
            print("Failed to initialize. Need to declare two, and only two, filters.")
            return
    
        #Now, for each of the bands, we want to get the 5 sigma depth. 
        #Also, add the band name to the metric name
        metricName += "_"
        self.ExgalM5_bands = list()
        self.lsstFilters = lsstFilters
        for lsstFilter in self.lsstFilters:
            self.ExgalM5_bands.append(metrics.ExgalM5(lsstFilter=lsstFilter))
            metricName += lsstFilter
       
        #Maximum extinction to be tolerated
        self.extinction_cut = extinction_cut
    
        #Function to apply to the color excess array.
        self.func = func
    
        #Add the function name to the metric name
        metricName += "_"+func.__name__
    
        #Initiate the metric.
        super(Exgalm5_color_with_cuts_AGN, self).__init__(
            col=m5Col, metricName=metricName, maps=self.ExgalM5_bands[0].maps, **kwargs)
        
    def run(self, dataSlice, slicePoint=None):
        
        # exclude areas with high extinction
        if slicePoint['ebv'] > self.extinction_cut:
            return self.badval
        
        # Get the m5 depths.
        mlim5_bands = list()
        for ExgalM5_band in self.ExgalM5_bands:
            mlim5_bands.append(ExgalM5_band.run(dataSlice, slicePoint))
            
        #Get the color of 5sigma limits.
        color_lim = mlim5_bands[0]-mlim5_bands[1]
        return color_lim

