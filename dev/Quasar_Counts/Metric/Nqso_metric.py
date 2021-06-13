import numpy as np
import astropy.units as u
import lsst.sims.maf.metrics as metrics
from lsst.sims.maf.metrics import BaseMetric

import os
import re
root_path = re.search("(.*/AGN_Photoz_LSST_OpSim)/*",os.path.realpath(__file__)).group(1)

import sys
sys.path.append(root_path+"/Quasar_Counts/")
from Fast_Nqso import Fast_Nqso

__all__ = ['NqsoMetric']

class NqsoMetric(BaseMetric):

    """
    This metric returns the number of quasars expected to be detected down to a 5sigma depth in a given band, as judged by tge ExgalM5 metric.

    """

    def __init__(self, lsstFilter, m5Col='fiveSigmaDepth', \
                 extinction_cut=1.0, filterCol='filter', units='mag',\
                 metricName='Nqso', zmin=0.3, zmax=6.7, m_bright=15.7, **kwargs):

        #Save the input parameters.
        self.lsstFilter = lsstFilter
        self.filterCol = filterCol
        self.zmin      = zmin
        self.zmax      = zmax
        self.m_bright  = m_bright

        #Add the LSST filter to the metric name.
        metricName += "_"+lsstFilter

        #Set up the ExgalM5 metric.
        self.EM5 = metrics.ExgalM5(m5Col=m5Col, units=units, lsstFilter=lsstFilter)

        #Create the Fast_Nqso object, which we will call qso_counter. Calculate the number of quasars per square degree.
        area = 1.0*u.deg**2
        self.qso_counter = Fast_Nqso('LSST'+lsstFilter, 'Shen20', 'A', area=area)

        #Maximum extinction to be tolerated
        self.extinction_cut = extinction_cut

        #Initiate the metric.
        super(NqsoMetric, self).__init__(
            col=[m5Col, filterCol], metricName=metricName, maps=self.EM5.maps, **kwargs)

    def run(self, dataSlice, slicePoint=None):

        # exclude areas with high extinction
        if slicePoint['ebv'] > self.extinction_cut:
            return self.badval

        #Get the 5sigma depth.
        dS = dataSlice[dataSlice[self.filterCol] == self.lsstFilter]
        if dS.shape[0]==0:
            return self.badval
        mlim5 = self.EM5.run(dS, slicePoint)
        if self.m_bright > mlim5:
            return 0

        #Get the number of quasars detected to that depth per unit magnitude.
        return self.qso_counter.Nqso(self.zmin, self.zmax, self.m_bright, mlim5)[0]
