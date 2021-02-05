# AGN_Photoz_LSST_OpSim
Metric development for the Photo-z workgroup of the LSST AGN Science Collaboration.

### Description

This repository has routines and metrics developed to assess how the different OpSims developed for LSST affect the estimation of AGN photometric redshifts. Due to the complicated nature of determining AGN photo-zs, instead of attempting a full simulation we developed simple metrics that should provide a clear idea of how different aspects of the observations would affect these estimates.

To all the estimates done here, one should add the tests for [DCR conducted by Gordon Richard's group](https://github.com/RichardsGroup/LSST_DCR) to get a full picture of LSST in the context of AGN photo-zs.

All calculations here use routines developed for the [OpSims tutorial by Gordon Richard's group](https://github.com/RichardsGroup/LSST_OpSim).

Below, this document describes the metrics developed and links to the results. 

### u-band depth

A key spectral feature of AGNs that helps determine the photometric redshift is the Lyman break. More specifically, the drop in flux shortwards of Lyman alpha. The two bands that are expected to be the shallowest are *u* and *y*, and hence the exact *u*-band depth critically determines the number of AGNs for which we will be able to detect this feature in the range z~1.5-3 . This redshift range covers the peak of luminous AGN activity. 

The assumption, then, is that to detect the Lyman break color with at least a 3 sigma significance, we need to detect *u* at the 3 sigma level, as *g* will always be much deeper. 

[This notebook](https://github.com/rjassef/AGN_Photoz_LSST_OpSim/blob/main/uband_depth/Lstar_depth_dust.ipynb) carries out the depth calculations, and is based on [this metric](https://github.com/rjassef/AGN_Photoz_LSST_OpSim/blob/main/uband_depth/ExgalM5_with_cuts_AGN.py). 

[This notebook](https://github.com/rjassef/AGN_Photoz_LSST_OpSim/blob/main/uband_depth/Visualize_Lstar_depth.ipynb) shows the results both as a the actual 3 sigma depth of *u* and the associated luminosity of the faintest AGN that would be detected as a function of L*. For comparison, it also shows the 5 sigma depth of *g* and the associated luminosity in terms of L*
