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

[This notebook](https://github.com/rjassef/AGN_Photoz_LSST_OpSim/blob/main/uband_depth/Visualize_Lstar_depth.ipynb) shows the results both as a the actual 3 sigma depth of *u* and the associated luminosity of the faintest AGN that would be detected as a function of L*. For comparison, it also shows the 5 sigma depth of *g* and the associated luminosity in terms of L*.

### Color Excess

A more general view of this is that accurate colors are needed to detect the spectral features of AGN in broad-band SEDs. The different depths of the bands will determine how accurately we can determine the colors of the faintest AGN at a given redshift. The metric then determines the color between the 5sigma depths of two consecutive bands, and in the visualization notebooks we compare those depths to the expected colors for an AGN as a function of redshift. 

Colors from two models are implemented. One uses the AGN SED template of [Assef et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...713..970A/abstract) with a standard amount of AGN absorption in their model. The other is based on the models (to be published very soon) by Matthew Temple for three *i*-band magnitudes: 21.5, 24.5 and 26.0. 

[This script](https://github.com/rjassef/AGN_Photoz_LSST_OpSim/blob/main/color_excess/script_Color_Excess.py) was used for carrying out the calculations, which are also implemented in [this notebook](https://github.com/rjassef/AGN_Photoz_LSST_OpSim/blob/main/color_excess/Color_Excess.ipynb). Both use [this metric] (https://github.com/rjassef/AGN_Photoz_LSST_OpSim/blob/main/color_excess/Exgalm5_color_with_cuts_AGN.py). 

[This notebook](https://github.com/rjassef/AGN_Photoz_LSST_OpSim/blob/main/color_excess/Visualize_Color_Excess.ipynb) shows the results for z=2. *Note that values closest to 0 are best, as that means that the 5sigma depths are in the exact same ratio as the expected AGN colors.*

[This notebook](https://github.com/rjassef/AGN_Photoz_LSST_OpSim/blob/main/color_excess/Redshift_Color_Excess.ipynb) shows the median color excess as a function of redshift for all the models.

### Mean time between Observations

AGN variability can have a major impact for photo-zs.
