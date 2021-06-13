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

**Note:** We only consider areas with foreground reddening E(B-V)<1.0 mag, which effectively cuts out the Galactic Plane and Galactic Center.

[This script](uband_depth_WFD/Script_Lstar_depth_NSIDE64.py) was used for carrying out the calculations, based on [this metric](uband_depth_WFD/ExgalM5_with_cuts_AGN.py).

[This notebook](uband_depth_WFD/Visualize_Lstar_depth_All_OpSims_WFD_only.ipynb) shows the results both as a the actual 3 sigma depth of *u* and the associated luminosity of the faintest AGN that would be detected as a function of L*. For comparison, it also shows the 5 sigma depth of *g* and the associated luminosity in terms of L*.

Hard copies of all the plots can found [here](uband_depth_WFD/plots_all_opsims_extremes_WFDonly_64).

### u-band Depth for DDFs

The same analysis described above for the WFD has also been done for the DDFs, for each one of them separately (COSMOS, ECDFS, EDFS, ELAIS-S1 and XMM-LSS), as well as for all of them together (AllDDFs).

[This script](uband_depth_DDFs/Script_DDF_uband_depth.py) was used for carrying out the calculations, based on [this metric](uband_depth_WFD/ExgalM5_with_cuts_AGN.py).

[This notebook](uband_depth_DDFs/Visualize_DDF_uband_depth_All_OpSims.ipynb) shows the results both as a the actual 3 sigma depth of *u* and the associated luminosity of the faintest AGN that would be detected as a function of L*. For comparison, it also shows the 5 sigma depth of *g* and the associated luminosity in terms of L*.

Hard copies of all the plots can found [here](uband_depth_DDFs/plots_all_opsims_extremes_64).


### Color Excess

A more general view of this is that accurate colors are needed to detect the spectral features of AGN in broad-band SEDs. The different depths of the bands will determine how accurately we can determine the colors of the faintest AGN at a given redshift. The metric then determines the color between the 5sigma depths of two consecutive bands, and in the visualization notebooks we compare those depths to the expected colors for an AGN as a function of redshift.

Colors from two models are implemented. One uses the AGN SED template of [Assef et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...713..970A/abstract) with a standard amount of AGN absorption in their model. The other is based on the models (to be published very soon) by Matthew Temple for three *i*-band magnitudes: 21.5, 24.5 and 26.0.

[This script](color_excess_WFD/script_Color_Excess_NSIDE_64.py) was used for carrying out the calculations using [this metric](color_excess_WFD/Exgalm5_color_with_cuts_AGN.py).

[This notebook](color_excess_WFD/Visualize_Color_Excess_All_OpSims_WFDonly.ipynb) shows the results for z=2. *Note that values closest to 0 are best, as that means that the 5sigma depths are in the exact same ratio as the expected AGN colors.*

Hard copies of all the plots can found [here](color_excess_WFD/plots_all_opsims_extremes_WFDonly_64).

[This notebook](color_excess_WFD/Redshift_Color_Excess_All_OpSims_WFDonly.ipynb) shows the median color excess as a function of redshift for all the models.

Hard copies of all the plots can found [here](color_excess_WFD/redshift_plots_all_opsims_extremes_WFDonly_64).

### Mean time between Observations

AGN variability can have a major impact for photo-zs by distorting the SED in between visits. This could be modeled as additional uncertainty on the observed colors. The first step taken here is to calculate the average time in between observations between two consecutive bands. Unlike the previous iteration, here we do not separate by seasons to avoid issues with the start/end of seasons (although the issue is unlikely very significant).

[This script](mean_time_between_observations_WFD/Script_Mean_Night_Separation.py) was used to carry out the calculations using [this metric](mean_time_between_observations_WFD/MeanNightFilterSeparationMetric.py), which consists on a modified version of a metric contributed by Phil Marshal for lensing observations.

[This notebook](mean_time_between_observations_WFD/Visualize_Mean_Night_Separation_All_OpSims_WFDonly.ipynb) shows the results.

Hard copies of all the plots can found [here](mean_time_between_observations_WFD/plots_all_opsims_extremes_WFDonly_64).

### Quasar Counts

Combining the survey depths and a quasar luminosity function, we can estimate the number of objects that will be detected in the 10 year survey depth of LSST. We follow a similar approach to that taken to create Table 10.2 of the [LSST Science Book](https://www.lsst.org/scientists/scibook). The main difference from that effort is that we we update the luminosity function to that of [Shen et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.3252S/abstract). 

The survey depths are calculated using [this script](Nqso_WFD_miniSurveys/Metric/EM5_script.py) based on [this metric](uband_depth_WFD/ExgalM5_with_cuts_AGN.py). Then, the quasar counts are estimated and displayed using [this notebook](Nqso_WFD_miniSurveys/Metric/EM5_filtered_Nqso_counter_all_bands.ipynb) in the WFD and [this notebook](Nqso_WFD_miniSurveys/Metric/EM5_Nqso_counter_all_bands.ipynb) for the combination of the WFD and the mini surveys. 

For a given depth, the estimation is done by interpolating [these tables](Nqso_WFD_miniSurveys/General) through [this function](Nqso_WFD_miniSurveys/Fast_Nqso.py). The tables which are calculated using [this function](Nqso_WFD_miniSurveys/Nqso_v2_MiLim.py), which uses [this integration scheme](Nqso_WFD_miniSurveys/Phi_Obs_v3_MiLim.py). See the comments in the python scripts for more details. 
