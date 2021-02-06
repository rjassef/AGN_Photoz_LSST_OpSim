Changes with respect to discussion in 2021-01-27 telecon
--------------------------------------------------------

* Loaded all the codes and jupyter notebooks into a [public github repository](https://github.com/rjassef/AGN_Photoz_LSST_OpSim). This repository will be regularly updated.  

* Improved the readability of the plots using as a guideline the plots and codes of Gordon Richard’s group for their RNAAS paper. Still some are hard to read, so will look into further improvements.

* Added FBS 1.5 Opsims, so that now all plots use FBS 1.5 and 1.6. FBS 1.7 needs to be added.  

* Implemented two physical models for the quasar colors instead of the broken power-law: the AGN template from Assef et al. (2010), and the semi-empirical colors from Matthew Temple’s work. (Maybe Matthew can describe them in depth in the next telecon?)

* Color excess’ are shown now for filter pairs independently instead of as a mean between all consecutive filter pairs.

* Color excess’ are shown now in two versions: 
  1. Area histograms for all OpSims and quasar color models at z=2, and
  2. Median of the distribution for each OpSim as a function of redshift for all quasar color models.

* Fixed a major bug with the color depth calculations. The calculations are based on the ExgalM5 metric that is part of MAF which takes the filter as an input. After some more testing, it became clear that the filter the metric takes as input is only used for applying the reddening, but is not applied to the data slices. The new codes now apply the filter selection to the data slices as well.

* Following a question asked by JT, I removed for the “mean time between bands metric” the requirement that only observations within a season are considered. This is not necessary, and considering the seasons could have two potential issues (granted both are extreme cases that we are very unlikely to find):

  1. that one band is taken right at the end of the season and the next band is taken right at the beginning of the next season, and
  2. that only one of the bands is taken during a season.

As expected, though, results are not significantly affected.
