# NoStRa: Non-Stationary Random field modeling

A doubly stochastic advection-diffusion-decay model (DSADM) on the circle

Michael Tsyrulnikov and Alexander Rakitko

The repository contains the R code of a doubly stochastic advection-diffusion-decay model (DSADM) defined on the 1D spatial domain (the circle). The DSADM is intended to be used in testing and developing data assimilation methodologies as well as in a broader context of non-stationary spatio-temporal random field modeling. The model is hierarchical: it is a stochastic (forced by the white noise) linear partial differential equation whose coefficients are transformed spatio-temporal random fields by themselves satisfying their own stochastic partial differential equations with constant coefficients. The solution to the model is a conditionally Gaussian random field that has complex spatio-temporal covariances with the tunable degree of non-stationarity. 

The repository also contains the R code of several hybrid ensemble filters (EnKF, EnVar, HBEF by Tsyrulnikov and Rakitko, Physica D, 2017, and the hybrid-HBEF or HHBEF that combines all the above filters) in an environment with the DSADM as the ``model of truth'', and R scripts used to compute the model’s and the filters’ spatio-temporal statistics.

Reference:
Michael Tsyrulnikov and Alexander Rakitko. Impact of non-stationarity of prior covariances on hybrid ensemble filters: a study with a doubly stochastic advection-diffusion-decay model of truth. To be submitted to Quart. J. Roy. Meteorol. Soc., September 2018.

