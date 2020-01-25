*****************************
Power spectral response method
*****************************

Introduction
============
The scripts in this repo establish the power spectral response method (PSRESP) for the analysis of the power spectral density (PSD) of light curves.
It is a Monte Carlo approach that takes the sampling directly into account and reveals the underlying PSD.
As model for the PSD, an unbroken power law is assumed.
For each model parameter (i.e. slope of the power law), a success fraction (SUF) is calculated that defines the inverse rejection level for this model.
The slope of the PSD model for the light curve is estimated as the mean of all slopes providing a significant SUF.
The corresponding error is given by the full width at half maximum (FWHM) of the SUF distribution.

Getting Started
===============
Input
-----
`psresp` takes a light curve in format time, flux and flux error.
For the PSD model, the trial slopes have to be forwarded via `slopes`.
The PSRESP method bins the light curve and the periodogram as defined by `dt` and `df`.
To determine the significant SUF, the percentile for the SUF distribution, `percentile`, needs to be given.
The number of simulations can be defined by `number_simulations`, it is 100 by default.
Additionally, the oversampling of the artificial light curves can be defined by `oversampling`.
`~plot_psresp` takes the output of `~gammapy.time.psresp` as input.

Output
------
`psresp` returns the mean slope and its error,
the success fraction as a function of model parameters (`slopes`, `dt`, `df`),
parameters `dt` and `df` providing a significant SUF
and the statistics used to calculate the mean slope and its error.

Example
=======
An example for estimating the slope of the power spectral density of an AGN light curve is shown in the figure below.
The light curve is from the gamma-ray source Mrk 421 observed with MAGIC at energies above 0.3 TeV in 2009 [1]_.
The PSRESP reveals the slope of the underlying PSD model to :math:`(1.5 \pm 1.3)` days
in agreement with :math:`(1.6 \pm 0.9)` [1]_.
Please note that in [1]_, only the half width at half maximum is used as an estimate for the slope error,
resulting in a less restricting estimate.

.. gp-extra-image:: time/example_PSRESP.png
   :width: 100 %

.. [1] MAGIC collaboration, The 2009 multiwavelength campaign on Mrk 421: Variability and correlation studies,
   `Link <https://arxiv.org/pdf/1502.02650.pdf>`_
