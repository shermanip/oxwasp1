# oxwasp1
# (OxWaSP Mini-Project 1)

##Inside-Out: Characterisation of CT noise in projection and image space with applications to 3D printing
X-ray computed tomography can be used to do quality control on 3D printed samples. However there are sources of error in the 3D printing, how the photons behave and in the X-ray detector. This project aims to find a relationship between the sample mean and sample variance grey values in images obtained from the X-ray detector by fitting linear regressions. In addition, latent variable models such as principle component analysis, factor analysis and the compound Poisson were attempted to be fitted to find sources of variance.

This repository contains MATLAB code for this project. *For the code to run, the location of the dataset must be modified.*

##/latex/
Contains files for the LaTeX report.

##/MATLAB/initial/histogram_greyValues.m
Script for obtaining a histogram of grey values of all pixels in the dataset. (Figure 3.2 p. 11)

##MATLAB/initial/time_series.m
Script for obtaining a time series and sample autocorrelation of the sample mean and standard error grey value for each image. (Figures 3.5 and 3.6 p. 14)

##MATLAB/initial/normalgof.m
Script for conduction $\chi^2$ goodness of fit test for fitting the Normal distribution to the grey values for each pixel. (Figure 3.4 p. 13)
