# DepInfeR-GP: Bayesian inference for identifying tumour-specific cancer dependencies through integration of ex-vivo drug response assays and drug-protein profiling
This repository hosts the implementation of DepInfeR-GP [link](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05682-0), a Bayesian model designed to identify and estimate specific molecular dependencies of individual cancers from their ex-vivo drug sensitivity profiles, and scripts that rerpocude the numerical simulations and figures in the paper. 

<p align="center"><img src="https://github.com/hwxing3259/depinfer_gp/blob/main/depinfer_linear_graphical-gp.drawio.pdf" alt="gears" width="900px" /></p>

my_helper_funcs.R --> Helper functions for e.g. MCMC sampler and cross validation 

CV_GDSC1.R, CV_beatAML.R, CV_EMBL.R  --> Applying the proposed model to three real datasets 

examles_and_figures_1,2.R --> generate figures presented in the paper
