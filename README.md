# DepInfeR-GP: Bayesian inference for identifying tumour-specific cancer dependencies through integration of ex-vivo drug response assays and drug-protein profiling
This repository hosts the implementation of DepInfeR-GP [(link)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05682-0), a Bayesian model designed to identify and estimate specific molecular dependencies of individual cancers from their ex-vivo drug sensitivity profiles, and scripts that rerpocude the numerical simulations and figures in the paper. 

<p align="center"><img src="https://github.com/hwxing3259/depinfer_gp/blob/main/depinfer_gp_graphical.png" alt="depinfer-gp" width="900px" /></p>

## Core API Interface
Here we use the BeatAML dataset from [Batzilla et al 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9436053/) as an example
```
source("depinfer_gp.R")
set.seed(31415926)
# load the beatAML dataset from Batzilla et al 2022.
load("beatAML_dataset.RData")

N <- 200 # number of MCMC steps, each step consists of thin times Gibbs update
burn_in<-50 # discard the first burn_in number of samples, results would have lenth = N-burn_in
thin<-2 # number of Gibbs updates in each MCMC step

GP_X <- tarMat_BeatAML
GP_X[GP_X==GP_X[1,1]] <- NA
GP_Y <- viabMat_BeatAML_raw_log

MCMC_CV <- CV_nu_par_0(N=N, burn_in=burn_in, thin=thin, X=GP_X, Y=GP_Y, 
                       vec_nu_1=seq(0.01, 0.3, length.out=5), vec_nu_2=seq(0.01, 0.3, length.out=6),
                       inclusion_prob=.1, back_fit_iter=20, k_fold=3, ncor=10)
```

## Reproducing numerical examples
Codes for reproducing the BeatAML ([Batzilla et al 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9436053/)) example: [Link](https://github.com/hwxing3259/depinfer_gp/blob/main/CV_BeatAML.R)
Codes for reproducing the EMBL ([Batzilla et al 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9436053/)) example: [Link](https://github.com/hwxing3259/depinfer_gp/blob/main/CV_EMBL.R)
Codes for reproducing the GDSC1 ([Batzilla et al 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9436053/)) example: [Link](https://github.com/hwxing3259/depinfer_gp/blob/main/CV_GDSC1.R)
Codes for reproducing figures in our paper: [Link1](https://github.com/hwxing3259), [Link2](https://github.com/hwxing3259)

## Cite us
```
@article{xing2024bayesian,
  title={Bayesian inference for identifying tumour-specific cancer dependencies through integration of ex-vivo drug response assays and drug-protein profiling},
  author={Xing, Hanwen and Yau, Christopher},
  journal={BMC bioinformatics},
  volume={25},
  number={1},
  pages={104},
  year={2024},
  publisher={Springer}
}
```

my_helper_funcs.R --> Helper functions for e.g. MCMC sampler and cross validation 

CV_GDSC1.R, CV_beatAML.R, CV_EMBL.R  --> Applying the proposed model to three real datasets 

examles_and_figures_1,2.R --> generate figures presented in the paper
