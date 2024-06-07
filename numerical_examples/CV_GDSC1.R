library(DepInfeR)
library(missForest)
library(tidyverse)
library(glmnet)
library(BiocParallel)
library(matrixStats)
library(randomForestSRC)
source('my_helper_funcs.R')
library(ggplot2)
library(ggtext)
library(doParallel)
library(foreach)

set.seed(31415926)
load('GDSC1_dataset.RData')
X_original <- tarMat_GDSC
X_new <- tarMat_GDSC
X_nan <- X_new[1,1]
X_new[X_new == X_nan] <- NaN
X_zeroed = tarMat_GDSC
X_zeroed[X_zeroed == X_nan] <- 0.

Y_original <- viabMat_GDSC
Y_new <- viabMat_GDSC_logit

P = dim(X_original)[2]


N<-50 # number of MCMC steps, each step consists of thin times Gibbs update
burn_in<-20 # discard the first burn_in number of samples, results would have lenth = N-burn_in
thin<-2 # number of Gibbs updates in each MCMC step
a0<-1. # controlling the prior on sigma2
b0<-1. # controlling the prior on sigma2
inclusion_prob<-0.1 # controlling the prior on z

vec_c0 <- seq(0.1, 2.5, length.out=30)



start_time = Sys.time()
# run this on the original dataset, compare model fit
MCMC_CV <- CV_par(N, burn_in, thin, X_original, Y_original, vec_c0, a0, b0, inclusion_prob, ncore=10, k_fold=3)
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_linear_original.RData')
# run this on the new Y dataset, compare the selected proteins with the original
MCMC_CV <- CV_par(N, burn_in, thin, X_original, Y_new, vec_c0, a0, b0, inclusion_prob, ncore=10, k_fold=3)
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_linear_logit_Y_with_NA.RData')
# run this on the new X dataset (set NAN=0), see if there is any difference
MCMC_CV <- CV_par(N, burn_in, thin, X_zeroed, Y_new, vec_c0, a0, b0, inclusion_prob, ncore=10, k_fold=3)
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_linear_X_zeroed_Y_logit.RData')




set.seed(31415926)
# full set of parameter
# run this on the original dataset, compare model fit
MCMC_CV <- CV_par(N, burn_in, thin, X_original, Y_original, vec_c0, a0, b0, inclusion_prob, ncore=10, k_fold=3, updated_z = FALSE, old_z = rep(1, P))
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_linear_original_all.RData')
# run this on the new Y dataset, compare the selected proteins with the original
MCMC_CV <- CV_par(N, burn_in, thin, X_original, Y_new, vec_c0, a0, b0, inclusion_prob, ncore=10, k_fold=3, updated_z = FALSE, old_z = rep(1, P))
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_linear_logit_Y_with_NA_all.RData')
# run this on the new X dataset (set NAN=0), see if there is any difference
MCMC_CV <- CV_par(N, burn_in, thin, X_zeroed, Y_new, vec_c0, a0, b0, inclusion_prob, ncore=10, k_fold=3, updated_z = FALSE, old_z = rep(1, P))
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_linear_X_zeroed_Y_logit_all.RData')
print(Sys.time()-start_time)



N <- 120 # number of MCMC steps, each step consists of thin times Gibbs update
burn_in<-20 # discard the first burn_in number of samples, results would have lenth = N-burn_in
thin<-2 # number of Gibbs updates in each MCMC step

GP_X <- tarMat_GDSC
GP_X[GP_X==GP_X[1,1]] <- NA
head(GP_X)
GP_Y <- viabMat_GDSC_logit
head(GP_Y)


set.seed(31415926)
start_time = Sys.time()
MCMC_CV <- CV_nu_par(N=N, burn_in=burn_in, thin=thin, X=GP_X, Y=GP_Y, 
                     vec_nu_1=seq(0.01, 0.1, length.out=5), vec_nu_2=c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3), inclusion_prob=.1, back_fit_iter=20, k_fold=3, ncor=10)
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_GP_X_NA_Y_logit.RData')
print(Sys.time()-start_time)


start_time = Sys.time()
MCMC_CV <- CV_nu_par(N=N, burn_in=burn_in, thin=thin, X=GP_X, Y=Y_original, 
                     vec_nu_1=seq(0.01, 0.1, length.out=5), vec_nu_2=c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3), inclusion_prob=.1, back_fit_iter=20, k_fold=3, ncor=10)
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_GP_X_NA_Y_original.RData')
print(Sys.time()-start_time)



start_time = Sys.time()
MCMC_CV <- CV_nu_par_0(N=N, burn_in=burn_in, thin=thin, X=GP_X, Y=GP_Y, 
                       vec_nu_1=seq(0.01, 0.3, length.out=5), vec_nu_2=seq(0.01, 0.3, length.out=6), inclusion_prob=.1, back_fit_iter=20, k_fold=3, ncor=10)
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_GP_0_X_NA_Y_logit.RData')
print(Sys.time()-start_time)


start_time = Sys.time()
MCMC_CV <- CV_nu_par_0(N=N, burn_in=burn_in, thin=thin, X=GP_X, Y=Y_original, 
                       vec_nu_1=seq(0.01, 0.3, length.out=6), vec_nu_2=seq(0.01, 0.3, length.out=6), inclusion_prob=.1, back_fit_iter=20, k_fold=3, ncor=12)
save(MCMC_CV, file = './GDSC1_cv_res/GDSC1_GP_0_X_NA_Y_original.RData')
print(Sys.time()-start_time)

# collected results are saved in GDSC1_CV.RData
