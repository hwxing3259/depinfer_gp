library(ggplot2)
library(gridExtra)
library(DepInfeR)
library(missForest)
library(tidyverse)
library(glmnet)
library(BiocParallel)
library(matrixStats)
library(randomForestSRC)
source('my_helper_funcs.R')
library(ggtext)



# reproducing the figs in main body


mean_fps <- function(protein, X, GP_est){
  protein_id = which(colnames(X)==protein)
  protein_mean = GP_est$curve_mean[,,protein_id,]
  dtfrm <- apply(X=protein_mean, MARGIN = c(2,3), FUN = mean, na.rm=TRUE)
  dtfrm <- as.data.frame(cbind(as.vector(t(dtfrm)), rep(seq(0, 1, length.out=10), times=dim(dtfrm)[1])))
  colnames(dtfrm) <- c('Contribution', 'Affinity')
  dtfrm$Sample <- rep(sapply(1:(dim(dtfrm)[1]/10), 
                             FUN=function(x){paste0('cell', x)}), each=10)
  dtfrm$Kinase <- protein
  return(dtfrm)
  
  # plot(seq(0, 1, length.out=10), colMeans(protein_mean[,1,], na.rm=TRUE), type='l', ylim=ylim, 
  #      ylab='Function value', xlab = 'Affinity score')
  # title(paste0('beatAML, f_ps with protein ', protein, ' of all sample cells'))
  # for (i in 2:dim(Y_new)[2]){
  #   lines(seq(0, 1, length.out=10), colMeans(protein_mean[,i,], na.rm=TRUE), type='l')
  # }
}



# beat_AML

# reproducing the figs in main body
load('beatAML_dataset.RData')
load('beatAML_CV.RData')

beatAML_X_original <- tarMat_BeatAML
beatAML_X_new <- tarMat_BeatAML
beatAML_X_nan <- beatAML_X_new[1,1]
beatAML_X_new[beatAML_X_new == beatAML_X_nan] <- NaN
beatAML_Y_new <- viabMat_BeatAML_raw_log
beatAML_Y_original <- viabMat_BeatAML

# run the spike-and-slab GP model, save the results
if (!file.exists('beatAML_GP0_inference.RData')){
  # visualizing GP curves
  set.seed(314159)
  GP_0_z <- my_gibbs_sampler_gp_0(N=120, burn_in=20, thin=1, X=beatAML_X_new, Y=beatAML_Y_new, 
                                  inclusion_prob=0.1, 
                                  nu=CV_gp_0_log$c0[[which.min(CV_gp_0_log$mse)]], 
                                  update_Z=TRUE)
  
  GP_0_para <- GP_para_0(burn_in=20, X=beatAML_X_new, Y=beatAML_Y_new, Z=GP_0_z$Z, 
                         sigma2=GP_0_z$sigma2, gamma2=GP_0_z$gamma2, 
                         nu=CV_gp_0_log$c0[[which.min(CV_gp_0_log$mse)]])
  save(GP_0_z, GP_0_para, file = 'beatAML_GP0_inference.RData')
} else{
  load('beatAML_GP0_inference.RData')
}


beatAML_CDK6 <- mean_fps('CDK6', beatAML_X_new, GP_0_para)
beatAML_PDXK <- mean_fps('PDXK', beatAML_X_new, GP_0_para)
beatAML_AKT1 <- mean_fps('AKT1', beatAML_X_new, GP_0_para)
beatAML_CDK17 <- mean_fps('CDK17', beatAML_X_new, GP_0_para)



set.seed(314159)
if (!file.exists('DepInfeR_beatAML_res.RData')){
  beatAML_DepInfeR_yhat <- MultiLasso(TargetMatrix=beatAML_X_original, ResponseMatrix=beatAML_Y_original, repeats = 100, BPPARAM = bpparam())
} else {load('DepInfeR_beatAML_res.RData')}
beatAML_DepInfeR_Y_hat <- t(t(beatAML_X_original %*% beatAML_DepInfeR_yhat$coefMat) + beatAML_DepInfeR_yhat$intercept)
beatAML_DepInfeR_res <- beatAML_Y_original - t(t(beatAML_X_original %*% beatAML_DepInfeR_yhat$coefMat) + beatAML_DepInfeR_yhat$intercept)

beatAML_GP0_res <- as.vector(beatAML_Y_new - apply(GP_0_para$y_hat, c(1,2), mean))
beatAML_GP0_Y_hat <- as.vector(apply(GP_0_para$y_hat, c(1,2), mean))
beatAML_obs_new <- as.vector(beatAML_Y_new)
beatAML_obs_old <- as.vector(beatAML_Y_original)


beatAML_fit <- data.frame(cbind(as.vector(beatAML_DepInfeR_Y_hat), 
                                as.vector(beatAML_DepInfeR_res), 
                                beatAML_GP0_res, beatAML_GP0_Y_hat, beatAML_obs_new, beatAML_obs_old))
colnames(beatAML_fit) <- c('DepInfeR_yhat', 'DepInfeR_res', 'GP0_res', 'GP0_yhat', 'Yobs_new', 'Yobs_old')

beatAML_fit_fig <- list(ggplot(data = beatAML_fit) + geom_point(aes(x=Yobs_old, y=DepInfeR_yhat)) + 
                          geom_abline(intercept=0, slope=1) + xlab('Observed values') + ylab('Fitted values') + ggtitle('beatAML, DepInfeR'),
                        ggplot(data = beatAML_fit) + geom_point(aes(x=Yobs_new, y=GP0_yhat)) + 
                          geom_abline(intercept=0, slope=1) + xlab('Observed values') + ylab('Fitted values') + ggtitle('beatAML, GP'),
                        ggplot(data = beatAML_fit) + geom_qq(aes(sample=DepInfeR_res)) + geom_qq_line(aes(sample=DepInfeR_res))+ 
                          ggtitle('beatAML, DepInfeR') + xlab('Theoretical quantiles') + ylab('Sample quantiles'),
                        ggplot(data = beatAML_fit) + geom_qq(aes(sample=GP0_res)) + geom_qq_line(aes(sample=GP0_res))+ 
                          ggtitle('beatAML, GP') + xlab('Theoretical quantiles') + ylab('Sample quantiles')
)



# EMBL:


load('EMBL_dataset.RData')
load('EMBL_CV.RData')

EMBL_X_new <- tarMat_EMBL
EMBL_X_nan <- EMBL_X_new[1,1]
EMBL_X_new[EMBL_X_new == EMBL_X_nan] <- NaN
EMBL_Y_new <- viabMat_EMBL_log
EMBL_X_original <- tarMat_EMBL
EMBL_Y_original <- viabMat_EMBL


if (!file.exists('EMBL_GP0_inference.RData')){
  set.seed(314159)
  GP_0_z <- my_gibbs_sampler_gp_0(N=120, burn_in=20, thin=1, X=EMBL_X_new, Y=EMBL_Y_new, 
                                  inclusion_prob=0.1, 
                                  nu=CV_gp_0_log$c0[[which.min(CV_gp_0_log$mse)]], 
                                  update_Z=TRUE)
  
  GP_0_para <- GP_para_0(burn_in=20, X=EMBL_X_new, Y=EMBL_Y_new, Z=GP_0_z$Z, 
                         sigma2=GP_0_z$sigma2, gamma2=GP_0_z$gamma2, 
                         nu=CV_gp_0_log$c0[[which.min(CV_gp_0_log$mse)]])
  
  save(GP_0_z, GP_0_para, file = 'EMBL_GP0_inference.RData')
} else{
  load('EMBL_GP0_inference.RData')
}


EMBL_ACAD11 <- mean_fps('ACAD11', EMBL_X_new, GP_0_para)
EMBL_STK26 <- mean_fps('STK26', EMBL_X_new, GP_0_para)
EMBL_CCNT1 <- mean_fps('CCNT1', EMBL_X_new, GP_0_para)



set.seed(314159)
if (!file.exists('DepInfeR_EMBL_res.RData')){
  EMBL_DepInfeR_yhat <- MultiLasso(TargetMatrix=EMBL_X_original, ResponseMatrix=EMBL_Y_original, repeats = 100, BPPARAM = bpparam())
} else {load('DepInfeR_EMBL_res.RData')}
EMBL_DepInfeR_Y_hat <- t(t(EMBL_X_original %*% EMBL_DepInfeR_yhat$coefMat) + EMBL_DepInfeR_yhat$intercept)
EMBL_DepInfeR_res <- EMBL_Y_original - t(t(EMBL_X_original %*% EMBL_DepInfeR_yhat$coefMat) + EMBL_DepInfeR_yhat$intercept)

EMBL_GP0_res <- as.vector(EMBL_Y_new - apply(GP_0_para$y_hat, c(1,2), mean))
EMBL_GP0_Y_hat <- as.vector(apply(GP_0_para$y_hat, c(1,2), mean))
EMBL_obs_new <- as.vector(EMBL_Y_new)
EMBL_obs_old <- as.vector(EMBL_Y_original)


EMBL_fit <- data.frame(cbind(as.vector(EMBL_DepInfeR_Y_hat), 
                             as.vector(EMBL_DepInfeR_res), 
                             EMBL_GP0_res, EMBL_GP0_Y_hat, EMBL_obs_new, EMBL_obs_old))
colnames(EMBL_fit) <- c('DepInfeR_yhat', 'DepInfeR_res', 'GP0_res', 'GP0_yhat', 'Yobs_new', 'Yobs_old')

EMBL_fit_fig <- list(ggplot(data = EMBL_fit) + geom_point(aes(x=Yobs_old, y=DepInfeR_yhat)) + 
                       geom_abline(intercept=0, slope=1) + xlab('Observed values') + ylab('Fitted values') + ggtitle('EMBL, DepInfeR'),
                     ggplot(data = EMBL_fit) + geom_point(aes(x=Yobs_new, y=GP0_yhat)) + 
                       geom_abline(intercept=0, slope=1) + xlab('Observed values') + ylab('Fitted values') + ggtitle('EMBL, GP'),
                     ggplot(data = EMBL_fit) + geom_qq(aes(sample=DepInfeR_res)) + geom_qq_line(aes(sample=DepInfeR_res))+ 
                       ggtitle('EMBL, DepInfeR') + xlab('Theoretical quantiles') + ylab('Sample quantiles'),
                     ggplot(data = EMBL_fit) + geom_qq(aes(sample=GP0_res)) + geom_qq_line(aes(sample=GP0_res))+ 
                       ggtitle('EMBL, GP') + xlab('Theoretical quantiles') + ylab('Sample quantiles')
)





# GDSC1:
load('GDSC1_GP0_inference.RData')
load('GDSC1_dataset.RData')

GDSC1_X_new <- tarMat_GDSC
GDSC1_X_nan <- GDSC1_X_new[1,1]
GDSC1_X_new[GDSC1_X_new == GDSC1_X_nan] <- NaN
GDSC1_Y_new <- viabMat_GDSC_logit

GDSC1_X_original <- tarMat_GDSC
GDSC1_Y_original <- viabMat_GDSC


if (!file.exists('GDSC1_GP0_inference.RData')){
  GP_0_z <- my_gibbs_sampler_gp_0(N=120, burn_in=20, thin=1, X=GDSC1_X_new, Y=GDSC1_Y_new, 
                                  inclusion_prob=0.1, 
                                  nu=CV_gp_0_logit$c0[[which.min(CV_gp_0_logit$mse)]], 
                                  update_Z=TRUE)
  
  GP_0_para <- GP_para_0(burn_in=20, X=GDSC1_X_new, Y=GDSC1_Y_new, Z=GP_0_z$Z, 
                         sigma2=GP_0_z$sigma2, gamma2=GP_0_z$gamma2, 
                         nu=CV_gp_0_logit$c0[[which.min(CV_gp_0_logit$mse)]])
  
  save(GP_0_z, GP_0_para, file = 'GDSC1_GP0_inference.RData')
} else{
  load('GDSC1_GP0_inference.RData')
}


GDSC1_CCNK <- mean_fps('CCNK', GDSC1_X_new, GP_0_para)
GDSC1_DYRK1A <- mean_fps('DYRK1A', GDSC1_X_new, GP_0_para)
GDSC1_AFF4 <- mean_fps('AFF4', GDSC1_X_new, GP_0_para)
GDSC1_MAP2K2 <- mean_fps('MAP2K2', GDSC1_X_new, GP_0_para)


plot_data <- list(GDSC1_CCNK, GDSC1_DYRK1A, GDSC1_AFF4, GDSC1_MAP2K2,
                  beatAML_CDK6, beatAML_PDXK, beatAML_AKT1, beatAML_CDK17,
                  EMBL_ACAD11, EMBL_STK26, EMBL_CCNT1)
kinases <- c("CCNK", "DYRK1A", "AFF4", "MAP2K2", 
             "CDK6", "PDXK", "AKT1", "CDK17",
             "ACAD11", "STK26", "CCNT1")
dataset_name <- c(rep("GDSC1", 4), rep("beatAML", 4), rep('EMBL', 3))
color_vec <- rep(c('black', 'black', 'navyblue', 'navyblue'), 4)

p <- lapply(1:11, function(i){ggplot(data = plot_data[[i]]) + 
    geom_line(aes(x=Affinity, y=Contribution, group=Sample), alpha=0.4, color = color_vec[i])+
    ggtitle(paste0(dataset_name[i], ', kinase=', kinases[i]))})

grid.arrange(grobs=p, ncol=4)

set.seed(314159)
if (!file.exists('DepInfeR_GDSC1_res.RData')){
  GDSC1_DepInfeR_yhat <- MultiLasso(TargetMatrix=GDSC1_X_original, ResponseMatrix=GDSC1_Y_original, repeats = 100, BPPARAM = bpparam())
} else {load('DepInfeR_GDSC1_res.RData')}
GDSC1_DepInfeR_Y_hat <- t(t(GDSC1_X_original %*% GDSC1_DepInfeR_yhat$coefMat) + GDSC1_DepInfeR_yhat$intercept)
GDSC1_DepInfeR_res <- GDSC1_Y_original - t(t(GDSC1_X_original %*% GDSC1_DepInfeR_yhat$coefMat) + GDSC1_DepInfeR_yhat$intercept)

GDSC1_GP0_res <- as.vector(GDSC1_Y_new - apply(GP_0_para$y_hat, c(1,2), mean))
GDSC1_GP0_Y_hat <- as.vector(apply(GP_0_para$y_hat, c(1,2), mean))
GDSC1_obs_new <- as.vector(GDSC1_Y_new)
GDSC1_obs_old <- as.vector(GDSC1_Y_original)


GDSC1_fit <- data.frame(cbind(as.vector(GDSC1_DepInfeR_Y_hat), 
                              as.vector(GDSC1_DepInfeR_res), 
                              GDSC1_GP0_res, GDSC1_GP0_Y_hat, GDSC1_obs_new, GDSC1_obs_old))
colnames(GDSC1_fit) <- c('DepInfeR_yhat', 'DepInfeR_res', 'GP0_res', 'GP0_yhat', 'Yobs_new', 'Yobs_old')

GDSC1_fit_fig <- list(ggplot(data = GDSC1_fit) + geom_point(aes(x=Yobs_old, y=DepInfeR_yhat)) + 
                        geom_abline(intercept=0, slope=1) + xlab('Observed values') + ylab('Fitted values') + ggtitle('GDSC1, DepInfeR'),
                      ggplot(data = subset(GDSC1_fit, !is.na(Yobs_new))) + geom_point(aes(x=Yobs_new, y=GP0_yhat)) + 
                        geom_abline(intercept=0, slope=1) + xlab('Observed values') + ylab('Fitted values') + ggtitle('GDSC1, GP'),
                      ggplot(data = GDSC1_fit) + geom_qq(aes(sample=DepInfeR_res)) + geom_qq_line(aes(sample=DepInfeR_res))+ 
                        ggtitle('GDSC1, DepInfeR') + xlab('Theoretical quantiles') + ylab('Sample quantiles'),
                      ggplot(data = subset(GDSC1_fit, !is.na(GP0_res))) + geom_qq(aes(sample=GP0_res)) + geom_qq_line(aes(sample=GP0_res))+ 
                        ggtitle('GDSC1, GP') + xlab('Theoretical quantiles') + ylab('Sample quantiles')
)




grid.arrange(grobs=c(GDSC1_fit_fig, beatAML_fit_fig, EMBL_fit_fig), ncol=4)





