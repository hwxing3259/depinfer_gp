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
load('GDSC1_dataset.RData')
load('GDSC1_CV.Rdata')

X_original <- tarMat_GDSC
X_new <- tarMat_GDSC
X_nan <- X_new[1,1]
X_new[X_new == X_nan] <- NaN

Y_original <- viabMat_GDSC
Y_new <- viabMat_GDSC_logit



k_fold <- 3
sample_size <- dim(Y_original)[1]
random_partition <- sample(1:sample_size, size=sample_size, replace = FALSE)
group_id <- split(random_partition, sort(1:sample_size %% k_fold))
response_mat <- rep(0, k_fold)
for (i in 1:k_fold){
  test_lasso <- MultiLasso(TargetMatrix=X_original[-group_id[[i]], ], ResponseMatrix=Y_original[-group_id[[i]], ], repeats = 100, BPPARAM = bpparam())
  response_mat[i] <- mean((Y_original[group_id[[i]],] - t(t(X_original[group_id[[i]], ] %*% test_lasso$coefMat) + test_lasso$intercept)
  )^2)
}

DepInfeR_MSE = mean(response_mat)
DepInfeR_nMSE = mean(response_mat)/var(as.vector(Y_original))

# how about a full non-linear treatment?
rf_test_100 <- rf_test(100, X_original, Y_original)
rf_test_200 <- rf_test(200, X_original, Y_original)
rf_test_300 <- rf_test(300, X_original, Y_original)
rf_test_400 <- rf_test(400, X_original, Y_original)
rf_test_500 <- rf_test(500, X_original, Y_original)

rf_MSE <- sapply(list(rf_test_100, rf_test_200, rf_test_300, rf_test_400, rf_test_500),
                 FUN = function(X){mean((X - Y_original)^2)})
rf_nMSE <- sapply(list(rf_test_100, rf_test_200, rf_test_300, rf_test_400, rf_test_500),
                  FUN = function(X){mean((X - Y_original)^2)}/var(as.vector(Y_original)))

# GDSC1 set
GDSC1_original <- c('RET', 'AURKA', 'AFF4', 'FLT3', 'AAK1', 'JAK1', 'EGFR', 'MAP2K2', 'BCR/ABL')

plot_data_GDSC1 <- data.frame(cbind(c(CV_gp_0_logit$mse/var(as.vector(Y_new), na.rm=TRUE), 
                                      CV_gp_0_original$mse/var(as.vector(Y_original), na.rm=TRUE)),
                                    c(int_over_union(CV_gp_0_logit$selected, GDSC1_original), 
                                      int_over_union(CV_gp_0_original$selected, GDSC1_original))))
colnames(plot_data_GDSC1) <- c('nMSE', 'IoU')
plot_data_GDSC1$Dataset <- c(rep('logit-transformed, incomplete', 30), rep('original, imputed', 30))
temp_data_1 <- data.frame(val=c(DepInfeR_nMSE, min(rf_nMSE)), Competing.model=c('DepInfeR','MultiRF'))

res_GDSC1 <- ggplot(data = plot_data_GDSC1, aes(x=nMSE, y=IoU)) + geom_point(aes(color=Dataset)) + 
  xlab('3-fold CV normalized MSE') + ylab('Similarity between sets of selected kinases') + 
  ggtitle('GDSC1 dataset') + geom_vline(aes(xintercept=val, linetype=Competing.model), temp_data_1)





protein_select <- function(cv_res){
  p <- cv_res$selected[[which.min(cv_res$mse)]]
  return(sort(p))
}

protein_list <- list()
protein_list[[1]] <- protein_select(CV_gp_0_logit)
protein_list[[2]] <- protein_select(CV_gp_0_original)
protein_list[[3]] <- GDSC1_original


protein_counter <- list()
for (i in 1:3){
  for (j in protein_list[[i]]){
    if (is.null(protein_counter[[j]])){
      protein_counter[[j]] <- 1
    }
    else{
      protein_counter[[j]] <- protein_counter[[j]] + 1
    }
  }
}

protein_res <- sort(unlist(protein_counter), decreasing = TRUE)
common_proteins <- sort(names(protein_res[protein_res > 1]))


all_protein <- unique(sort(unlist(protein_list)))
id_mtrx <- matrix(NA, ncol=length(all_protein), nrow=length(protein_list))


my_data <- as.matrix(expand_grid(1:length(protein_list), 1:length(all_protein)))
idct <- rep(NA, length(protein_list)*length(all_protein))

for (i in 1:length(idct)){
  if (all_protein[my_data[i,2]] %in% common_proteins){
    idct[i] <- 2* (all_protein[my_data[i,2]] %in% protein_list[[my_data[i,1]]])
  }
  else{
    idct[i] <- all_protein[my_data[i,2]] %in% protein_list[[my_data[i,1]]]
  }
}

my_face = all_protein
for (i in 1:length(all_protein)){
  if (all_protein[i] %in% common_proteins){
    my_face[i] = 'bold'
  }
  else{
    my_face[i] = 'plain'
  }
}


my_data <- as.data.frame(cbind(my_data, idct))
colnames(my_data) <- c('Group_id', 'Protein_id', 'Label')
my_data$Selection <- as.factor(sapply(my_data$Label, FUN=function(i){c('Not selected', 'Selected', 'Consistent')[i+1]}))

GDSC1_selected <- ggplot(my_data, aes(Protein_id, Group_id, fill= Selection)) + geom_tile(color='white') +
  # scale_color_continuous(breaks=c(0,1,2), labels=c('Not selected', 'Selected', 'Consistent')) + 
  scale_y_continuous(breaks = 1:3, labels = c('GP_0_logit', 'GP_0_original','DepInfeR')) + 
  scale_x_continuous(breaks = 1:length(all_protein), labels = all_protein)+ xlab("Kinases") + ylab("") + 
  # guides(fill = guide_colourbar(label = FALSE,
  #                               ticks = FALSE))+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, face=my_face),
        axis.line.y = element_blank(),
        axis.line.x = element_blank()) + ggtitle('Selected Kinases in GDSC1') + 
  scale_fill_manual(breaks = c('Not selected', 'Selected', 'Consistent'), values=c('#122b44', '#326a97', '#54b1f7'))



# Res with the linear model 
plot(CV_gp_logit$mse/var(as.vector(Y_new), na.rm=TRUE), int_over_union(CV_gp_logit$selected, GDSC1_original),
     xlab='3-fold CV normalized MSE', ylab='Similarity between sets of selected proteins', ylim=c(0, 0.75), xlim=c(0.7, 2.0), col='purple')
points(CV_gp_original$mse/var(as.vector(Y_original), na.rm=TRUE), int_over_union(CV_gp_original$selected, GDSC1_original), col='blue')
points(CV_gp_0_logit$mse/var(as.vector(Y_new), na.rm=TRUE), int_over_union(CV_gp_0_logit$selected, GDSC1_original), col='grey')
points(CV_gp_0_original$mse/var(as.vector(Y_original), na.rm=TRUE), int_over_union(CV_gp_0_original$selected, GDSC1_original), col='red')

abline(v=DepInfeR_nMSE, lty=1, col='black')
abline(v=min(rf_nMSE), lty=2, col='black')

points(CV_linear_original$rmse, int_over_union(CV_linear_original$selected, GDSC1_original), col='brown')
abline(v=min(CV_linear_original$rmse), lty=2, col='brown')
# points(CV_linear_original_all$rmse, int_over_union(CV_linear_original_all$selected, GDSC1_original), col='black')
points(CV_linear_logit$rmse, int_over_union(CV_linear_logit$selected, GDSC1_original), col='green')
abline(v=min(CV_linear_logit$rmse), lty=2, col='green')
points(CV_linear_zeroed_logit$rmse, int_over_union(CV_linear_zeroed_logit$selected, GDSC1_original), col='pink')
abline(v=min(CV_linear_zeroed_logit$rmse), lty=2, col='pink')
legend(1.43, 0.7, legend=c('GP with NaN in X and logit Y', 'GP with NaN in X and z-score Y', 
                           'GP0 with NaN in X and logit Y', 'GP0 with NaN in X and z-score Y',
                           "DepInfeR with the original dataset", "MultiRF with the original dataset", 
                           'linear with the original dataset',
                           'linear with the original dataset, all proteins',
                           'linear with original X + logit Y',
                           'linear with original X + logit Y, all proteins',
                           'linear with NaN in X + logit Y',
                           'linear with NaN in X + logit Y, all proteins'),
       col=c('purple', 'blue', 'grey', 'red', 'black', 'black', 'brown','brown', 'green','green', 'pink', 'pink'), 
       lty=c(NA, NA, NA, NA,1,2,NA,2,NA,2,NA,2), 
       pch=c(1,1,1,1,NA,NA, 1, NA, 1, NA, 1, NA), cex=0.8)
title('GDSC1 dataset')









load('beatAML_dataset.RData')
load('beatAML_CV.Rdata')
X_original <- tarMat_BeatAML
X_new <- tarMat_BeatAML
X_nan <- X_new[1,1]
X_new[X_new == X_nan] <- NaN

Y_original <- viabMat_BeatAML
Y_new <- viabMat_BeatAML_raw_log




set.seed(314159)

k_fold <- 3
sample_size <- dim(Y_original)[1]
random_partition <- sample(1:sample_size, size=sample_size, replace = FALSE)
group_id <- split(random_partition, sort(1:sample_size %% k_fold))
response_mat <- rep(0, k_fold)
for (i in 1:k_fold){
  test_lasso <- MultiLasso(TargetMatrix=X_original[-group_id[[i]], ], ResponseMatrix=Y_original[-group_id[[i]], ], repeats = 100, BPPARAM = bpparam())
  response_mat[i] <- mean((Y_original[group_id[[i]],] - t(t(X_original[group_id[[i]], ] %*% test_lasso$coefMat) + test_lasso$intercept)
  )^2)
}

DepInfeR_MSE = mean(response_mat)
DepInfeR_nMSE = mean(response_mat)/var(as.vector(Y_original))


# how about a full non-linear treatment?
rf_test_100 <- rf_test(100, X_original, Y_original)
rf_test_200 <- rf_test(200, X_original, Y_original)
rf_test_300 <- rf_test(300, X_original, Y_original)
rf_test_400 <- rf_test(400, X_original, Y_original)
rf_test_500 <- rf_test(500, X_original, Y_original)

rf_MSE <- sapply(list(rf_test_100, rf_test_200, rf_test_300, rf_test_400, rf_test_500),
                 FUN = function(X){mean((X - Y_original)^2)})
rf_nMSE <- sapply(list(rf_test_100, rf_test_200, rf_test_300, rf_test_400, rf_test_500),
                  FUN = function(X){mean((X - Y_original)^2)}/var(as.vector(Y_original)))


plot_data_beatAML <- data.frame(cbind(c(CV_gp_0_log$mse/var(as.vector(Y_new), na.rm=TRUE), 
                                        CV_gp_0_original$mse/var(as.vector(Y_original), na.rm=TRUE)),
                                      c(int_over_union(CV_gp_0_log$selected, beatAML_original), 
                                        int_over_union(CV_gp_0_original$selected, beatAML_original))))
colnames(plot_data_beatAML) <- c('nMSE', 'IoU')
plot_data_beatAML$Dataset <- c(rep('log-transformed, incomplete', 30), rep('original, imputed', 30))
temp_data_1 <- data.frame(val=c(DepInfeR_nMSE, min(rf_nMSE)), Competing.model=c('DepInfeR','MultiRF'))

res_beatAML <- ggplot(data = plot_data_beatAML, aes(x=nMSE, y=IoU)) + geom_point(aes(color=Dataset)) + 
  xlab('3-fold CV normalized MSE') + ylab('Similarity between sets of selected kinases') + 
  ggtitle('beatAML dataset') + geom_vline(aes(xintercept=val, linetype=Competing.model), temp_data_1)




protein_select <- function(cv_res){
  p <- cv_res$selected[[which.min(cv_res$mse)]]
  return(sort(p))
}

protein_list <- list()
protein_list[[1]] <- protein_select(CV_gp_0_log)
protein_list[[2]] <- protein_select(CV_gp_0_original)
protein_list[[3]] <- beatAML_original


protein_counter <- list()
for (i in 1:3){
  for (j in protein_list[[i]]){
    if (is.null(protein_counter[[j]])){
      protein_counter[[j]] <- 1
    }
    else{
      protein_counter[[j]] <- protein_counter[[j]] + 1
    }
  }
}

protein_res <- sort(unlist(protein_counter), decreasing = TRUE)
common_proteins <- sort(names(protein_res[protein_res > 1]))



all_protein <- unique(sort(unlist(protein_list)))
id_mtrx <- matrix(NA, ncol=length(all_protein), nrow=length(protein_list))


my_data <- as.matrix(expand_grid(1:length(protein_list), 1:length(all_protein)))
idct <- rep(NA, length(protein_list)*length(all_protein))

for (i in 1:length(idct)){
  if (all_protein[my_data[i,2]] %in% common_proteins){
    idct[i] <- 2* (all_protein[my_data[i,2]] %in% protein_list[[my_data[i,1]]])
  }
  else{
    idct[i] <- all_protein[my_data[i,2]] %in% protein_list[[my_data[i,1]]]
  }
}


my_face = all_protein
for (i in 1:length(all_protein)){
  if (all_protein[i] %in% common_proteins){
    my_face[i] = 'bold'
  }
  else{
    my_face[i] = 'plain'
  }
}

my_data <- as.data.frame(cbind(my_data, idct))
colnames(my_data) <- c('Group_id', 'Protein_id', 'Label')
my_data$Selection <- as.factor(sapply(my_data$Label, FUN=function(i){c('Not selected', 'Selected', 'Consistent')[i+1]}))

beatAML_selected <- ggplot(my_data, aes(Protein_id, Group_id, fill= Selection)) + geom_tile(color='white') +
  # scale_color_continuous(breaks=c(0,1,2), labels=c('Not selected', 'Selected', 'Consistent')) + 
  scale_y_continuous(breaks = 1:3, labels = c('GP_0_log', 'GP_0_original','DepInfeR')) + 
  scale_x_continuous(breaks = 1:length(all_protein), labels = all_protein)+ xlab("Kinases") + ylab("") + 
  # guides(fill = guide_colourbar(label = FALSE,
  #                               ticks = FALSE))+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, face=my_face),
        axis.line.y = element_blank(),
        axis.line.x = element_blank()) + ggtitle('Selected Kinases in beatAML') + 
  scale_fill_manual(breaks = c('Not selected', 'Selected', 'Consistent'), values=c('#122b44', '#326a97', '#54b1f7'))
# 




# Res with the linear model
beatAML_original <- c('AKT1', 'AURKA', 'BCR', 'CDK17', 'FLT3', 'INPPL1',
                      'LCK', 'MAP2K2', 'MAP4K5', 'NQO2', 'NTRK1', 'PDGFRB', 'RET', 'RIPK2', 'TNK2')
plot(CV_gp_log$mse/var(as.vector(Y_new), na.rm=TRUE), int_over_union(CV_gp_log$selected, beatAML_original),
     xlab='3-fold CV normalized MSE', ylab='Similarity between sets of selected proteins', ylim=c(0, 0.75), xlim=c(0.7, 2.0), col='purple')
points(CV_gp_original$mse/var(as.vector(Y_original), na.rm=TRUE), int_over_union(CV_gp_original$selected, beatAML_original), col='blue')
points(CV_gp_0_log$mse/var(as.vector(Y_new), na.rm=TRUE), int_over_union(CV_gp_0_log$selected, beatAML_original), col='grey')
points(CV_gp_0_original$mse/var(as.vector(Y_original), na.rm=TRUE), int_over_union(CV_gp_0_original$selected, beatAML_original), col='red')

abline(v=DepInfeR_nMSE, lty=1, col='black')
abline(v=min(rf_nMSE), lty=2, col='black')

points(CV_linear_original$rmse, int_over_union(CV_linear_original$selected, beatAML_original), col='brown')
abline(v=min(CV_linear_original$rmse), lty=2, col='brown')
# points(CV_linear_original_all$rmse, int_over_union(CV_linear_original_all$selected, beatAML_original), col='black')
points(CV_linear_log$rmse, int_over_union(CV_linear_log$selected, beatAML_original), col='green')
abline(v=min(CV_linear_log$rmse), lty=2, col='green')
points(CV_linear_zeroed_log$rmse, int_over_union(CV_linear_zeroed_log$selected, beatAML_original), col='pink')
abline(v=min(CV_linear_zeroed_log$rmse), lty=2, col='pink')
legend(1.43, 0.7, legend=c('GP with NaN in X and log Y', 'GP with NaN in X and z-score Y', 
                           'GP0 with NaN in X and log Y', 'GP0 with NaN in X and z-score Y',
                           "DepInfeR with the original dataset", "MultiRF with the original dataset", 
                           'linear with the original dataset',
                           'linear with the original dataset, all proteins',
                           'linear with original X + log Y',
                           'linear with original X + log Y, all proteins',
                           'linear with NaN in X + log Y',
                           'linear with NaN in X + log Y, all proteins'),
       col=c('purple', 'blue', 'grey', 'red', 'black', 'black', 'brown','brown', 'green','green', 'pink', 'pink'), 
       lty=c(NA, NA, NA, NA,1,2,NA,2,NA,2,NA,2), 
       pch=c(1,1,1,1,NA,NA, 1, NA, 1, NA, 1, NA), cex=0.8)
title('beatAML dataset')














# reproducing the figs in main body
load('EMBL_dataset.RData')
load('EMBL_CV.Rdata')

X_original <- tarMat_EMBL
X_new <- tarMat_EMBL
X_nan <- X_new[1,1]
X_new[X_new == X_nan] <- NaN

Y_original <- viabMat_EMBL
Y_new <- viabMat_EMBL_log


set.seed(314159)
# how about a full non-linear treatment?
rf_test_100 <- rf_test(100, X_original, Y_original)
rf_test_200 <- rf_test(200, X_original, Y_original)
rf_test_300 <- rf_test(300, X_original, Y_original)
rf_test_400 <- rf_test(400, X_original, Y_original)
rf_test_500 <- rf_test(500, X_original, Y_original)

rf_MSE <- sapply(list(rf_test_100, rf_test_200, rf_test_300, rf_test_400, rf_test_500),
                 FUN = function(X){mean((X - Y_original)^2)})
rf_nMSE <- sapply(list(rf_test_100, rf_test_200, rf_test_300, rf_test_400, rf_test_500),
                  FUN = function(X){mean((X - Y_original)^2)}/var(as.vector(Y_original)))


k_fold <- 3
sample_size <- dim(Y_original)[1]
random_partition <- sample(1:sample_size, size=sample_size, replace = FALSE)
group_id <- split(random_partition, sort(1:sample_size %% k_fold))
response_mat <- rep(0, k_fold)
for (i in 1:k_fold){
  test_lasso <- MultiLasso(TargetMatrix=X_original[-group_id[[i]], ], ResponseMatrix=Y_original[-group_id[[i]], ], repeats = 100, BPPARAM = bpparam())
  response_mat[i] <- mean((Y_original[group_id[[i]],] - t(t(X_original[group_id[[i]], ] %*% test_lasso$coefMat) + test_lasso$intercept)
  )^2)
}

DepInfeR_MSE = mean(response_mat)
DepInfeR_nMSE = mean(response_mat)/var(as.vector(Y_original))
# EMBL set
EMBL_original <- c('AFF1', 'AURKA', 'BTK', 'C2CD5', 'CCNT1', 'CDK17', 'CDK6', 'CDK7',
                   'CHEK1', 'DDR1', 'FLT3', 'GAK', 'INPPL1', 'LATS1', 'LIMK1', 'MAP2K2',
                   'MAP4K2', 'MAP4K5', 'PRKAB2', 'PTK2', 'RIPK2', 'RIPK3', 'RPL4', 'SIK2')

plot_data_EMBL <- data.frame(cbind(c(CV_gp_0_log$mse/var(as.vector(Y_new), na.rm=TRUE), 
                                     CV_gp_0_original$mse/var(as.vector(Y_original), na.rm=TRUE)),
                                   c(int_over_union(CV_gp_0_log$selected, EMBL_original), 
                                     int_over_union(CV_gp_0_original$selected, EMBL_original))))
colnames(plot_data_EMBL) <- c('nMSE', 'IoU')
plot_data_EMBL$Dataset <- c(rep('log-transformed, incomplete', 30), rep('original, imputed', 30))
temp_data_1 <- data.frame(val=c(DepInfeR_nMSE, min(rf_nMSE)), Competing.model=c('DepInfeR','MultiRF'))

res_EMBL <- ggplot(data = plot_data_EMBL, aes(x=nMSE, y=IoU)) + geom_point(aes(color=Dataset)) + 
  xlab('3-fold CV normalized MSE') + ylab('Similarity between sets of selected kinases') + 
  ggtitle('EMBL dataset') + geom_vline(aes(xintercept=val, linetype=Competing.model), temp_data_1)







protein_select <- function(cv_res){
  p <- cv_res$selected[[which.min(cv_res$mse)]]
  return(sort(p))
}

protein_list <- list()
protein_list[[1]] <- protein_select(CV_gp_0_log)
protein_list[[2]] <- protein_select(CV_gp_0_original)
protein_list[[3]] <- EMBL_original


protein_counter <- list()
for (i in 1:3){
  for (j in protein_list[[i]]){
    if (is.null(protein_counter[[j]])){
      protein_counter[[j]] <- 1
    }
    else{
      protein_counter[[j]] <- protein_counter[[j]] + 1
    }
  }
}

protein_res <- sort(unlist(protein_counter), decreasing = TRUE)
common_proteins <- sort(names(protein_res[protein_res > 1]))

all_protein <- unique(sort(unlist(protein_list)))

id_mtrx <- matrix(NA, ncol=length(all_protein), nrow=length(protein_list))


my_data <- as.matrix(expand_grid(1:length(protein_list), 1:length(all_protein)))
idct <- rep(NA, length(protein_list)*length(all_protein))

for (i in 1:length(idct)){
  if (all_protein[my_data[i,2]] %in% common_proteins){
    idct[i] <- 2* (all_protein[my_data[i,2]] %in% protein_list[[my_data[i,1]]])
  }
  else{
    idct[i] <- all_protein[my_data[i,2]] %in% protein_list[[my_data[i,1]]]
  }
}

my_face = all_protein
for (i in 1:length(all_protein)){
  if (all_protein[i] %in% common_proteins){
    my_face[i] = 'bold'
  }
  else{
    my_face[i] = 'plain'
  }
}

my_data <- as.data.frame(cbind(my_data, idct))
colnames(my_data) <- c('Group_id', 'Protein_id', 'Label')
my_data$Selection <- as.factor(sapply(my_data$Label, FUN=function(i){c('Not selected', 'Selected', 'Consistent')[i+1]}))

EMBL_selected <- ggplot(my_data, aes(Protein_id, Group_id, fill= Selection)) + geom_tile(color='white') +
  # scale_color_continuous(breaks=c(0,1,2), labels=c('Not selected', 'Selected', 'Consistent')) + 
  scale_y_continuous(breaks = 1:3, labels = c('GP_0_log', 'GP_0_original','DepInfeR')) + 
  scale_x_continuous(breaks = 1:length(all_protein), labels = all_protein)+ xlab("Kinases") + ylab("") + 
  # guides(fill = guide_colourbar(label = FALSE,
  #                               ticks = FALSE))+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, face=my_face),
        axis.line.y = element_blank(),
        axis.line.x = element_blank()) + ggtitle('Selected Kinases in EMBL') + 
  scale_fill_manual(breaks = c('Not selected', 'Selected', 'Consistent'), values=c('#122b44', '#326a97', '#54b1f7'))
# 




# Res with linear models
EMBL_original <- c('AFF1', 'AURKA', 'BTK', 'C2CD5', 'CCNT1', 'CDK17', 'CDK6', 'CDK7',
                   'CHEK1', 'DDR1', 'FLT3', 'GAK', 'INPPL1', 'LATS1', 'LIMK1', 'MAP2K2',
                   'MAP4K2', 'MAP4K5', 'PRKAB2', 'PTK2', 'RIPK2', 'RIPK3', 'RPL4', 'SIK2')
plot(CV_gp_log$mse/var(as.vector(Y_new), na.rm=TRUE), int_over_union(CV_gp_log$selected, EMBL_original),
     xlab='3-fold CV normalized MSE', ylab='Similarity between sets of selected proteins', ylim=c(0, 0.8), xlim=c(0.6, 2.1), col='purple')
points(CV_gp_original$mse/var(as.vector(Y_original), na.rm=TRUE), int_over_union(CV_gp_original$selected, EMBL_original), col='blue')
points(CV_gp_0_log$mse/var(as.vector(Y_new), na.rm=TRUE), int_over_union(CV_gp_0_log$selected, EMBL_original), col='grey')
points(CV_gp_0_original$mse/var(as.vector(Y_original), na.rm=TRUE), int_over_union(CV_gp_0_original$selected, EMBL_original), col='red')

abline(v=DepInfeR_nMSE, lty=1, col='black')
abline(v=min(rf_MSE_norm), lty=2, col='black')

points(CV_linear_log$rmse, int_over_union(CV_linear_log$selected, EMBL_original), col='green')
abline(v=min(CV_linear_log$rmse), lty=2, col='green')
points(CV_linear_zeroed_log$rmse, int_over_union(CV_linear_zeroed_log$selected, EMBL_original), col='pink')
abline(v=min(CV_linear_zeroed_log$rmse), lty=2, col='pink')
legend(1.45, 0.78, legend=c('GP with NaN in X and log Y', 'GP with NaN in X and z-score Y', 
                            'GP0 with NaN in X and log Y', 'GP0 with NaN in X and z-score Y',
                            "DepInfeR with the original dataset", "MultiRF with the original dataset", 
                            # 'linear with the original dataset',
                            # 'linear with the original dataset, all proteins',
                            'linear with original X and log Y',
                            'linear with original X and log Y, all proteins',
                            'linear with NaN in X + log Y',
                            'linear with NaN in X + log Y, all proteins'),
       col=c('purple', 'blue', 'grey', 'red', 'black', 'black', # 'brown','brown', 
             'green','green', 'pink', 'pink'), 
       lty=c(NA, NA, NA, NA,1,2,NA,2,NA,2), # ,NA,2), 
       pch=c(1,1,1,1,NA,NA, 1, NA, 1, NA), # 1, NA), 
       cex=0.8)
title('EMBL dataset')



selected_kinases_list <- list(GDSC1_selected, beatAML_selected, EMBL_selected)
grid.arrange(grobs=selected_kinases_list, ncol=1)



res_list <- list(res_GDSC1, res_beatAML, res_EMBL)
grid.arrange(grobs=res_list, ncol=1)







