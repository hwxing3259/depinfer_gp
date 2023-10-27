library(LaplacesDemon)

# here we define functions used our collapsed Gibbs sampler
log_prior_z <- function(z, inclusion_prob){
  # length of z = total number of proteins, 0=off, 1=on for each protein
  # z \sim Bernoulli(inclusion_prob)
  return(sum(z)*log(inclusion_prob) + sum(1-z)*log(1-inclusion_prob))
}

log_prior_hyper <- function(a0, b0, c0, inclusion_prob){
  # beta 2, 10 prior on inclusion prob, half Normal prior on a0, b0, c0
  ans <- dbeta(inclusion_prob, shape1 = 2, shape2 = 10)
  return(ans + sum(dhalfcauchy(c(a0, b0, c0), log=TRUE)))
}


log_lkd_y <- function(z, X, Y, a0, b0, c0){
  # marginal likelihood of the observation y, will integrate out all missing parts
  # X is a D*P matrix, Y is a D*S matrix, z is a length P binary vector
  # a0, b0 are prior parameters on sigma^2, c0 is the prior var of betas
  D <- dim(X)[1]
  P <- 1 + dim(X)[2] # remember the intercept
  S <- dim(Y)[2]
  my_X <- cbind(1, X)
  my_z <- c(1, z) # always include intercept
  non_NA_id <- !is.na(Y) # also has dimension D*S
  log_lkd <- -0.5*sum(non_NA_id)*log(2*pi) + lgamma(a0 + 0.5*sum(non_NA_id)) - lgamma(a0) + a0 * log(b0)
  # easy part of the marginal likelihood
  
  base_cov <- diag(D) + (my_X %*% diag(my_z) %*% t(my_X))/c0
  base_eigen_decomp <- eigen(base_cov) # VAV^T
  quadratic_form <- 0
  log_det <- 0
  for (s in 1:S){
    if (all(!non_NA_id[, s])){
      print('o boy u just got a column of NA')
    }
    else if (all(non_NA_id[, s])){ # when there is no missing values
      quadratic_form <- quadratic_form + sum(1/base_eigen_decomp$values * (t(base_eigen_decomp$vectors) %*% Y[,s])^2)
      log_det <- log_det - 0.5 * sum(log(base_eigen_decomp$values))
    }
    else { # when we do have missing values at column s of Y
      missing_pattern <- non_NA_id[, s]
      sub_y <- Y[missing_pattern,s]
      sub_cov <- base_cov[missing_pattern, missing_pattern]
      sub_cov_eigen <- eigen(sub_cov)
      quadratic_form <- quadratic_form + sum(1/sub_cov_eigen$values * (t(sub_cov_eigen$vectors) %*% sub_y)^2)
      log_det <- log_det - 0.5 * sum(log(sub_cov_eigen$values))
    }
  }
  log_lkd <- log_lkd + log_det - (a0 + 0.5*sum(non_NA_id))*log(b0 + 0.5*quadratic_form)
  post_IG_parameter <- c(a0 + 0.5*sum(non_NA_id), b0 + 0.5*quadratic_form)
  return(list(log_lkd, post_IG_parameter))
}


reg_par_post <- function(X, z, Y, post_IG_parameter, c0){
  # posterior sampling for beta and sigma^2 given Z are in closed form due to conjugacy
  my_sigma2 <- rgamma(1, shape=post_IG_parameter[1], rate=post_IG_parameter[2])
  # now work out the (block-wise) posterior mean and variance of the regression coefficient beta
  D <- dim(X)[1]
  P <- 1 + dim(X)[2] # remember the intercept
  S <- dim(Y)[2]
  my_X <- cbind(1, X)
  my_z <- c(1, z)
  non_NA_id <- !is.na(Y) # also has dimension D*S
  
  beta <- matrix(NA, nrow = P, ncol = S)
  base_xz <- t(t(my_X) * my_z) # D*(P+1), all columns with z_p=0 are now 0
  base_xtx <- t(base_xz) %*% base_xz
  base_precision <- (base_xtx + c0*diag(P))/my_sigma2 # inverse of the covariance matrix
  base_weight <- solve(base_xtx + c0*diag(P)) %*% t(base_xz)
  
  for (s in 1:S){
    if (all(!non_NA_id[, s])){
      print('o boy u just got a column of NA')
    }
    else if (all(non_NA_id[, s])){ # when there is no missing values
      beta[, s] <- as.vector(rmvnp(n=1, mu=as.vector(base_weight%*%Y[, s]), Omega=base_precision))
    }
    else { # when we do have missing values at column s of Y
      missing_pattern <- non_NA_id[, s]
      sub_y <- Y[missing_pattern,s]
      sub_xz <- base_xz[missing_pattern, ]
      sub_xtx <- t(sub_xz) %*% sub_xz
      sub_precision <- (sub_xtx + c0*diag(P))/my_sigma2
      sub_weight <- solve(sub_xtx + c0*diag(P)) %*% t(sub_xz)
      sub_mean <- sub_weight %*% sub_y
      beta[, s] <- as.vector(rmvnp(n=1, mu=as.vector(sub_mean), Omega=sub_precision))
    }
  }
  
  return(list(my_sigma2, beta))
}


my_gibbs_marginal <- function(N, burn_in, thin, X, Y, c0=1.){
  old_a0 <- 1.
  old_b0 <- 1.
  old_c0 <- 1.
  old_inclusion_prob <- 0.1
  
  old_z <- rbinom(dim(X)[2], 1, old_inclusion_prob) # random init from prior
  old_prior <- log_prior_z(old_z, old_inclusion_prob)
  old_hyper_prior <- log_prior_hyper(old_a0, old_b0, old_c0, old_inclusion_prob)
  old_list <- log_lkd_y(old_z, X, Y, old_a0, old_b0, old_c0)
  old_lkd <- old_list[[1]]
  old_IG <- old_list[[2]]
  
  post_density_record <- rep(NA, N)
  my_Z <- matrix(NA, nrow = N-burn_in, ncol = dim(X)[2])
  my_sig2 <- rep(NA, N-burn_in)
  my_beta <- array(NA, dim = c(N-burn_in, dim(X)[2]+1, dim(Y)[2]))
  my_hyper <- matrix(NA, nrow = N-burn_in, ncol=4)
  P <- dim(X)[2]
  
  for (iter in 1:N){
    for (sub_iter in 1:thin){
      # Gibbs update on Z
      for (id in 1:P){
        # chosen_id <- sample(1:P, 1)
        proposed_z <- old_z
        proposed_z[id] = 1-old_z[id] # flip one bit
        proposed_list <- log_lkd_y(proposed_z, X, Y, old_a0, old_b0, old_c0)
        proposed_lkd <- proposed_list[[1]]
        proposed_prior <- log_prior_z(proposed_z, old_inclusion_prob)
        log_mh_ratio <- proposed_prior + proposed_lkd - old_prior - old_lkd
        if (log(runif(1)) < log_mh_ratio){
          old_z <- proposed_z
          old_lkd <- proposed_lkd
          old_prior <- proposed_prior
          old_IG <- proposed_list[[2]]
        }
      }
      # Gibbs update on hyper parameters
      u <- runif(1, 0.8, 1.25)
      proposed_a0 <- old_a0*u
      proposed_list <- log_lkd_y(old_z, X, Y, proposed_a0, old_b0, old_c0)
      proposed_lkd <- proposed_list[[1]]
      proposed_hyper_prior <- log_prior_hyper(proposed_a0, old_b0, old_c0, old_inclusion_prob)
      log_mh_ratio <- proposed_hyper_prior + proposed_lkd - old_hyper_prior - old_lkd - log(u)
      if (log(runif(1)) < log_mh_ratio){
        old_a0 <- proposed_a0
        old_lkd <- proposed_lkd
        old_hyper_prior <- proposed_hyper_prior
        old_IG <- proposed_list[[2]]
      }
      
      u <- runif(1, 0.8, 1.25)
      proposed_b0 <- old_b0*u
      proposed_list <- log_lkd_y(old_z, X, Y, old_a0, proposed_b0, old_c0)
      proposed_lkd <- proposed_list[[1]]
      proposed_hyper_prior <- log_prior_hyper(old_a0, proposed_b0, old_c0, old_inclusion_prob)
      log_mh_ratio <- proposed_hyper_prior + proposed_lkd - old_hyper_prior - old_lkd - log(u)
      if (log(runif(1)) < log_mh_ratio){
        old_b0 <- proposed_b0
        old_lkd <- proposed_lkd
        old_hyper_prior <- proposed_hyper_prior
        old_IG <- proposed_list[[2]]
      }
      
      u <- runif(1, 0.8, 1.25)
      proposed_c0 <- old_c0*u
      proposed_list <- log_lkd_y(old_z, X, Y, proposed_a0, old_b0, proposed_c0)
      proposed_lkd <- proposed_list[[1]]
      proposed_hyper_prior <- log_prior_hyper(old_a0, old_b0, proposed_c0, old_inclusion_prob)
      log_mh_ratio <- proposed_hyper_prior + proposed_lkd - old_hyper_prior - old_lkd - log(u)
      if (log(runif(1)) < log_mh_ratio){
        old_c0 <- proposed_c0
        old_lkd <- proposed_lkd
        old_hyper_prior <- proposed_hyper_prior
        old_IG <- proposed_list[[2]]
      }
      
      u <- runif(1, 0.8, 1.25)
      proposed_inclusion_prob <- old_inclusion_prob*u
      if (proposed_inclusion_prob>1.){}
      else{
        proposed_prior <- log_prior_z(old_z, proposed_inclusion_prob)
        proposed_hyper_prior <- log_prior_hyper(old_a0, old_b0, old_c0, proposed_inclusion_prob)
        log_mh_ratio <- proposed_hyper_prior + proposed_prior - old_hyper_prior - old_prior - log(u)
        if (log(runif(1)) < log_mh_ratio){
          old_inclusion_prob <- proposed_inclusion_prob
          old_prior <- proposed_prior
          old_hyper_prior <- proposed_hyper_prior
        }
      }
    }
    # after the Gibbs update, record samples if passed burn_in
    post_density_record[iter] <- old_prior + old_lkd + old_hyper_prior
    print(paste('iter ', iter, 'lkd: ', post_density_record[iter]))
    if (iter > burn_in){
      my_Z[(iter-burn_in), ] <- old_z
      # sample and recordregression parameters from the conjugate posterior
      my_reg_para <- reg_par_post(X, old_z, Y, old_IG, old_c0)
      my_sig2[(iter-burn_in)] <- my_reg_para[[1]]
      my_beta[(iter-burn_in),,] <- my_reg_para[[2]]
      my_hyper[(iter-burn_in), ] <- c(old_a0, old_b0, old_c0, old_inclusion_prob)
      print(my_reg_para[[1]])
    }
  }
  return(list(my_Z, my_sig2, my_beta, my_hyper, post_density_record))
}




my_gibbs <- function(N, burn_in, thin, X, Y, a0, b0, c0, inclusion_prob, updated_Z=TRUE, old_z=NULL){
  if (is.null(old_z)){
    old_z <- rbinom(dim(X)[2], 1, inclusion_prob) # random init from prior
  }
  old_prior <- log_prior_z(old_z, inclusion_prob)
  old_list <- log_lkd_y(old_z, X, Y, a0, b0, c0)
  old_lkd <- old_list[[1]]
  old_IG <- old_list[[2]]
  
  post_density_record <- rep(NA, N)
  my_Z <- matrix(NA, nrow = N-burn_in, ncol = dim(X)[2])
  my_sig2 <- rep(NA, N-burn_in)
  my_beta <- array(NA, dim = c(N-burn_in, dim(X)[2]+1, dim(Y)[2]))
  P <- dim(X)[2]
  
  for (iter in 1:N){
    for (sub_iter in 1:thin){
      # Gibbs update on Z
      if (updated_Z){
        for (id in 1:P){
          # chosen_id <- sample(1:P, 1)
          proposed_z <- old_z
          proposed_z[id] = 1-old_z[id] # flip one bit
          proposed_list <- log_lkd_y(proposed_z, X, Y, a0, b0, c0)
          proposed_lkd <- proposed_list[[1]]
          proposed_prior <- log_prior_z(proposed_z, inclusion_prob)
          log_mh_ratio <- proposed_prior + proposed_lkd - old_prior - old_lkd
          if (log(runif(1)) < log_mh_ratio){
            old_z <- proposed_z
            old_lkd <- proposed_lkd
            old_prior <- proposed_prior
            old_IG <- proposed_list[[2]]
          }
        }        
      }
    }
    # after the Gibbs update, record samples if passed burn_in
    post_density_record[iter] <- old_prior + old_lkd
    if (iter%%10 == 0){
      print(paste('iter ', iter, 'log post: ', post_density_record[iter]))
    }
    
    if (iter > burn_in){
      my_Z[(iter-burn_in), ] <- old_z
      # sample and record regression parameters from the conjugate posterior
      my_reg_para <- reg_par_post(X, old_z, Y, old_IG, c0)
      my_sig2[(iter-burn_in)] <- my_reg_para[[1]]
      my_beta[(iter-burn_in),,] <- my_reg_para[[2]]
    }
  }
  return(list(my_Z, my_sig2, my_beta, post_density_record))
}


choose_c0 <- function(N, burn_in, thin, X, Y, vec_c0, a0=1., b0=1., inclusion_prob=.1, k_fold=3, updated_z=TRUE, old_z=NULL){
  cv_mse <- rep(NA, length(vec_c0))
  cv_selected_protein <- lapply(1:length(vec_c0), function(x){NA})
  store_WAIC <- rep(NA, length(vec_c0))
  store_R2 <- rep(NA, length(vec_c0))
  for (j in 1:length(vec_c0)){
    c0 <- vec_c0[j]
    sample_size <- dim(X)[1]
    random_partition <- sample(1:sample_size, size=sample_size, replace = FALSE)
    group_id <- split(random_partition, sort(1:sample_size %% k_fold))
    # group_id <- list(random_partition[1:(sample_size%/%3)],
    #                  random_partition[(1+sample_size%/%3):(2*(sample_size%/%3))],
    #                  random_partition[(1+2*(sample_size%/%3)):sample_size])
    
    record_mse <- matrix(NA, nrow = dim(X)[1], ncol = dim(Y)[2])
    for (fold in 1:k_fold){
      my_group_id <- group_id[[fold]]
      test_MCMC <- my_gibbs(N, burn_in, thin, X[-my_group_id, ], Y[-my_group_id, ], a0, b0, c0, inclusion_prob, updated_z, old_z)
      
      my_mse <- array(NA, dim = c(dim(test_MCMC[[3]])[1], length(my_group_id), dim(Y)[2]))
      
      for (i in 1:dim(my_mse)[1]){
        my_fit <- cbind(1, X[my_group_id, ]) %*% diag(c(1, test_MCMC[[1]][i, ])) %*% test_MCMC[[3]][i,,]
        my_mse[i,,] <- (Y[my_group_id, ] - my_fit)^2
      }
      record_mse[my_group_id, ] <- apply(my_mse, c(2,3), mean, na.rm=TRUE)
    }
    cv_mse[j] <- mean(record_mse, na.rm=TRUE)
    
    test_full_MCMC <- my_gibbs(N+30, burn_in, thin, X, Y, a0, b0, c0, inclusion_prob, updated_z, old_z)
    inclusion_mean <- colMeans(test_full_MCMC[[1]])
    cv_selected_protein[[j]] <- colnames(X)[inclusion_mean>0.5]
    # how about working out the WAIC and see if it is inline with 3-fold CV?
    my_beta <- test_full_MCMC[[3]]
    my_Z <- test_full_MCMC[[1]]
    my_sigma <- test_full_MCMC[[2]]
    my_lkd <- matrix(NA, nrow = dim(Y)[1]*dim(Y)[2], ncol = length(my_sigma))
    my_R2 <- rep(NA, length(my_sigma))
    for (i in 1:length(my_sigma)){
      my_fit <- cbind(1, X) %*% diag(c(1, my_Z[i, ])) %*% my_beta[i,,]
      my_lkd[, i] <- as.vector(dnorm(Y, mean = my_fit, sd = sqrt(my_sigma[i]), log = TRUE))
      my_R2[i] <- cor(Y[!is.na(Y)], my_fit[!is.na(Y)])^2
    }
    store_WAIC[j] <- WAIC(my_lkd[!is.na(my_lkd[,1]), ])$WAIC
    store_R2[j] <- mean(my_R2)
    
  }
  return(list('c0'=c0, 'mse'=cv_mse, 'WAIC'=store_WAIC, 
              'R2'=store_R2, 'selected'=cv_selected_protein))
}


CV_par <- function(N, burn_in, thin, X, Y, vec_c0, a0=1., b0=1., inclusion_prob=.1, ncore=2, k_fold=3, updated_z=TRUE, old_z=NULL){
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  res <- foreach(c0=vec_c0, .packages = 'LaplacesDemon', 
                 .export = c('choose_c0', 'my_gibbs', 'log_prior_z', 'log_prior_hyper', 
                             'log_lkd_y', 'reg_par_post')) %dopar% 
    choose_c0(N, burn_in, thin, X, Y, c0, a0, b0, inclusion_prob, k_fold, updated_z, old_z)
  stopCluster(cl)
  
  pool_c0 <- vec_c0
  pool_mse <- rep(NA, length(pool_c0))
  pool_WAIC <- rep(NA, length(pool_c0))
  pool_R2 <- rep(NA, length(pool_c0))
  pool_cv_selected_protein <- lapply(1:length(vec_c0), function(x){NA})
  for (i in 1:length(vec_c0)){
    pool_mse[i] <- res[[i]]$'mse'
    pool_WAIC[i] <- res[[i]]$'WAIC'
    pool_R2[i] <- res[[i]]$'R2'
    pool_cv_selected_protein[[i]] <- res[[i]]$'selected'[[1]]
  }
  return(list('c0'=pool_c0, 'mse'=pool_mse, 'WAIC'=pool_WAIC, 
              'R2'=pool_R2, 'selected'=pool_cv_selected_protein))
}



log_prior_z_gp <- function(z, inclusion_prob){
  # length of z = total number of proteins, 0=off, 1=on for each protein
  # z \sim Bernoulli(inclusion_prob)
  return(sum(z)*log(inclusion_prob) + sum(1-z)*log(1-inclusion_prob))
}

log_prior_sigma2_gp <- function(sigma2){
  return(dhalft(sigma2, log=TRUE))
}

log_prior_gamma2_gp <- function(gamma2){
  return(dhalfnorm(gamma2, log=TRUE))
}

log_prior_nu <- function(nu){
  return(sum(dhalfnorm(nu, log=TRUE)))
}

kernel_matrix_gp <- function(position, nu){ # 1D Gaussian RBF kernel
  return(nu[1]*exp(sapply(position, FUN = function(x){-((position-x)^2)/(2*nu[2])})))
}

my_kernels_gp <- function(X, nu){ # how one compute the kernel matrices when we ignore the NA terms
  # also try other possibilities e.g. fixing f(NA)=0 or treat f(NA) as an unknown parameter?
  # work out the P kernel matrices
  kernels <- array(NA, dim=c(dim(X)[2], dim(X)[1], dim(X)[1]))
  for (p in 1:dim(X)[2]){
    non_NA_X_id <- !is.na(X[,p])
    if (any(non_NA_X_id)){
      non_NA_X <- X[non_NA_X_id, p]
      K_star <- kernel_matrix_gp(non_NA_X, nu)
      K <- matrix(0, ncol=dim(X)[1], nrow=dim(X)[1])
      K[non_NA_X_id, non_NA_X_id] <- K_star
      kernels[p,,] <- K      
    }
    else{
      kernels[p,,] <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[1])
    }
    
  }
  return(kernels)
}


y_lkd_gp <- function(X, Y, Z, sigma2, gamma2, kernels){
  base_cov <- apply(kernels[as.logical(Z),,,drop=FALSE], c(2,3), sum) + gamma2 + sigma2*diag(dim(X)[1])
  base_eigen_decomp <- eigen(base_cov) # VAV^T
  quadratic_form <- 0
  log_det <- 0
  for (s in 1:dim(Y)[2]){
    y_sub_obs <- Y[,s]
    non_na_id <- !is.na(y_sub_obs)
    if (all(!non_na_id)){
      print('o boy u just got a column of NA')
    }
    else if (all(non_na_id)){ # when there is no missing values
      quadratic_form <- quadratic_form + sum(1/base_eigen_decomp$values * (t(base_eigen_decomp$vectors) %*% Y[,s])^2)
      log_det <- log_det - 0.5 * sum(log(base_eigen_decomp$values))
    }
    else{ # when we do have missing values at column s of Y
      missing_pattern <- non_na_id
      sub_y <- Y[missing_pattern,s]
      sub_cov <- base_cov[missing_pattern, missing_pattern]
      sub_cov_eigen <- eigen(sub_cov)
      quadratic_form <- quadratic_form + sum(1/sub_cov_eigen$values * (t(sub_cov_eigen$vectors) %*% sub_y)^2)
      log_det <- log_det - 0.5 * sum(log(sub_cov_eigen$values))
    }
    
    # rest should be similar to what we have earlier, copy and paste?
  }
  return(-0.5*sum(!is.na(Y))*log(2*pi) + log_det - 0.5 * quadratic_form )
}


my_gibbs_sampler_gp <- function(N, burn_in, thin, X, Y, inclusion_prob, nu, update_Z=TRUE){
  old_sigma2 <- rgamma(1, 1, 1)
  old_gamma2 <- rgamma(1, 1, 1)
  old_z <- rbinom(dim(X)[2], 1, inclusion_prob)
  old_kernel_matrices <- my_kernels_gp(X, nu)
  old_lkd <- y_lkd_gp(X, Y, old_z, old_sigma2, old_gamma2, old_kernel_matrices)
  
  post_density_record <- rep(NA, N)
  lkd_record <- rep(NA, N)
  my_Z <- matrix(NA, nrow = N-burn_in, ncol = dim(X)[2])
  my_sigma2 <- rep(NA, N-burn_in)
  my_gamma2 <- rep(NA, N-burn_in)
  P <- dim(X)[2]
  
  for (iter in 1:N){
    for (sub_iter in 1:thin){
      if (update_Z){
        # Gibbs update on Z
        for (id in 1:P){
          # chosen_id <- sample(1:P, 1)
          proposed_z <- old_z
          proposed_z[id] = 1-old_z[id] # flip one bit
          proposed_lkd <- y_lkd_gp(X, Y, proposed_z, old_sigma2, old_gamma2, old_kernel_matrices)
          proposed_prior <- log_prior_z_gp(proposed_z, inclusion_prob)
          log_mh_ratio <- log_prior_z_gp(proposed_z, inclusion_prob) + proposed_lkd - 
            log_prior_z_gp(old_z, inclusion_prob) - old_lkd
          if (log(runif(1)) < log_mh_ratio){
            old_z <- proposed_z
            old_lkd <- proposed_lkd
          }
        }        
      }
      
      for (xx in 1:7){
        # now update sigma2 and gamma2
        u <- runif(1, 0.8, 1.25)
        proposed_sigma2 <- old_sigma2 * u
        proposed_lkd <- y_lkd_gp(X, Y, old_z, proposed_sigma2, old_gamma2, old_kernel_matrices)
        log_mh_ratio <- log_prior_sigma2_gp(proposed_sigma2) + proposed_lkd - 
          log_prior_sigma2_gp(old_sigma2) - old_lkd - log(u)
        if (log(runif(1)) < log_mh_ratio){
          old_sigma2 <- proposed_sigma2
          old_lkd <- proposed_lkd
        }
        
        u <- runif(1, 0.8, 1.25)
        proposed_gamma2 <- old_gamma2 * u
        proposed_lkd <- y_lkd_gp(X, Y, old_z, old_sigma2, proposed_gamma2, old_kernel_matrices)
        log_mh_ratio <- log_prior_gamma2_gp(proposed_gamma2) + proposed_lkd - 
          log_prior_gamma2_gp(old_gamma2) - old_lkd - log(u)
        if (log(runif(1)) < log_mh_ratio){
          old_gamma2 <- proposed_gamma2
          old_lkd <- proposed_lkd
        }    
      }
    }
    
    # after the Gibbs update, record samples if passed burn_in
    lkd_record[iter] <- old_lkd
    post_density_record[iter] <- old_lkd + log_prior_gamma2_gp(old_gamma2) + 
      log_prior_sigma2_gp(proposed_sigma2) + log_prior_z_gp(old_z, inclusion_prob)
    if (iter%%10 == 0){
      print(paste('iter ', iter, 'log post: ', post_density_record[iter]))
    }
    
    if (iter > burn_in){
      my_Z[(iter-burn_in), ] <- old_z
      my_sigma2[(iter-burn_in)] <- old_sigma2
      my_gamma2[(iter-burn_in)] <- old_gamma2
    }
  }
  return(list('Z'=my_Z, 'sigma2'=my_sigma2, 'gamma2'=my_gamma2, 
              'log_post_density'=post_density_record, 'loglkd'=lkd_record))
}

# still need a function inferring the GP function values!
GP_para <- function(burn_in, X, Y, Z, sigma2, gamma2, nu, test_loc=NULL, X_pred=NULL){
  old_kernel_matrices <- my_kernels_gp(X, nu)
  D <- dim(X)[1]
  P <- dim(X)[2]
  S <- dim(Y)[2]
  R <- length(sigma2)
  if (is.null(test_loc)){
    test_loc <- seq(0.01,0.99,length.out=10)
  }
  curve_matrix <- array(NA, dim = c(R, S, P, length(test_loc)))
  curve_var_matrix <- array(NA, dim = c(R, S, P, length(test_loc)))
  intercept_matrix <- matrix(NA, nrow = R, ncol = S)
  
  y_hat <- array(NA, dim=c(D, S, R))
  if (!is.null(X_pred)){
    y_pred_hat <- array(0, dim=c(dim(X_pred)[1], S, R)) 
    y_pred_var <- array(0, dim=c(dim(X_pred)[1], S, R)) 
  }
  else{
    y_pred_hat <- NULL
    y_pred_var <- NULL
  }
  
  
  for (r in 1:R){ # run analysis for each (Z, sigma2, gamma2) pair
    para_store <- array(0, dim=c(D, P, S))
    intercept_store <- rep(0, S)
    for (s in 1:S){ # we essentially repeat all analysis for each S
      Y_s_non_na <- !is.na(Y[,s])
      for (sub_iter in 1:burn_in){ # running some kind of Bayesian backfitting
        for (p in 1:P){ # checking out all proteins
          if (Z[r, p] != 0){ # if the pth protein is relevant/chosen
            X_p_non_na <- !is.na(X[,p]) # identify proteins with non-na affinity  
            total_non_na_id <- as.logical(Y_s_non_na * X_p_non_na)
            if (any(total_non_na_id)){
              y_sub_p <- (Y[,s] - intercept_store[s] - rowSums(para_store[, -p, s]))[total_non_na_id]
              X_p_relative_position <- X[total_non_na_id, p]
              K_sp_obs_kernel <- as.matrix(kernel_matrix_gp(X_p_relative_position, nu))
              K_star_star <- as.matrix(kernel_matrix_gp(X[X_p_non_na, p], nu))
              K_star <- nu[1]*exp(sapply(X[X_p_non_na, p], 
                                         FUN = function(x){-((X_p_relative_position-x)^2)/(2*nu[2])}))
              if (is.vector(K_star)){
                K_star <- as.matrix(K_star)
              }
              else{
                K_star <- t(K_star)
              }
              C <- solve(K_sp_obs_kernel + sigma2[r]*diag(dim(K_sp_obs_kernel)[1]))
              C <- round(C, 9)
              GP_mean_p <- K_star %*% C %*% y_sub_p
              GP_cov_p <- (K_star_star - round(K_star %*% C %*% t(K_star), 9))
              GP_cov_p <- 0.5*(GP_cov_p + t(GP_cov_p)) + 1e-6*diag(dim(K_star_star)[1]) # for numerical stability
              para_store[X_p_non_na, p, s] <- as.vector(rmvn(n = 1, mu = as.vector(GP_mean_p), Sigma = GP_cov_p))
            }
          }
        }
        # now update the sth intercept
        y_sub_a_s <- (Y[,s] - rowSums(para_store[, , s]))[Y_s_non_na]
        intercept_store[s] <- rnorm(n=1, mean=mean(y_sub_a_s)*gamma2[r]/(gamma2[r]+sigma2[r]/sum(Y_s_non_na)),
                                    sd=sqrt(1/(1/gamma2[r] + sum(Y_s_non_na)/sigma2[r])))
      }
      
      # after the burn-in in period of back fitting, we assume the parameter estimation has converged,
      # record the posterior mean curve of each GP whose p is chosen by Z[r,]!
      
      for (p in 1:P){ # gibbs sampler targeting f_sp(curve), f_sp(obs points) | everything else, 
        #turns out to be normal and the f_sp(obs points) can be integrated out nicely
        if (Z[r, p] != 0){
          X_p_non_na <- !is.na(X[,p]) # identify proteins with non-na affinity
          total_non_na_id <- as.logical(Y_s_non_na * X_p_non_na)
          if (any(total_non_na_id)){
            y_sub_p <- (Y[,s] - intercept_store[s] - rowSums(para_store[, -p, s]))[total_non_na_id]
            X_p_relative_position <- X[total_non_na_id, p]
            K_sp_obs_kernel <- as.matrix(kernel_matrix_gp(X_p_relative_position, nu))
            K_star_star <- as.matrix(kernel_matrix_gp(test_loc, nu))
            K_star <- nu[1]*exp(sapply(test_loc, FUN = function(x){-((X_p_relative_position-x)^2)/(2*nu[2])}))
            if (is.vector(K_star)){
              K_star <- as.matrix(K_star)
            }
            else{
              K_star <- t(K_star)
            }
            C <- solve(K_sp_obs_kernel + sigma2[r]*diag(dim(K_sp_obs_kernel)[1]))
            GP_mean_p <- K_star %*% C %*% y_sub_p
            GP_cov_p <- K_star_star - K_star %*% C %*% t(K_star)
            curve_matrix[r, s, p, ] <- GP_mean_p
            curve_var_matrix[r, s, p, ] <- diag(GP_cov_p)
          }
          
          # now let's worry about the prediction, given a set of new X,
          # compute the predicted value of f_{sp}(X_dp)
          # still, one can integrate out f_sp samples in para_store
          # adding contribution form each selected protein
          if (!is.null(X_pred)){
            X_pred_p_non_na <- !is.na(X_pred[,p]) # identify proteins with non-na affinity in X_pred
            if (any(total_non_na_id) && any(X_pred_p_non_na)){ # if we have observed some data from the training set and it also appears in 
              # the testing set. If either is false, the contribution of this p-s pair to the test point should be 0
              X_pred_p_relative_position <- X_pred[X_pred_p_non_na, p]
              K_pred_star_star <- as.matrix(kernel_matrix_gp(X_pred_p_relative_position, nu))
              K_pred_star <- nu[1]*exp(sapply(X_pred_p_relative_position,
                                              FUN = function(x){-((X_p_relative_position-x)^2)/(2*nu[2])}))
              if (is.vector(K_pred_star)){
                K_pred_star <- as.matrix(K_pred_star)
              }
              else{
                K_pred_star <- t(K_pred_star)
              }
              
              pred_GP_mean_p <- K_pred_star %*% C %*% y_sub_p
              pred_GP_var_p <- diag(K_pred_star_star - K_pred_star %*% C %*% t(K_pred_star))
              
              y_pred_hat[X_pred_p_non_na, s, r] <- y_pred_hat[X_pred_p_non_na, s, r] + pred_GP_mean_p
              y_pred_var[X_pred_p_non_na, s, r] <- y_pred_var[X_pred_p_non_na, s, r] + pred_GP_var_p
            }
          }
        }
      }
      intercept_matrix[r, s] <- intercept_store[s]
      # y_pred_hat_non_zero <- y_pred_hat[, s, r] != 0
      # y_pred_hat[y_pred_hat_non_zero, s, r] <- y_pred_hat[y_pred_hat_non_zero, s, r] + intercept_store[s]
      y_pred_hat[, s, r] <- y_pred_hat[, s, r] + intercept_store[s]
    }
    y_hat[,,r] <- t(t(apply(para_store, MARGIN = c(1,3), sum)) + intercept_store)
  }
  return(list('curve_mean'=curve_matrix, 
              'curve_sd'=curve_var_matrix, 
              'intercept'=intercept_matrix,
              'y_hat'=y_hat,
              'y_pred_hat'=y_pred_hat,
              'y_pred_var'=y_pred_var))
}


# also need a function taking y_hat as input and return WAIC
my_WAIC_gp <- function(Y, yhat, my_sigma){
  non_NA_id <- !is.na(Y)
  my_lkd_gp_mtrx <- matrix(NA, ncol = length(my_sigma), nrow = sum(non_NA_id))
  for (r in 1:length(my_sigma)){
    my_yhat <- yhat[,,r]
    my_lkd_gp_mtrx[, r] <- dnorm(Y[non_NA_id], mean = my_yhat[non_NA_id], sd=sqrt(my_sigma[r]), log = TRUE)
  }
  return(WAIC(my_lkd_gp_mtrx))
}


# No lets stick with 3-fold CV!
choose_nu <- function(N, burn_in, thin, X, Y, nu, inclusion_prob=.1, back_fit_iter=20, k_fold=3){
  vec_nu_1 <- nu[1]
  vec_nu_2 <- nu[2]
  cv_mse <- rep(NA, length(vec_nu_1)*length(vec_nu_2))
  cv_selected_protein <- lapply(1:length(vec_nu_1)*length(vec_nu_2), function(x){NA})
  store_WAIC <- rep(NA, length(vec_nu_1)*length(vec_nu_2))
  store_R2 <- rep(NA, length(vec_nu_1)*length(vec_nu_2))
  nu_vec <- cbind(rep(vec_nu_1, each=length(vec_nu_2)), 
                  rep(vec_nu_2, times=length(vec_nu_1))) # each row is a para pair
  
  for (j in 1:(dim(nu_vec)[1])){
    nu <- nu_vec[j, ]
    sample_size <- dim(X)[1]
    random_partition <- sample(1:sample_size, size=sample_size, replace = FALSE)
    group_id <- split(random_partition, sort(1:sample_size %% k_fold))
    
    record_mse <- matrix(NA, nrow = dim(Y)[1], ncol = dim(Y)[2])
    for (fold in 1:k_fold){
      my_group_id <- group_id[[fold]]
      my_gibbs <- my_gibbs_sampler_gp(N, burn_in, thin, X[-my_group_id, ], Y[-my_group_id, ], inclusion_prob, nu)
      my_para <- GP_para(back_fit_iter, X[-my_group_id, ], Y[-my_group_id, ], 
                         my_gibbs$Z, my_gibbs$sigma2, my_gibbs$gamma2, nu, X_pred = X[my_group_id, ])
      
      my_mse <- array(NA, dim = c(length(my_gibbs$'sigma2'), length(my_group_id), dim(Y)[2])) # this is RD'S
      for (i in 1:dim(my_mse)[1]){
        my_fit <- my_para$'y_pred_hat'[,,i] # recall it has dim D'SR
        my_mse[i,,] <- (Y[my_group_id, ] - my_fit)^2
        test_na <- !is.na(as.vector(Y[my_group_id, ]))
        # plot(as.vector(Y[my_group_id, ])[test_na], as.vector(my_fit)[test_na])
        # abline(a=0, b=1)
      }
      record_mse[my_group_id, ] <- apply(my_mse, c(2,3), mean, na.rm=TRUE)
    }
    cv_mse[j] <- mean(record_mse, na.rm=TRUE)
    print(record_mse)
    
    test_full_MCMC <- my_gibbs_sampler_gp(N, burn_in, thin, X, Y, inclusion_prob, nu)
    test_full_para <- GP_para(back_fit_iter, X, Y, 
                              test_full_MCMC$Z, test_full_MCMC$sigma2, test_full_MCMC$gamma2, nu)  
    inclusion_mean <- colMeans(test_full_MCMC$'Z')
    cv_selected_protein[[j]] <- colnames(X)[inclusion_mean>0.5]
    
    my_sigma <- test_full_MCMC$'sigma2'
    print(my_sigma)
    my_lkd_gp <- matrix(NA, nrow = dim(Y)[1]*dim(Y)[2], ncol = length(my_sigma))
    my_R2 <- rep(NA, length(my_sigma))
    for (i in 1:length(my_sigma)){
      my_lkd_gp[, i] <- as.vector(dnorm(Y, mean = test_full_para$'y_hat'[,,i], sd = sqrt(my_sigma[i]), log = TRUE))
      my_R2[i] <- cor(Y[!is.na(Y)], (test_full_para$'y_hat'[,,i])[!is.na(Y)])^2
    }    
    store_WAIC[j] <- WAIC(my_lkd_gp[!is.na(my_lkd_gp[,1]), ])$WAIC
    store_R2[j] <- mean(my_R2) 
  }
  return(list('nu'=nu_vec, 'mse'=cv_mse, 'WAIC'=store_WAIC, 
              'R2'=store_R2, 'selected'=cv_selected_protein))
}


CV_nu_par <- function(N, burn_in, thin, X, Y, vec_nu_1, vec_nu_2, 
                      inclusion_prob=.1, back_fit_iter=20, k_fold=3, ncor=2){
  nu_vec <- cbind(rep(vec_nu_1, each=length(vec_nu_2)), 
                  rep(vec_nu_2, times=length(vec_nu_1))) # each row is a para pair
  nu_list <- lapply(seq_len(nrow(nu_vec)), function(i) nu_vec[i,])
  cl <- makeCluster(ncor)
  registerDoParallel(cl)
  # start from here!
  res <- foreach(nu=nu_list, .packages = 'LaplacesDemon', 
                 .export = c('choose_nu', 'GP_para', 'my_gibbs_sampler_gp', 
                             'y_lkd_gp', 'my_kernels_gp', 'kernel_matrix_gp', 
                             'log_prior_z_gp', 'log_prior_gamma2_gp', 
                             'log_prior_sigma2_gp')) %dopar% 
    choose_nu(N=N, burn_in=burn_in, thin=thin, X=X, Y=Y, nu=nu, inclusion_prob=inclusion_prob, 
              back_fit_iter=back_fit_iter, k_fold=k_fold)
  stopCluster(cl)
  
  
  pool_nu <- nu_vec
  pool_mse <- rep(NA, length(nu_list))
  pool_WAIC <- rep(NA, length(nu_list))
  pool_R2 <- rep(NA, length(nu_list))
  pool_cv_selected_protein <- lapply(1:length(nu_list), function(x){NA})
  for (i in 1:length(nu_list)){
    pool_mse[i] <- res[[i]]$'mse'
    pool_WAIC[i] <- res[[i]]$'WAIC'
    pool_R2[i] <- res[[i]]$'R2'
    pool_cv_selected_protein[[i]] <- res[[i]]$'selected'[[1]]
  }
  return(list('c0'=nu_list, 'mse'=pool_mse, 'WAIC'=pool_WAIC, 
              'R2'=pool_R2, 'selected'=pool_cv_selected_protein))
  
}



log_prior_z_gp_0 <- function(z, inclusion_prob){
  # length of z = total number of proteins, 0=off, 1=on for each protein
  # z \sim Bernoulli(inclusion_prob)
  return(sum(z)*log(inclusion_prob) + sum(1-z)*log(1-inclusion_prob))
}

log_prior_sigma2_gp_0 <- function(sigma2){
  return(dhalft(sigma2, log=TRUE))
}

log_prior_gamma2_gp_0 <- function(gamma2){
  return(dhalfnorm(gamma2, log=TRUE))
}

log_prior_nu_0 <- function(nu){
  return(sum(dhalfnorm(nu, log=TRUE)))
}

kernel_matrix_gp_0 <- function(position, nu){ # 1D Gaussian RBF kernel
  return(nu[1]*exp(sapply(position, FUN = function(x){-((position-x)^2)/(2*nu[2])})))
  # K0 <- nu[1]*exp(-(position^2)/(2*nu[2]))
  # return(res - K0 %*% t(K0)/nu[1])
}

# my_kernels_gp <- function(X, nu){ # how one compute the kernel matrices when we ignore the NA terms
#   # also try other possibilities e.g. fixing f(NA)=0 or treat f(NA) as an unknown parameter?
#   # work out the P kernel matrices
#   kernels <- array(NA, dim=c(dim(X)[2], dim(X)[1], dim(X)[1]))
#   for (p in 1:dim(X)[2]){
#     non_NA_X_id <- !is.na(X[,p])
#     if (any(non_NA_X_id)){
#       non_NA_X <- X[non_NA_X_id, p]
#       K_star <- kernel_matrix_gp_0(non_NA_X, nu)
#       K <- matrix(0, ncol=dim(X)[1], nrow=dim(X)[1])
#       K[non_NA_X_id, non_NA_X_id] <- K_star
#       kernels[p,,] <- K      
#     }
#     else{
#       kernels[p,,] <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[1])
#     }
#     
#   }
#   return(kernels)
# }

my_kernels_gp_0 <- function(X, nu){ # how one compute the kernel matrices when we ignore the NA terms
  # also try other possibilities e.g. fixing f(NA)=0 or treat f(NA) as an unknown parameter?
  # work out the P kernel matrices
  kernels <- array(NA, dim=c(dim(X)[2], dim(X)[1], dim(X)[1]))
  for (p in 1:dim(X)[2]){
    non_NA_X_id <- !is.na(X[,p])
    if (any(non_NA_X_id)){
      non_NA_X <- X[non_NA_X_id, p]
      K_star <- kernel_matrix_gp_0(non_NA_X, nu)
      K0 <- nu[1]*exp(-(non_NA_X^2)/(2*nu[2]))
      K_star <- K_star - K0 %*% t(K0)/nu[1]
      
      K <- matrix(0, ncol=dim(X)[1], nrow=dim(X)[1])
      K[non_NA_X_id, non_NA_X_id] <- K_star
      kernels[p,,] <- K      
    }
    else{
      kernels[p,,] <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[1])
    }
    
  }
  return(kernels)
}



y_lkd_gp_0 <- function(X, Y, Z, sigma2, gamma2, kernels){
  base_cov <- apply(kernels[as.logical(Z),,,drop=FALSE], c(2,3), sum) + gamma2 + sigma2*diag(dim(X)[1])
  base_eigen_decomp <- eigen(base_cov) # VAV^T
  quadratic_form <- 0
  log_det <- 0
  for (s in 1:dim(Y)[2]){
    y_sub_obs <- Y[,s]
    non_na_id <- !is.na(y_sub_obs)
    if (all(!non_na_id)){
      print('o boy u just got a column of NA')
    }
    else if (all(non_na_id)){ # when there is no missing values
      quadratic_form <- quadratic_form + sum(1/base_eigen_decomp$values * (t(base_eigen_decomp$vectors) %*% Y[,s])^2)
      log_det <- log_det - 0.5 * sum(log(base_eigen_decomp$values))
    }
    else{ # when we do have missing values at column s of Y
      missing_pattern <- non_na_id
      sub_y <- Y[missing_pattern,s]
      sub_cov <- base_cov[missing_pattern, missing_pattern]
      sub_cov_eigen <- eigen(sub_cov)
      quadratic_form <- quadratic_form + sum(1/sub_cov_eigen$values * (t(sub_cov_eigen$vectors) %*% sub_y)^2)
      log_det <- log_det - 0.5 * sum(log(sub_cov_eigen$values))
    }
    
    # rest should be similar to what we have earlier, copy and paste?
  }
  return(-0.5*sum(!is.na(Y))*log(2*pi) + log_det - 0.5 * quadratic_form )
}


my_gibbs_sampler_gp_0 <- function(N, burn_in, thin, X, Y, inclusion_prob, nu, update_Z=TRUE){
  old_sigma2 <- rgamma(1, 1, 1)
  old_gamma2 <- rgamma(1, 1, 1)
  old_z <- rbinom(dim(X)[2], 1, inclusion_prob)
  old_kernel_matrices <- my_kernels_gp_0(X, nu)
  old_lkd <- y_lkd_gp_0(X, Y, old_z, old_sigma2, old_gamma2, old_kernel_matrices)
  
  post_density_record <- rep(NA, N)
  lkd_record <- rep(NA, N)
  my_Z <- matrix(NA, nrow = N-burn_in, ncol = dim(X)[2])
  my_sigma2 <- rep(NA, N-burn_in)
  my_gamma2 <- rep(NA, N-burn_in)
  P <- dim(X)[2]
  
  for (iter in 1:N){
    for (sub_iter in 1:thin){
      if (update_Z){
        # Gibbs update on Z
        for (id in 1:P){
          # chosen_id <- sample(1:P, 1)
          proposed_z <- old_z
          proposed_z[id] = 1-old_z[id] # flip one bit
          proposed_lkd <- y_lkd_gp_0(X, Y, proposed_z, old_sigma2, old_gamma2, old_kernel_matrices)
          proposed_prior <- log_prior_z_gp_0(proposed_z, inclusion_prob)
          log_mh_ratio <- log_prior_z_gp_0(proposed_z, inclusion_prob) + proposed_lkd - 
            log_prior_z_gp_0(old_z, inclusion_prob) - old_lkd
          if (log(runif(1)) < log_mh_ratio){
            old_z <- proposed_z
            old_lkd <- proposed_lkd
          }
        }        
      }
      
      for (xx in 1:7){
        # now update sigma2 and gamma2
        u <- runif(1, 0.8, 1.25)
        proposed_sigma2 <- old_sigma2 * u
        proposed_lkd <- y_lkd_gp_0(X, Y, old_z, proposed_sigma2, old_gamma2, old_kernel_matrices)
        log_mh_ratio <- log_prior_sigma2_gp_0(proposed_sigma2) + proposed_lkd - 
          log_prior_sigma2_gp_0(old_sigma2) - old_lkd - log(u)
        if (log(runif(1)) < log_mh_ratio){
          old_sigma2 <- proposed_sigma2
          old_lkd <- proposed_lkd
        }
        
        u <- runif(1, 0.8, 1.25)
        proposed_gamma2 <- old_gamma2 * u
        proposed_lkd <- y_lkd_gp_0(X, Y, old_z, old_sigma2, proposed_gamma2, old_kernel_matrices)
        log_mh_ratio <- log_prior_gamma2_gp_0(proposed_gamma2) + proposed_lkd - 
          log_prior_gamma2_gp_0(old_gamma2) - old_lkd - log(u)
        if (log(runif(1)) < log_mh_ratio){
          old_gamma2 <- proposed_gamma2
          old_lkd <- proposed_lkd
        }    
      }
    }
    
    # after the Gibbs update, record samples if passed burn_in
    lkd_record[iter] <- old_lkd
    post_density_record[iter] <- old_lkd + log_prior_gamma2_gp_0(old_gamma2) + 
      log_prior_sigma2_gp_0(proposed_sigma2) + log_prior_z_gp_0(old_z, inclusion_prob)
    if (iter%%10 == 0){
      print(paste('iter ', iter, 'log post: ', post_density_record[iter]))
    }
    
    if (iter > burn_in){
      my_Z[(iter-burn_in), ] <- old_z
      my_sigma2[(iter-burn_in)] <- old_sigma2
      my_gamma2[(iter-burn_in)] <- old_gamma2
    }
  }
  return(list('Z'=my_Z, 'sigma2'=my_sigma2, 'gamma2'=my_gamma2, 
              'log_post_density'=post_density_record, 'loglkd'=lkd_record))
}

# still need a function inferring the GP function values!
GP_para_0 <- function(burn_in, X, Y, Z, sigma2, gamma2, nu, test_loc=NULL, X_pred=NULL){
  D <- dim(X)[1]
  P <- dim(X)[2]
  S <- dim(Y)[2]
  R <- length(sigma2)
  if (is.null(test_loc)){
    test_loc <- seq(0.01,0.99,length.out=10)
  }
  curve_matrix <- array(NA, dim = c(R, S, P, length(test_loc)))
  curve_var_matrix <- array(NA, dim = c(R, S, P, length(test_loc)))
  intercept_matrix <- matrix(NA, nrow = R, ncol = S)
  
  y_hat <- array(NA, dim=c(D, S, R))
  if (!is.null(X_pred)){
    y_pred_hat <- array(0, dim=c(dim(X_pred)[1], S, R)) 
    y_pred_var <- array(0, dim=c(dim(X_pred)[1], S, R)) 
  }
  else{
    y_pred_hat <- NULL
    y_pred_var <- NULL
  }
  
  
  for (r in 1:R){ # run analysis for each (Z, sigma2, gamma2) pair
    para_store <- array(0, dim=c(D, P, S))
    intercept_store <- rep(0, S)
    for (s in 1:S){ # we essentially repeat all analysis for each S
      Y_s_non_na <- !is.na(Y[,s])
      for (sub_iter in 1:burn_in){ # running some kind of Bayesian backfitting
        for (p in 1:P){ # checking out all proteins
          if (Z[r, p] != 0){ # if the pth protein is relevant/chosen
            X_p_non_na <- !is.na(X[,p]) # identify proteins with non-na affinity  
            total_non_na_id <- as.logical(Y_s_non_na * X_p_non_na)
            if (any(total_non_na_id)){
              y_sub_p <- c(0, (Y[,s] - intercept_store[s] - rowSums(para_store[, -p, s]))[total_non_na_id]) # append f(0)=0
              X_p_relative_position <- c(0, X[total_non_na_id, p]) # append noiseless f(0)=0
              K_sp_obs_kernel <- as.matrix(kernel_matrix_gp_0(X_p_relative_position, nu))
              K_star_star <- as.matrix(kernel_matrix_gp_0(X[X_p_non_na, p], nu))
              K_star <- nu[1]*exp(sapply(X[X_p_non_na, p], 
                                         FUN = function(x){-((X_p_relative_position-x)^2)/(2*nu[2])}))
              if (is.vector(K_star)){
                K_star <- as.matrix(K_star)
              }
              else{
                K_star <- t(K_star)
              }
              C <- solve(K_sp_obs_kernel + diag(c(0, rep(sigma2[r], dim(K_sp_obs_kernel)[1]-1)))) # recall that f(0)=0 is noiseless
              # sigma2[r]*diag(dim(K_sp_obs_kernel)[1]))
              C <- round(C, 9)
              GP_mean_p <- K_star %*% C %*% y_sub_p
              GP_cov_p <- (K_star_star - round(K_star %*% C %*% t(K_star), 9))
              GP_cov_p <- 0.5*(GP_cov_p + t(GP_cov_p)) + 1e-6*diag(dim(K_star_star)[1]) # for numerical stability
              para_store[X_p_non_na, p, s] <- as.vector(rmvn(n = 1, mu = as.vector(GP_mean_p), Sigma = GP_cov_p))
            }
          }
        }
        # now update the sth intercept
        y_sub_a_s <- (Y[,s] - rowSums(para_store[, , s]))[Y_s_non_na]
        intercept_store[s] <- rnorm(n=1, mean=mean(y_sub_a_s)*gamma2[r]/(gamma2[r]+sigma2[r]/sum(Y_s_non_na)),
                                    sd=sqrt(1/(1/gamma2[r] + sum(Y_s_non_na)/sigma2[r])))
      }
      
      # after the burn-in in period of back fitting, we assume the parameter estimation has converged,
      # record the posterior mean curve of each GP whose p is chosen by Z[r,]!
      
      for (p in 1:P){ # gibbs sampler targeting f_sp(curve), f_sp(obs points) | everything else, 
        #turns out to be normal and the f_sp(obs points) can be integrated out nicely
        if (Z[r, p] != 0){
          X_p_non_na <- !is.na(X[,p]) # identify proteins with non-na affinity
          total_non_na_id <- as.logical(Y_s_non_na * X_p_non_na)
          if (any(total_non_na_id)){
            y_sub_p <- c(0, (Y[,s] - intercept_store[s] - rowSums(para_store[, -p, s]))[total_non_na_id]) # adding noiseless f0)=0
            X_p_relative_position <- c(0, X[total_non_na_id, p]) # adding noiseless f(0)=0
            K_sp_obs_kernel <- as.matrix(kernel_matrix_gp_0(X_p_relative_position, nu))
            K_star_star <- as.matrix(kernel_matrix_gp_0(test_loc, nu))
            K_star <- nu[1]*exp(sapply(test_loc, FUN = function(x){-((X_p_relative_position-x)^2)/(2*nu[2])}))
            if (is.vector(K_star)){
              K_star <- as.matrix(K_star)
            }
            else{
              K_star <- t(K_star)
            }
            C <- solve(K_sp_obs_kernel + diag(c(0, rep(sigma2[r], dim(K_sp_obs_kernel)[1]-1)))) # recall that f(0)=0 is noiseless
            # + sigma2[r]*diag(dim(K_sp_obs_kernel)[1]))
            GP_mean_p <- K_star %*% C %*% y_sub_p
            GP_cov_p <- K_star_star - K_star %*% C %*% t(K_star)
            curve_matrix[r, s, p, ] <- GP_mean_p
            curve_var_matrix[r, s, p, ] <- diag(GP_cov_p)
          }
          
          # now let's worry about the prediction, given a set of new X,
          # compute the predicted value of f_{sp}(X_dp)
          # still, one can integrate out f_sp samples in para_store
          # adding contribution form each selected protein
          if (!is.null(X_pred)){
            X_pred_p_non_na <- !is.na(X_pred[,p]) # identify proteins with non-na affinity in X_pred
            if (any(total_non_na_id) && any(X_pred_p_non_na)){ # if we have observed some data from the training set and it also appears in 
              # the testing set. If either is false, the contribution of this p-s pair to the test point should be 0
              X_pred_p_relative_position <- X_pred[X_pred_p_non_na, p]
              K_pred_star_star <- as.matrix(kernel_matrix_gp_0(X_pred_p_relative_position, nu))
              K_pred_star <- nu[1]*exp(sapply(X_pred_p_relative_position,
                                              FUN = function(x){-((X_p_relative_position-x)^2)/(2*nu[2])}))
              if (is.vector(K_pred_star)){
                K_pred_star <- as.matrix(K_pred_star)
              }
              else{
                K_pred_star <- t(K_pred_star)
              }
              
              pred_GP_mean_p <- K_pred_star %*% C %*% y_sub_p
              pred_GP_var_p <- diag(K_pred_star_star - K_pred_star %*% C %*% t(K_pred_star))
              
              y_pred_hat[X_pred_p_non_na, s, r] <- y_pred_hat[X_pred_p_non_na, s, r] + pred_GP_mean_p
              y_pred_var[X_pred_p_non_na, s, r] <- y_pred_var[X_pred_p_non_na, s, r] + pred_GP_var_p
            }
          }
        }
      }
      intercept_matrix[r, s] <- intercept_store[s]
      # y_pred_hat_non_zero <- y_pred_hat[, s, r] != 0
      # y_pred_hat[y_pred_hat_non_zero, s, r] <- y_pred_hat[y_pred_hat_non_zero, s, r] + intercept_store[s]
      y_pred_hat[, s, r] <- y_pred_hat[, s, r] + intercept_store[s]
    }
    y_hat[,,r] <- t(t(apply(para_store, MARGIN = c(1,3), sum)) + intercept_store)
  }
  return(list('curve_mean'=curve_matrix, 
              'curve_sd'=curve_var_matrix, 
              'intercept'=intercept_matrix,
              'y_hat'=y_hat,
              'y_pred_hat'=y_pred_hat,
              'y_pred_var'=y_pred_var))
}


# also need a function taking y_hat as input and return WAIC
my_WAIC_gp_0 <- function(Y, yhat, my_sigma){
  non_NA_id <- !is.na(Y)
  my_lkd_gp_0_mtrx <- matrix(NA, ncol = length(my_sigma), nrow = sum(non_NA_id))
  for (r in 1:length(my_sigma)){
    my_yhat <- yhat[,,r]
    my_lkd_gp_0_mtrx[, r] <- dnorm(Y[non_NA_id], mean = my_yhat[non_NA_id], sd=sqrt(my_sigma[r]), log = TRUE)
  }
  return(WAIC(my_lkd_gp_0_mtrx))
}


# No lets stick with 3-fold CV!
choose_nu_0 <- function(N, burn_in, thin, X, Y, nu, inclusion_prob=.1, back_fit_iter=20, k_fold=3){
  vec_nu_1 <- nu[1]
  vec_nu_2 <- nu[2]
  cv_mse <- rep(NA, length(vec_nu_1)*length(vec_nu_2))
  cv_selected_protein <- lapply(1:length(vec_nu_1)*length(vec_nu_2), function(x){NA})
  store_WAIC <- rep(NA, length(vec_nu_1)*length(vec_nu_2))
  store_R2 <- rep(NA, length(vec_nu_1)*length(vec_nu_2))
  nu_vec <- cbind(rep(vec_nu_1, each=length(vec_nu_2)), 
                  rep(vec_nu_2, times=length(vec_nu_1))) # each row is a para pair
  
  for (j in 1:(dim(nu_vec)[1])){
    nu <- nu_vec[j, ]
    sample_size <- dim(X)[1]
    random_partition <- sample(1:sample_size, size=sample_size, replace = FALSE)
    group_id <- split(random_partition, sort(1:sample_size %% k_fold))
    
    record_mse <- matrix(NA, nrow = dim(Y)[1], ncol = dim(Y)[2])
    for (fold in 1:k_fold){
      my_group_id <- group_id[[fold]]
      my_gibbs <- my_gibbs_sampler_gp_0(N, burn_in, thin, X[-my_group_id, ], Y[-my_group_id, ], inclusion_prob, nu)
      my_para <- GP_para_0(back_fit_iter, X[-my_group_id, ], Y[-my_group_id, ], 
                           my_gibbs$Z, my_gibbs$sigma2, my_gibbs$gamma2, nu, X_pred = X[my_group_id, ])
      
      my_mse <- array(NA, dim = c(length(my_gibbs$'sigma2'), length(my_group_id), dim(Y)[2])) # this is RD'S
      for (i in 1:dim(my_mse)[1]){
        my_fit <- my_para$'y_pred_hat'[,,i] # recall it has dim D'SR
        my_mse[i,,] <- (Y[my_group_id, ] - my_fit)^2
        test_na <- !is.na(as.vector(Y[my_group_id, ]))
        # plot(as.vector(Y[my_group_id, ])[test_na], as.vector(my_fit)[test_na])
        # abline(a=0, b=1)
      }
      record_mse[my_group_id, ] <- apply(my_mse, c(2,3), mean, na.rm=TRUE)
    }
    cv_mse[j] <- mean(record_mse, na.rm=TRUE)
    print(record_mse)
    
    test_full_MCMC <- my_gibbs_sampler_gp_0(N, burn_in, thin, X, Y, inclusion_prob, nu)
    test_full_para <- GP_para_0(back_fit_iter, X, Y, 
                                test_full_MCMC$Z, test_full_MCMC$sigma2, test_full_MCMC$gamma2, nu)  
    inclusion_mean <- colMeans(test_full_MCMC$'Z')
    cv_selected_protein[[j]] <- colnames(X)[inclusion_mean>0.5]
    
    my_sigma <- test_full_MCMC$'sigma2'
    print(my_sigma)
    my_lkd_gp_0 <- matrix(NA, nrow = dim(Y)[1]*dim(Y)[2], ncol = length(my_sigma))
    my_R2 <- rep(NA, length(my_sigma))
    for (i in 1:length(my_sigma)){
      my_lkd_gp_0[, i] <- as.vector(dnorm(Y, mean = test_full_para$'y_hat'[,,i], sd = sqrt(my_sigma[i]), log = TRUE))
      my_R2[i] <- cor(Y[!is.na(Y)], (test_full_para$'y_hat'[,,i])[!is.na(Y)])^2
    }    
    store_WAIC[j] <- WAIC(my_lkd_gp_0[!is.na(my_lkd_gp_0[,1]), ])$WAIC
    store_R2[j] <- mean(my_R2) 
  }
  return(list('nu'=nu_vec, 'mse'=cv_mse, 'WAIC'=store_WAIC, 
              'R2'=store_R2, 'selected'=cv_selected_protein))
}


CV_nu_par_0 <- function(N, burn_in, thin, X, Y, vec_nu_1, vec_nu_2, 
                        inclusion_prob=.1, back_fit_iter=20, k_fold=3, ncor=2){
  nu_vec <- cbind(rep(vec_nu_1, each=length(vec_nu_2)), 
                  rep(vec_nu_2, times=length(vec_nu_1))) # each row is a para pair
  nu_list <- lapply(seq_len(nrow(nu_vec)), function(i) nu_vec[i,])
  cl <- makeCluster(ncor)
  registerDoParallel(cl)
  # start from here!
  res <- foreach(nu=nu_list, .packages = 'LaplacesDemon', 
                 .export = c('choose_nu_0', 'GP_para_0', 'my_gibbs_sampler_gp_0', 
                             'y_lkd_gp_0', 'my_kernels_gp_0', 'kernel_matrix_gp_0', 
                             'log_prior_z_gp_0', 'log_prior_gamma2_gp_0', 
                             'log_prior_sigma2_gp_0')) %dopar% 
    choose_nu_0(N=N, burn_in=burn_in, thin=thin, X=X, Y=Y, nu=nu, inclusion_prob=inclusion_prob, 
                back_fit_iter=back_fit_iter, k_fold=k_fold)
  stopCluster(cl)
  
  
  pool_nu <- nu_vec
  pool_mse <- rep(NA, length(nu_list))
  pool_WAIC <- rep(NA, length(nu_list))
  pool_R2 <- rep(NA, length(nu_list))
  pool_cv_selected_protein <- lapply(1:length(nu_list), function(x){NA})
  for (i in 1:length(nu_list)){
    pool_mse[i] <- res[[i]]$'mse'
    pool_WAIC[i] <- res[[i]]$'WAIC'
    pool_R2[i] <- res[[i]]$'R2'
    pool_cv_selected_protein[[i]] <- res[[i]]$'selected'[[1]]
  }
  return(list('c0'=nu_list, 'mse'=pool_mse, 'WAIC'=pool_WAIC, 
              'R2'=pool_R2, 'selected'=pool_cv_selected_protein))
  
}



# this chuck of code is taken directly from DepInfeR
MultiLasso <- function(TargetMatrix, ResponseMatrix, repeats = 100, BPPARAM = bpparam(), lambda = "lambda.min"){
  runGlm.multi <- function(i, X, y, folds = 3, lambda = "lambda.min", 
                           standardize = FALSE) {
    res <- cv.glmnet(X, y, family = "mgaussian", nfolds = folds, 
                     alpha = 1, standardize = standardize)
    res
  }
  results <- bplapply(seq(repeats), runGlm.multi, TargetMatrix, 
                      ResponseMatrix, BPPARAM = BPPARAM)
  modelList <- list()
  lambdaList <- rep(NA, length(results))
  varExplain.all <- rep(NA, length(results))
  varExplain.cv <- rep(NA, length(results))
  coefTab <- data.frame(a = rep(colnames(TargetMatrix), times = ncol(ResponseMatrix)), 
                        b = rep(colnames(ResponseMatrix), each = ncol(TargetMatrix)))
  intercept_store <- matrix(NA, nrow = length(results), ncol=dim(ResponseMatrix)[2])
  
  for (i in seq(length(results))) {
    res <- results[[i]]
    lambdaList[i] <- res[[lambda]]
    coefModel_original <- glmnet::coef.glmnet(res, s = lambda)
    coefModel <- as.vector(Reduce(cbind, coefModel_original)[-1, ])
    coefTab[[paste0("r", i)]] <- coefModel
    y.pred <- predict(res, s = lambda, newx = TargetMatrix)
    intercept_store[i,] <- Reduce(cbind, coefModel_original)[1, ]
    varExp <- stats::cor(as.vector(ResponseMatrix), as.vector(y.pred))^2
    varExplain.all[i] <- varExp
  }
  
  resMat <- -as.matrix(coefTab[, !colnames(coefTab) %in% c("a", 
                                                           "b")])
  coefTab[["freq"]] <- rowSums(!resMat == 0)/ncol(resMat)
  coefTab[["med"]] <- rowMedians(resMat)
  freqMat <- coefMat <- matrix(data = NA, nrow = ncol(TargetMatrix), ncol = ncol(ResponseMatrix), 
                               dimnames = list(colnames(TargetMatrix), colnames(ResponseMatrix)))
  for (eachCol in colnames(freqMat)) {
    freqMat[, eachCol] <- coefTab[coefTab$b %in% eachCol, 
    ]$freq
    coefMat[, eachCol] <- coefTab[coefTab$b %in% eachCol, 
    ]$med
  }
  return(list(coefMat = coefMat, freqMat = freqMat, lambdaList = lambdaList, 
       varExplain.all = varExplain.all, inputX = TargetMatrix, inputY = ResponseMatrix,
       intercept = colMedians(intercept_store)))
}

int_over_union <- function(s1, s2){
  dist <- rep(NA, length(s1))
  for (i in 1:length(s1)){
    dist[i] <- length(intersect(s1[[i]], s2))/length(union(s1[[i]], s2))
  }
  return(dist)
}


library(randomForestSRC)
rf_test <- function(n_tree, X, Y){
  k_fold <- 3
  sample_size <- dim(Y)[1]
  random_partition <- sample(1:sample_size, size=sample_size, replace = FALSE)
  group_id <- split(random_partition, sort(1:sample_size %% k_fold))
  response_mat <- matrix(0, ncol=dim(Y)[2], nrow=dim(Y)[1])  
  for (k in 1:k_fold){
    test <- rfsrc(get.mv.formula(colnames(Y)), 
                  data=as.data.frame(cbind(Y, X)[-group_id[[k]], ]), 
                  importance = 'permute')
    test_pred <- predict.rfsrc(test, newdata = as.data.frame(X[group_id[[k]] ,]))
    response_mat[group_id[[k]], ] <- get.mv.predicted(test_pred)
  }
  return(response_mat)
}