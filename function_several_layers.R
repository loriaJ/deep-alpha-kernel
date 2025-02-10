source('several_layers_metropolish.R')
source('sequential_anisotropic.R')
stable_KP <- function(x0, y0, x_new, n_iters=3000, alpha2=1,sigma_offset=1,
                      sigma_known = FALSE,n_layers = 2,sd_offset = NA,
                      anisotropic=TRUE,
                      delta=1,save_sigmas_indic=F){
  alpha0 <- alpha2/2
  n <- nrow(x0)
  m <- nrow(x_new)
  # print(m)
  p <- ncol(x0)
  x_all <- rbind(x0,
                 x_new)
  old_scales <- vector(mode = 'list',length = n_layers)
  scales_orig <- ralpha_stable_true0(n=p,alpha = alpha0,beta = 1,scale = 1)
  if(save_sigmas_indic){
    save_sigmas <- array(NA,dim = c(n_iters,n_layers+1,n+m,n+m))
  }
  
  A0 <- t(scales_orig * t(x0)) %*% t(x0)/p

  old_scales <- ralpha_stable_true0(n=n_layers,alpha = alpha0,beta = 1,scale = 1)
  
  Sigma_list <- list()
  if(anisotropic){
    Sigma_list <- update_sigma_list(A_0 = A0,
                                    delta = delta,
                                    s_ell = old_scales)
      # print(length(Sigma_list))
      # print(n_layers+1)
      # print(length(old_scales))
    # seq_kernels_anisotrop()
  } else {
    Sigma_list[[1]] <- compute_sigma_from_A(A0,K = k_0)
    for(j in 2:(n_layers+1)){
      Sigma_list[[j]] <- compute_sigma_from_A(old_scales[j-1] * Sigma_list[[j-1]],K = K)
    }
  }
  # n_layers <- 4
  if(sigma_known){
    sigma_offset_l <- rep(sigma_offset,n_iters+1)
  } else {
    if(is.na(sd_offset)){
      sd_offset <- 1
    }
    sigma_offset_l <- rep(NA,n_iters+1)
    sigma_offset_l[1] <- abs(rcauchy(1,0,sd_offset))
  }
  # print('fails here')
  diag(Sigma_list[[n_layers+1]]) <- diag(Sigma_list[[n_layers+1]]) + sigma_offset_l[1]
  
  # print('started')
  pb <- txtProgressBar(max = n_iters)
  scales_0 <- matrix(NA,nrow = n_iters+1,ncol = p)
  scales_0[1,] <- scales_orig
  scales_layers <- matrix(NA,nrow = n_iters+1,ncol = n_layers) 
  scales_layers[1,] <- old_scales
  log_ps <- rep(NA,n_iters)
  y_star <- matrix(NA,nrow = n_iters,ncol = m)
  for(t1 in 1:n_iters){
    setTxtProgressBar(pb,t1)
    # print(scales_0[t1,])
    # print(sigma_offset_l[t1])
    samp1 <- mcmc_sample(scales_0 = scales_0[t1,],
                         layer_scales = scales_layers[t1,], 
                         alpha_2 = alpha0, x=x0,y=y0,
                         A_0=A0,
                         Sigma_list = Sigma_list,
                         k_0=k_0,
                         K=K,
                         sigma_offset=sigma_offset_l[t1],
                         anisotropic=anisotropic,
                         delta=delta)
    # print('sample')
    scales_0[t1+1,] <- samp1[[1]]
    scales_layers[t1+1,] <- samp1[[2]]
    Sigma_list <- samp1[[3]]# Sigma_list
    log_ps[t1] <- samp1[[4]]
    if(!sigma_known){
      #print(sigma_known)
      sigma_prop <- # 1/rgamma(1,1,1)
        abs(rcauchy(1,0,sd_offset))
      ### TODO: check how numerically unstable this is
      #print('')
      Sigma_L_prop <- Sigma_list[[n_layers+1]] + diag(sigma_prop - sigma_offset_l[t1],n)
      log_prop <- mvtnorm::dmvnorm(y0,sigma = Sigma_L_prop,
                                   log = T,checkSymmetry = F)
      log_u0 <- log(runif(1))
      
      if((!is.na(log_prop - log_ps[t1])) & (log_u0 < (log_prop - log_ps[t1]))){
        ## accept
        Sigma_list[[n_layers+1]] <- Sigma_L_prop
        sigma_offset_l[t1+1] <- sigma_prop
        log_ps[t1] <- log_prop
      } else {
        ## reject
        sigma_offset_l[t1+1] <- sigma_offset_l[t1]
        
      }
    }
    if(anisotropic){
      # tempA0 <- 
      sigma0s <- compute_sigmas_anisotropic(x_all = x_all,scales_0= scales_0[t1+1,],s_ell = scales_layers[t1+1,],
                                            delta = delta,p=p)
      diag(sigma0s[[n_layers+1]]) <- diag(sigma0s[[n_layers+1]]) + sigma_offset_l[t1+1]
    } else {
      
      sigma0s <- compute_sigmas_all(x_all = x_all,k_0 = k_0,K = K,
                                    scales0 = scales_0[t1+1,],
                                    scales_layers = scales_layers[t1+1,],
                                    sigma_offset = sigma_offset_l[t1+1],
                                    n_layers = n_layers,
                                    p=p)
    }
    sigma0 <- sigma0s[[n_layers+1]]
    sigma_nm <- sigma0[1:n,n+(1:m)]
    sigma_nn <- sigma0[1:n,1:n]
    sigma_mm <- sigma0[n+(1:m),n + (1:m)]
    mu <- t(sigma_nm) %*%solve(sigma_nn,y0)
    v1 <- sigma_mm - t(sigma_nm) %*% solve(sigma_nn,(sigma_nm))
    # print('pre-check')
    if(any(is.na(v1)|is.infinite(v1))){
      print('stop')
    }
    # ev1 <- eigen(v1)$values
    # if(any(is.na(ev1)) | any(ev1<1)){
    #   print('stop')
    # }
    y_star[t1,] <- mvtnorm::rmvnorm(n=1,mean = mu,
                                    sigma = v1,
                                    checkSymmetry = F)
    if(save_sigmas_indic){
      
      for(ell in 1:(n_layers+1)){
        save_sigmas[t1,ell,,] <- sigma0s[[ell]]
      }
    }
  }
  if(save_sigmas_indic){
    return(list(scales0 = scales_0,
                scales_layers = scales_layers,
                log_likelihood = log_ps,
                kernel_list = Sigma_list,
                offset = sigma_offset_l,
                y_new = y_star,
                save_sigmas = save_sigmas))
  } else{
    return(list(scales0 = scales_0,
                scales_layers = scales_layers,
                log_likelihood = log_ps,
                kernel_list = Sigma_list,
                offset = sigma_offset_l,
                y_new = y_star))
  }
  
  
}

compute_sigmas_all <- function(x_all,k_0,K,scales0,
                               scales_layers,
                               sigma_offset,n_layers,p){
  # all(t((p:1) * t(x)) %*% t(x) == (x) %*% diag(p:1) %*% t(x))
  A0 <- t(scales0 * t(x_all)) %*% t(x_all)/p
    # x_all %*% diag(scales0)%*% t(x_all)/p
  Sigma_list <- list()
  Sigma_list[[1]] <- compute_sigma_from_A(A0,K = k_0)
  # n_layers <- 4
  for(j in 2:(n_layers+1)){
    Sigma_list[[j]] <- compute_sigma_from_A(
      scales_layers[j-1] * Sigma_list[[j-1]],
      K = K)
  }
  
  diag(Sigma_list[[n_layers+1]]) <- diag(Sigma_list[[n_layers+1]])+ sigma_offset
  
  return(Sigma_list)
}




