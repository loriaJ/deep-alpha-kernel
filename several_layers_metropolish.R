

ralpha_stable_true0 <- function(n,alpha,beta=1,scale=1){
  phi_0 <- runif(n,-pi/2,pi/2) # ? 
  E_i <- rexp(n)
  if(alpha == 1){
    scale * ( 2/pi * 
                ( (pi/2 + beta * phi_0)*tan(phi_0) - 
                    beta*log(E_i*cos(phi_0)/(pi/2 + beta*phi_0) ) )  )  + 
      (2/pi) * scale * beta *log(scale)
    
  } else if(beta == 1) { 
    b_0 <- atan(beta*tanpi(alpha/2))/alpha
    exp(-(1/alpha)*log(cos(phi_0)) + log(sin(alpha*phi_0+alpha*b_0)) +
          ((1-alpha)/alpha)*(log(cos(phi_0 -alpha*phi_0 -alpha*b_0 ))-log(E_i)) + 
          log(scale)/alpha)
    
  } else {
    b_0 <- atan(beta*tanpi(alpha/2))/alpha
    sin(alpha*phi_0+alpha*b_0)/cos(phi_0)^(1/alpha) *
      (cos(phi_0 -alpha*phi_0 -alpha*b_0 )/E_i)^((1-alpha)/alpha) *
      scale^(1/alpha)
  }
  
}

compute_sigma_from_A <- function(A_L,K=function(x)exp(-x)){
  n <- nrow(A_L)
  sigma0 <- matrix(NA,n,n)
  for(i in 1:n){
    for(j in i:n){
      R_entry <- A_L[i,i] + A_L[j,j] - 2 * A_L[i,j]
      if(R_entry < 0){
        print('problems in R entries')
      }
      sigma0[i,j] <- K(R_entry)
      if(sigma0[i,j] < 0){
        print('problems in sigma0 entries')
      }
      sigma0[j,i] <- sigma0[i,j]
    }
  }
  return(sigma0)
}


update_sigma_ells <- function(curr_new_sigma_ell,
                              old_scales_list,
                              curr_ell,Sigma_list,
                              K=function(s){exp(-s)},
                              sigma_offset=1){
  L <- length(old_scales_list)
  new_Sigma_List <- Sigma_list
  # new_Sigma_List <- sigma_offset
  new_Sigma_List[[curr_ell]] <- curr_new_sigma_ell
  for(ell in (curr_ell+1):(L+1)){
    # a_ell <- 
    new_Sigma_List[[ell]] <- compute_sigma_from_A(
      old_scales_list[[ell-1]] * new_Sigma_List[[ell-1]],
      K)
  }
  diag(new_Sigma_List[[L+1]]) <- diag(new_Sigma_List[[L+1]]) + sigma_offset
  return(new_Sigma_List)
}

# list 1 and 2 are lists with the same length: L
# want to put all the entries after (inclusive) index_from of list 2
# zip_lists <- function(list_1,list_2,index_from){
#   
# }



mcmc_sample <- function(scales_0,layer_scales,alpha_2, x,y,
                        A_0,Sigma_list,
                        k_0=function(s){1/(1+s^2)},
                        K=function(s){exp(-s)},
                        sigma_offset=1,
                        anisotropic = FALSE,
                        delta=1){
  n <- length(y)
  p <- ncol(x)
  L <- length(layer_scales)
  prop_0 <- ralpha_stable_true0(p,alpha=alpha_2,beta = 1)
  # ell == 1
  prop_temp <- scales_0
  Sigma_L_old <- Sigma_list[[L+1]]
  log_p_old <- mvtnorm::dmvnorm(y,sigma = Sigma_L_old,log = T,
                                checkSymmetry = F)
  curr_old_A_ell <- A_0
  log_proposals0 <- rep(NA,p)
  # first layer
  for(j in 1:p){
    # print(j)
    prop_temp[j] <- prop_0[j]
    # all(t((p:1) * t(x)) %*% t(x) == (x) %*% diag(p:1) %*% t(x))
    curr_new_A_ell <- t(prop_temp * t(x)) %*% t(x)/p
      # (x) %*% diag(prop_temp,p) %*%t(x)/p
    if(anisotropic){
      updated_sigma_ells <- update_sigma_list(
        A_0 = curr_new_A_ell,
        delta=delta,
        s_ell = layer_scales)
      diag(updated_sigma_ells[[L+1]]) <- 
        diag(updated_sigma_ells[[L+1]]) + sigma_offset
    } else {
      
      curr_new_sigma_ell <- compute_sigma_from_A(curr_new_A_ell,K = k_0)
      updated_sigma_ells <- update_sigma_ells(Sigma_list = Sigma_list,
                                              old_scales_list = layer_scales,
                                              curr_ell = 1,K=K,
                                              sigma_offset = sigma_offset,
                                              curr_new_sigma_ell = curr_new_sigma_ell)
    }
    Sigma_L_star <- updated_sigma_ells[[L+1]]
    log_p_new <- mvtnorm::dmvnorm(y,sigma = Sigma_L_star,
                                  log = T,
                                  checkSymmetry = F)
    log_u0 <- log(runif(1))
    log_proposals0[j] <- log_p_new
    # print(log_p_new)
    # print(log_p_old)
    
    if((!is.na(log_p_new - log_p_old)) & (log_u0 < (log_p_new - log_p_old))){
      ## accept
      Sigma_L_old <- Sigma_L_star
      Sigma_list <- updated_sigma_ells
      curr_old_A_ell <- curr_new_A_ell
      log_p_old <- log_p_new
    } else {
      ## reject
      # print('reject')
      prop_temp[j] <- scales_0[j]
    }
  }
  # old_scales_list[[1]] <- prop_temp 
  scales_new <- layer_scales
  for(ell in 1:(L)){
    # print(ell)
    prop_ell <- ralpha_stable_true0(1,alpha=alpha_2,beta = 1)
    scales_new[ell] <- prop_ell
    
    if(anisotropic){
      new_sigmas <- update_anisotropic(new_sigma_ell = Sigma_list[[ell]],
                                       scales_list = scales_new,
                                       curr_ell=ell,
                                       delta = delta,
                                       Sigma_list = Sigma_list)
      diag(new_sigmas[[L+1]]) <- diag(new_sigmas[[L+1]]) + sigma_offset
      
    } else {
      
      new_sigmas <- update_sigma_ells(curr_new_sigma_ell = Sigma_list[[ell]],
                                      old_scales_list = scales_new,
                                      sigma_offset = sigma_offset,
                                      curr_ell = ell,
                                      Sigma_list = Sigma_list,
                                      K = K)
    }
    Sigma_L_star <- new_sigmas[[L+1]]
    log_p_new <- mvtnorm::dmvnorm(y,sigma = Sigma_L_star,
                                  checkSymmetry = F,log = T)
    log_u0 <- log(runif(1))
    if(is.na(log_p_new - log_p_old)){
      print('got NA')
    }
    # print((log_p_new - log_p_old))
    if((!is.na(log_p_new - log_p_old)) & (log_u0 < (log_p_new - log_p_old))){
      ## accept
      # print(c('accept',ell))
      Sigma_L_old <- Sigma_L_star
      Sigma_list <- new_sigmas
      # curr_old_A_ell <- curr_new_A_ell
      log_p_old <- log_p_new
      # 
    } else {
      ## reject
      # print(c('reject',ell))
      scales_new[ell] <- layer_scales[ell]
    }
  }
  
  return(list(prop_temp,
              scales_new,# scales, sigmas, log_p
              Sigma_list,log_p_old,log_proposals0
              ))
}

k_0 <- function(s){1/(1+s^2)}
K <- function(s){exp(-s)}
