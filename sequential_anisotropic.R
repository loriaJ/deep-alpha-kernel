
J_delta <- function(delta,theta){
  # u0 <- 
  if(delta == 0){
    pi-theta
  } else if(delta == 1){
    sin(theta) + (pi-theta)*cos(theta)
  } else {
    theta <- max(theta,1e-4) # to handle the limit
    # if(theta>0.999999){
    # print(theta)
    # }
    ctheta <- cos(theta) ## should be equal to x * y/(||x|| * ||y||)
    st <- sin(theta)
    f0 <- function(u){
      cu <- cos(u)
      (st^2)^delta * st*cu^delta/((1-cu*ctheta)^(delta+1))
      # exp(
      # delta * log(cu) - 
      #   (delta + 1) * log(1 - cu * ctheta))
    }
    i0 <- integrate(f0,lower = 0,upper = pi/2,rel.tol = 1e-5)
    
    gamma(delta+1)*(((st)^(2))^delta) *st * #sign(st)*
      i0$value # * pi/2 
    ## need to multiply by pi/2 because of sampling
  }
}


seq_kernels_anisotrop <- function(delta,Sigma_l,s_l){
  n <- dim(Sigma_l)[1]
  angle_matrix <- matrix(NA,n,n)
  k_matrix <- matrix(NA,n,n)
  for(i in 1:n){
    for(j in i:n){
      # adding the +1 to the sigma[i,j]  acts as if we had a 
      # bias term in each layer, and makes computations
      # much more numerically stable
      
      pre_ratio <- sqrt((1+s_l*Sigma_l[i,i])*(1+s_l*Sigma_l[j,j]))
      ratio <- (1+s_l*Sigma_l[i,j])/pre_ratio
      # the min is used to prevent small round-off issues
      angle_matrix[i,j] <- acos(min(1,ratio)) 
      J0 <- J_delta(delta = delta,
                    theta = angle_matrix[i,j])
      if(is.na(angle_matrix[i,j]) & (i!= j)){
        print(c('help',J0))
      } 
      if(i == j){
        k_matrix[i,j] <-  (1+s_l*Sigma_l[i,i])^(delta) #* J0/pi
      } else {
        
      k_matrix[i,j] <- 1/pi * (pre_ratio)^(delta) * J0
      k_matrix[j,i] <- k_matrix[i,j]
      }
    }
  }
  return(list(k_matrix,angle_matrix))
}

update_sigma_list <- function(A_0,delta,s_ell){
  L <- length(s_ell)
  sigma_list0 <- vector(mode='list',length = L+1)
  sigma_list0[[1]] <- A_0
  for(ell in 1:L){
    sigma_list0[[ell+1]] <- seq_kernels_anisotrop(
      delta=delta,
      sigma_list0[[ell]],
      s_l = s_ell[ell]
    )[[1]]
  }
  return(sigma_list0)
}

update_anisotropic <- function(new_sigma_ell, #= Sigma_list[[ell]],
                   scales_list,#  = scales_new,
                   curr_ell,#=ell,
                   delta,# = delta,
                   Sigma_list){#  = Sigma_list
  L <- length(scales_list)
  new_Sigma_list <- Sigma_list
  new_Sigma_list[[curr_ell]] <- new_sigma_ell
  for(ell in (curr_ell+1):(L+1)){
    # print(scales_list[ell-1])
    new_Sigma_list[[ell]] <- seq_kernels_anisotrop(delta = delta,
                          Sigma_l = new_Sigma_list[[ell-1]],
                          s_l = scales_list[ell-1])[[1]]
  }
  
  return(new_Sigma_list)
} 

compute_sigmas_anisotropic <- function(x_all,scales_0,
                                       s_ell,delta,p){
  A0 <- t(scales_0 * t(x_all)) %*% t(x_all)/p
  sigmas <- update_sigma_list(A_0 = A0,delta = delta,s_ell = s_ell)
  return(sigmas)
}

