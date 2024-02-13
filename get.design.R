
get.design <- function(a = seq(0.5, 1, 0.1), 
                       lam = seq(0.95, 0.99, 0.001), nsim, t, pp, arm_size,
                       target_fwer = 0.10, n_ia){
  A <- length(arm_size)
  N <- matrix(rep(arm_size[-1] + arm_size[1], length(t)), nrow = length(t))
  n_trial <- rep(NA, nsim)
  
  if(length(t)>1){
    ## parameters defining the stopping boundaries
    mat <- expand.grid(lam, a)
    names(mat) <- c("lam", "a")
    fwer <- pow <- NULL
    for(k in 1:nrow(mat)){
      rej <- NULL
      for(i in 1:nsim){
        n_ctr <- matrix(rep(n_ia[, i, 1], A-1), ncol = A-1)
        n_trt <- n_ia[, i, -1]
        ## decision cutoff value
        ct <- as.numeric(mat$lam[k])*((n_ctr+n_trt)/N)^as.numeric(mat$a[k])
        #ct <- as.numeric(mat$lam[k])*(t/N)^as.numeric(mat$a[k])
        rej <- rbind(rej, apply(pp[, i, -1]>ct, 2, all))
        
        ## go/no decisions:
        dec <- (pp[, i, -1]<=ct)
        idx <- nstop <- NULL
        for(j in 1:(A-1)){
          if(any(dec[, j])){
            idx[j] <- min(which(dec[, j]))
          }else{
            idx[j] <- length(t)
          }
          nstop[j] <- n_trt[idx[j], j]
        }
        n_trial[i] <- sum(c(n_ctr[max(idx), 1], nstop))
      }
      fwer <- c(fwer,mean(rowSums(rej)>0))
      pow <- rbind(pow, colMeans(rej))
    }
  }
  
  if(length(t)==1){
    fwer <- pow <- NULL
    for(k in 1:length(lam)){
      rej <- NULL
      for(i in 1:nsim){
        rej <- rbind(rej, pp[i, -1]>lam[k])
      }
      fwer <- c(fwer, mean(rowSums(rej)>0))
      pow <- rbind(pow, colMeans(rej))
    }
  }
  
  fwer2 <- ifelse(fwer <= target_fwer, fwer-target_fwer, 1000)
  idx <- which.min(abs(fwer2))
  
  if(length(t)>1){
    return(c(mat[idx, ], fwer[idx], pow[idx, ], mean(n_trial)))
  }
  if(length(t)==1){
    return(c(lam[idx], NA, fwer[idx], pow[idx, ]))
  }
}

