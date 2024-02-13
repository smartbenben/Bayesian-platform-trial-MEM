require("partitions")

# R: the number of trials to consider
# max_cl: maximum number of blocks
get.part <- function(R, max_cl){
  ## generate all the possible partitions 
  ## and store them in a matrix 
  part_mat <- t(setparts(R))
  part_mat <- part_mat[apply(part_mat, 1, function(x){
    length(unique(x))<=max_cl
  }), ]
  part_mat <- data.frame(part_mat)
  names(part_mat) <- LETTERS[1:R]
  return(part_mat)
}

#enumerate all possible partitions of studies
#and calculate the uncertainties related to each partition

##x: vector of responses in each study
##n: vector of sample sizes in each study
##beta prior for response rate: a0, b0 [default: a0 = b0 = 1]
update.part.bin <- function(x, n, prior_part, part, a0 = 1, b0 = 1){
  R <- length(x)
  K <- nrow(part)
  
  p <- foreach(k = 1:K,.combine = "c")%do%{
    grp <- unlist(part[k, ])
    S <- aggregate(x, by=list(grp), sum)$x
    N <- aggregate(n, by=list(grp), sum)$x
    #calculate marginal probs m(s_j)
    prod((beta(a0+S, b0+N-S)/beta(a0,b0)))*prior_part[k]
  }
  post_part <- p/sum(p)
  
  idx <- which.max(post_part)
  
  sim_mat <- part[, 1:(R-1)]==part[, R]
  
  if(is.null(nrow(sim_mat))){
    post_sim <- sum(sim_mat*post_part)
  }else{
    post_sim <- colSums(sim_mat*post_part)
  }
  
  return(list("part_hat" = part[idx, ], "phat" = post_part[idx], 
              "post_part" = post_part,
              "post_sim" = post_sim))
}

## partitions of studies
## assume the last study is the current study
set.mem.prior <- function(num_study, delta){
  part <- get.part(R = num_study, max_cl = num_study)
  K <- nrow(part)
  # number of blocks/unique response rate in each partition
  n_bk <- apply(part, 1, function(x){
    length(unique(x))
  })
  prior <- n_bk^(delta)/sum(n_bk^(delta))
  
  sim_mat <- part[, 1:(num_study-1)]==part[, num_study]
  if(is.null(nrow(sim_mat))){
    prior_sim <- sum(sim_mat*prior)
  }else{
    prior_sim <- colSums(sim_mat*prior)
  }
  list("part" = part, "prior" = prior, "prior_sim" = prior_sim)
}
  

maxEN <- function(n_c, hdat, prior){
  H <- nrow(hdat)
  fit0 <- set.mem.prior(num_study = H+1, delta = prior)
  #enumerate possible outcomes for current data
  cdat <- expand.grid(0:n_c, n_c)
  names(cdat) <- c("x", "n")
  out <- foreach(k = 1:nrow(cdat),.combine="rbind")%do%{
    fit <- update.part.bin(x = c(hdat$x, cdat$x[k]), 
                           n = c(hdat$n, cdat$n[k]), 
                           prior_part = fit0$prior, 
                           part = fit0$part)
    wt <- fit$post_sim
    #wt <- ifelse(unlist(fit$part_hat[1:H])==unlist(fit$part_hat[1+H]), 1, 0)*fit$phat
    astar <- 1 + sum(wt*hdat$x) + cdat$x[k]
    bstar <- 1 + sum(wt*hdat$y) + cdat$n[k]-cdat$x[k]
    
    c(astar, bstar, astar+bstar)
  }
  out <- data.frame(out)
  rownames(out) <- NULL
  names(out) <- c("alpha", "beta", "N")
  
  max(out$N)
}

find.mem.prior <- function(n_c, hdat, target_N){
  f <- function(x){
    maxEN(n_c = n_c, hdat = hdat, prior = x) - target_N
  }
  uniroot(f, c(-100, 100))$root
}

