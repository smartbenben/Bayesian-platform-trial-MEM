rm(list=ls(all=TRUE))
library(foreach)
library(RBesT)

setwd("C:\\Users\\willi\\OneDrive - Yale University\\Desktop\\research\\hybrid platform design\\resubmission\\code")

source("mem functions.R")
source("cal posterior.R")
source("get.design.R")

v_H <- 0.025

load(paste("simulation v=", v_H, ".RData", sep =""))
##----------------------------------------------------
##  analysis:
##----------------------------------------------------
# time of interim analysis
#t_analysis <- seq(90, sum(arm_size), 30)
t_analysis <- c(90, sum(arm_size))
##----------------------------------------------------
##  run analysese through scenarios
##----------------------------------------------------
J <- length(t_analysis)

pp_fixed <- array(NA, dim = c(S, length(t_analysis), 2*nsim, A))
pp_map <- array(NA, dim = c(S, length(t_analysis), 2*nsim, A))
pp_rob <- array(NA, dim = c(S, length(t_analysis), 2*nsim, A))
pp_mem <- array(NA, dim = c(S, length(t_analysis), 2*nsim, A))

n_array <- array(NA, dim = c(S, length(t_analysis), 2*nsim, A))

beta_a <- beta_b <- 1

for(s in 1:S){
  fit0 <- set.mem.prior(num_study = H+1, delta = mem_prior)
  
  for(i in 1:(2*nsim)){
    response <- D[[i]]$resp
    arm <- D[[i]]$arm_code
    
    for(j in 1:J){
      ##interim analysis data
      ia_dta <- response[1:t_analysis[j]]
      ia_arm <- arm[1:t_analysis[j]]
      
      #number of responses in each arm:
      xvec <- aggregate(ia_dta, by = list(ia_arm), sum)$x
      nvec <- aggregate(ia_dta, by = list(ia_arm), length)$x
      yvec <- nvec - xvec
      #save interim sample sizes for treatment arms
      n_array[s, j, i, ] <- nvec
      
      map_post <- data.frame(t(postmix(map_mix[[s]], n = nvec[1], r = xvec[1])[, ]))
      names(map_post) <- c("p", "alpha", "beta")
      
      rob_post <- data.frame(t(postmix(map_rob[[s]], n = nvec[1], r = xvec[1])[, ]))
      names(rob_post) <- c("p", "alpha", "beta")
      
      fit <- update.part.bin(x = c(trial_H[[s]]$x, xvec[1]), 
                             n = c(trial_H[[s]]$n, nvec[1]), 
                             prior_part = fit0$prior, 
                             part = fit0$part)
      #rstar <- ifelse(fit$part_hat[1,1:H]==fit$part_hat[1,1+H], 1, 0)*fit$phat
      
      rstar <- fit$post_sim
      astar <- beta_a + sum(rstar*trial_H[[s]]$x) + xvec[1]
      bstar <- beta_b + sum(rstar*trial_H[[s]]$y) + yvec[1]
      
      for(a in 2:A){
        ## fixed prior, no borrowing
        pp_fixed[s, j, i, a] <- beta.ineq(a = beta_a+xvec[a], b = beta_b+yvec[a], 
                          c = beta_a+xvec[1], d = beta_b+yvec[1], 
                          delta = 0)
        #map
        pp_map[s, j, i, a] <- sum(apply(map_post, 1, function(x){
          x[1]*beta.ineq(a = beta_a+xvec[a], b = beta_b+yvec[a], 
                         c = x[2], d = x[3], delta = 0)
        }))
        
        #rmap
        pp_rob[s, j, i, a] <- sum(apply(rob_post, 1, function(x){
          x[1]*beta.ineq(a = beta_a+xvec[a], b = beta_b+yvec[a], 
                         c = x[2], d = x[3], delta = 0)
        }))
        #mem
        pp_mem[s, j, i, a] <- beta.ineq(a = beta_a+xvec[a], b = beta_b+yvec[a], 
                                 c = astar, d = bstar, 
                                 delta = 0)
      }
    }
  }
}

pp0_list <- list()
pp1_list <- list()
pp0_list[[1]] <- pp_fixed[ , ,1:nsim, ]
pp0_list[[2]] <- pp_map[ , ,1:nsim, ]
pp0_list[[3]] <- pp_rob[ , ,1:nsim, ]
pp0_list[[4]] <- pp_mem[ , ,1:nsim, ]
pp1_list[[1]] <- pp_fixed[ , ,(1+nsim):(2*nsim), ]
pp1_list[[2]] <- pp_map[ , ,(1+nsim):(2*nsim), ]
pp1_list[[3]] <- pp_rob[ , ,(1+nsim):(2*nsim), ]
pp1_list[[4]] <- pp_mem[ , ,(1+nsim):(2*nsim), ]

save.image(paste("analysis v=", v_H, ".RData", sep =""))
