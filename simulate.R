rm(list=ls(all=TRUE))
library(foreach)
library(RBesT)
require(randomizr)
require(extraDistr)

setwd("C:\\Users\\willi\\OneDrive - Yale University\\Desktop\\research\\hybrid platform design\\resubmission\\code")
source("mem functions.R")

# sd of historical controls
v_H <- 0.025
##----------------------------------------------------
##  current trial:
##----------------------------------------------------
#sample size of each arm
arm_size <- c(30, 60, 60) 
# number of arms
A <- length(arm_size)
# total sample size
N <- sum(arm_size)
# response rate with the 1st element represent control
arm_rate1 <- c(0.15, 0.30, 0.40)
arm_rate0 <- c(0.15, 0.15, 0.15)
arm_name <- c("Ctr", paste("T", 1:(A-1), sep=""))

##----------------------------------------------------
##  simulate internal data:
##----------------------------------------------------
nsim <- 1000

dat <- dat0 <- list()
#mat <- matrix(NA, nrow = 3, ncol = nsim)
## simulation
for(i in 1:nsim){
  arm_code <- as.vector(complete_ra(N = N, m_each = arm_size,
                                    conditions = 1:A))
  pvec <- arm_rate1[arm_code]
  set.seed(i)
  resp <- rbinom(N, size = 1, prob = pvec)
  
  dat[[i]] <- data.frame(arm_code, resp)
  
  #p1 <- mean(resp[arm_code==1])
  #p2 <- mean(resp[arm_code==2])
  #p3 <- mean(resp[arm_code==3])
  #mat[, i] <- c(p1, p2, p3)
}
#rowMeans(mat)
#apply(mat, 1, sd)

for(i in 1:nsim){
  set.seed(i)
  arm_code <- as.vector(complete_ra(N = N, m_each = arm_size,
                                    conditions = 1:A))
  pvec <- arm_rate0[arm_code]
  resp <- rbinom(N, size = 1, prob = pvec)
  
  dat0[[i]] <- data.frame(arm_code, resp)
}

D <- c(dat0, dat)

##----------------------------------------------------
##  simulate historical data:
##----------------------------------------------------
# number of historical trials
H <- 3
# sample size of historical trials
n_H <- rep(120, H)

## scenarios of historical trials
rate_diff <- c(-0.10, -0.05, 0, 0.05, 0.1, 0.15, 0.2)
## number of scenarios
S <- length(rate_diff)

hdat <- matrix(NA, nrow = S, ncol = H)
trial_H <- map_mix <- map_rob <- list()

for(s in 1:S){
  if(v_H>0){
    rate_H <- rtnorm(H, mean = arm_rate1[1] + rate_diff[s], sd = v_H, a = 0, b = 1)
  }
  if(v_H==0){
    rate_H <- arm_rate1[1] + rate_diff[s]
  }
  
  ## simulate historical trial data
  hdat[s, ] <- rbinom(H, size = n_H, prob = rate_H)
  
  htrial <- data.frame("trial" = 1:H, "x" = hdat[s, ], "n" = n_H)
  htrial$y <- htrial$n - htrial$x
  
  options(RBesT.MC.control=list(adapt_delta=0.999))
  set.seed(5)
  map_mc <- gMAP(cbind(htrial$x, htrial$n) ~ 1 | htrial$trial, 
                 family = binomial, tau.dist = "HalfNormal", 
                 tau.prior = 1, beta.prior = 2)
  map_mix[[s]] <- automixfit(map_mc, Nc = 2)
  map_rob[[s]] <- robustify(map_mix[[s]], weight = 0.5, mean = 1/2, n = 1)
  
  #target_N <- arm_size[1] + min(arm_size[2] - arm_size[1], ess(map_mix[[s]]))
  target_N <- arm_size[1] + min(arm_size[2] - arm_size[1], 30)
  print(target_N)

  trial_H[[s]] <- htrial
}

mem_prior <- find.mem.prior(n_c = arm_size[1], hdat = trial_H[[which(rate_diff==0)]], 
                               target_N = 60)

maxEN(n_c = arm_size[1], hdat =trial_H[[3]], prior = mem_prior)
save.image(paste("simulation v=", v_H, ".RData", sep =""))

