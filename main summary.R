rm(list=ls(all=TRUE))

load("analysis v=0.025.RData")

##----------------------------------------------------
##  HY design:
##----------------------------------------------------
opt_par <- NULL

for(m in 1:4){
  design <- get.design(lam = seq(0.95, 0.99, 0.001), nsim = nsim, t = t_analysis, 
                       pp = pp0_list[[m]][which(rate_diff==0), , , ], 
                       arm_size = arm_size, target_fwer = 0.051,
                       n_ia = n_array[which(rate_diff==0), ,1:nsim, ])
  
  opt_par <- rbind(opt_par, design[1:5])
}
opt_par <- data.frame(opt_par)
names(opt_par) <- c("lam", "a", "fwp", "powA", "powB")
opt_par$met <- c("Fixed", "MAP", "rMAP", "MEM")

tab0 <- tab1 <- NULL

for(s in 1:S){
  for(m in 1:4){
    design <- get.design(lam = opt_par$lam[m], a = opt_par$a[m], nsim = nsim, 
                         t = t_analysis, pp = pp0_list[[m]][s, , , ], 
                         arm_size = arm_size, n_ia = n_array[s, ,1:nsim, ])
    
    tab0 <- rbind(tab0, c(m, rate_diff[s], unlist(design[3:6])))
  }
}

for(s in 1:S){
  for(m in 1:4){
    design <- get.design(lam = opt_par$lam[m], a = opt_par$a[m], nsim = nsim, 
                         t = t_analysis, pp = pp1_list[[m]][s, , , ], 
                         arm_size = arm_size, n_ia = n_array[s, ,(1+nsim):(2*nsim), ])
    
    tab1 <- rbind(tab1, c(m, rate_diff[s], unlist(design[3:6])))
  }
}

tab1 <- as.data.frame(tab1)
tab0 <- as.data.frame(tab0)
names(tab0) <- c("met", "p0", "fwp", "powA", "powB")
names(tab1) <- c("met", "p0", "fwp", "powA", "powB")

##----------------------------------------------------
##  power figs:
##----------------------------------------------------
setEPS()
postscript(paste("power curve_v=", v_H, ".eps", sep = ""), width = 10, height = 10)
par(mfrow = c(1, 1), mar = c(5, 5, 1, 1), oma = c(1, 1, 1, 1))
plot(0, 0, type = "n", xaxt = "n", yaxt = "n",  
     xlim = c(-0.1, 0.15), ylim = c(0, 1),
     xlab = "", ylab = "", main = "")
axis(1, seq(-0.1, 0.15, 0.05), seq(-0.1, 0.15, 0.05), cex.axis = 1.25, las = 1)
axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1), cex.axis = 1.25, las = 1)
abline ( h =0.1, lty = 2, col = "grey")
for(m in 1:4){
  lines(rate_diff[1:6], tab0$fwp[tab0$met == m][1:6], type = "b", lwd = 3, lty = 1, col = m)
  lines(rate_diff[1:6], tab1$fwp[tab1$met == m][1:6], type = "b", lwd = 3, lty = 1, col = m)
}
title(xlab="Difference from IC", 
      ylab="FWP/FWER", cex.lab = 1.5, line = 3.2)

legend(0.1, 0.5, bty = "n", pch = 21, col = c(1:4), 
       lty = c(1, 1, 1, 1), cex = 1.2, lwd = 3, 
       c("IC only", "HY-MAP", "HY-rMAP", "HY-MEM"))
dev.off()
##----------------------------------------------------
##  EN figs:
##----------------------------------------------------
setEPS()
postscript(paste("en curve_v=", v_H, ".eps", sep = ""), width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(5, 5, 1, 1), oma = c(1, 1, 1, 1))
plot(NA, NA, type = "n", xaxt = "n", yaxt = "n",  
     xlim = c(-0.1, 0.15), ylim = c(t_analysis[1], 180+5),
     xlab = "", ylab = "", main = "")
axis(1, seq(-0.1, 0.2, 0.05), seq(-0.1, 0.2, 0.05), cex.axis = 1.25, las = 1)
axis(2, seq(t_analysis[1], 180, 10), seq(t_analysis[1], 180, 10), 
     cex.axis = 1.25, las = 1)
abline( h = 150, lwd = 2, col = "grey", lty = 2)
for(m in 1:4){
  lines(rate_diff[1:6], tab0[tab0$met == m, 6][1:6], type = "b", lwd = 3, col = m)
}
title(xlab="Difference from IC", 
      ylab="Expected sample size under H0", cex.lab = 1.5, line = 3.2)
legend(-0.1, 190, bty = "n", pch = 21, col = 1:4, cex = 1.2, lwd = 3,
       c("IC only", "HY-MAP", "HY-rMAP", "HY-MEM"))

plot(NA, NA, type = "n", xaxt = "n", yaxt = "n",  
     xlim = c(-0.1, 0.15), ylim = c(90, 180),
     xlab = "", ylab = "", main = "")
axis(1, seq(-0.1, 0.2, 0.05), seq(-0.1, 0.2, 0.05), cex.axis = 1.25, las = 1)
axis(2, seq(t_analysis[1], 180, 10), seq(t_analysis[1], 180, 10), 
     cex.axis = 1.25, las = 1)
abline( h = 150, lwd = 2, col = "grey", lty = 2)
for(m in 1:4){
  lines(rate_diff[1:6], tab1[tab1$met == m, 6][1:6], type = "b", lwd = 3, col = m)
}
title(xlab="Difference from IC", 
      ylab="Expected sample size under H1", cex.lab = 1.5, line = 3.2)
dev.off()

ftab0 <- ftab1 <- NULL
for(m in 1:4){
  ftab0 <- cbind(ftab0, tab0$fwp[tab0$met == m])
  ftab1 <- cbind(ftab1, tab1$fwp[tab1$met == m])
}
ftab <- cbind(ftab0, ftab1)
colnames(ftab) <- c(paste(rep(c("H0_", "H1_"), each = 4), 
                          rep(c("IC", "MAP", "rMAP", "MEM"), 2), 
                          sep=""))
rownames(ftab) <- rate_diff

library(xtable)
xtable(ftab, digits = 3)
xtable(opt_par[, c(6, 1:3)], digits = 3)
ftab

names(tab0)[6] <- names(tab1)[6] <- "en"
tab0[tab0$met==1, c(1,2, 4, 5)]

tab <- data.frame(met = tab0$met, p=tab0$p0, 
                  fw = paste(tab1$fwp,"(" ,tab0$fwp, ")", sep = ""), en = tab0$en)
write.csv(tab, paste("tab_v=", v_H, ".csv", sep =""))


tab3 <- read.csv("tab_v=0.csv")
tab4 <- read.csv("tab_v=0.025.csv")

tab5 <- cbind(tab3[tab3$p%in%c(-0.10, 0, 0.10), 2:4], tab4[tab4$p%in%c(-0.10, 0, 0.10), 4])
xtable(tab5)

