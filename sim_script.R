
# ------------------------------------------------------------------------------------------
# This is a script for one simulation study. 
# A for-loop is used for repeated simulations. 
# In each simulation run, 3 steps are included: 
#  1) data generation
#  2) two-phase design implementation
#  3) estimation and inference
# ------------------------------------------------------------------------------------------

library(MASS)
source("/utility/design.r") # functions for phase II designs
source("/utility/design_resid_interval.r") # functions for phase II designs
Rcpp::sourceCpp("/utility/loglik_interval.cpp") # functions for estimation via max likelihood
source("/utility/est_pwc_interval.R") # functions for estimation via max likelihood

####### Test Data configuration #######

# phase I sample size 
n=8000 
# log of baseline hazard; and true regression coefs. for T|X,Z=1
log.inv.alp=5; be_vec=c(log(1), log(3.5), log(1)); 
# true parameters for X2|X1
et_vec = c(log(0.44), log(3))
# true regression coefs. for Z|X
ga_vec = c(log(0.01), log(2.8), log(7), log(1))
# true parameters for R1|X1
nu_vec=c(log(0.1/0.9), 0)

# phase I stratifying factor: prevalent and incident event status, and phase I covariate x1
phI.design.factor= c("d1", "d2", "x1")
# size of the phase II sample
n2_samp=1000; 

w=0.5; phII.design= paste0("BRSD-eff-seq", w) # bivariate residual-dependent sampling, 0.5*n2_samp were selected in each wave
#phII.design="bal"; # balanced design used in phase II

###### A sample simulation run ######

# a simulated complete data sample 
dat=read.csv(file="sample_dat.csv", header=1)

# Phase I stratification
phI_strat = phaesI_strat.f(pcuts=c(1/3,2/3), dat, design.factor=phI.design.factor)

# Phase II selection 
para0 <- est.H0_wei.f(dat)
if(phII.design %in% c("RSD-Smu1", "RSD-Smu2", paste0("BRSD-eff-seq",w))){
  
  # residual-dependent sampling 
  phII_sel = designPhII_resid(phaseI_strat=phI_strat, n2samp=n2_samp, 
                              design=phII.design, w, 
                              para0$par)$sel
  
}else{
  # Phase II design
  phII_sel = design_strat.f(n2_samp, phII.design, phI_strat)
}    

dt.phII = dat
dt.phII$x2[-phII_sel$s_id] <- NA
dt.phII$r <- 1
dt.phII$r[-phII_sel$s_id] <- 0

# estimation
est <- est.obs_pwc.f(dt.phII, interac.ind, brks)

# inference 
ase <- ase_pwc.f(grad=1e-06, est, dt.phII, interac.ind, brks)

rds.list[[iter]] <- list(est=est, ase=ase, dt.phII=dt.phII, phII_sel=phII_sel)

  
# save the results
saveRDS(rds.list, file=paste(n2_samp, "pwc-", paste0(brks, collapse = "-"), "intrc", paste0(interac.ind, collapse = ""), phII.design, phI.strat.style, 
                              ".RDS", sep=""))
  



