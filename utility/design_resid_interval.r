
logLn_cal1_mu_wei.f <- function(mu1, mu2, d1, d2, t1, t2, alp){
  
  expcovZ= exp(mu1)
  expcovT=exp(mu2)
  Hft1 = (alp[1]*t1)^(alp[2])*expcovT
  Sft1 = exp(-Hft1)
  Hft2 = (alp[1]*t2)^(alp[2])*expcovT
  Sft2 = exp(-Hft2)
  piZ = expcovZ/(1+expcovZ)
  
  res <- log(piZ)
  res[d1>1 & d2==1] <- (log((1-Sft2)*(1-piZ)+piZ))[d1>1 & d2==1]
  res[d1>1 & d2==0] <- rep(0, sum(d1>1 & d2==0))
  res[d2>1] <- rep(0, sum(d2>1))
  res[d1==0 & d2==1] <- (log( Sft1-Sft2 )+log(1-piZ))[d1==0 & d2==1]
  res[d1==0 & d2==0] <- (log(1-piZ)-Hft2)[d1==0 & d2==0]
  
  return(res)
}

Smu_mu_wei.f <- function(mu1, mu2, d1, d2, t1, t2, alp){
  
  expcovZ= exp(mu1)
  expcovT=exp(mu2)
  Hft1 = (alp[1]*t1)^(alp[2])*expcovT
  Sft1 = exp(-Hft1)
  Hft2 = (alp[1]*t2)^(alp[2])*expcovT
  Sft2 = exp(-Hft2)
  piZ = expcovZ/(1+expcovZ)
  
  # res <- log(piZ)
  # res[d1>1 & d2==1] <- (log((1-Sft2)*(1-piZ)+piZ))[d1>1 & d2==1]
  # res[d1>1 & d2==0] <- rep(0, sum(d1>1 & d2==0))
  # res[d2>1] <- rep(0, sum(d2>1))
  # res[d1==0 & d2==1] <- (log( Sft1-Sft2 )+log(1-piZ))[d1==0 & d2==1]
  # res[d1==0 & d2==0] <- (log(1-piZ)-Hft2)[d1==0 & d2==0]
  
  res1 <- (1-piZ)
  res1[d1>1 & d2==1] <- (Sft2*piZ*(1-piZ)/((1-Sft2)*(1-piZ)+piZ))[d1>1 & d2==1]#( (1-Sft2)*(-piZ)*(1-piZ) + piZ*(1-piZ) )/((1-Sft2)*(1-piZ)+piZ)
  res1[d1>1 & d2==0] <- rep(0, sum(d1>1 & d2==0))
  res1[d2>1] <- rep(0, sum(d2>1))
  res1[d1==0 & d2==1] <- -piZ[d1==0 & d2==1]#(log( Sft1-Sft2 )+log(1-piZ))[d1==0 & d2==1]
  res1[d1==0 & d2==0] <- -piZ[d1==0 & d2==0]
  
  res2 <- rep(0,length(mu1))
  res2[d1>1 & d2==1] <- ((1-piZ)*(Sft2*Hft2)/((1-Sft2)*(1-piZ)+piZ))[d1>1 & d2==1]#(log((1-Sft2)*(1-piZ)+piZ))[d1>1 & d2==1]
  res2[d1>1 & d2==0] <- rep(0, sum(d1>1 & d2==0))
  res2[d2>1] <- rep(0, sum(d2>1))
  res2[d1==0 & d2==1] <- ((Sft2*Hft2-Sft1*Hft1)/(Sft1-Sft2))[d1==0 & d2==1]#(log( Sft1-Sft2 )+log(1-piZ))[d1==0 & d2==1]
  res2[d1==0 & d2==0] <- (-Hft2)[d1==0 & d2==0]#(log(1-piZ)-Hft2)[d1==0 & d2==0]
  
  res=list(Smu1=res1, Smu2=res2)
  
  return(res)
}

Smu_num_vec.f <- function(grad, d1, d2, a1, a2, x1, para.H0){
  
  alp=c(exp(para.H0[1]), exp(para.H0[2])); beta=para.H0[3]; gam=para.H0[4:5]
  mu1=gam[1]+gam[2]*x1
  mu2=beta[1]*x1
  
  l0 = logLn_cal1_mu_wei.f(mu1, mu2, d1, d2, a1, a2, alp)
  
  lp_mu2 = logLn_cal1_mu_wei.f(mu1, mu2+grad, d1, d2, a1, a2, alp)
  
  lp_mu1 = logLn_cal1_mu_wei.f(mu1+grad, mu2, d1, d2, a1, a2, alp)
  
  res <- list(Smu1=(lp_mu1-l0)/grad, Smu2=(lp_mu2-l0)/grad)
  
  return(res)
}

est.H0_wei.f <- function(dat){
  
  obj.f <- function(para){
    
    alp = exp(para[1:2])
    beta = c(para[3],0,0)
    gamma = c(para[4:5],0,0)
    
    #res0 <- log.calLn.f(dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, alp, beta, gamma)
    res0 <- logLn_cal1_wei(dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gamma, alp, beta)
    
    res <- -sum(res0)
    
    # print(para)
    # print(res)
    
    return(res)
    
  }
  
  est <- optim(par=rep(0.01, 5), obj.f, method=("L-BFGS-B"), hessian=TRUE)
  #est <- nlm(f=obj.f, p=rep(0.01, p1+p2+p3+4), hessian=TRUE)
  return(est)
}

designPhII_resid <- function(phaseI_strat, n2samp, design="RSD-Smu1", w=NULL,
                             para.H0){
  
  # qzi0 <- Smu_vec.f(phaseI_strat$dt_ext$d1, phaseI_strat$dt_ext$d2,
  #                  phaseI_strat$dt_ext$a1, phaseI_strat$dt_ext$a2,
  #                  phaseI_strat$dt_ext$x1, para.H0)
  
  alp=c(exp(para.H0[1]), exp(para.H0[2])); beta=para.H0[3]; gam=para.H0[4:5]
  mu1=gam[1]+gam[2]*phaseI_strat$dt_ext$x1
  mu2=beta[1]*phaseI_strat$dt_ext$x1
  
  qzi <- Smu_mu_wei.f(mu1, mu2,
                      phaseI_strat$dt_ext$d1, phaseI_strat$dt_ext$d2,
                      phaseI_strat$dt_ext$a1, phaseI_strat$dt_ext$a2, alp)

  
  k= round(min(min(n2samp/2, sum(phaseI_strat$dt_ext$x1==0)/2, sum(phaseI_strat$dt_ext$x1==1)/2)*0.15), 10)
  #k=round(n2samp/2*0.05)
  if(design %in% c("RSD-Smu1")){
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[1]], phaseI_strat) 
  }else if(design %in% c("RSD-Smu2")){
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[2]], phaseI_strat)
    
  }else if(design %in% paste0("BRSD-eff-seq",w)){
    
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    
    covS=cov(mat)
    Seff1=mat[,1]#-covS[2,1]/covS[2,2]*mat[,1]
    Seff2=mat[,2]#-covS[1,2]/covS[1,1]*mat[,2]
    
    if(w==0 | w==1){
      qzi_use=w*Seff1+(1-w)*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a=round(n2samp*w)
      #ka=round(n2a/2*0.15)
      ka= round(min(min(n2a/2, sum(phaseI_strat$dt_ext$x1==0)/2, sum(phaseI_strat$dt_ext$x1==1)/2)*0.15, 10))
      res01 <- rkd_ext_opt.f(ka, n2a, Seff1, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% res01$s_id), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a
      #kb=round(n2b/2*0.15)
      kb= round(min(min(n2b/2, sum(dt_b$x1==0)/2, sum(dt_b$x1==1)/2)*0.15, 10))
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, res01$s_id)
      
    }
  }

  res = list(sel=res0, qzi=qzi)
  
  return(res)
}
