
# max. comp loglik
est.cal_pwc.f <- function(dat, interac.ind=c(0,0), brks){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  obj.f <- function(para){
    
    #eta = para[1:2]
    alp = exp(para[1:len])
    beta = c(para[len+1:len.beta], 0)
    gam = c(para[len+len.beta+1:len.gam],0)
    
    #res0 <- logLn_cal_pwc( dat$rz, dat$z.o, dat$d, dat$a, dat$x1, dat$x2, gam, alp, beta, brks)
    res0 <- logLn_cal1_pwc(dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gam, alp, beta, brks)
    
    res <- -sum(res0)
    
    # print(para)
    # print(res)
    
    return(res)
    
  }
  
  est <- optim(par=rep(0.01, len+len.beta+len.gam), obj.f, method=("L-BFGS-B"), hessian=TRUE)
  #est <- nlm(f=obj.f, p=rep(0.01, p1+p2+p3+4), hessian=TRUE)
  return(est)
}


# max. obs loglik
est.obs_pwc.f <- function(dat, interac.ind=c(0,0), brks){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  obj.f <- function(para){
    
    # alp = exp(para[1:len]) #c(exp(para[1]), 1)
    # beta = para[len+(1:2)]
    # eta = para[len+2+(1:3)]
    # gam = para[len+6]
    
    #eta0=exp(para[1:3])
    #eta = eta0/(1+eta0) #para[1:2]
    eta = para[1:2]
    p1 = expit.f(para[3])
    alp = exp(para[3+1:len])
    beta = c(para[3+len+1:len.beta], 0)
    gam = c(para[3+len+len.beta+1:len.gam],0)
    
    #res0 <- logLn_obs_pwc(dat$r, dat$rz, dat$z.o, dat$d, dat$a, dat$x1, dat$x2, gam, alp, beta, eta, brks)
    res0 <- logLn_obs_pwc(dat$r, dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gam, alp, beta, eta, p1, brks)
    
    res <- -sum(res0)
    
    # print(para)
    # print(res)
    
    return(res)
    
  }
  
  est <- optim(par=rep(0.01,len+3+len.beta+len.gam), obj.f, method=("L-BFGS-B"), hessian=TRUE)
  #est <- nlm(f=obj.f, p=rep(0.01, p1+p2+p3+4), hessian=TRUE)
  return(est)
}


logL_pwc.f <- function(para, dat, interac.ind=c(0,0), brks){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  
  #res0 <- logLn_obs_pwc(dat$r, dat$rz, dat$z.o, dat$d, dat$a, dat$x1, dat$x2, gam, alp, beta, eta, brks)
  res0 <- logLn_obs_pwc(dat$r, dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gam, alp, beta, eta, p1, brks)
  
  return(res0)
}

score_pwc.f <- function(grad=1e-06, para, dat, interac.ind=c(0,0), brks){
  
  logL0 <- logL_pwc.f(para, dat, interac.ind, brks)
  sc.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    logL.p <- logL_pwc.f(para.p, dat, interac.ind, brks)
    (logL0-logL.p)/grad 
  })
  
  res <- do.call("cbind", sc.list)
  return(res)
}

hess_pwc.f <- function(grad=1e-06, para, dat, interac.ind, brks){
  
  s0 <- score_pwc.f(grad, para, dat, interac.ind, brks)
  h.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    s.p <- score_pwc.f(grad, para.p, dat, interac.ind, brks)
    (colSums(s.p)-colSums(s0))/grad 
  })
  
  res <- do.call("cbind", h.list)
  return(res)
  
}

ase_pwc.f <- function(grad=1e-06, est, dat, interac.ind, brks){
  # score vector
  score.mat = score_pwc.f(grad, est$par, dat, interac.ind, brks)
  #score.mat = score_trueh.f(grad, est$par, del, time, v.mat, pos.cov.T, pos.cov.Z, pos.cov.X)
  
  #hess.mat = hess_pwc.f(grad, est, dat, interac.ind, brks)
  hess.mat = est$hess 
  #hessian(lik.f, est, x1=dt.cov[,1], v=dt.cov[,-1, drop=FALSE], del1, del2, b1, b2, cov_dist) # hessian matrix by numeircal approximation
  B = t(score.mat) %*% score.mat
  A = (hess.mat)
  #A = (hesmat)
  avar = ginv(A)%*% B%*% ginv(A)
  ase = sqrt(diag(avar))
  res = list(ase=ase, avar=avar, A=A, B=B, score=colSums(score.mat))
  return(res)
}


cumrisk_g_x1_x2_pwc.f <- function(t, x1, x2, para, brks, interac.ind){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  alp = exp(para[1:len])
  beta = c(para[len+1:len.beta], 0)
  gam = c(para[len+len.beta+1:len.gam],0)
  
  
  #res0 <- logLn_cal1_pwc(dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gam, alp, beta, brks)
  
  res0 <- logL_cal1_pwc(2, 1, 0, t, x1, x2, gam, alp, beta, brks)
  
  return(exp(res0))
}

ase.cumrisk_g_x1_x2_pwc.f <- function(grad, para.varmat, t, x1, x2, para, brks, interac.ind){
  
  res0 <- cumrisk_pwc.f(t, x1, x2, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- cumriskg_x1_x2_pwc.f(t, x1, x2, para.p, brks, interac.ind)
    (res0-res0.p)/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=res0, partial=dres, lb=res0-1.96*ase, ub=res0+1.96*ase))
}


cumrisk_g_x2_pwc.f <- function(t, x2, para, brks, interac.ind){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  res1 <- exp(logL_cal_pwc(2, 1, 0, t, 1, x2, gam, alp, beta, eta, p1, brks))#*p1
  res0 <- exp(logL_cal_pwc(2, 1, 0, t, 0, x2, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  p11 = expit.f(sum(eta))*p1
  p10 = expit.f(eta[1])*(1-p1)
  
  denom=x2*(p11+p10)+(1-x2)*( 1-p11-p10 )
  
  return((res0+res1)/denom)
}

ase.cumrisk_g_x2_pwc.f <- function(grad, para.varmat, t, x2, para, brks, interac.ind){
  
  res0 <- cumrisk_g_x2_pwc.f(t, x2, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- cumrisk_g_x2_pwc.f(t, x2, para.p, brks, interac.ind)
    (res0-res0.p)/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=res0, partial=dres, ase=ase, lb=res0-1.96*ase, ub=res0+1.96*ase))
}

cumrisk_g_x1_pwc.f <- function(t, x1, para, brks, interac.ind){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  res1 <- exp(logL_cal_pwc(2, 1, 0, t, x1, 1, gam, alp, beta, eta, p1, brks))#*p1
  res0 <- exp(logL_cal_pwc(2, 1, 0, t, x1, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  denom=x1*p1+(1-x1)*( 1-p1)
  
  return((res0+res1)/denom)
}

ase.cumrisk_g_x1_pwc.f <- function(grad, para.varmat, t, x1, para, brks, interac.ind){
  
  res0 <- cumrisk_g_x1_pwc.f(t, x1, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- cumrisk_g_x1_pwc.f(t, x1, para.p, brks, interac.ind)
    (res0-res0.p)/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=res0, ase=ase, lb=res0-1.96*ase, ub=res0+1.96*ase))
}

pre1_risk_g_x2_pwc.f <- function(x2, para, brks, interac.ind){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  #alp = exp(para[2+1:len])
  #beta = c(para[2+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  p11=expit.f(sum(eta))*p1
  p10=expit.f(eta[1])*(1-p1)
  
  p_x2_1 = (x2*p11+(1-x2)*(p1-p11))
  p_x2_0 = (x2*p10+(1-x2)*(1-p1-p10))
  
  res1 <- p_x2_1*expit.f(gam[1]+gam[2]+gam[3]*x2)#exp(logL_cal_pwc(0, NA, 1, t, 1, x2, gam, alp, beta, eta, brks))*p1
  res0 <- p_x2_0*expit.f(gam[1]+gam[3]*x2)#exp(logL_cal_pwc(0, NA, 1, t, 0, x2, gam, alp, beta, eta, brks))*(1-p1)
  
  denom=x2*(p11+p10)+(1-x2)*( 1-p11-p10 )
  
  return((res0+res1)/denom)
}

pre1_risk_g_x1_pwc.f <- function(x1, para, brks, interac.ind){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  #alp = exp(para[2+1:len])
  #beta = c(para[2+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  
  
  p1_x1=expit.f(eta[1]+eta[2]*x1)#*p1
  #p_x1_0=expit.f(eta[1])#*(1-p1)
  
  #p_x1_1 = (x1*p11+(1-x1)*p10)
  #p_x1_0 = (x1*(p1-p11)+(1-x1)*(1-p1-p10))
  
  res1 <- p1_x1*expit.f(gam[1]+gam[2]*x1+gam[3]+gam[4]*x1)#exp(logL_cal_pwc(0, NA, 1, t, 1, x2, gam, alp, beta, eta, brks))*p1
  res0 <- (1-p1_x1)*expit.f(gam[1]+gam[2]*x1)#exp(logL_cal_pwc(0, NA, 1, t, 0, x2, gam, alp, beta, eta, brks))*(1-p1)
  
  #denom=x1*p1+(1-x1)*( 1-p1)
  
  return((res0+res1))
}

ase.pre1_risk_g_x1_pwc.f <- function(grad, para.varmat, x1, para, brks, interac.ind){
  
  res0 <- pre1_risk_g_x1_pwc.f(x1, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- pre1_risk_g_x1_pwc.f(x1, para.p, brks, interac.ind)
    (res0-res0.p)/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=res0, lb=res0-1.96*ase, ub=res0+1.96*ase))
}

ase.pre1_risk_g_x2_pwc.f <- function(grad, para.varmat, x2, para, brks, interac.ind){
  
  res0 <- pre1_risk_g_x2_pwc.f(x2, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- pre1_risk_g_x2_pwc.f(x2, para.p, brks, interac.ind)
    (res0-res0.p)/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=res0, lb=res0-1.96*ase, ub=res0+1.96*ase))
}

inc1_risk_g_x1_pwc.f <- function(t, x1, para, brks, interac.ind){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  #p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  #gam = c(para[2+len+len.beta+1:len.gam],0)
  
  p1_x1=expit.f(eta[1]+eta[2]*x1)#*p1
  #p10=expit.f(eta[1])*(1-p1)
  
  #p_x1_1 = (x1*p11+(1-x1)*p10)
  #p_x1_0 = (x1*(p1-p11)+(1-x1)*(1-p1-p10))
  
  levs_x1_1 = alp*exp(beta[1]*x1+beta[2]+beta[3]*x1); 
  levs_x1_0 = alp*exp(beta[1]*x1); 
  surv_x1_1 = ppwc(t, brks, levs_x1_1, lower=FALSE,0)
  surv_x1_0 = ppwc(t, brks, levs_x1_0, lower=FALSE,0)
  
  res1 <- p1_x1 * (1-surv_x1_1) #exp(logL_cal_pwc(0, NA, 1, t, x1, 1, gam, alp, beta, eta, brks))*p1
  res0 <- (1-p1_x1) * (1-surv_x1_0) #exp(logL_cal_pwc(0, NA, 1, t, x1, 0, gam, alp, beta, eta, brks))*(1-p1)
  
  #denom=x1*p1+(1-x1)*( 1-p1)
  
  #return((res0+res1)/denom)
  return(res0+res1)
}

inc1_risk_g_x2_pwc.f <- function(t, x2, para, brks, interac.ind){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  #gam = c(para[2+len+len.beta+1:len.gam],0)
  
  p11=expit.f(sum(eta))*p1
  p10=expit.f(eta[1])*(1-p1)
  
  p_x2_1 = (x2*p11+(1-x2)*(p1-p11))
  p_x2_0 = (x2*p10+(1-x2)*(1-p1-p10))
  
  levs_x2_1 = alp*exp(beta[1]+beta[2]*x2+beta[3]*x2); 
  levs_x2_0 = alp*exp(beta[2]*x2); 
  surv_x2_1 = ppwc(t, brks, levs_x2_1, lower=FALSE,0)
  surv_x2_0 = ppwc(t, brks, levs_x2_0, lower=FALSE,0)
  
  res1 <- p_x2_1 * (1-surv_x2_1) #exp(logL_cal_pwc(0, NA, 1, t, x1, 1, gam, alp, beta, eta, brks))*p1
  res0 <- p_x2_0 * (1-surv_x2_0) #exp(logL_cal_pwc(0, NA, 1, t, x1, 0, gam, alp, beta, eta, brks))*(1-p1)
  
  denom=x2*(p11+p10)+(1-x2)*( 1-p11-p10 )
  
  return((res0+res1)/denom)
}

wgt1_risk_g_x1_pwc.f <- function(t, x1, para, brks, interac.ind){
  
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  p11=expit.f(sum(eta))*p1
  p10=expit.f(eta[1])*(1-p1)
  
  p_x1_1 = (x1*p11+(1-x1)*p10)
  p_x1_0 = (x1*(p1-p11)+(1-x1)*(1-p1-p10))
  
  
  levs_x1_1 = alp*exp(beta[1]*x1+beta[2]+beta[3]*x1); 
  levs_x1_0 = alp*exp(beta[1]*x1); 
  surv_x1_1 = ppwc(t, brks, levs_x1_1, lower=FALSE,0)
  surv_x1_0 = ppwc(t, brks, levs_x1_0, lower=FALSE,0)
  
  pz_x1_1 = expit.f(gam[1]+gam[2]*x1+gam[3]+gam[4]*x1)
  pz_x1_0 = expit.f(gam[1]+gam[2]*x1)
  
  res1 <- p_x1_1*pz_x1_1/(1-(1-pz_x1_1)*surv_x1_1)  #exp(logL_cal_pwc(0, NA, 1, t, 1, x2, gam, alp, beta, eta, brks))*p1
  res0 <- p_x1_0*pz_x1_0/(1-(1-pz_x1_0)*surv_x1_0)#exp(logL_cal_pwc(0, NA, 1, t, 0, x2, gam, alp, beta, eta, brks))*(1-p1)
  
  denom=x1*p1+(1-x1)*( 1-p1)
  
  return((res0+res1)/denom)
  
  #return(exp(res0))
  
}

wgt2_risk_g_x1_pwc.f <- function(t, x1, para, brks, interac.ind){
  
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  p11=expit.f(sum(eta))*p1
  p10=expit.f(eta[1])*(1-p1)
  
  p_x1_1 = (x1*p11+(1-x1)*p10)
  p_x1_0 = (x1*(p1-p11)+(1-x1)*(1-p1-p10))
  
  
  levs_x1_1 = alp*exp(beta[1]*x1+beta[2]+beta[3]*x1); 
  levs_x1_0 = alp*exp(beta[1]*x1); 
  surv_x1_1 = ppwc(t, brks, levs_x1_1, lower=FALSE,0)
  surv_x1_0 = ppwc(t, brks, levs_x1_0, lower=FALSE,0)
  
  pz_x1_1 = expit.f(gam[1]+gam[2]*x1+gam[3]+gam[4]*x1)
  pz_x1_0 = expit.f(gam[1]+gam[2]*x1)
  
  res1 <- p_x1_1*(1-pz_x1_1)*(1-surv_x1_1)/(1-(1-pz_x1_1)*surv_x1_1)  #exp(logL_cal_pwc(0, NA, 1, t, 1, x2, gam, alp, beta, eta, brks))*p1
  res0 <- p_x1_0*(1-pz_x1_0)*(1-surv_x1_0)/(1-(1-pz_x1_0)*surv_x1_0)#exp(logL_cal_pwc(0, NA, 1, t, 0, x2, gam, alp, beta, eta, brks))*(1-p1)
  
  denom=x1*p1+(1-x1)*( 1-p1)
  
  return((res0+res1)/denom)
  
  #return(exp(res0))
  
}


# TPR=P(X1=1|D0=1); FPR=P(X1=1|D0=0)
tpr.fpr.auc_of_x1.f <- function(para, brks, interac.ind){
  
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  
  # P(x1, x2)
  p11 = expit.f( sum(eta) )*p1
  p10 = p1 - p11
  p01 = expit.f( eta[1] )*(1-p1)
  p00 = (1-p1)-p01
  
  # numerator
  res11 = expit.f( cbind(1, 1, 1, 1) %*% gam )*p11 # P(D0=1|x=(1,1))*P(x=(1,1))
  res10 = expit.f( cbind(1, 1, 0, 0) %*% gam )*p10 # P(D0=1|x=(1,0))*P(x=(1,0))
  res01 = expit.f( cbind(1, 0, 1, 0) %*% gam )*p01 # P(D0=1|x=(0,1))*P(x=(0,1))
  res00 = expit.f( cbind(1, 0, 0, 0) %*% gam )*p00 # P(D0=1|x=(0,0))*P(x=(0,0))
  denom1 = res11 + res10 + res01 + res00 # P(D0=1)
  
  tpr= (res11+res10)/denom1
  fpr= (p11+p10-(res11+res10))/(1-denom1)
  
  
  auc=tpr*(1-fpr)
  auc.ties=0.5*(tpr+(1-fpr))
  
  return(list(tpr=tpr, fpr=fpr, #auc=auc, 
              auc.ties=auc.ties))
}

# TPR=P(X2=1|D0=1); FPR=P(X2=1|D0=0)
tpr.fpr.auc_of_x2.f <- function(para, brks, interac.ind){

  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  
  # P(x1, x2)
  p11 = expit.f( sum(eta) )*p1
  p10 = p1 - p11
  p01 = expit.f( eta[1] )*(1-p1)
  p00 = (1-p1)-p01
  
  # numerator
  res11 = expit.f( cbind(1, 1, 1, 1) %*% gam )*p11 # P(D0=1|x=(1,1))*P(x=(1,1))
  res10 = expit.f( cbind(1, 1, 0, 0) %*% gam )*p10 # P(D0=1|x=(1,0))*P(x=(1,0))
  res01 = expit.f( cbind(1, 0, 1, 0) %*% gam )*p01 # P(D0=1|x=(0,1))*P(x=(0,1))
  res00 = expit.f( cbind(1, 0, 0, 0) %*% gam )*p00 # P(D0=1|x=(0,0))*P(x=(0,0))
  denom1 = res11 + res10 + res01 + res00 # P(D0=1)
  
  tpr= (res11+res01)/denom1
  fpr= (p11+p01-(res11+res01))/(1-denom1)
  
  auc=tpr*(1-fpr)
  auc.ties=0.5*(tpr+(1-fpr))
  
  return(list(tpr=tpr, fpr=fpr, #auc=auc, 
              auc.ties=auc.ties))
}

# TPR=P(X2=1|D(t)=1); FPR=P(X2=1|D(t)=0)
dyn.tpr.fpr.auc_of_x2.f <- function(t, para, brks, interac.ind){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  #res1 <- exp(logL_cal_pwc(2, 1, 0, t, 1, x2, gam, alp, beta, eta, p1, brks))#*p1
  #res0 <- exp(logL_cal_pwc(2, 1, 0, t, 0, x2, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  res11 <- exp(logL_cal_pwc(2, 1, 0, t, 1, 1, gam, alp, beta, eta, p1, brks))#*p1
  res01 <- exp(logL_cal_pwc(2, 1, 0, t, 0, 1, gam, alp, beta, eta, p1, brks))#*(1-p1)
  res10 <- exp(logL_cal_pwc(2, 1, 0, t, 1, 0, gam, alp, beta, eta, p1, brks))#*p1
  res00 <- exp(logL_cal_pwc(2, 1, 0, t, 0, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  denom1 = res11 + res01 + res10 + res00
  
  p11 = expit.f(sum(eta))*p1
  p10 = expit.f(eta[1])*(1-p1)
  
  tpr= (res11+res01)/denom1
  fpr= (p11+p10-(res11+res01))/(1-denom1)
  
  auc.ties=0.5*(tpr+(1-fpr))
  
  return(list(tpr=tpr, fpr=fpr, #auc=auc, 
              auc.ties=auc.ties))
}

# TPR=P(X2=1|D(t)=1); FPR=P(X2=1|D(t)=0)
dyn.tpr.fpr.auc_of_x1.f <- function(t, para, brks, interac.ind){
  
  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  #res1 <- exp(logL_cal_pwc(2, 1, 0, t, 1, x2, gam, alp, beta, eta, p1, brks))#*p1
  #res0 <- exp(logL_cal_pwc(2, 1, 0, t, 0, x2, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  res11 <- exp(logL_cal_pwc(2, 1, 0, t, 1, 1, gam, alp, beta, eta, p1, brks))#*p1
  res01 <- exp(logL_cal_pwc(2, 1, 0, t, 1, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  res10 <- exp(logL_cal_pwc(2, 1, 0, t, 0, 1, gam, alp, beta, eta, p1, brks))#*p1
  res00 <- exp(logL_cal_pwc(2, 1, 0, t, 0, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  denom1 = res11 + res01 + res10 + res00
  
  tpr= (res11+res01)/denom1
  fpr= (p1-(res11+res01))/(1-denom1)
  
  auc.ties=0.5*(tpr+(1-fpr))
  
  return(list(tpr=tpr, fpr=fpr, #auc=auc, 
              auc.ties=auc.ties))
}

ase.tpr.fpr.auc_of_x2.f <- function(grad, para.varmat, para, brks, interac.ind){
  
  res0 <- tpr.fpr.auc_of_x2.f(para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- tpr.fpr.auc_of_x2.f(para.p, brks, interac.ind)
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}

ase.tpr.fpr.auc_of_x1.f <- function(grad, para.varmat, para, brks, interac.ind){
  
  res0 <- tpr.fpr.auc_of_x1.f(para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- tpr.fpr.auc_of_x1.f(para.p, brks, interac.ind)
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}

ase.dyn.tpr.fpr.auc_of_x2.f <- function(t, grad, para.varmat, para, brks, interac.ind){
  
  res0 <- dyn.tpr.fpr.auc_of_x2.f(t, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- dyn.tpr.fpr.auc_of_x2.f(t, para.p, brks, interac.ind)
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}

ase.dyn.tpr.fpr.auc_of_x1.f <- function(t, grad, para.varmat, para, brks, interac.ind){
  
  res0 <- dyn.tpr.fpr.auc_of_x1.f(t, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- dyn.tpr.fpr.auc_of_x1.f(t, para.p, brks, interac.ind)
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  
  
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}


ase.dyn.auc_of_x2.f <- function(t, grad, para.varmat, para, brks, interac.ind){
  
  res0 <- dyn.tpr.fpr.auc_of_x2.f(t, para, brks, interac.ind)$auc.ties
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- dyn.tpr.fpr.auc_of_x2.f(t, para.p, brks, interac.ind)$auc.ties
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}

ase.dyn.auc_of_x1.f <- function(t, grad, para.varmat, para, brks, interac.ind){
  
  res0 <- dyn.tpr.fpr.auc_of_x1.f(t, para, brks, interac.ind)$auc.ties
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- dyn.tpr.fpr.auc_of_x1.f(t, para.p, brks, interac.ind)$auc.ties
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}

# wgt1_risk_g_x1_pwc.f(1, 1, para=c(test2_pwc_cov$par,p1), brks=c(1,2), interac.ind=c(0,0))
# 
# sapply(seq(0,5, by=0.1), function(t){
#   wgt2_risk_g_x1_pwc.f(t, 1, para=c(test2_pwc_cov$par,p1), brks=c(1,2), interac.ind=c(0,0))
# })
# 
# sapply(seq(0,5, by=0.1), function(t){
#   wgt2_risk_g_x1_pwc.f(t, 0, para=c(test2_pwc_cov$par,p1), brks=c(1,2), interac.ind=c(0,0))
# })

# cumrisk_pwc.f(1, 1, 1, para=test2_pwc_cov$par[-c(1,2)], brks=1, interac.ind=c(0,0))
# cumrisk_pwc.f(1, 1, 0, para=test2_pwc_cov$par[-c(1,2)], brks=1, interac.ind=c(0,0))
# cumrisk_pwc.f(1, 0, 1, para=test2_pwc_cov$par[-c(1,2)], brks=1, interac.ind=c(0,0))
# cumrisk_pwc.f(1, 0, 0, para=test2_pwc_cov$par[-c(1,2)], brks=1, interac.ind=c(0,0))

# res11=ase.cumrisk_pwc.f(grad=1e-06, para.varmat=pwc_ase$avar[-(1:2), -(1:2)], 1,1,1, para=test2_pwc_cov$par[-(1:2)], brks=1, interac.ind=c(0,0))
# res10=ase.cumrisk_pwc.f(grad=1e-06, para.varmat=pwc_ase$avar[-(1:2), -(1:2)], 1,1,0, para=test2_pwc_cov$par[-(1:2)], brks=1, interac.ind=c(0,0))
# res01=ase.cumrisk_pwc.f(grad=1e-06, para.varmat=pwc_ase$avar[-(1:2), -(1:2)], 1,0,1, para=test2_pwc_cov$par[-(1:2)], brks=1, interac.ind=c(0,0))
# res00=ase.cumrisk_pwc.f(grad=1e-06, para.varmat=pwc_ase$avar[-(1:2), -(1:2)], 1,0,0, para=test2_pwc_cov$par[-(1:2)], brks=1, interac.ind=c(0,0))





