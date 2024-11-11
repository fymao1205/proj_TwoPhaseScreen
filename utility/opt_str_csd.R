
# try optimal stratified designs for Current Status Data

logit <- function(x){log(x/(1-x))}
# joint probab of (Y, X) ; Y=(D0, D1) if rz=1 and Y=D1 if rz=0
P_Y_X_g_A_rz_csd_pwc.f <- function(rz, z, d, t, x1, x2, para, brks){
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:2])
  beta = c(para[3+2+1:2], 0)
  gam = c(para[3+2+2+1:3],0)
  
  #res0 <- logLn_csd_cal_wei(rz, z, d, t, x1, x2, gam, alp, beta, eta, p1)
  d1=ifelse(rz==1, z, ifelse(d==0, 0, 2))
  res0 <- logLn_cal_pwc(d1, d, rep(0, length(rz)), t, x1, x2, gam, alp, beta, eta, p1, brks)
  res = exp(res0)
  
  return(res)
}

# joint prob of (A, rz)
# A ~ weibull(sc0, sh0)
# rz|x1 ~ bin(p(x1)); indexed by para.rz


# logL_csd_obs_wei(int r, int rz, int z, int d, double t, double x1, double x2, 
  # NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta)

logL_csd_obs_pwc.f <- function(para, brks, r, rz, z, d, t, x1, x2){
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:2])
  beta = c(para[3+2+1:2], 0)
  gam = c(para[3+2+2+1:3],0)
  
  #res = logL_csd_obs_wei(r, rz, z, d, t, x1, x2, gam, alp, beta, eta, p1)
  d1=ifelse(rz==1, z, ifelse(d==0, 0, 2))
  res = logLn_obs_pwc(r, d1, d, rep(0, length(r)), t, x1, x2, gam, alp, beta, eta, p1, brks)
  
  return(res)
}

score_csd_pwc_ele.f <- function(grad=1e-06, pos, para, brks, r, rz, z, d, t, x1, x2){
  
  para.p = para
  para.p[pos] = para[pos]+grad
  
  logL0 <- logL_csd_obs_pwc.f(para, brks, r, rz, z, d, t, x1, x2)
  logL.p <- logL_csd_obs_pwc.f(para.p, brks, r, rz, z, d, t, x1, x2)
  
  res <- (logL0-logL.p)/grad 
  
  return(res)
}

SS_csd_pwc_ele.f <- function(grad=1e-06, pos, para, brks, r, rz, z, d, t, x1, x2){
  
  #t=ifelse(t>10,10,t)
  
  s0 <- score_csd_pwc_ele.f(grad, pos, para, brks,r, rz, z, d, t, x1, x2)
  
  s0^2
}

S1S2_csd_pwc_ele.f <- function(grad=1e-06, pos1, pos2, para, brks,r, rz, z, d, t, x1, x2){
  
  #t=ifelse(t>10,10,t)
  
  s1 <- score_csd_pwc_ele.f(grad, pos1, para, brks,r, rz, z, d, t, x1, x2)
  s2 <- score_csd_pwc_ele.f(grad, pos2, para, brks,r, rz, z, d, t, x1, x2)
  
  s1*s2
}


# stratification on (rz, z, d, t, x1); can be collapsible - rely on specifying sp.list 
calSS_ele_approx.f <- function(grad, pos, sp, #sp11, sp12, sp13, sp01, sp02, sp03, 
                               para, brks,sc0, sh0, para.rz){
  
  # p1.rz
  p1.rz_g_x1_0 = expit.f(para.rz[1])
  p1.rz_g_x1_1 = expit.f(sum(para.rz))
  
  
  
  # rz = 1: need to 
  int1.f <- function(sp, a){
    
    # joint pdf of P_Y_X_g_A_rz_csd_pwc.f(rz, z, d, t, x1, x2, para, brks) 
    p_z_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 1, 1, para, brks)*p1.rz_g_x1_1
    p_z_0_d_0_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 1, 1, para, brks)*p1.rz_g_x1_1
    p_z_0_d_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 1, 1, para, brks)*p1.rz_g_x1_1
    
    p_z_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 1, 0, para, brks)*p1.rz_g_x1_1
    p_z_0_d_0_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 1, 0, para, brks)*p1.rz_g_x1_1
    p_z_0_d_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 1, 0, para, brks)*p1.rz_g_x1_1
    
    p_z_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 0, 1, para, brks)*p1.rz_g_x1_0
    p_z_0_d_0_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 0, 1, para, brks)*p1.rz_g_x1_0
    p_z_0_d_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 0, 1, para, brks)*p1.rz_g_x1_0
    
    p_z_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 0, 0, para, brks)*p1.rz_g_x1_0
    p_z_0_d_0_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 0, 0, para, brks)*p1.rz_g_x1_0
    p_z_0_d_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 0, 0, para, brks)*p1.rz_g_x1_0
    
    jp_mat = cbind(p_z_1_x_11_rz, p_z_0_d_0_x_11_rz, p_z_0_d_1_x_11_rz,
               p_z_1_x_10_rz, p_z_0_d_0_x_10_rz, p_z_0_d_1_x_10_rz,
               p_z_1_x_01_rz, p_z_0_d_0_x_01_rz, p_z_0_d_1_x_01_rz,
               p_z_1_x_00_rz, p_z_0_d_0_x_00_rz, p_z_0_d_1_x_00_rz)
    
    
    # score^2: SS_csd_pwc_ele.f(grad=1e-06, pos, para, r, rz, z, d, t, x1, x2)
    
    # rx = 1; selected in phase II
    ss_rx_1_z_1_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 1, 1, a, 1, 1)*sp[7]
    ss_rx_1_z_0_d_0_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 0, a, 1, 1)*sp[5]
    ss_rx_1_z_0_d_1_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 1, a, 1, 1)*sp[6]
    
    ss_rx_1_z_1_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 1, 1, a, 1, 0)*sp[7]
    ss_rx_1_z_0_d_0_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 0, a, 1, 0)*sp[5]
    ss_rx_1_z_0_d_1_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 1, a, 1, 0)*sp[6]
    
    ss_rx_1_z_1_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 1, 1, a, 0, 1)*sp[3]
    ss_rx_1_z_0_d_0_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 0, a, 0, 1)*sp[1]
    ss_rx_1_z_0_d_1_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 1, a, 0, 1)*sp[2]
    
    ss_rx_1_z_1_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 1, 1, a, 0, 0)*sp[3]
    ss_rx_1_z_0_d_0_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 0, a, 0, 0)*sp[1]
    ss_rx_1_z_0_d_1_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 1, a, 0, 0)*sp[2]
    
    ss_mat_rx_1 = cbind(ss_rx_1_z_1_x_11, ss_rx_1_z_0_d_0_x_11, ss_rx_1_z_0_d_1_x_11,
               ss_rx_1_z_1_x_10, ss_rx_1_z_0_d_0_x_10, ss_rx_1_z_0_d_1_x_10,
               ss_rx_1_z_1_x_01, ss_rx_1_z_0_d_0_x_01, ss_rx_1_z_0_d_1_x_01,
               ss_rx_1_z_1_x_00, ss_rx_1_z_0_d_0_x_00, ss_rx_1_z_0_d_1_x_00)
    
    # rx = 0; not selected in phase II
    ss_rx_0_z_1_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 1, 1, a, 1, NA)*(1-sp[7])
    ss_rx_0_z_0_d_0_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 0, 0, a, 1, NA)*(1-sp[5])
    ss_rx_0_z_0_d_1_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 0, 1, a, 1, NA)*(1-sp[6])
    
    ss_rx_0_z_1_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 1, 1, a, 0, NA)*(1-sp[3])
    ss_rx_0_z_0_d_0_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 0, 0, a, 0, NA)*(1-sp[1])
    ss_rx_0_z_0_d_1_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 0, 1, a, 0, NA)*(1-sp[2])
    
    ss_mat_rx_0 = cbind( ss_rx_0_z_1_x1_1, ss_rx_0_z_0_d_0_x1_1, ss_rx_0_z_0_d_1_x1_1,
                     ss_rx_0_z_1_x1_0, ss_rx_0_z_0_d_0_x1_0, ss_rx_0_z_0_d_1_x1_0 )
    
    res0 <- sum( (drop( jp_mat %*% t(ss_mat_rx_1) ) + drop( (jp_mat[,c(1:3, 7:9)] + jp_mat[,3+c(1:3, 7:9)]) %*% t(ss_mat_rx_0) ) )*dgamma(a, sh0, sc0) )
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }

  # rz = 0 
  int0.f <- function(sp, a){
    
    # joint pdf of P_Y_X_g_A_rz_csd_pwc.f(rz, z, d, t, x1, x2, para, brks) 
    p_z_0_d_0_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 1, 1, para, brks)*(1-p1.rz_g_x1_1)
    p_d_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 1, 1, para, brks)*(1-p1.rz_g_x1_1)
    
    p_z_0_d_0_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 1, 0, para, brks)*(1-p1.rz_g_x1_1)
    p_d_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 1, 0, para, brks)*(1-p1.rz_g_x1_1)
    
    p_z_0_d_0_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 0, 1, para, brks)*(1-p1.rz_g_x1_0)
    p_d_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 0, 1, para, brks)*(1-p1.rz_g_x1_0)
    
    p_z_0_d_0_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 0, 0, para, brks)*(1-p1.rz_g_x1_0)
    p_d_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 0, 0, para, brks)*(1-p1.rz_g_x1_0)
    
    jp_mat = cbind(p_z_0_d_0_x_11_rz, p_d_1_x_11_rz,
                   p_z_0_d_0_x_10_rz, p_d_1_x_10_rz,
                   p_z_0_d_0_x_01_rz, p_d_1_x_01_rz,
                   p_z_0_d_0_x_00_rz, p_d_1_x_00_rz)
    
    
    # score^2: SS_csd_pwc_ele.f(grad=1e-06, pos, para, r, rz, z, d, t, x1, x2)
    
    # rx = 1; selected in phase II
    ss_rx_1_z_0_d_0_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, 0, 0, a, 1, 1)*sp[5]
    ss_rx_1_d_1_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, NA, 1, a, 1, 1)*sp[8]
    
    ss_rx_1_z_0_d_0_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, 0, 0, a, 1, 0)*sp[5]
    ss_rx_1_d_1_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, NA, 1, a, 1, 0)*sp[8]
    
    ss_rx_1_z_0_d_0_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, 0, 0, a, 0, 1)*sp[1]
    ss_rx_1_d_1_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, NA, 1, a, 0, 1)*sp[4]
    
    ss_rx_1_z_0_d_0_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, 0, 0, a, 0, 0)*sp[1]
    ss_rx_1_d_1_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, NA, 1, a, 0, 0)*sp[4]
    
    ss_mat_rx_1 = cbind(ss_rx_1_z_0_d_0_x_11, ss_rx_1_d_1_x_11,
                        ss_rx_1_z_0_d_0_x_10, ss_rx_1_d_1_x_10,
                        ss_rx_1_z_0_d_0_x_01, ss_rx_1_d_1_x_01,
                        ss_rx_1_z_0_d_0_x_00, ss_rx_1_d_1_x_00)
    
    # rx = 0; not selected in phase II
    ss_rx_0_z_0_d_0_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 0, 0, 0, a, 1, NA)*(1-sp[5])
    ss_rx_0_d_1_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 0, NA, 1, a, 1, NA)*(1-sp[8])
    
    ss_rx_0_z_0_d_0_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 0, 0, 0, a, 0, NA)*(1-sp[1])
    ss_rx_0_d_1_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 0, NA, 1, a, 0, NA)*(1-sp[4])
    
    ss_mat_rx_0 = cbind( ss_rx_0_z_0_d_0_x1_1, ss_rx_0_d_1_x1_1,
                         ss_rx_0_z_0_d_0_x1_0, ss_rx_0_d_1_x1_0 )
    
    res0 <- sum( (drop( jp_mat %*% t(ss_mat_rx_1) ) + drop( (jp_mat[,c(1:2, 5:6)] + jp_mat[,2+c(1:2, 5:6)]) %*% t(ss_mat_rx_0) )  )*dgamma(a, sh0, sc0) )
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  
  # # a.cuts: 0-1
  # res1 <- (pcubature(int1.f, lower = 0, upper=cuts[1],sp=sp.list[[1]])$integral + 
  #            pcubature(int0.f, lower = 0, upper=cuts[1], sp=sp.list[[1]])$integral )
  # # a.cuts: 1-2
  # res2 <- (pcubature(int1.f, lower = cuts[1], upper=cuts[2], sp=sp.list[[2]])$integral + 
  #            pcubature(int0.f, lower = cuts[1], upper=cuts[2], sp=sp.list[[2]])$integral)
  # # a.cuts: 2-3
  # res3 <- (pcubature(int1.f, lower = cuts[2], upper=Inf, sp=sp.list[[3]])$integral + 
  #            pcubature(int0.f, lower = cuts[2], upper=Inf, sp=sp.list[[3]])$integral ) #+ int0_2.f(sp03)
  # 
  # res = res1+res2+res3#+res4
  
  res = (pcubature(int1.f, lower = 0, upper=100,sp=sp)$integral + 
                  pcubature(int0.f, lower = 0, upper=100, sp=sp)$integral )
  
  return(res)
}


calS1S2_ele_approx.f <- function(grad, pos1, pos2, sp, #sp11, sp12, sp13, sp01, sp02, sp03, 
                               para, brks, sc0, sh0, para.rz){
  
  # p1.rz
  p1.rz_g_x1_0 = expit.f(para.rz[1])
  p1.rz_g_x1_1 = expit.f(sum(para.rz))
  
  # rz = 1: need to 
  int1.f <- function(sp, a){
    
    # joint pdf of P_Y_X_g_A_rz_csd_pwc.f(rz, z, d, t, x1, x2, para, brks) 
    p_z_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 1, 1, para, brks) *p1.rz_g_x1_1
    p_z_0_d_0_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 1, 1, para, brks) *p1.rz_g_x1_1
    p_z_0_d_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 1, 1, para, brks) *p1.rz_g_x1_1
    
    p_z_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 1, 0, para, brks) *p1.rz_g_x1_1
    p_z_0_d_0_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 1, 0, para, brks) *p1.rz_g_x1_1
    p_z_0_d_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 1, 0, para, brks) *p1.rz_g_x1_1
    
    p_z_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 0, 1, para, brks) *p1.rz_g_x1_0
    p_z_0_d_0_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 0, 1, para, brks) *p1.rz_g_x1_0
    p_z_0_d_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 0, 1, para, brks) *p1.rz_g_x1_0
    
    p_z_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 0, 0, para, brks) *p1.rz_g_x1_0
    p_z_0_d_0_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 0, 0, para, brks) *p1.rz_g_x1_0
    p_z_0_d_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 0, 0, para, brks) *p1.rz_g_x1_0
    
    jp_mat = cbind(p_z_1_x_11_rz, p_z_0_d_0_x_11_rz, p_z_0_d_1_x_11_rz,
                   p_z_1_x_10_rz, p_z_0_d_0_x_10_rz, p_z_0_d_1_x_10_rz,
                   p_z_1_x_01_rz, p_z_0_d_0_x_01_rz, p_z_0_d_1_x_01_rz,
                   p_z_1_x_00_rz, p_z_0_d_0_x_00_rz, p_z_0_d_1_x_00_rz)
    
    
    # score^2: SS_csd_pwc_ele.f(grad=1e-06, pos, para, brks, r, rz, z, d, t, x1, x2)
    
    # rx = 1; selected in phase II
    ss_rx_1_z_1_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 1, 1, a, 1, 1)*sp[7]
    ss_rx_1_z_0_d_0_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 0, a, 1, 1)*sp[5]
    ss_rx_1_z_0_d_1_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 1, a, 1, 1)*sp[6]
    
    ss_rx_1_z_1_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 1, 1, a, 1, 0)*sp[7]
    ss_rx_1_z_0_d_0_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 0, a, 1, 0)*sp[5]
    ss_rx_1_z_0_d_1_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 1, a, 1, 0)*sp[6]
    
    ss_rx_1_z_1_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 1, 1, a, 0, 1)*sp[3]
    ss_rx_1_z_0_d_0_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 0, a, 0, 1)*sp[1]
    ss_rx_1_z_0_d_1_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 1, a, 0, 1)*sp[2]
    
    ss_rx_1_z_1_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 1, 1, a, 0, 0)*sp[3]
    ss_rx_1_z_0_d_0_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 0, a, 0, 0)*sp[1]
    ss_rx_1_z_0_d_1_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 1, a, 0, 0)*sp[2]
    
    ss_mat_rx_1 = cbind(ss_rx_1_z_1_x_11, ss_rx_1_z_0_d_0_x_11, ss_rx_1_z_0_d_1_x_11,
                        ss_rx_1_z_1_x_10, ss_rx_1_z_0_d_0_x_10, ss_rx_1_z_0_d_1_x_10,
                        ss_rx_1_z_1_x_01, ss_rx_1_z_0_d_0_x_01, ss_rx_1_z_0_d_1_x_01,
                        ss_rx_1_z_1_x_00, ss_rx_1_z_0_d_0_x_00, ss_rx_1_z_0_d_1_x_00)
    
    # rx = 0; not selected in phase II
    ss_rx_0_z_1_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 1, 1, a, 1, NA)*(1-sp[7])
    ss_rx_0_z_0_d_0_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 0, 0, a, 1, NA)*(1-sp[5])
    ss_rx_0_z_0_d_1_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 0, 1, a, 1, NA)*(1-sp[6])
    
    ss_rx_0_z_1_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 1, 1, a, 0, NA)*(1-sp[3])
    ss_rx_0_z_0_d_0_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 0, 0, a, 0, NA)*(1-sp[1])
    ss_rx_0_z_0_d_1_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 0, 1, a, 0, NA)*(1-sp[2])
    
    ss_mat_rx_0 = cbind( ss_rx_0_z_1_x1_1, ss_rx_0_z_0_d_0_x1_1, ss_rx_0_z_0_d_1_x1_1,
                         ss_rx_0_z_1_x1_0, ss_rx_0_z_0_d_0_x1_0, ss_rx_0_z_0_d_1_x1_0 )
    
    res0 <- sum( (drop( jp_mat %*% t(ss_mat_rx_1) ) + drop( (jp_mat[,c(1:3, 7:9)] + jp_mat[,3+c(1:3, 7:9)]) %*% t(ss_mat_rx_0) ) )*dgamma(a, sh0, sc0) )
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  # rz = 0 
  int0.f <- function(sp, a){
    
    # joint pdf of P_Y_X_g_A_rz_csd_pwc.f(rz, z, d, t, x1, x2, para, brks) 
    p_z_0_d_0_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 1, 1, para, brks) *(1-p1.rz_g_x1_1)
    p_d_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 1, 1, para, brks) *(1-p1.rz_g_x1_1)
    
    p_z_0_d_0_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 1, 0, para, brks) *(1-p1.rz_g_x1_1)
    p_d_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 1, 0, para, brks) *(1-p1.rz_g_x1_1)
    
    p_z_0_d_0_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 0, 1, para, brks) *(1-p1.rz_g_x1_0)
    p_d_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 0, 1, para, brks) *(1-p1.rz_g_x1_0)
    
    p_z_0_d_0_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 0, 0, para, brks) *(1-p1.rz_g_x1_0)
    p_d_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 0, 0, para, brks) *(1-p1.rz_g_x1_0)
    
    jp_mat = cbind(p_z_0_d_0_x_11_rz, p_d_1_x_11_rz,
                   p_z_0_d_0_x_10_rz, p_d_1_x_10_rz,
                   p_z_0_d_0_x_01_rz, p_d_1_x_01_rz,
                   p_z_0_d_0_x_00_rz, p_d_1_x_00_rz)
    
    
    # score^2: SS_csd_pwc_ele.f(grad=1e-06, pos, para, brks, r, rz, z, d, t, x1, x2)
    
    # rx = 1; selected in phase II
    ss_rx_1_z_0_d_0_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, 0, 0, a, 1, 1)*sp[5]
    ss_rx_1_d_1_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, NA, 1, a, 1, 1)*sp[8]
    
    ss_rx_1_z_0_d_0_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, 0, 0, a, 1, 0)*sp[5]
    ss_rx_1_d_1_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, NA, 1, a, 1, 0)*sp[8]
    
    ss_rx_1_z_0_d_0_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, 0, 0, a, 0, 1)*sp[1]
    ss_rx_1_d_1_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, NA, 1, a, 0, 1)*sp[4]
    
    ss_rx_1_z_0_d_0_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, 0, 0, a, 0, 0)*sp[1]
    ss_rx_1_d_1_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, NA, 1, a, 0, 0)*sp[4]
    
    ss_mat_rx_1 = cbind(ss_rx_1_z_0_d_0_x_11, ss_rx_1_d_1_x_11,
                        ss_rx_1_z_0_d_0_x_10, ss_rx_1_d_1_x_10,
                        ss_rx_1_z_0_d_0_x_01, ss_rx_1_d_1_x_01,
                        ss_rx_1_z_0_d_0_x_00, ss_rx_1_d_1_x_00)
    
    # rx = 0; not selected in phase II
    ss_rx_0_z_0_d_0_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 0, 0, 0, a, 1, NA)*(1-sp[5])
    ss_rx_0_d_1_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 0, NA, 1, a, 1, NA)*(1-sp[8])
    
    ss_rx_0_z_0_d_0_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 0, 0, 0, a, 0, NA)*(1-sp[1])
    ss_rx_0_d_1_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 0, NA, 1, a, 0, NA)*(1-sp[4])
    
    ss_mat_rx_0 = cbind( ss_rx_0_z_0_d_0_x1_1, ss_rx_0_d_1_x1_1,
                         ss_rx_0_z_0_d_0_x1_0, ss_rx_0_d_1_x1_0 )
    
    res0 <- sum( (drop( jp_mat %*% t(ss_mat_rx_1) ) + drop( (jp_mat[,c(1:2, 5:6)] + jp_mat[,2+c(1:2, 5:6)]) %*% t(ss_mat_rx_0) )  )*dgamma(a, sh0, sc0) )
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  
  # # a.cuts: 0-1
  # res1 <- (pcubature(int1.f, lower = 0, upper=cuts[1],sp=sp.list[[1]]$r1_1)$integral + 
  #            pcubature(int0.f, lower = 0, upper=cuts[1], sp=sp.list[[1]]$r1_0)$integral )
  # # a.cuts: 1-2
  # res2 <- (pcubature(int1.f, lower = cuts[1], upper=cuts[2], sp=sp.list[[2]]$r1_1)$integral + 
  #            pcubature(int0.f, lower = cuts[1], upper=cuts[2], sp=sp.list[[2]]$r1_0)$integral)
  # # a.cuts: 2-3
  # res3 <- (pcubature(int1.f, lower = cuts[2], upper=Inf, sp=sp.list[[3]]$r1_1)$integral + 
  #            pcubature(int0.f, lower = cuts[2], upper=Inf, sp=sp.list[[3]]$r1_0)$integral ) #+ int0_2.f(sp03)
  # res = res1+res2+res3
  
  res = (pcubature(int1.f, lower = 0, upper=100,sp=sp)$integral + 
           pcubature(int0.f, lower = 0, upper=100, sp=sp)$integral )
  
  return(res)
}

# gaussian legendre 
# prepare gaussleg.f(): reference of Jooyoung's code
gaussleg.f <- function(n, x1, x2) {
  EPS <- 3e-14
  
  m <- (n + 1)/2
  xm <- 0.5 * (x2 + x1)
  xl <- 0.5 * (x2 - x1)
  
  x <- rep(0, n)
  w <- rep(0, n)
  for (i in 1:m) {
    z <- cos(pi * (i - 0.25)/(n + 0.5))
    
    tol <- 9999
    while (tol > EPS) {
      p1 <- 1
      p2 <- 0
      for (j in 1:n) {
        p3 <- p2
        p2 <- p1
        p1 <- ((2 * j - 1) * z * p2 - (j - 1) * p3)/j
      }
      
      pp <- n * (z * p1 - p2)/(z * z - 1)
      z1 <- z
      z <- z1 - p1/pp
      
      tol <- abs(z - z1)
      if (tol <= EPS) {
        break
      }
    }
    
    x[i] = xm - xl * z
    x[n + 1 - i] = xm + xl * z
    w[i] = (2 * xl)/((1 - z * z) * pp * pp)
    w[n + 1 - i] = w[i]
  }
  
  return(as.matrix(data.frame(location = x, weight = w)))
}


calSS_ele_approx2.f <- function(grad, pos, sp, #sp11, sp12, sp13, sp01, sp02, sp03, 
                                para, brks,sc0, sh0, para.rz){
  
  # p1.rz
  p1.rz_g_x1_0 = expit.f(para.rz[1])
  p1.rz_g_x1_1 = expit.f(sum(para.rz))
  
  
  
  # rz = 1: need to 
  int1.f <- function(sp, a){
    
    #a <- ifelse(a>10, 10, a)
    
    # joint pdf of P_Y_X_g_A_rz_csd_pwc.f(rz, z, d, t, x1, x2, para, brks) 
    p_z_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 1, 1, para, brks)*p1.rz_g_x1_1
    p_z_0_d_0_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 1, 1, para, brks)*p1.rz_g_x1_1
    p_z_0_d_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 1, 1, para, brks)*p1.rz_g_x1_1
    
    p_z_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 1, 0, para, brks)*p1.rz_g_x1_1
    p_z_0_d_0_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 1, 0, para, brks)*p1.rz_g_x1_1
    p_z_0_d_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 1, 0, para, brks)*p1.rz_g_x1_1
    
    p_z_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 0, 1, para, brks)*p1.rz_g_x1_0
    p_z_0_d_0_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 0, 1, para, brks)*p1.rz_g_x1_0
    p_z_0_d_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 0, 1, para, brks)*p1.rz_g_x1_0
    
    p_z_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 0, 0, para, brks)*p1.rz_g_x1_0
    p_z_0_d_0_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 0, 0, para, brks)*p1.rz_g_x1_0
    p_z_0_d_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 0, 0, para, brks)*p1.rz_g_x1_0
    
    jp_mat = cbind(p_z_1_x_11_rz, p_z_0_d_0_x_11_rz, p_z_0_d_1_x_11_rz,
                   p_z_1_x_10_rz, p_z_0_d_0_x_10_rz, p_z_0_d_1_x_10_rz,
                   p_z_1_x_01_rz, p_z_0_d_0_x_01_rz, p_z_0_d_1_x_01_rz,
                   p_z_1_x_00_rz, p_z_0_d_0_x_00_rz, p_z_0_d_1_x_00_rz)
    
    
    # score^2: SS_csd_pwc_ele.f(grad=1e-06, pos, para, r, rz, z, d, t, x1, x2)
    
    # rx = 1; selected in phase II
    ss_rx_1_z_1_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 1, 1, a, 1, 1)*sp[6] # 111
    ss_rx_1_z_0_d_0_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 0, a, 1, 1)*sp[2] #001
    ss_rx_1_z_0_d_1_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 1, a, 1, 1)*sp[4] #011
    
    ss_rx_1_z_1_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 1, 1, a, 1, 0)*sp[6]
    ss_rx_1_z_0_d_0_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 0, a, 1, 0)*sp[2]
    ss_rx_1_z_0_d_1_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 1, a, 1, 0)*sp[4]
    
    ss_rx_1_z_1_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 1, 1, a, 0, 1)*sp[5] # 110
    ss_rx_1_z_0_d_0_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 0, a, 0, 1)*sp[1] #000
    ss_rx_1_z_0_d_1_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 1, a, 0, 1)*sp[3] #010
    
    ss_rx_1_z_1_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 1, 1, a, 0, 0)*sp[5]
    ss_rx_1_z_0_d_0_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 0, a, 0, 0)*sp[1]
    ss_rx_1_z_0_d_1_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 1, 0, 1, a, 0, 0)*sp[3]
    
    ss_mat_rx_1 = cbind(ss_rx_1_z_1_x_11, ss_rx_1_z_0_d_0_x_11, ss_rx_1_z_0_d_1_x_11,
                        ss_rx_1_z_1_x_10, ss_rx_1_z_0_d_0_x_10, ss_rx_1_z_0_d_1_x_10,
                        ss_rx_1_z_1_x_01, ss_rx_1_z_0_d_0_x_01, ss_rx_1_z_0_d_1_x_01,
                        ss_rx_1_z_1_x_00, ss_rx_1_z_0_d_0_x_00, ss_rx_1_z_0_d_1_x_00)
    
    # rx = 0; not selected in phase II
    ss_rx_0_z_1_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 1, 1, a, 1, NA)*(1-sp[6])
    ss_rx_0_z_0_d_0_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 0, 0, a, 1, NA)*(1-sp[2])
    ss_rx_0_z_0_d_1_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 0, 1, a, 1, NA)*(1-sp[4])
    
    ss_rx_0_z_1_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 1, 1, a, 0, NA)*(1-sp[5])
    ss_rx_0_z_0_d_0_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 0, 0, a, 0, NA)*(1-sp[1])
    ss_rx_0_z_0_d_1_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 1, 0, 1, a, 0, NA)*(1-sp[3])
    
    ss_mat_rx_0 = cbind( ss_rx_0_z_1_x1_1, ss_rx_0_z_0_d_0_x1_1, ss_rx_0_z_0_d_1_x1_1,
                         ss_rx_0_z_1_x1_0, ss_rx_0_z_0_d_0_x1_0, ss_rx_0_z_0_d_1_x1_0 )
    
    res0 <- sum( ( drop( jp_mat %*% t(ss_mat_rx_1) ) + drop( (jp_mat[,c(1:3, 7:9)] + jp_mat[,3+c(1:3, 7:9)]) %*% t(ss_mat_rx_0) ) )*dgamma(a, sh0, sc0) )
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  # rz = 0 
  int0.f <- function(sp, a){
    
    #a <- ifelse(a>10, 10, a)
    
    # joint pdf of P_Y_X_g_A_rz_csd_pwc.f(rz, z, d, t, x1, x2, para, brks) 
    p_z_0_d_0_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 1, 1, para, brks)*(1-p1.rz_g_x1_1)
    p_d_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 1, 1, para, brks)*(1-p1.rz_g_x1_1)
    
    p_z_0_d_0_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 1, 0, para, brks)*(1-p1.rz_g_x1_1)
    p_d_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 1, 0, para, brks)*(1-p1.rz_g_x1_1)
    
    p_z_0_d_0_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 0, 1, para, brks)*(1-p1.rz_g_x1_0)
    p_d_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 0, 1, para, brks)*(1-p1.rz_g_x1_0)
    
    p_z_0_d_0_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 0, 0, para, brks)*(1-p1.rz_g_x1_0)
    p_d_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 0, 0, para, brks)*(1-p1.rz_g_x1_0)
    
    jp_mat = cbind(p_z_0_d_0_x_11_rz, p_d_1_x_11_rz,
                   p_z_0_d_0_x_10_rz, p_d_1_x_10_rz,
                   p_z_0_d_0_x_01_rz, p_d_1_x_01_rz,
                   p_z_0_d_0_x_00_rz, p_d_1_x_00_rz)
    
    
    # score^2: SS_csd_pwc_ele.f(grad=1e-06, pos, para, r, rz, z, d, t, x1, x2)
    
    # rx = 1; selected in phase II
    ss_rx_1_z_0_d_0_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, 0, 0, a, 1, 1)*sp[2] #001
    ss_rx_1_d_1_x_11 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, NA, 1, a, 1, 1)*sp[8] #211
    
    ss_rx_1_z_0_d_0_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, 0, 0, a, 1, 0)*sp[2]
    ss_rx_1_d_1_x_10 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, NA, 1, a, 1, 0)*sp[8]
    
    ss_rx_1_z_0_d_0_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, 0, 0, a, 0, 1)*sp[1] #000
    ss_rx_1_d_1_x_01 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, NA, 1, a, 0, 1)*sp[7] #210
    
    ss_rx_1_z_0_d_0_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, 0, 0, a, 0, 0)*sp[1]
    ss_rx_1_d_1_x_00 = SS_csd_pwc_ele.f(grad, pos, para, brks, 1, 0, NA, 1, a, 0, 0)*sp[7]
    
    ss_mat_rx_1 = cbind(ss_rx_1_z_0_d_0_x_11, ss_rx_1_d_1_x_11,
                        ss_rx_1_z_0_d_0_x_10, ss_rx_1_d_1_x_10,
                        ss_rx_1_z_0_d_0_x_01, ss_rx_1_d_1_x_01,
                        ss_rx_1_z_0_d_0_x_00, ss_rx_1_d_1_x_00)
    
    # rx = 0; not selected in phase II
    ss_rx_0_z_0_d_0_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 0, 0, 0, a, 1, NA)*(1-sp[2])
    ss_rx_0_d_1_x1_1 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 0, NA, 1, a, 1, NA)*(1-sp[8])
    
    ss_rx_0_z_0_d_0_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 0, 0, 0, a, 0, NA)*(1-sp[1])
    ss_rx_0_d_1_x1_0 = SS_csd_pwc_ele.f(grad, pos, para, brks, 0, 0, NA, 1, a, 0, NA)*(1-sp[7])
    
    ss_mat_rx_0 = cbind( ss_rx_0_z_0_d_0_x1_1, ss_rx_0_d_1_x1_1,
                         ss_rx_0_z_0_d_0_x1_0, ss_rx_0_d_1_x1_0 )
    
    res0 <-  sum(  ( drop( jp_mat %*% t(ss_mat_rx_1) ) + drop( (jp_mat[,c(1:2, 5:6)] + jp_mat[,2+c(1:2, 5:6)]) %*% t(ss_mat_rx_0) )  ) *dgamma(a, sh0, sc0) 
                 )
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  
  # # a.cuts: 0-1
  # res1 <- (pcubature(int1.f, lower = 0, upper=cuts[1],sp=sp.list[[1]])$integral + 
  #            pcubature(int0.f, lower = 0, upper=cuts[1], sp=sp.list[[1]])$integral )
  # # a.cuts: 1-2
  # res2 <- (pcubature(int1.f, lower = cuts[1], upper=cuts[2], sp=sp.list[[2]])$integral + 
  #            pcubature(int0.f, lower = cuts[1], upper=cuts[2], sp=sp.list[[2]])$integral)
  # # a.cuts: 2-3
  # res3 <- (pcubature(int1.f, lower = cuts[2], upper=Inf, sp=sp.list[[3]])$integral + 
  #            pcubature(int0.f, lower = cuts[2], upper=Inf, sp=sp.list[[3]])$integral ) #+ int0_2.f(sp03)
  # 
  # res = res1+res2+res3#+res4
  
   #vres1 = pcubature(int1.f, lower = 0, upper=100,sp=sp)$integral
   
   #vres0 = pcubature(int0.f, lower = 0, upper=100, sp=sp)$integral 
  
  uu=gaussleg.f(20,0,1)
  
  res11=sapply(1:20, function(i){
    int1.f(sp, uu[i,1])*uu[i,2]
  })#sum(int1.f(sp, uu[,1])*uu[,2])
  res10=sapply(1:20, function(i){
    int0.f(sp, uu[i,1])*uu[i,2]
  })#sum(int0.f(sp, uu[,1])*uu[,2])
  
  res21=sapply(1:20, function(i){
    int1.f(sp,1/uu[i,1])*(uu[i,1])^(-2)*uu[i,2]
  })#sum(int1.f(sp,1/uu[,1])*(uu[,1])^(-2)*uu[,2])
  res20=sapply(1:20, function(i){
    int0.f(sp,1/uu[i,1])*(uu[i,1])^(-2)*uu[i,2]
  })#sum(int0.f(sp,1/uu[,1])*(uu[,1])^(-2)*uu[,2])
  
  res=sum(res11+res10+res21+res20)
  
  return(res)
}


calS1S2_ele_approx2.f <- function(grad, pos1, pos2, sp, #sp11, sp12, sp13, sp01, sp02, sp03, 
                                  para, brks, sc0, sh0, para.rz){
  
  # p1.rz
  p1.rz_g_x1_0 = expit.f(para.rz[1])
  p1.rz_g_x1_1 = expit.f(sum(para.rz))
  
  # rz = 1: need to 
  int1.f <- function(sp, a){
    
    #a <- ifelse(a>10, 10, a)
    
    # joint pdf of P_Y_X_g_A_rz_csd_pwc.f(rz, z, d, t, x1, x2, para, brks) 
    p_z_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 1, 1, para, brks) *p1.rz_g_x1_1
    p_z_0_d_0_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 1, 1, para, brks) *p1.rz_g_x1_1
    p_z_0_d_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 1, 1, para, brks) *p1.rz_g_x1_1
    
    p_z_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 1, 0, para, brks) *p1.rz_g_x1_1
    p_z_0_d_0_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 1, 0, para, brks) *p1.rz_g_x1_1
    p_z_0_d_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 1, 0, para, brks) *p1.rz_g_x1_1
    
    p_z_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 0, 1, para, brks) *p1.rz_g_x1_0
    p_z_0_d_0_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 0, 1, para, brks) *p1.rz_g_x1_0
    p_z_0_d_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 0, 1, para, brks) *p1.rz_g_x1_0
    
    p_z_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 1, 1, a, 0, 0, para, brks) *p1.rz_g_x1_0
    p_z_0_d_0_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 0, a, 0, 0, para, brks) *p1.rz_g_x1_0
    p_z_0_d_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(1, 0, 1, a, 0, 0, para, brks) *p1.rz_g_x1_0
    
    jp_mat = cbind(p_z_1_x_11_rz, p_z_0_d_0_x_11_rz, p_z_0_d_1_x_11_rz,
                   p_z_1_x_10_rz, p_z_0_d_0_x_10_rz, p_z_0_d_1_x_10_rz,
                   p_z_1_x_01_rz, p_z_0_d_0_x_01_rz, p_z_0_d_1_x_01_rz,
                   p_z_1_x_00_rz, p_z_0_d_0_x_00_rz, p_z_0_d_1_x_00_rz)
    
    
    # score^2: SS_csd_pwc_ele.f(grad=1e-06, pos, para, brks, r, rz, z, d, t, x1, x2)
    
    # rx = 1; selected in phase II
    ss_rx_1_z_1_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 1, 1, a, 1, 1)*sp[6] #111
    ss_rx_1_z_0_d_0_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 0, a, 1, 1)*sp[2] #001
    ss_rx_1_z_0_d_1_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 1, a, 1, 1)*sp[4] #011
    
    ss_rx_1_z_1_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 1, 1, a, 1, 0)*sp[6]
    ss_rx_1_z_0_d_0_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 0, a, 1, 0)*sp[2]
    ss_rx_1_z_0_d_1_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 1, a, 1, 0)*sp[4]
    
    ss_rx_1_z_1_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 1, 1, a, 0, 1)*sp[5]
    ss_rx_1_z_0_d_0_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 0, a, 0, 1)*sp[1]
    ss_rx_1_z_0_d_1_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 1, a, 0, 1)*sp[3]
    
    ss_rx_1_z_1_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 1, 1, a, 0, 0)*sp[5]
    ss_rx_1_z_0_d_0_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 0, a, 0, 0)*sp[1]
    ss_rx_1_z_0_d_1_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 1, 0, 1, a, 0, 0)*sp[3]
    
    ss_mat_rx_1 = cbind(ss_rx_1_z_1_x_11, ss_rx_1_z_0_d_0_x_11, ss_rx_1_z_0_d_1_x_11,
                        ss_rx_1_z_1_x_10, ss_rx_1_z_0_d_0_x_10, ss_rx_1_z_0_d_1_x_10,
                        ss_rx_1_z_1_x_01, ss_rx_1_z_0_d_0_x_01, ss_rx_1_z_0_d_1_x_01,
                        ss_rx_1_z_1_x_00, ss_rx_1_z_0_d_0_x_00, ss_rx_1_z_0_d_1_x_00)
    
    # rx = 0; not selected in phase II
    ss_rx_0_z_1_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 1, 1, a, 1, NA)*(1-sp[6])
    ss_rx_0_z_0_d_0_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 0, 0, a, 1, NA)*(1-sp[2])
    ss_rx_0_z_0_d_1_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 0, 1, a, 1, NA)*(1-sp[4])
    
    ss_rx_0_z_1_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 1, 1, a, 0, NA)*(1-sp[5])
    ss_rx_0_z_0_d_0_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 0, 0, a, 0, NA)*(1-sp[1])
    ss_rx_0_z_0_d_1_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 1, 0, 1, a, 0, NA)*(1-sp[3])
    
    ss_mat_rx_0 = cbind( ss_rx_0_z_1_x1_1, ss_rx_0_z_0_d_0_x1_1, ss_rx_0_z_0_d_1_x1_1,
                         ss_rx_0_z_1_x1_0, ss_rx_0_z_0_d_0_x1_0, ss_rx_0_z_0_d_1_x1_0 )
    
    res0 <- sum( (drop( jp_mat %*% t(ss_mat_rx_1) ) + drop( (jp_mat[,c(1:3, 7:9)] + jp_mat[,3+c(1:3, 7:9)]) %*% t(ss_mat_rx_0) ) )*dgamma(a, sh0, sc0) )
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  # rz = 0 
  int0.f <- function(sp, a){
    
    #a <- ifelse(a>10, 10, a)
    
    # joint pdf of P_Y_X_g_A_rz_csd_pwc.f(rz, z, d, t, x1, x2, para, brks) 
    p_z_0_d_0_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 1, 1, para, brks) *(1-p1.rz_g_x1_1)
    p_d_1_x_11_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 1, 1, para, brks) *(1-p1.rz_g_x1_1)
    
    p_z_0_d_0_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 1, 0, para, brks) *(1-p1.rz_g_x1_1)
    p_d_1_x_10_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 1, 0, para, brks) *(1-p1.rz_g_x1_1)
    
    p_z_0_d_0_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 0, 1, para, brks) *(1-p1.rz_g_x1_0)
    p_d_1_x_01_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 0, 1, para, brks) *(1-p1.rz_g_x1_0)
    
    p_z_0_d_0_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(0, 0, 0, a, 0, 0, para, brks) *(1-p1.rz_g_x1_0)
    p_d_1_x_00_rz = P_Y_X_g_A_rz_csd_pwc.f(0, NA, 1, a, 0, 0, para, brks) *(1-p1.rz_g_x1_0)
    
    jp_mat = cbind(p_z_0_d_0_x_11_rz, p_d_1_x_11_rz,
                   p_z_0_d_0_x_10_rz, p_d_1_x_10_rz,
                   p_z_0_d_0_x_01_rz, p_d_1_x_01_rz,
                   p_z_0_d_0_x_00_rz, p_d_1_x_00_rz)
    
    
    # score^2: SS_csd_pwc_ele.f(grad=1e-06, pos, para, brks, r, rz, z, d, t, x1, x2)
    
    # rx = 1; selected in phase II
    ss_rx_1_z_0_d_0_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, 0, 0, a, 1, 1)*sp[2] #001
    ss_rx_1_d_1_x_11 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, NA, 1, a, 1, 1)*sp[8] #211
    
    ss_rx_1_z_0_d_0_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, 0, 0, a, 1, 0)*sp[2]
    ss_rx_1_d_1_x_10 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, NA, 1, a, 1, 0)*sp[8]
    
    ss_rx_1_z_0_d_0_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, 0, 0, a, 0, 1)*sp[1]
    ss_rx_1_d_1_x_01 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, NA, 1, a, 0, 1)*sp[7]
    
    ss_rx_1_z_0_d_0_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, 0, 0, a, 0, 0)*sp[1]
    ss_rx_1_d_1_x_00 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 1, 0, NA, 1, a, 0, 0)*sp[7]
    
    ss_mat_rx_1 = cbind(ss_rx_1_z_0_d_0_x_11, ss_rx_1_d_1_x_11,
                        ss_rx_1_z_0_d_0_x_10, ss_rx_1_d_1_x_10,
                        ss_rx_1_z_0_d_0_x_01, ss_rx_1_d_1_x_01,
                        ss_rx_1_z_0_d_0_x_00, ss_rx_1_d_1_x_00)
    
    # rx = 0; not selected in phase II
    ss_rx_0_z_0_d_0_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 0, 0, 0, a, 1, NA)*(1-sp[2])
    ss_rx_0_d_1_x1_1 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 0, NA, 1, a, 1, NA)*(1-sp[8])
    
    ss_rx_0_z_0_d_0_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 0, 0, 0, a, 0, NA)*(1-sp[1])
    ss_rx_0_d_1_x1_0 = S1S2_csd_pwc_ele.f(grad, pos1, pos2, para, brks, 0, 0, NA, 1, a, 0, NA)*(1-sp[7])
    
    ss_mat_rx_0 = cbind( ss_rx_0_z_0_d_0_x1_1, ss_rx_0_d_1_x1_1,
                         ss_rx_0_z_0_d_0_x1_0, ss_rx_0_d_1_x1_0 )
    
    res0 <- sum( (drop( jp_mat %*% t(ss_mat_rx_1) ) + drop( (jp_mat[,c(1:2, 5:6)] + jp_mat[,2+c(1:2, 5:6)]) %*% t(ss_mat_rx_0) )  )*dgamma(a, sh0, sc0) )
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  
  # # a.cuts: 0-1
  # res1 <- (pcubature(int1.f, lower = 0, upper=cuts[1],sp=sp.list[[1]]$r1_1)$integral + 
  #            pcubature(int0.f, lower = 0, upper=cuts[1], sp=sp.list[[1]]$r1_0)$integral )
  # # a.cuts: 1-2
  # res2 <- (pcubature(int1.f, lower = cuts[1], upper=cuts[2], sp=sp.list[[2]]$r1_1)$integral + 
  #            pcubature(int0.f, lower = cuts[1], upper=cuts[2], sp=sp.list[[2]]$r1_0)$integral)
  # # a.cuts: 2-3
  # res3 <- (pcubature(int1.f, lower = cuts[2], upper=Inf, sp=sp.list[[3]]$r1_1)$integral + 
  #            pcubature(int0.f, lower = cuts[2], upper=Inf, sp=sp.list[[3]]$r1_0)$integral ) #+ int0_2.f(sp03)
  # res = res1+res2+res3
  
  # res = (pcubature(int1.f, lower = 0, upper=100,sp=sp)$integral +
  #         pcubature(int0.f, lower = 0, upper=100, sp=sp)$integral )
  
  uu=gaussleg.f(20,0,1)
  
  res11=sapply(1:20, function(i){
    int1.f(sp, uu[i,1])*uu[i,2]
  })#sum(int1.f(sp, uu[,1])*uu[,2])
  res10=sapply(1:20, function(i){
    int0.f(sp, uu[i,1])*uu[i,2]
  })#sum(int0.f(sp, uu[,1])*uu[,2])
  
  res21=sapply(1:20, function(i){
    int1.f(sp,1/uu[i,1])*(uu[i,1])^(-2)*uu[i,2]
  })#sum(int1.f(sp,1/uu[,1])*(uu[,1])^(-2)*uu[,2])
  res20=sapply(1:20, function(i){
    int0.f(sp,1/uu[i,1])*(uu[i,1])^(-2)*uu[i,2]
  })#sum(int0.f(sp,1/uu[,1])*(uu[,1])^(-2)*uu[,2])
  
  res=sum(res11+res10+res21+res20)
  
  return(res)
}

info_mat.f <- function(pos_vec, grad, sp, para, brks, sc0, sh0, para.rz){
  
  
  p = length(pos_vec)
  res <- matrix(NA, p, p)
  
  run0=expand.grid(pos_vec, pos_vec)
  run1=run0[run0[,1]<run0[,2],]
  run2=run0[run0[,1]==run0[,2],]
  
  vec1 = mclapply(1:nrow(run1), function(i){

    x=as.vector(unlist(run1[i,]))
    calS1S2_ele_approx2.f(grad, x[1], x[2], sp, para, brks, sc0, sh0, para.rz)
  
  }, mc.cores=8)
  
  vec2 = mclapply(1:nrow(run2), function(i){
    # if(x[1]==x[2]){
    #   calSS_ele_approx.f(grad, x[1], cuts, sp11, sp12, sp13, sp01, sp02, sp03, para, brks, sc0, sh0, pv)
    # }else{
    x=as.vector(unlist(run2[i,]))
    calSS_ele_approx2.f(grad, x[1], sp, para, brks, sc0, sh0, para.rz)
    #}
  }, mc.cores=8)
  
  diag(res) <- as.vector(unlist(vec2))
  res[upper.tri(res)] <-as.vector(unlist(vec1))
  tres <- t(res)
  tres[upper.tri(tres)] <- as.vector(unlist(vec1))
  
  return(tres)
}


# r.z=1, z/d (3); r.z=0, z/d (2); each with x1 (2)
# list length of 3 crsp to str on A
opt_one_ODS_asym.f <- function(target_pos, grad, n_sample, phaseI_res, ini, para, brks, sc0, sh0, para.rz){
  
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  #len=length(ini)/3
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
  
    print(pi_mat)
    
    #info = E_DR.Abar_AY(pi_vec=pi_mat, cuts=phaseI_res$Abar_cut, beta=beta, eta=eta, P2=P2, param.A=param.A)
    #var_beta1 = sqrt(diag(ginv(info))/n)[2]
    
    info = info_mat.f(pos_vec=1:length(para),grad=grad,pi_mat, para, brks, sc0, sh0, para.rz)
    
    se_eta1 = sqrt(diag(ginv(info))/n)[target_pos]
    #se_eta1 = sqrt(1/var_eta1/n)
    #se_eta1 <- ifelse(is.nan(se_eta1), 100, se_eta1)
    
    print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=logit(ini), objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-05, tolFun = 1e-05)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}

optA_ODS_asym.f <- function(target_pos_vec, grad, n_sample, phaseI_res, ini, para, brks, sc0, sh0, para.rz){
  
  pos_vec = 1:length(para)
  
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  #len=length(ini)/3
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    print(pi_mat)
    
    #info = E_DR.Abar_AY(pi_vec=pi_mat, cuts=phaseI_res$Abar_cut, beta=beta, eta=eta, P2=P2, param.A=param.A)
    #var_beta1 = sqrt(diag(ginv(info))/n)[2]
    
    # info = info_mat.f(pos_vec=pos_vec,grad=grad, cuts=phaseI_res$cuts, 
    #                   pi_mat, 
    #                   #list(pi_mat[1:len], pi_mat[len+(1:len)], pi_mat[2*len+(1:len)]),
    #                   para, brks, sc0, sh0, para.rz)
    
    info = info_mat.f(pos_vec=1:length(para),grad=grad,pi_mat, para, brks, sc0, sh0, para.rz)
    
    se_eta1 = sum(diag(ginv(info*n))[target_pos_vec])
    #se_eta1 = sqrt(1/var_eta1/n)
    #se_eta1 <- ifelse(is.nan(se_eta1), 100, se_eta1)
    
    print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-05, tolFun = 1e-05)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}

optD_ODS_asym.f <- function(target_pos_vec, grad, n_sample, phaseI_res, ini, para, brks, sc0, sh0, para.rz){
  
  pos_vec = 1:length(para)
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  #len=length(ini)/3
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    print(pi_mat)
    
    #info = E_DR.Abar_AY(pi_vec=pi_mat, cuts=phaseI_res$Abar_cut, beta=beta, eta=eta, P2=P2, param.A=param.A)
    #var_beta1 = sqrt(diag(ginv(info))/n)[2]
    
    # info = info_mat.f(pos_vec=pos_vec,grad=grad, cuts=phaseI_res$cuts, 
    #                   pi_mat, 
    #                   #list(pi_mat[1:len], pi_mat[len+(1:len)], pi_mat[2*len+(1:len)]),
    #                   #sp11=pi_mat[4], sp12=pi_mat[5], sp13=pi_mat[6], sp01=pi_mat[1], sp02=pi_mat[2], sp03=pi_mat[3], 
    #                   para, brks, sc0, sh0, para.rz)
    info = info_mat.f(pos_vec=1:length(para),grad=grad,pi_mat, para, brks, sc0, sh0, para.rz)
    
    se_eta1 = det(ginv(info*n)[target_pos_vec, target_pos_vec])
    #se_eta1 = sqrt(1/var_eta1/n)
    #se_eta1 <- ifelse(is.nan(se_eta1), 100, se_eta1)
    
    print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-05, tolFun = 1e-05)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}


opt_cumt_x2_asym.f <- function(t, x2, grad, n_sample, phaseI_res, ini, para, brks, sc0, sh0, para.rz){
  
  pos_vec = 1:length(para)
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  #len=length(ini)/3
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    print(pi_mat)
    
    #info = E_DR.Abar_AY(pi_vec=pi_mat, cuts=phaseI_res$Abar_cut, beta=beta, eta=eta, P2=P2, param.A=param.A)
    #var_beta1 = sqrt(diag(ginv(info))/n)[2]
    
    info = info_mat.f(pos_vec=pos_vec,grad=grad, 
                      pi_mat, 
                      para, brks, sc0, sh0, para.rz)
    
    #se_eta1 = det(ginv(info)[target_pos_vec, target_pos_vec])
    
    se_eta1 = ase.cumrisk_g_x2_pwc.f(grad, ginv(info)/n, t, x2, para, brks, c(0,0))$ase*100
    
    #print(logit.pi_mat)
    #print(se_eta1)
    
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-05, tolFun = 1e-05)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}

opt_cum0_x2_1_asym.f <- function(grad, n_sample, phaseI_res, ini, para, brks, sc0, sh0, para.rz){
  
  pos_vec = 1:length(para)
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  #len=length(ini)/3
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    print(pi_mat)
    
    #info = E_DR.Abar_AY(pi_vec=pi_mat, cuts=phaseI_res$Abar_cut, beta=beta, eta=eta, P2=P2, param.A=param.A)
    #var_beta1 = sqrt(diag(ginv(info))/n)[2]
    
    info = info_mat.f(pos_vec=pos_vec,grad=grad, 
                      pi_mat, 
                      para, brks, sc0, sh0, para.rz)

    #se_eta1 = det(ginv(info)[target_pos_vec, target_pos_vec])
    
    se_eta1 = ase.cumrisk_g_x2_pwc.f(grad, ginv(info)/n, 0, 1, para, brks, c(0,0))$ase
    
    #se_eta1 = sqrt(1/var_eta1/n)
    #se_eta1 <- ifelse(is.nan(se_eta1), 100, se_eta1)
    
    print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-07, tolFun = 1e-07)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}

opt_cum3_x2_1_asym.f <- function(grad, n_sample, phaseI_res, ini, para, brks, sc0, sh0, para.rz){
  
  pos_vec = 1:length(para)
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  #len=length(ini)/3
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    print(pi_mat)
    
    #info = E_DR.Abar_AY(pi_vec=pi_mat, cuts=phaseI_res$Abar_cut, beta=beta, eta=eta, P2=P2, param.A=param.A)
    #var_beta1 = sqrt(diag(ginv(info))/n)[2]
    
    info = info_mat.f(pos_vec=pos_vec,grad=grad, 
                      pi_mat, 
                      para, brks, sc0, sh0, para.rz)
    
    #se_eta1 = det(ginv(info)[target_pos_vec, target_pos_vec])
    
    se_eta1 = ase.cumrisk_g_x2_pwc.f(grad, ginv(info)/n, 3, 1, para, brks, c(0,0))$ase
    
    #se_eta1 = sqrt(1/var_eta1/n)
    #se_eta1 <- ifelse(is.nan(se_eta1), 100, se_eta1)
    
    print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-07, tolFun = 1e-07)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}

opt_cum0_x2_0_asym.f <- function(grad, n_sample, phaseI_res, ini, para, brks, sc0, sh0, para.rz){
  
  pos_vec = 1:length(para)
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  #len=length(ini)/3
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    print(pi_mat)
    
    #info = E_DR.Abar_AY(pi_vec=pi_mat, cuts=phaseI_res$Abar_cut, beta=beta, eta=eta, P2=P2, param.A=param.A)
    #var_beta1 = sqrt(diag(ginv(info))/n)[2]
    
    info = info_mat.f(pos_vec=pos_vec,grad=grad, 
                      pi_mat, 
                      para, brks, sc0, sh0, para.rz)
    
    #se_eta1 = det(ginv(info)[target_pos_vec, target_pos_vec])
    
    se_eta1 = ase.cumrisk_g_x2_pwc.f(grad, ginv(info)/n, 0, 0, para, brks, c(0,0))$ase
    
    #se_eta1 = sqrt(1/var_eta1/n)
    #se_eta1 <- ifelse(is.nan(se_eta1), 100, se_eta1)
    
    print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-07, tolFun = 1e-07)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}


opt_cum3_x2_0_asym.f <- function(grad, n_sample, phaseI_res, ini, para, brks, sc0, sh0, para.rz){
  
  pos_vec = 1:length(para)
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  #len=length(ini)/3
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    print(pi_mat)
    
    #info = E_DR.Abar_AY(pi_vec=pi_mat, cuts=phaseI_res$Abar_cut, beta=beta, eta=eta, P2=P2, param.A=param.A)
    #var_beta1 = sqrt(diag(ginv(info))/n)[2]
    
    info = info_mat.f(pos_vec=pos_vec,grad=grad, 
                      pi_mat, 
                      para, brks, sc0, sh0, para.rz)
    
    #se_eta1 = det(ginv(info)[target_pos_vec, target_pos_vec])
    
    se_eta1 = ase.cumrisk_g_x2_pwc.f(grad, ginv(info)/n, 3, 0, para, brks, c(0,0))$ase
    
    #se_eta1 = sqrt(1/var_eta1/n)
    #se_eta1 <- ifelse(is.nan(se_eta1), 100, se_eta1)
    
    print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-07, tolFun = 1e-07)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}

# cum0 given x2=0 or 1
#opt_cum0_x2_ODS_asym.f
# cum3 given x2=0 or 1
#opt_cum3_x2_ODS_asym.f
#optA_cum0_ODS_asym.f
#optA_cum3_ODS_asym.f

# test
# para=c(et_vec, log(0.6/0.4), c(-5, -5), be_vec[1:2], ga_vec[1:3])
# #sp.list = lapply(1:3, function(x){list(r1_1=rep(1,6), r1_0=rep(1,4))})
# calSS_ele_approx.f(1e-06, 1, rep(n2_samp/n,8), para, brks, sc0, sh0, para.rz)
# calS1S2_ele_approx.f(1e-06, 1, 2, rep(n2_samp/n,8), para, brks, sc0, sh0, para.rz)
# calSS_ele_approx2.f(1e-06, 1, rep(n2_samp/n,8), para, brks, sc0, sh0, para.rz)
# calS1S2_ele_approx2.f(1e-06, 1, 2, rep(n2_samp/n,8), para, brks, sc0, sh0, para.rz)
# info=info_mat.f( 1:length(para),  1e-06, rep(1,8), para, 1, 1, 2, c(1,2) )
# # sqrt(diag(ginv(infotest))/n)



