//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>

#include "commonf.h"

using namespace Rcpp;
using namespace arma;


// Obtain environment containing function
//Rcpp::Environment package_env("package:mnormt"); 

// Make function callable from C++
//Rcpp::Function bipmvnrm = package_env["pmnorm"];    

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export()]]
double logL_csd_cal_wei(int rz, int z, int d, double t, double x1, double x2,
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); //exp(eta[0]);//
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  double Hf = pow(alp[0]*t, alp[1])*exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2);
  double Sf = exp(-Hf); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(rz==1){
    
    if(d==1){
      logL1=z*log(piZ) + (1-z)*(log(1-Sf)+log(1-piZ));
    }else{
      logL1=(-Hf+log(1-piZ));
    }
    //logL1 = d*z*log(piZ) + d*(1-z)*(log(1-Sf)+log(1-piZ)) + (1-d)*(-Hf+log(1-piZ));
  }else{
    //d*(p.z1+(1-surv.a)*p.z0)+(1-d)*surv.a*p.z0
    logL1 = d*log(piZ+(1-Sf)*(1-piZ))+(1-d)*(-Hf+log(1-piZ));
  }
  
  //double logL2 = x2*(eta[0]+eta[1]*x1) - log(1+exp(eta[0]+eta[1]*x1)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  double logL2 = x1*(x2*(eta[0]+eta[1]) - log(1+exp(eta[0]+eta[1])) + log(p1) ) + (1-x1)*( x2*(eta[0]) - log(1+exp(eta[0])) + log(1-p1) );
  
  double res=logL1+logL2;
  
  return res;
}

// [[Rcpp::export()]]
double logL_csd_obs_wei(int r, int rz, int z, int d, double t, double x1, double x2, 
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  double res=0.0;
  if(r==1){
    res = logL_csd_cal_wei(rz, z, d, t, x1, x2, gam, alp, beta, eta, p1);
  }else{
    double res0 = logL_csd_cal_wei(rz, z, d, t, x1, 0.0, gam, alp, beta, eta, p1);
    double res1 = logL_csd_cal_wei(rz, z, d, t, x1, 1.0, gam, alp, beta, eta, p1);
    res = log(exp(res0)+exp(res1));
  }
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_csd_cal_wei(IntegerVector rz, IntegerVector z, IntegerVector d, NumericVector t, NumericVector x1, NumericVector x2, 
                            NumericVector gam, NumericVector alp, NumericVector beta, 
                            NumericVector eta, double p1){
  
  NumericVector res=no_init(d.size());
  
  for(int i=0; i<d.size(); i++){
    res[i] = logL_csd_cal_wei(rz[i], z[i], d[i], t[i], x1[i], x2[i], gam, alp, beta, eta, p1); //wgt[i]*(temp);
  }
  
  return res;
}

// [[Rcpp::export()]]
double logL_cal_wei(int d1, int d2, double t1, double t2, double x1, double x2,
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); //exp(eta[0]);//
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  double expcovT=exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2);
  double Hft1 = pow(alp[0]*t1, alp[1])*expcovT;
  double Sft1 = exp(-Hft1); 
  double Hft2 = pow(alp[0]*t2, alp[1])*expcovT;
  double Sft2 = exp(-Hft2); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(d1==1){
    logL1=log(piZ);
  }else{
    
    if(d1>1){
      //logL1=d2*log((1-Sft2)*(1-piZ)+piZ)+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log((1-Sft2)*(1-piZ)+piZ);
      }else{
        //logL1=(log(1-piZ)-Hft2);
        logL1=0;
      }
    }else{//d1==0 or 1
      //logL1=d2*(log( Sft1-Sft2 )+log(1-piZ))+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log( Sft1-Sft2 )+log(1-piZ);
      }else{
        logL1=(log(1-piZ)-Hft2);
      }
    }
  }
  
  double logL2 = x1*(x2*(eta[0]+eta[1]) - log(1+exp(eta[0]+eta[1])) + log(p1) ) + (1-x1)*( x2*(eta[0]) - log(1+exp(eta[0])) + log(1-p1) ); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  //double logL2 = x1*x2*log(eta[0]) + x1*(1-x2)*log(eta[1]) + (1-x1)*x2*log(eta[2]) + (1-x1)*(1-x2)*log(1-eta[0]-eta[1]-eta[2]);
  
  double res=logL1+logL2;
  
  return res;
}

// [[Rcpp::export()]]
double logL_cal1_wei(int d1, int d2, double t1, double t2, double x1, double x2,
                    NumericVector gam, NumericVector alp, NumericVector beta){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); //exp(eta[0]);//
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  double expcovT=exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2);
  double Hft1 = pow(alp[0]*t1, alp[1])*expcovT;
  double Sft1 = exp(-Hft1); 
  double Hft2 = pow(alp[0]*t2, alp[1])*expcovT;
  double Sft2 = exp(-Hft2); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(d1==1){
    
    logL1=log(piZ);
  }else{
    
    if(d1>1){
      //logL1=d2*log((1-Sft2)*(1-piZ)+piZ)+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log((1-Sft2)*(1-piZ)+piZ);
      }else{
        logL1=0;//(log(1-piZ)-Hft2);
      }
    }else{//d1==0 or 1
      //logL1=d2*(log( Sft1-Sft2 )+log(1-piZ))+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log( Sft1-Sft2 )+log(1-piZ);
      }else{
        logL1=(log(1-piZ)-Hft2);
      }
    }
  }
  
  //double logL2 = x2*(eta[0]+eta[1]*x1) - log(1+exp(eta[0]+eta[1]*x1)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  
  double res=logL1;
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_cal1_wei(IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, NumericVector x1, NumericVector x2, 
                            NumericVector gam, NumericVector alp, NumericVector beta){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_cal1_wei(d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta); //wgt[i]*(temp);
  }
  
  return res;
}

double logL_obs_wei(int r, int d1, int d2, double t1, double t2, double x1, double x2, 
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  double res=0.0;
  if(r==1){
    res = logL_cal_wei(d1, d2, t1, t2, x1, x2, gam, alp, beta, eta, p1);
  }else{
    double res0 = logL_cal_wei(d1, d2, t1, t2, x1, 0.0, gam, alp, beta, eta, p1);
    double res1 = logL_cal_wei(d1, d2, t1, t2, x1, 1.0, gam, alp, beta, eta, p1);
    res = log(exp(res0)+exp(res1));
  }
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_obs_wei(IntegerVector r, IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, NumericVector x1, NumericVector x2, 
                            NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_obs_wei(r[i], d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta, eta, p1); //wgt[i]*(temp);
  }
  
  return res;
}

// [[Rcpp::export()]]
double logL_cal_pwc(double d1, double d2, double t1, double t2, double x1, double x2,
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1,  
                    NumericVector brks){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); 
  NumericVector levs = alp*exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2); 
  double Hft1 = Hpwc_double(t1, brks, levs, 0);
  double Sft1 = exp(-Hft1); 
  double Hft2 = Hpwc_double(t2, brks, levs, 0);
  double Sft2 = exp(-Hft2); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(d1==1){
    
    logL1=log(piZ);
  }else{
    
    if(d1>1){
      //logL1=d2*log((1-Sft2)*(1-piZ)+piZ)+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log((1-Sft2)*(1-piZ)+piZ);
      }else{
        logL1=0;//(log(1-piZ)-Hft2);
      }
    }else{//d1==0 or 1
      //logL1=d2*(log( Sft1-Sft2 )+log(1-piZ))+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log( Sft1-Sft2 )+log(1-piZ);
      }else{
        logL1=(log(1-piZ)-Hft2);
      }
    }
  }
  
  //double logL2 = x2*(eta[0]+eta[1]*x1) - log(1+exp(eta[0]+eta[1]*x1)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  //double logL2 = x1*x2*log(eta[0]) + x1*(1-x2)*log(eta[1]) + (1-x1)*x2*log(eta[2]) + (1-x1)*(1-x2)*log(1-eta[0]-eta[1]-eta[2]);
  double logL2 = x1*(x2*(eta[0]+eta[1]) - log(1+exp(eta[0]+eta[1])) + log(p1) ) + (1-x1)*( x2*eta[0] - log(1+exp(eta[0])) + log(1-p1) );
  
  double res=logL1+logL2;
  
  return res;
}

// [[Rcpp::export()]]
double logL_cal1_pwc(int d1, int d2, double t1, double t2, double x1, double x2,
                     NumericVector gam, NumericVector alp, NumericVector beta, 
                     NumericVector brks){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); //exp(eta[0]);//
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  //double Hf = pow(alp[0]*t, alp[1])*exp(beta[0]*x1+beta[1]*x2);
  NumericVector levs = alp*exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2); 
  //double hf = hpwc_double(t, brks, levs, 0); 
  double Hft1 = Hpwc_double(t1, brks, levs, 0);
  double Sft1 = exp(-Hft1); 
  double Hft2 = Hpwc_double(t2, brks, levs, 0);
  double Sft2 = exp(-Hft2); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(d1==1){
    
    logL1=log(piZ);
  }else{
    
    if(d1>1){
      //logL1=d2*log((1-Sft2)*(1-piZ)+piZ)+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log((1-Sft2)*(1-piZ)+piZ);
      }else{
        logL1=0;//(log(1-piZ)-Hft2);
      }
    }else{//d1==0 or 1
      //logL1=d2*(log( Sft1-Sft2 )+log(1-piZ))+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log( Sft1-Sft2 )+log(1-piZ);
      }else{
        logL1=(log(1-piZ)-Hft2);
      }
    }
  }
  
  //double logL2 = x2*(eta[0]+eta[1]*x1) - log(1+exp(eta[0]+eta[1]*x1)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  
  double res=logL1;//+logL2;
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_cal1_pwc(IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, NumericVector x1, NumericVector x2, 
                            NumericVector gam, NumericVector alp, NumericVector beta, 
                            NumericVector brks){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_cal1_pwc(d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta, brks); //wgt[i]*(temp);
  }
  
  return res;
}


// [[Rcpp::export()]]
NumericVector logLn_cal_pwc(IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, NumericVector x1, NumericVector x2, 
                             NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1,  
                             NumericVector brks){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_cal_pwc(d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta, eta, p1, brks); //wgt[i]*(temp);
  }
  
  return res;
}


double logL_obs_pwc(int r, int d1, int d2, double t1, double t2, double x1, double x2, 
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1, 
                    NumericVector brks){
  double res=0.0;
  if(r==1){
    res = logL_cal_pwc(d1, d2, t1, t2, x1, x2, gam, alp, beta, eta, p1, brks);
  }else{
    double res0 = logL_cal_pwc(d1, d2, t1, t2, x1, 0.0, gam, alp, beta, eta, p1, brks);
    double res1 = logL_cal_pwc(d1, d2, t1, t2, x1, 1.0, gam, alp, beta, eta, p1, brks);
    res = log(exp(res0)+exp(res1));
  }
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_obs_pwc(IntegerVector r, IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, 
                            NumericVector x1, NumericVector x2, 
                            NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1, 
                            NumericVector brks){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_obs_pwc(r[i], d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta, eta, p1, brks); //wgt[i]*(temp);
  }
  
  return res;
}


// [[Rcpp::export()]]
NumericVector vHpwc_H1(NumericVector t, NumericVector x1, NumericVector x2, 
                       NumericVector alp, NumericVector beta, NumericVector brks){
  
  NumericVector res=no_init(t.size());
  for(int i=0; i<t.size(); i++){
    
    //double piX = expit.f(eta[0]+eta[1]*x2[i]);
    NumericVector levs = alp*exp(beta[0]*x1[i]+beta[1]*x2[i]);
    res[i] = Hpwc_double(t[i], brks, levs, 0); 
  }
  
  return res;
  
}

// [[Rcpp::export()]]
NumericVector vhpwc_H1(NumericVector t, NumericVector x1, NumericVector x2, 
                       NumericVector alp, NumericVector beta, NumericVector brks){
  
  NumericVector res=no_init(t.size());
  for(int i=0; i<t.size(); i++){
    
    //double piX = expit.f(eta[0]+eta[1]*x2[i]);
    NumericVector levs = alp*exp(beta[0]*x1[i]+beta[1]*x2[i]);
    res[i] = hpwc_double(t[i], brks, levs, 0); 
  }
  
  return res;
  
}


// [[Rcpp::export()]]
NumericVector vppwc_H1(NumericVector t, NumericVector x1, NumericVector x2, 
                       NumericVector alp, NumericVector beta, NumericVector brks){
  
  NumericVector res=no_init(t.size());
  for(int i=0; i<t.size(); i++){
    
    //double piX = expit.f(eta[0]+eta[1]*x2[i]);
    NumericVector levs = alp*exp(beta[0]*x1[i]+beta[1]*x2[i]);
    res[i] = ppwc_double(t[i], brks, levs, 1,0); 
  }
  
  return res;
  
}

