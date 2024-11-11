## Introduction
This repository houses the R code used for simulation studies within the project titled “Two-Phase Designs for Evaluating Tests for Cancer Screening and Management.” Inside, a folder labeled utility contains all utility functions. Additionally, there’s a running script named sim_script.R designed for executing simulations.
utility

This folder contains:
-	*design.r* and *design_resid_interval.r* contain utitlity functions for phase II designs including stratified and residual-dependent sampling schemes;
-	*commonf.cpp*, *commonf.h*, *loglik_interval.cpp* and *est_pwc_interval.R* contains functions for estimation and inference via maximum likelihood method; and
-	*opt_str_csd.r* contains functions for optimal stratified designs, which are not feasible in practice but implementable in simulation studies.
  
### sim_script.R
This script is designed to execute an example simulation study. By modifying the data configurations, it can replicate all the simulation results outlined in the manuscript. The script employs a for-loop to conduct repeated simulations. Within each simulation iteration, three key steps are executed: 1) data generation; 2) implementation of the two-phase design; and 3) estimation and inference.

### dat_interval.csv
A toy data set showing the data structure as an input. 
