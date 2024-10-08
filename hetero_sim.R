if(!dir.exists('simulation_results')) dir.create('simulation_results')

n <- as.numeric(Sys.getenv("E_N"))
Q <- as.numeric(Sys.getenv("E_Q"))
B <- as.numeric(Sys.getenv("E_B"))
C <- as.numeric(Sys.getenv("E_C"))
trial_number <- as.numeric(Sys.getenv("E_TRIAL"))


# Remove the decimal point in Q
Q_formatted <- gsub("\\.", "", as.character(Q))

result_filename <- sprintf("results_n%s_Q%s_B%s_C%s_trial%s.rds", n, Q_formatted, B, C, trial_number)
if(file.exists(file.path("simulation_results", result_filename))) quit()

print(c(n,Q,B,C,trial_number))

library(foreach)
library(doParallel)
library(dplyr)
library(SuperLearner)
library(hal9001)
#library(tidyverse)
library(caret)
library(data.table)
library(ggplot2)
library(earth)

source('functions.R')

#library <- c('SL.mean', 'SL.glm', 'SL.caretMP', 'SL.caretRF')
library <- c('SL.mean', 'SL.glm', 'SL.gam',
	     'SL.xgboost',#'SL.bartMachine', 
	     'SL.randomForest', 'SL.earth')


fw = function(n) runif(n, min = 0, max = 1)

fs = function(w){
  rbinom(length(w), size = 1, prob = clamp(w, Q, 1-Q))
}

fa = function(w, s){
  rbinom(length(w), prob = plogis(s+w+s*w), size = 1)
}

fm = function(w, s, a){
  rbinom(length(w), prob = (a+w+B*s)/3, size = 1)
}

fy = function(w, s, a, m){
  a + C*(1+s)*w*a + m*a + rnorm(length(w), sd = 1)
}

set.seed(1e9 * (Q) + 1e7 * (B + 1) + 1e6 * (C + 1) + n + trial_number)


data <- generate(fw, fs, fa, fm, fy, n)

out <- estima(data, library)
vardecomp <- fun_decomp(out)
bindecomp <- binary_decomp(out)

sim_dcm_result <- sim_dcm(fw, fs, fa, fm, fy)
sim_deh_result <- sim_deh(fw, fs, fa, fm, fy)
sim_dem_result <- sim_dem(fw, fs, fa, fm, fy)
sim_dmv_result <- sim_dmv(fw, fs, fa, fm, fy)

results <- list(out=out, vardecomp=vardecomp, bindecomp=bindecomp,
  sim_dcm=sim_dcm_result, sim_deh=sim_deh_result, sim_dem=sim_dem_result, sim_dmv=sim_dmv_result)


saveRDS(results, file.path("simulation_results", result_filename))

