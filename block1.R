library(tidyverse)
library(lme4)
library(bartCause)
library(stan4bart)

source('load_ihdp.R')
source('get_balance.R')
source('models.R')
source('observational_random_intercepts.R')

get.seed <- Sys.getenv('SLURM_ARRAY_TASK_ID')

taus <- c(0, .1, .2, .5, .8)

A <- lapply(1:length(taus), function(i) observational_random_intercepts(tau = taus[i], type = 'A', seed = get.seed))
B <- lapply(1:length(taus), function(i) observational_random_intercepts(tau = taus[i], type = 'B', seed = get.seed))
C <- lapply(1:length(taus), function(i) observational_random_intercepts(tau = taus[i], type = 'C', seed = get.seed))

results <- loo::nlist(A, B, C)

out_file <- paste0('iteration_', iteration, '.csv')
write_rds(results, out_file)
