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

A <- lapply(1:5, function(i){observational_random_intercepts(tau = taus[i], type = 'A', seed = get.seed)})

out_file <- paste0('results/A/iteration_', get.seed, '.csv')
write_rds(A, out_file)
