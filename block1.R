library(tidyverse)
library(lme4)
library(bartCause)
library(stan4bart)

source('load_ihdp.R')
source('get_balance.R')
source('models.R')
source('observational_random_intercepts.R')

iteration <- Sys.getenv('SLURM_ARRAY_TASK_ID')

taus <- c(0, .1, .2, .5, .8)

results <- map(taus, function(i){
  observational_random_intercepts(tau = taus[i], seed = iteration)
})

out_file <- paste0('iteration_', iteration, '.csv')
write_rds(results, out_file)
