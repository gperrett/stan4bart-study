source('observational_random_slopes.R')
source('models.R')
source('load_ihdp.R')
source('get_balance.R')

get.seed <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# g1
A <- observational_random_slopes(tau = .2, type = 'A',  seed = get.seed)
out_file <- paste0('results_random_slopes/A/g1_iteration_', get.seed, '.rds')
write_rds(A, out_file)
rm(A)
gc()

B <- observational_random_slopes(tau = .2, type = 'B',  seed = get.seed)
out_file <- paste0('results_random_slopes/B/g1_iteration_', get.seed, '.rds')
write_rds(B, out_file)
rm(B)
gc()

C <- observational_random_slopes(tau = .2, type = 'C',  seed = get.seed)
out_file <- paste0('results_random_slopes/C/g1_iteration_', get.seed, '.rds')
write_rds(C, out_file)
rm(C)
gc()


#g2


A <- observational_random_slopes(tau = .2, type = 'A',  seed = get.seed, group = 'g2')
out_file <- paste0('results_random_slopes/A/g2_iteration_', get.seed, '.rds')
write_rds(A, out_file)
rm(A)
gc()

B <- observational_random_slopes(tau = .2, type = 'B',  seed = get.seed, group = 'g2')
out_file <- paste0('results_random_slopes/B/g2_iteration_', get.seed, '.rds')
write_rds(B, out_file)
rm(B)
gc()

C <- observational_random_slopes(tau = .2, type = 'C',  seed = get.seed, group = 'g2')
out_file <- paste0('results_random_slopes/C/g2_iteration_', get.seed, '.rds')
write_rds(C, out_file)
rm(C)
gc()
