source('observational_random_intercepts.R')
source('models.R')
source('load_ihdp.R')
source('get_balance.R')

get.seed <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# g1
A <- observational_random_intercepts(tau = .2, type = 'A',  seed = get.seed)
out_file <- paste0('results/A/g1_iteration_', get.seed, '.rds')
write_rds(A, out_file)
rm(A)
gc()

B <- observational_random_intercepts(tau = .2, type = 'B',  seed = get.seed)
out_file <- paste0('results/B/g1_iteration_', get.seed, '.rds')
write_rds(B, out_file)
rm(B)
gc()

C <- observational_random_intercepts(tau = .2, type = 'C',  seed = get.seed)
out_file <- paste0('results/C/g1_iteration_', get.seed, '.rds')
write_rds(C, out_file)
rm(C)
gc()


#g2


A <- observational_random_intercepts(tau = .2, type = 'A',  seed = get.seed, group = 'g2')
out_file <- paste0('results/A/g2_iteration_', get.seed, '.rds')
write_rds(A, out_file)
rm(A)
gc()

B <- observational_random_intercepts(tau = .2, type = 'B',  seed = get.seed, group = 'g2')
out_file <- paste0('results/B/g2_iteration_', get.seed, '.rds')
write_rds(B, out_file)
rm(B)
gc()

C <- observational_random_intercepts(tau = .2, type = 'C',  seed = get.seed, group = 'g2')
out_file <- paste0('results/C/g2_iteration_', get.seed, '.rds')
write_rds(C, out_file)
rm(C)
gc()

