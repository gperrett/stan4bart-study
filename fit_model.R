source('functions.R')
source('dgp.R')

library(stan4bart)

# function from dgp.R
draw <- draw_nested_random_intercept(tau = .1, sigma_districts = .2, sigma_schools = .2, parametric_type = 'B', seed = 19)
draw$balance # check balance
df <- draw$sim.dat # get df

# fit
fit_stan4bart <- stan4bart(y ~ bart(. - schoolid - district) + 
                             (1|district/schoolid), 
                           data = df, treatment = z, 
                           cores = 1, 
                           chains = 4, 
                           seed = 0, 
                           verbose = T)
