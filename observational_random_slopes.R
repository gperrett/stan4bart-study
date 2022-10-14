library(tidyverse)
library(bartCause)
library(stan4bart)
library(rstanarm)
library(bcf)
source('load_ihdp.R')
source('models.R')

observational_random_slopes <- function(tau, type, seed = NULL, group = 'g1', .rho = .2){
  
  tau <<- tau
  ihdp <- load_ihdp()
  group <<- group
  
  treat <<- ihdp$treat
  
  set.seed(seed)
  
  `%not_in%` <- Negate(`%in%`)
  covs.cont <- c("bw", "b.head", "preterm", "birth.o", "nnhealth")
  covs.cat  <- c("sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                 "cig", "first", "booze", "drugs", "work.dur", "prenatal")
  
  p=length(c(covs.cont,covs.cat))
  X <- ihdp[,c(covs.cont, covs.cat)]
  treat <- ihdp$treat
  
  # scale continuous variables 
  X[, covs.cont] <- scale(X[, covs.cont])
  
  
  N = nrow(X)
  dimx = ncol(X)
  Xmat = as.matrix(X)
  
  g <<- switch (group,
                g1 = ihdp$g1, 
                g2 = ihdp$g2
  )
  
  n.g <- length(unique(g))
  
  # liner treatment version A
  if(type == 'A'){
    betaA = sample(c(0:4),dimx+1,replace=TRUE,prob=c(.5,.2,.15,.1,.05))
    y0hat = cbind(rep(1, N), Xmat) %*% betaA
    y1hat = y0hat+(tau*sd(y0hat)) 

  }
  
  if(type == 'B'){
    betaB = c(sample(c(.0,.1,.2,.3, .4),(dimx+1),replace=TRUE,prob=c(.6,.1,.1,.1, .1))) 
    y0hat = exp((cbind(rep(1, N), (Xmat + .5)) %*% betaB))
    y1hat = cbind(rep(1, N), (Xmat + .5)) %*% betaB 
    offset = c(mean(y1hat[ihdp$treat==1] - y0hat[ihdp$treat==1])) - (tau *sd(y0hat))
    y1hat = cbind(rep(1, N), (Xmat + .5)) %*% betaB -offset
    
  }
  
  if(type == 'C'){
    # get model matrix
    ytmp=rnorm(N)
    mod.bal <- glm(formula=ytmp~(bw+b.head+preterm+birth.o+nnhealth+sex+twin+b.marr+mom.lths+mom.hs+mom.scoll+cig+first+booze+drugs+work.dur+prenatal)^2 + I(bw^2) + I(b.head^2) 
                   + I(preterm^2) + I(birth.o^2) + I(nnhealth^2),x=T,data=cbind.data.frame(Xmat))
    coefs <- mod.bal$coef
    XX <- mod.bal$x
    XX <- XX[,!is.na(coefs)]
    
    # create y 
    betaC.m0 = sample(c(0,1,2),p+1,replace=T,prob=c(.5,.4,.1))
    betaC.m1 = sample(c(0,1,2),p+1,replace=T,prob=c(.5,.4,.1))
    # quadratic coefficients
    betaC.q0 = sample(c(0,.5,1),ncol(XX)-(p+1),replace=TRUE,prob=c(.7,.25,.05))
    betaC.q1 = sample(c(0,.5,1),ncol(XX)-(p+1),replace=TRUE,prob=c(.7,.25,.05))
    #
    betaC0 = c(betaC.m0,betaC.q0)
    betaC1 = c(betaC.m1,betaC.q1)
    y0hat = (XX) %*% betaC0
    y1hat = (XX) %*% betaC1 
    offset = c(mean(y1hat[ihdp$treat==1] - y0hat[ihdp$treat==1])) - (tau *sd(y0hat))
    y1hat = (XX) %*% betaC1 - offset

  }
  
  
  Mu <- c(0, 0)
  rho = .rho
  Rho <- matrix(c(1, rho, rho, 1), nrow = 2)
  sigmas <- c(1, .2)
  Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)
  random_effects <- MASS::mvrnorm(n.g, Mu, Sigma)
  
  y0hat <- y0hat + random_effects[, 1][g]
  y1hat <- y1hat + random_effects[, 1][g] + random_effects[, 2][g]
  y1 <- rnorm(N, y1hat, 1)
  y0 <- rnorm(N, y0hat, 1)
  
  y <- if_else(treat == 1, y1, y0)
  
  y <<- y
  dat <<- cbind(y, treat, X, g)
  
  # get causal stats
  average.truth <<- mean(y1[treat ==1] - y0[treat ==1])
  
  icate.truth <<- y1hat[treat == 1] - y0hat[treat == 1]
  
  g.truth <<-  tibble(g = dat$g, y1, y0, z = treat) %>% 
    filter(z ==1) %>% 
    group_by(g) %>% 
    summarise(g.truth = mean(y1 - y0)) %>% 
    select(g.truth) %>% 
    as_vector()
  
  
  ########## fit models #########
  results <- list()
  g.results <- list()
  
  # linear regression with no group info
  lin.reg <- lm(y ~ 0 + . -g, data = dat)
  results[[1]] <- linear.regression(lin.reg, .model = c('linear regression'))
  rm(lin.reg)
  gc()
  
  
  # linear regresson with fixed effects
  lin.reg.fix <- lm(y ~ 0 + ., data = dat)
  results[[length(results) + 1]] <- linear.regression(lin.reg.fix, .model = 'linear regression + fixed effects')
  rm(lin.reg.fix )
  gc()
  
  group.lin.reg <- lm(y ~ 0 + . + treat:g, data = dat)
  results[[length(results)]][5:8] <- linear.regression(group.lin.reg, .model = 'linear regression + fixed effects')[5:8]
  
  g.results[[length(g.results) + 1]] <- interval.linear.regression(group.lin.reg, .model = 'linear regression + fixed effects')
  
  rm(group.lin.reg)
  gc()
  
  # linear regression partial pooling
  partial_pool_random_intercepts <- stan_lmer(y ~ treat+bw+b.head+preterm+birth.o+nnhealth+sex+twin+b.marr+mom.lths+mom.hs+mom.scoll+cig+first+booze+drugs+work.dur+prenatal + (1|g), 
                                              data = dat, cores = 1, chains = 4, seed = 0)
  
  results[[length(results) + 1]] <- extract.rstanarm(partial_pool_random_intercepts, 'partial pooling random intercepts')
  g.results[[length(g.results) + 1]] <- interval.extract.rstanarm(partial_pool_random_intercepts, .model = 'partial pooling random intercepts')
  rm(partial_pool_random_intercepts)
  gc()
  
  # partial pool with random slopes
  partial_pool_random_slopes <- stan_lmer(y ~ treat+bw+b.head+preterm+birth.o+nnhealth+sex+twin+b.marr+mom.lths+mom.hs+mom.scoll+cig+first+booze+drugs+work.dur+prenatal+(treat|g), 
                                          data = dat, cores = 1, chains = 4, seed = 0)
  results[[length(results) + 1]] <- extract.rstanarm(partial_pool_random_slopes, 'partial pooling random slopes')
  g.results[[length(g.results) + 1]] <- interval.extract.rstanarm(partial_pool_random_slopes, .model = 'partial pooling random slopes')
  rm(partial_pool_random_slopes)
  gc()
  
  # vanilla bart
  bart <- bartc(y, treat, . -g, data = dat, 
                estimand = 'att', 
                seed = 0,
                n.threads = 1)
  
  results[[length(results) + 1]] <- extract.bart(bart, 'vanilla bart')
  g.results[[length(g.results) + 1]] <- interval.extract.bart(bart, 'vanilla bart')
  
  # bart with fixed effects
  f.bart <- bartc(y, treat, ., data = dat, estimand = 'att', seed = 0, n.threads = 1)
  results[[length(results) + 1]] <- extract.bart(f.bart, 'bart with fixed effects')
  g.results[[length(g.results) + 1]] <- interval.extract.bart(f.bart, 'bart with fixed effects')
  rm(f.bart)
  gc()
  
  # rbart for one of the two groups
  r.bart <- bartc(y, treat, . -g,
                  group.by = g,
                  data = dat,
                  group.effects = TRUE, 
                  use.ranef = TRUE,
                  estimand = 'att',
                  seed = 0,
                  n.threads = 1)
  
  results[[length(results) + 1]] <- extract.bart(r.bart, 'rbart')
  g.results[[length(g.results) + 1]] <- interval.extract.bart(r.bart, 'rbart')
  rm(r.bart)
  gc()
  
  
  # refit but now with p.scores
  dat$p.score <- bart$p.score
  
  # vanillia
  s4b_random_intercepts <- stan4bart(y ~ bart(. -g) + (1|g),
                                     data = dat,
                                     treatment = treat,
                                     cores = 1,
                                     chains = 10,
                                     iter = 4000, 
                                     seed = 0)
  
  results[[length(results) + 1]] <- extract.stan4bart(s4b_random_intercepts, .model = 'stan4bart random intercepts')
  g.results[[length(g.results) + 1]] <- interval.extract.stan4bart(s4b_random_intercepts, 'stan4bart random intercepts')
  rm(s4b_random_intercepts)
  gc()
  
  s4b_random_slopes <- stan4bart(y ~ bart(. -g) + (treat|g),
                                 data = dat,
                                 treatment = treat,
                                 cores = 1,
                                 chains = 10,
                                 iter = 4000, 
                                 seed = 0)
  
  
  results[[length(results) + 1]] <- extract.stan4bart(s4b_random_slopes, .model = 'stan4bart random slopes')
  g.results[[length(g.results) + 1]] <- interval.extract.stan4bart(s4b_random_slopes, 'stan4bart random slopes')
  rm(s4b_random_slopes)
  gc()
  
  
  # vanilla bcf 
  bcf_fit <- bcf(y , treat, as.matrix(X), as.matrix(X), bart$p.score, 2000, 2000)
  results[[length(results) + 1]] <- extract.bcf(bcf_fit, .model = 'vanilla bcf')
  g.results[[length(g.results) + 1]] <- interval.extract.bcf(bcf_fit, .model = 'vanilla bcf')
  
  rm(bcf_fit)
  gc()
  
  groups <- matrix(nrow = nrow(dat), ncol = length(unique(g)))
  for (i in 1:length(unique(g))) {
    groups[, i] <- ifelse(g == unique(g)[order(unique(g))][i], 1, 0)
  }
  
  names(groups) <- paste0('group_', unique(g)[order(unique(g))])
  X_mat <- as.matrix(cbind(X, groups))
  bcf_fit <- bcf(y , treat, X_mat, X_mat, bart$p.score, 2000, 2000)
  results[[length(results) + 1]] <- extract.bcf(bcf_fit, .model = 'bcf with groups')
  g.results[[length(g.results) + 1]] <- interval.extract.bcf(bcf_fit, .model = 'bcf with groups')
  rm(bcf_fit)
  gc()
  
  results <- bind_rows(results)
  rownames(results) <- 1:nrow(results)
  group.results <- bind_rows(g.results)
  
  out <-  loo::nlist(
    results,
    group.results 
  )
  
  return(out)
}
