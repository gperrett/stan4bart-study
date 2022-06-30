library(tidyverse)
library(lme4)
library(bartCause)
library(stan4bart)
source('load_ihdp.R')
source('get_balance.R')
source('models.R')

observational_random_slopes <- function(tau, type, g1.ICC = .1, g2.ICC = .2, .rho =.2, seed = NULL){
  tau <<- tau
  ihdp <- load_ihdp()
  
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
  g1 <- ihdp$g1
  g2 <- ihdp$g2
  n.g1 <- length(unique(g1))
  n.g2 <- length(unique(g2))
  
  # liner treatment version A
  if(type == 'A'){
    betaA = sample(c(0:4),dimx+1,replace=TRUE,prob=c(.5,.2,.15,.1,.05))
    y0hat = cbind(rep(1, N), Xmat) %*% betaA
    y0 = rnorm(N, y0hat, 1)
    y1 = rnorm(N, y0hat+(tau*sd(y0hat)), 1)
  }
  
  if(type == 'B'){
    betaB = c(sample(c(.0,.1,.2,.3, .4),(dimx+1),replace=TRUE,prob=c(.6,.1,.1,.1, .1))) 
    y0hat = exp((cbind(rep(1, N), (Xmat + .5)) %*% betaB))
    y1hat = cbind(rep(1, N), (Xmat + .5)) %*% betaB 
    offset = c(mean(y1hat[ihdp$treat==1] - y0hat[ihdp$treat==1])) - (tau *sd(y0hat))
    y1hat = cbind(rep(1, N), (Xmat + .5)) %*% betaB -offset
    
    y0 = rnorm(N, y0hat, 1)
    y1 = rnorm(N, y1hat, 1)
  }
  
  if(type == 'C'){
    # get model matrix
    ytmp=rnorm(N)
    mod.bal <- glm(formula=ytmp~(bw+b.head+preterm+birth.o+nnhealth+sex+twin+b.marr+mom.lths+mom.hs+mom.scoll+cig+first+booze+drugs+work.dur+prenatal)^2 + I(bw^2) + I(b.head^2) + I(preterm^2) + I(birth.o^2) + I(nnhealth^2),x=T,data=cbind.data.frame(Xmat))
    coefs <- mod.bal$coef
    XX <- mod.bal$x
    XX <- XX[,!is.na(coefs)]
    
    # create y 
    betaC.m0 = sample(c(0,1,2),p+1,replace=T,prob=c(.6,.3,.1))
    betaC.m1 = sample(c(0,1,2),p+1,replace=T,prob=c(.6,.3,.1))
    # quadratic coefficients
    betaC.q0 = sample(c(0,.5,1),ncol(XX)-(p+1),replace=TRUE,prob=c(.8,.15,.05))
    betaC.q1 = sample(c(0,.5,1),ncol(XX)-(p+1),replace=TRUE,prob=c(.8,.15,.05))
    #
    betaC0 = c(betaC.m0,betaC.q0)
    betaC1 = c(betaC.m1,betaC.q1)
    y0hat = (XX) %*% betaC0
    y1hat = (XX) %*% betaC1 
    offset = c(mean(y1hat[ihdp$treat==1] - y0hat[ihdp$treat==1])) - (tau *sd(y0hat))
    y1hat = (XX) %*% betaC1 - offset
    y0 = rnorm(N, y0hat, 1)
    y1 = rnorm(N, y1hat, 1)
  }
  
  y <- if_else(treat == 1, y1, y0)
  
  # add group effects 
  # create random effects 
  Mu <- c(0, 0)
  rho = .rho
  Rho <- matrix(c(1, rho, rho, 1), nrow = 2)
  sigmas <- c(sqrt(g1.ICC), .1)
  Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)
  g1.mat <- MASS::mvrnorm(n.g1, Mu, Sigma, empirical = TRUE)
  
  sigmas <- c(sqrt(g2.ICC), .1)
  Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)
  g2.mat <- MASS::mvrnorm(n.g1, Mu, Sigma, empirical = TRUE)

  # add random intercepts 
  y <- y + g1.mat[, 1][g1] + g2.mat[, 1][g2]
  
  # add random slopes
  y[treat ==1] <- y[treat ==1] + g1.mat[, 2][g1][treat==1] + g2.mat[, 2][g2][treat == 1]
  
  y <<- y
  dat <<- cbind(y, treat, X, g1, g2)
  
  # get causal stats
  balance <- get_balance(dat[, c(covs.cat, covs.cont)], dat$treat , estimand = 'ATT')
  
  average.truth <<- mean(y1[treat ==1] - y0[treat ==1])
  
  icate.truth <<- y1[treat ==1] - y0[treat ==1]
  
  g1.truth <<-  tibble(g1 = dat$g1, y1, y0, z = treat) %>% 
    filter(z ==1) %>% 
    group_by(g1) %>% 
    summarise(g1.truth = mean(y1 - y0)) %>% 
    select(g1.truth) %>% 
    as_vector()
  
  g2.truth <<- tibble(g2 = dat$g2, y1, y0, z = treat) %>% 
    filter(z ==1) %>% 
    group_by(g2) %>% 
    summarise(g2.truth = mean(y1 - y0)) %>% 
    select(g2.truth) %>% 
    as_vector()
  
  
  ########## fit models #########
  results <- list()
  # linear regression no pooling 
  no_pool <- lm(y ~ . -g1 -g2, data = dat)
  results[[1]] <- linear.regression(no_pool, .model = c('no pooling'))
  rm(no_pool)
  gc()
  # linear regresson full pooling
  full_pool <- lm(y ~ ., data = dat)
  results[[length(results) + 1]] <- linear.regression(full_pool, .model = 'full pooling')
  rm(full_pool)
  gc()
  
  # linear regression partial pooling i.e. lme4
  partial_pool <- lmer(y ~ . - g1 - g2 + (z|g1) + (z|g2), data = dat)
  results[[length(results) + 1]] <- extract.lme4(partial_pool)

  # vanilla bart
  bart <- bartc(y, treat, . -g1 -g2, data = dat, estimand = 'att', seed = 0)
  results[[length(results) + 1]] <- extract.bart(bart, 'vanilla bart no pooling')
  rm(bart)
  gc()
  
  # bart with fixed effects
  f.bart <- bartc(y, treat, ., data = dat, estimand = 'att', seed = 0)
  results[[length(results) + 1]] <- extract.bart(f.bart, 'bart with fixed effects')
  rm(f.bart)
  gc()
  
  # rbart for one of the two groups
  r.bart <- bartc(y, treat, . -g1, 
                  group.by = g1, 
                  data = dat,
                  use.ranef = TRUE, 
                  estimand = 'att', 
                  seed = 0)
  
  results[[length(results) + 1]] <- extract.bart(r.bart, 'rbart')
  rm(r.bart)
  gc()
  
  # stan4bart 
  
  # vanillia 
  s4b <- stan4bart(y ~ bart(. -g1 -g2) + (z|g1) + (z|g2), 
                   data = dat, 
                   treatment = treat, 
                   cores = 1, 
                   chains = 4, 
                   seed = 0)
  
  results[[length(results) + 1]] <- extract.stan4bart(s4b, .model = 'vanilla stan4bart')
  rm(s4b)
  gc()
  # with a fixed z
  s4b.z <- stan4bart(y ~ treat + bart(. -g1 -g2) + (z|g1) + (z|g2), 
                     data = dat, 
                     treatment = treat, 
                     cores = 1, 
                     chains = 4, 
                     seed = 0)
  
  results[[length(results) + 1]] <- extract.stan4bart(s4b.z, .model = 'stan4bart with linear z')
  rm(s4b.z)
  gc()
  
  # with many trees 
  s4b.trees <- stan4bart(y ~ bart(. -g1 -g2) + (z|g1) + (z|g2), 
                         data = dat, 
                         treatment = treat, 
                         cores = 1, 
                         chains = 4, 
                         bart_args = list(n.trees = 200), 
                         seed = 0)
  
  results[[length(results) + 1]] <- extract.stan4bart(s4b.trees, .model = 'stan4bart with many trees')
  rm(s4b.trees)
  gc()
  
  # refit but now with p.scores 
  dat$p.score <- bart$p.score
  
  # vanillia 
  ps.s4b <- stan4bart(y ~ bart(. -g1 -g2) + (z|g1) + (z|g2), 
                      data = dat, 
                      treatment = treat, 
                      cores = 1, 
                      chains = 4, 
                      seed = 0)
  
  results[[length(results) + 1]] <- extract.stan4bart(ps.s4b, .model = 'vanilla stan4bart with p.score')
  rm(ps.s4b)
  gc()
  
  # with a fixed z
  ps.s4b.z <- stan4bart(y ~ treat + bart(. -g1 -g2) + (z|g1) + (z|g2), 
                        data = dat, 
                        treatment = treat, 
                        cores = 1, 
                        chains = 4, 
                        seed = 0)
  
  results[[length(results) + 1]] <- extract.stan4bart(ps.s4b.z, .model = 'stan4bart with linear z with p.score')
  rm(ps.s4b.z)
  gc()
  # with many trees 
  ps.s4b.trees <- stan4bart(y ~ bart(. -g1 -g2) + (z|g1) + (z|g2), 
                            data = dat, 
                            treatment = treat, 
                            cores = 1, 
                            chains = 4, 
                            bart_args = list(n.trees = 200), 
                            seed = 0)
  
  results[[length(results) + 1]] <- extract.stan4bart(ps.s4b.trees, .model = 'stan4bart with many trees with p.score')
  rm(ps.s4b.trees)
  gc()
  
  singular <- isSingular(partial_pool)
  results <- results %>% bind_rows()
  rownames(results) <- 1:nrow(results)
  
  out <-  loo::nlist(
    results,
    singular,
    dat,
    balance,
    average.truth,
    icate.truth,
    g1.truth,
    g2.truth
  )
  
  return(out)
  cat('iteraction complete!')
}
