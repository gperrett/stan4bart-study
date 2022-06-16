source('functions.R')
library(tidyverse)

draw_nested_random_intercept <- function(tau, sigma_districts,sigma_schools, confounded = T, parametric_type, seed = seed){
  dat <- loadData()
  
  # set seed for response surface
  set.seed(seed)
  
  # isolate covariates
  covs.cont <- c(
    "z_ses_indonesian_base",
    "z_ses_math_base", 
    "z_ses_natscience_base", 
    "z_ses_socscience_base",
    "z_effort_base", 
    "z_belonging_base", 
    "z_selfesteem_base", 
    "z_perfgoal_base", 
    "z_perseverance_base", 
    "z_growthm_base", 
    "z_learngoal_base")
  
  covs.cat <- c("sd10a_base", "sd10b_base", "sd10c_base", "sd10d_base", 
                "sd10e_base", "sd10f_base", "sd10g_base", "sd10h_base", "girl", "quiet_study")
  
  conf.school <- c("urban", 
                   "teach_growthm_base", 
                   "teach_fail_enhance_base")
  
  
  dat <- dat %>% 
    count(district) %>% 
    rename(district_count = n) %>% 
    mutate(k = row_number()) %>% 
    select(district, k) %>% 
    right_join(dat)
  
  dat <- dat %>% 
    count(district, schoolid) %>% 
    rename(school_count = n) %>% 
    mutate(j = row_number()) %>% 
    select(district, schoolid, j) %>% 
    right_join(dat)
  
  # Set random effects for jth school and kth district 
  # number of schools and districts and index table for each
  num_districts <- length(unique(dat$district))
  num_schools <- length(unique(dat$schoolid))
  
  # apply assignment mechanism for z
  if(isTRUE(confounded)){
    #set.seed(2)
    # create indices for district and school random effects 
    asn <- dat %>%
      group_by(j) %>%
      summarise(z = mean(z)) %>%
      left_join(dat %>%
                  select(-c(district, schoolid, k,z, z_nes_school_mean_base, teach_self_efficacy_base))  %>%
                  group_by(j) %>%
                  summarise_all(.funs = list(
                    mean = mean
                    #,
                    # variance = var,
                    # min = min,
                    # max = max
                  )) %>%
                  ungroup()) 
    
    X <- asn[,3:length(asn)]
    #X <-  cbind(X, X[,"z_belonging_base_mean"] * X[,"sd10a_base_mean"])
    XX <- cbind(asn$z_perseverance_base_mean * asn$girl_mean*asn$z_perfgoal_base_mean)
    Xmat <- as.matrix(cbind(X, XX))
    
    #beta <- c(sample(c(0, .1, .25, .3), ncol(X), prob = c(.3, .2, .3, .2), replace = TRUE), .3, .2)
    beta <- c(1, #sd10a_base                
              0, #sd10b_base                
              0, #sd10c_base                
              0, #sd10d_base               
              -1, #sd10e_base               
              0, #sd10f_base                
              0, #sd10g_base               
              -.5, #sd10h_base               
              1, #girl                    
              0, #quiet_study
              -1, #z_ses_indonesian_base
              1, #z_ses_math_base
              0, #z_ses_natscience_base
              -1, #z_ses_socscience_base 
              -1, #z_effort_base
              -.5, #z_belonging_base
              .5, #z_selfesteem_base 
              -1, #z_perfgoal_base 
              .5, #z_perseverance_base
              0, # z_growthm_base 
              -1, # z_learngoal_base 
              1, # urban 
              #0, #teach_self_efficacy_base
              -1, #teach_growthm_base
              -1 ,#teach_fail_enhance_base
              1)
    
    which(colnames(X) == "z_perseverance_base_mean")
    p.score <- as.matrix(cbind(rep(1, nrow(Xmat)), Xmat)) %*% as.matrix(c(-.5, beta))
    hist(p.score)
    p.score <- exp(p.score)/(1 + exp(p.score))
    hist(p.score)
    #p.score <- pnorm(p.score)
    dat$p.score <- p.score[dat$j]
    z <- rbinom(nrow(dat), 1, dat$p.score)
    dat$z <- z
    dat <-dat %>% filter((z ==1 & p.score < .9 & p.score > .1) | z ==0)
    p <- plotBart::plot_overlap_pScores(dat, 'z',pscores = dat$p.score) + theme_classic() + theme(legend.position = 'top')
    print(p)
    dat %>% group_by(z) %>% summarise(range(p.score))

  }
  



  # random intercepts
  district_intercept <- rnorm(num_districts, 0, sigma_districts)[dat$k]
  school_intercept <- rnorm(num_schools, 0, sigma_schools)[dat$j] 
  
  
  #Create fixed effects: 
  # matrix of fixed effect predictors 
  covs <- c(covs.cat, covs.cont,conf.school)
  p=length(c(covs.cat, covs.cont,conf.school))
  X <- dat[, covs]
  N = nrow(X)
  dimx = ncol(X)
  Xmat = as.matrix(X)

  
  # linear 
  if(parametric_type == 'A'){
    beta = sample(c(0, 1, 2, 3, 4),
                  dimx + 1,
                  replace = TRUE,
                  prob = c(.5, .2, .15, .1, .05))
    
    # simulate outcome
    y0hat = cbind(rep(1, N), Xmat) %*% beta 
    y1hat = (cbind(rep(1, N), Xmat) %*% beta) + (tau * sd(y0hat + school_intercept + district_intercept))
  }
  
  if(parametric_type == 'B'){
    # individual level coefficents 
    beta <- c(sample(c(0,.1,.15, .2, .25),(dimx+1),replace=TRUE,prob=c(.6,.1,.1, .1, .1))) 
    y0hat = exp((cbind(rep(1, N), (Xmat + .5)) %*% beta))
    y1hat = (cbind(rep(1, N), (Xmat + .5)) %*% beta)
    # create offset: set outcome to standard deviation scale of y0|z = 1
    if(isTRUE(confounded)){
      offset <- (mean(y1hat[dat$z ==1] - y0hat[dat$z ==1])) -  (tau * sd(y0hat[dat$z ==1] + school_intercept[dat$z ==1] + district_intercept[dat$z ==1]))
    }else{
      offset <- (mean(y1hat- y0hat)) -  (tau * sd(y0hat + school_intercept + district_intercept))
    }
    
    y1hat = (cbind(rep(1, N), (Xmat + .5)) %*% beta) -offset
  }
  
  # non-linear 3 way interactions 
  if(parametric_type == 'C'){
    covs.cont <- c(covs.cont, conf.school)
    formula <- as.formula(paste0("y.tmp ~ (", paste0(c(covs.cont, covs.cat), collapse = " + "), ")^2 + ", paste0(paste0("I(", covs.cont, "^2)"), collapse = " + ")))
    y.tmp <- rnorm(N) 
    temp <- glm(formula, as.data.frame(Xmat), x = TRUE, family = gaussian)
    XXmat <- temp$x[,!is.na(temp$coef)]
    
    # main effects coefficients
    beta.m0 = sample(c(0,.05, .1, .15, .2),p+1,replace=T,prob=c(.4,.15, .15, .15, .15))
    beta.m1 = sample(c(0,.05, .1, .15, .2),p+1,replace=T,prob=c(.4,.15, .15, .15, .15))
    # quadratic coefficients
    #these we make pretty rare since they really represent 3-way interactions
    beta.q0 = sample(c(0,.1,.2),ncol(XXmat)-(p+1),replace=TRUE,prob=c(.8,.15,.05))
    beta.q1 = sample(c(0,.1,.2),ncol(XXmat)-(p+1),replace=TRUE,prob=c(.8,.15,.05))
    #
    beta0 = c(beta.m0,beta.q0)
    beta1 = c(beta.m1,beta.q1)
    y0hat = (XXmat %*% beta0) 
    y1hat = (XXmat %*% beta1) 
    # means <- as.vector(tapply(as.vector(y1hat - y0hat), as.factor(dat$j), mean)[dat$j])
    # offset <- means -  (tau * sd(y0hat[dat$z == 1] + school_intercept[dat$z == 1] + district_intercept[dat$z == 1]))
    if(isTRUE(confounded)){
      offset <- (mean(y1hat - y0hat)) -  (tau * sd(y0hat[dat$z == 1] + school_intercept[dat$z == 1] + district_intercept[dat$z == 1]))
    }else{
      offset <- (mean(y1hat - y0hat)) -  (tau * sd(y0hat + school_intercept + district_intercept))
    }
    y1hat = (XXmat %*% beta1)  - offset
    
  }
  
  
  # add random effects create counteractions and realize y
  y0 <- rnorm(N, y0hat, 1) + district_intercept + school_intercept
  y1 <- rnorm(N, y1hat, 1) + district_intercept + school_intercept
  
  
  y <- ifelse(dat$z ==1, y1, y0)
  
  p <- data.frame(y, z = dat$z) %>%
    ggplot(aes(y, fill = as.factor(z))) +
    geom_density(alpha = .5) +
    theme_classic()
  print(p)
  
  
  # calculate truth for results 
  if(isTRUE(confounded)){
    average.truth <- mean(y1[dat$z ==1] - y0[dat$z ==1])
    
    icate.truth <- y1[dat$z ==1] - y0[dat$z ==1]
    
    district.cate.truth <-  tibble(district = dat$district, y1, y0, z = dat$z) %>% 
      filter(z ==1) %>% 
      group_by(district) %>% 
      summarise(district.cate.truth = mean(y1 - y0)) %>% 
      select(district, district.cate.truth)
    
    school.cate.truth <- tibble(schoolid = dat$schoolid, y1, y0, z = dat$z) %>% 
      filter(z ==1) %>% 
      group_by(schoolid) %>% 
      summarise(school.cate.truth = mean(y1 - y0)) %>% 
      select(schoolid, school.cate.truth) 
  }else{
    average.truth <- mean(y1 - y0)
    
    icate.truth <- y1 - y0
    
    district.cate.truth <-  tibble(district = dat$district, y1, y0, z = dat$z) %>% 
      group_by(district) %>% 
      summarise(district.cate.truth = mean(y1 - y0)) %>% 
      select(district, district.cate.truth)
    
    school.cate.truth <- tibble(schoolid = dat$schoolid, y1, y0, z = dat$z) %>% 
      group_by(schoolid) %>% 
      summarise(school.cate.truth = mean(y1 - y0)) %>% 
      select(schoolid, school.cate.truth) 
  }
  
  
  
  # create simulation df
  dat$district <- as.factor(dat$district)
  dat$schoolid <- as.factor(dat$schoolid)
  sim.dat <- cbind.data.frame(y, 
                              z = dat$z, 
                              p.score = dat$p.score,
                              dat[, covs], 
                              district = dat$district, 
                              schoolid = dat$schoolid)
  
  # calculate balance 
  if(isTRUE(confounded)) est <- 'ATT' else est <- 'ATE'
  balance <- get_balance(sim.dat[, c(covs.cat, covs.cont, conf.school)], sim.dat$z , estimand = est)
  
  
  
  # return a list with data (sim.dat) for simulation [[1]], balance statistics [[2]], and truth for results [[3]]
  data <-
    loo::nlist(sim.dat,
               balance,
               average.truth,
               icate.truth,
               district.cate.truth,
               school.cate.truth
    )
  return(data)
}
