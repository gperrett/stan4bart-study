linear.regression <- function(fit, vary = FALSE, .model){
  
  average.pred <- fit$coefficients[1]
  se <- sqrt(diag(vcov(fit)))[1]
  interval <- c(average.pred - 1.96*se, average.pred + 1.96*se)
  
  z1 <- fit$model
  z1$treat <- 1
  
  z0 <- fit$model
  z0$treat <- 0
  
  
  icate.pred <- predict(fit, newdata = z1) - predict(fit, newdata = z0)
  
  g.pred <- tibble(icate.pred, g, treat) %>%
    filter(treat == 1) %>% 
    group_by(g) %>% 
    summarise(g.pred = mean(icate.pred)) %>% 
    select(g.pred) %>% 
    as_vector()
  
  icate.pred <- icate.pred[treat == 1]
  
  
  out <- data.frame(
    u_bias = average.pred - average.truth,
    s_bias = (average.pred - average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  average.truth &&  average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - icate.truth)^2)) / sd(y),
    pegste = sqrt(mean((g.pred - g.truth)^2)) / sd(y),
    model = .model, 
    tau = tau
  )
  
  out <- data.frame(out)
  
  return(out)
  
}

extract.bart <- function(fit, .model){
  estimates <- summary(fit, target = 'sate')$estimates
  average.pred <- estimates[nrow(estimates),1]
  interval <- c(estimates[nrow(estimates),3], estimates[nrow(estimates),4])
  samples.icate <- bartCause::extract(fit, 'icate')
  icate.pred <- apply(samples.icate, 2, mean)
  
  # get group level predictions 
  g <- dat$g[dat$treat == 1]
  g.pred <- unlist(lapply(split(icate.pred, g), mean))
  
  out <- data.frame(
    u_bias = average.pred - average.truth,
    s_bias = (average.pred - average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  average.truth &&  average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - icate.truth)^2)) / sd(y),
    pegste = sqrt(mean((g.pred - g.truth)^2)) / sd(y),
    model = .model, 
    tau = tau
  )
  
  
  return(out)
}

extract.rstanarm <- function(fit, .model){
  
  average.pred <- fit$coefficients[2]
  se <- sqrt(diag(vcov(fit)))[2]
  interval <- c(average.pred - 1.96*se, average.pred + 1.96*se)
  
  z1 <- fit$data
  z1$treat <- 1
  
  z0 <- fit$data
  z0$treat <- 0
  
  
  icate.pred <- apply(posterior_epred(fit, newdata = z1), 2, mean) - apply(posterior_epred(fit, newdata = z0), 2, mean)
  
  g.pred <- tibble(icate.pred, g, treat) %>%
    filter(treat == 1) %>% 
    group_by(g) %>% 
    summarise(g.pred = mean(icate.pred)) %>% 
    select(g.pred) %>% 
    as_vector()
  
  icate.pred <- icate.pred[treat == 1]
  
  
  out <- data.frame(
    u_bias = average.pred - average.truth,
    s_bias = (average.pred - average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  average.truth &&  average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - icate.truth)^2)) / sd(y),
    pegste = sqrt(mean((g.pred - g.truth)^2)) / sd(y),
    model = .model, 
    tau = tau
  )
  
  
  return(out)
  
}



extract.stan4bart <- function(fit, .model){
  
  samples.ppd.test <- extract(fit, type = "ppd", sample = "test")
  
  # Individual sample treatment effects
  samples.ite <- (dat$y - samples.ppd.test) * (2 * dat$treat - 1)
  samples.ite <- samples.ite[dat$treat == 1,]
  
  # Sample average treatment effect
  samples.sate <- apply(samples.ite, 2, mean)
  sate <- mean(samples.sate)
  sate.int <- quantile(samples.sate, c(.025, .975))
  
  # icate
  samples.mu.train <- stan4bart:::extract.stan4bartFit(fit)
  samples.mu.test  <- stan4bart:::extract.stan4bartFit(fit, sample = "test")
  samples.icate <- (samples.mu.train  - samples.mu.test) * (2 * fit$frame[[fit$treatment]] - 1)
  
  # results
  average.pred <- sate
  interval <- sate.int
  icate.pred <- apply(samples.icate, 1, mean)
  icate.pred <- icate.pred[dat$treat == 1]
  
  # get group level predictions 
  g <- dat$g[dat$treat == 1]
  g.pred <- unlist(lapply(split(icate.pred, g), mean))
  
  out <- data.frame(
    u_bias = average.pred - average.truth,
    s_bias = (average.pred - average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  average.truth &&  average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - icate.truth)^2)) / sd(y),
    pegste = sqrt(mean((g.pred - g.truth)^2)) / sd(y),
    model = .model, 
    tau = tau
  )
  
  return(out)
  
}

# extract.bcf <- function(fit, .model){
#   tau_post <- fit$tau
#   icate.pred <- colMeans(tau_post)[treat == 1]
#   average.pred <- mean(icate.pred)
#   interval <- f
#   bcf_fit
#   g.pred <- unlist(lapply(split(icate.pred, dat$g[treat == 1]), mean))
#   
# }

interval.linear.regression <- function(fit, .model){
  truth <- g.truth
  if(group == 'g1'){
    point <- coef(fit)[c('treat', paste0('treat:g', 16:40))] + coef(fit)[paste0('g', 15:40)]
    se <- sqrt(diag(vcov(fit)))[c('treat', paste0('treat:g', 16:40))]
  }else{
    point <- coef(fit)[c('treat', paste0('treat:g', 2:8))] + coef(fit)[paste0('g', 1:8)]
    se <- sqrt(diag(vcov(fit)))[c('treat', paste0('treat:g', 2:8))]
  }

  lci <- (point - 1.96*se)
  uci <- (point + 1.96*se)
  int.length <- (uci - lci)/sd(y)
  cover <- truth >= lci & truth <= uci
  point + 1.96*se
  u_bias <- truth - point
  s_bias <- u_bias /sd(y)
  model = .model
  tau = tau
  group <- paste0(group, '-', unique(g)[order(unique(g))])

  
  group_est <- cbind.data.frame(truth, cover, int.length, u_bias, s_bias, group, model, tau)
  
}


interval.extract.rstanarm <- function(fit, .model){
  # individual effects 
  para_name <- colnames(as.matrix(fit))
  para_name
  
  # Obtain school-level varying intercept a_j
  # draws for overall mean
  mu_a_sims <- as.matrix(fit, 
                         pars = "(Intercept)")
  u_sims <- as.matrix(fit, 
                      regex_pars = "b\\[\\(Intercept\\) g\\:")
  
  if(nrow(fit$stan_summary) == 79) u_sims <- u_sims + as.matrix(fit, regex_pars = "b\\[treat g\\:")
  
  a_sims <- as.numeric(mu_a_sims) + u_sims  
  
  # Posterior mean and SD of each alpha
  a_mean <- apply(X = a_sims,     # posterior mean
                  MARGIN = 2,
                  FUN = mean)
  a_sd <- apply(X = a_sims,       # posterior SD
                MARGIN = 2,
                FUN = sd)
  
  # Posterior median and 95% credible interval
  a_quant <- apply(X = a_sims, 
                   MARGIN = 2, 
                   FUN = quantile, 
                   probs = c(0.025, 0.975))
  
  estimates <- data.frame(t(a_quant))
  
  estimates$cover <- estimates$X2.5. <= g.truth & estimates$X97.5. >= g.truth
  estimates$int.length <- with(estimates, X97.5. - X2.5.) /sd(y)
  estimates$u_bias <- a_mean - g.truth
  estimates$s_bias <- estimates$u_bias/sd(y)
  estimates$group <- paste0(group, '-', unique(g)[order(unique(g))])
  estimates$model <- .model 
  estimates$tau = tau
  estimates$truth <- g.truth
  estimates <- estimates[, c('truth', 'cover', 'int.length', 'u_bias', 's_bias', 'group', 'model', 'tau')]
  
}

interval.extract.bart <- function(fit, .model){
  samples.icate <- bartCause::extract(fit, 'icate')
  
  group.split <- split(as.data.frame(t(samples.icate)), dat$g[dat$treat == 1])
  samples.gcatt<- lapply(group.split, function(i) apply(i, 2L, mean))
  
  n.obs <- lapply(group.split, nrow)
  est <- sapply(samples.gcatt,  mean)
  sd <- sapply(names(samples.gcatt), function(i) sqrt(var(samples.gcatt[[i]]) + mean((extract(fit, 'sigma')^2)) / n.obs[[i]]))
  ucl <- est + 1.96*sd
  lcl <- est - 1.96*sd
  cover <- lcl <= g.truth & ucl>=g.truth
  truth <- g.truth
  int.length <- (ucl - lcl) / sd(y)
  u_bias <- truth - est
  s_bias <- u_bias/sd(y)
  group <- paste0(group, '-', unique(g)[order(unique(g))])
  model <- .model 
  tau <- tau 
  
  out <- cbind.data.frame(truth, cover, int.length, u_bias, s_bias, group, model, tau)
  
  return(out)
  
}



interval.extract.stan4bart <- function(fit, .model){
  samples.ppd.test <- extract(fit, type = "ppd", sample = "test")
  
  # Individual sample treatment effects
  samples.ite <- (dat$y - samples.ppd.test) * (2 * dat$treat - 1)
  samples.ite <- samples.ite[dat$treat == 1,]
  
  # Sample average treatment effect
  samples.satt <- apply(samples.ite, 2, mean)
  satt <- mean(samples.satt)
  satt.int <- quantile(samples.satt, c(.025, .975))
  
  # get group level intervals 
  group.splits <- split(as.data.frame(samples.ite), dat$g[dat$treat == 1])
  group.splits <- lapply(names(group.splits), function(i) apply(group.splits[[i]], 2, mean))
  g.pred <- unlist(lapply(group.splits, mean))
  g.intervals <- lapply(1:length(group.splits), function(i)quantile(group.splits[[i]], c(.025, .975)))
  g.intervals <- bind_rows(g.intervals)
  g.intervals <- cbind(g.intervals,g.pred ,g.truth)
  g.intervals$truth <- g.truth
  g.intervals$cover <- with(g.intervals, truth >= `2.5%` & truth <= `97.5%`)
  g.intervals$u_bias <- g.pred - g.truth 
  g.intervals$s_bias <- (g.pred - g.truth)/sd(y)
  g.intervals$int.length <- (g.intervals$`97.5%` - g.intervals$`2.5%`)/sd(y)
  g.intervals$model <- .model
  g.intervals$group <- paste0(group, '-', unique(g)[order(unique(g))])
  g.intervals$tau <- tau
  g.intervals <- g.intervals[, c('truth', 'cover', 'int.length', 'u_bias', 's_bias', 'group', 'model', 'tau')]
  
  return(g.intervals)  
  
}

extract.bcf <- function(fit, .model){
  mu.cf <- t(fit$yhat) - (2 * treat - 1) * t(fit$tau)
  mu.obs <- t(fit$yhat)
  samples.icate <- (2 * treat - 1) * (mu.obs - mu.cf)
  y.cf <- mu.cf + rnorm(nrow(mu.cf) * ncol(mu.cf), 0, rep(fit$sigma, each = nrow(mu.cf)))
  samples.ite <- (2 * treat - 1) * (y - y.cf)
  samples.satt <- colMeans(samples.ite[treat == 1,])
  
  # calculate stats
  icate.est <- apply(samples.icate, 1L, mean)
  average.est <- mean(samples.satt)
  interval <- quantile(samples.satt, c(.025, .975))
  cover <- average.truth >= interval[1] & average.truth <= interval[2]
  
  est.gsatt <- sapply(levels(g), function(j) {
    if(!is.null(dim(samples.ite[treat == 1 & g == j,]))){
      mean(apply(samples.ite[treat == 1 & g == j,], 2L, mean))
    }else{mean(samples.ite[treat == 1 & g == j,])}
      })
  
  out <- data.frame(
    u_bias = average.est - average.truth,
    s_bias = (average.est - average.truth) /sd(y),
    ci_len =(interval[2] - interval[1])/sd(y), 
    cover = as.double(interval[1] <=  average.truth &&  average.truth <= interval[2]),
    pehe = sqrt(mean((icate.est[treat == 1] - icate.truth)^2)) / sd(y),
    pegste = sqrt(mean((est.gsatt - g.truth)^2)) / sd(y),
    model = .model, 
    tau = tau
  )
  
  return(out)
  
}

interval.extract.bcf <- function(fit, .model){
  mu.cf <- t(fit$yhat) - (2 * treat - 1) * t(fit$tau)
  y.cf <- mu.cf + rnorm(nrow(mu.cf) * ncol(mu.cf), 0, rep(fit$sigma, each = nrow(mu.cf)))
  samples.ite <- (2 * treat - 1) * (y - y.cf)
  samples.satt <- colMeans(samples.ite[treat == 1,])
  gsatt <- sapply(levels(g), function(j) {
    if(!is.null(dim(samples.ite[treat == 1 & g == j,]))){
      apply(samples.ite[treat == 1 & g == j,], 2L, mean)
    }else{samples.ite[treat == 1 & g == j,]}
  })
  
  means <- apply(gsatt, 2, mean)
  estimates <- as.data.frame(t(apply(gsatt, 2, function(i) quantile(i, c(.025, .975)))))
  estimates$cover <- g.truth >= estimates[, 1] & g.truth <= estimates[, 2]
  estimates$int.length <- (estimates[, 2] - estimates[, 1]) /sd(y)
  estimates$truth <- g.truth
  estimates$u_bias <- g.truth - means
  estimates$s_bias <- estimates$u_bias/sd(y)
  estimates$model = .model
  estimates$tau = tau
  estimates$group <- paste0(group, '-', unique(g)[order(unique(g))])
  estimates <- estimates[, c('truth', 'cover', 'int.length', 'u_bias', 's_bias', 'group', 'model', 'tau')]
  
  return(estimates)

}
