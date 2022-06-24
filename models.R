# loadData <- function(){
#   dat <- readr::read_csv('data/wb_sample.csv', show_col_types = F)
#   return(dat)
# }

linear.regression <- function(fit, .model){
  average.pred <- fit$coefficients[2]
  se <- sqrt(diag(vcov(fit)))[2]
  interval <- c(average.pred - 1.96*se, average.pred + 1.96*se)
  icate.pred <- average.pred
  g1.pred <- average.pred
  g2.pred <- average.pred
  
  out <- data.frame(
    u_bias = average.pred - average.truth,
    s_bias = (average.pred - average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  average.truth &&  average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - icate.truth)^2)) / sd(y),
    pehe.g1 = sqrt(mean((g1.pred - g1.truth)^2)) / sd(y),
    pehe.g2 = sqrt(mean((g2.pred - g2.truth)^2)) / sd(y), 
    model = .model, 
    tau = tau
  )
  
  out <- data.frame(out)
  
  return(out)

}


extract.lme4 <- function(fit){
  average.pred <- fit@beta[2]
  se <- sqrt(diag(vcov(fit)))[2]
  interval <- c(average.pred - 1.96*se, average.pred + 1.96*se)
  icate.pred <- average.pred
  g1.pred <- average.pred
  g2.pred <- average.pred
  
  out <- data.frame(
    u_bias = average.pred - average.truth,
    s_bias = (average.pred - average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  average.truth &&  average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - icate.truth)^2)) / sd(y),
    pehe.g1 = sqrt(mean((g1.pred - g1.truth)^2)) / sd(y),
    pehe.g2 = sqrt(mean((g2.pred - g2.truth)^2)) / sd(y), 
    model = 'partial pooling - lme4', 
    tau = tau
  )
  
  return(out)
  
}

extract.bart <- function(fit, .model){
  y.pred <- apply(bartCause::extract(fit, "mu.obs"), 2, mean)
  estimates <- summary(fit)$estimates
  average.pred <- estimates$estimate
  interval <- c(estimates$ci.lower, estimates$ci.upper)
  samples.icate <- bartCause::extract(fit, 'icate')
  icate.pred <- apply(samples.icate, 2, mean)
  
  g1.pred <- unlist(lapply(split(samples.icate, dat$g1[dat$treat == 1]), mean))
  g1.pred <- g1.pred[!is.na(g1.pred)]
  
  g2.pred <- unlist(lapply(split(samples.icate, dat$g2[dat$treat == 1]), mean))
  g2.pred <- g2.pred[!is.na(g2.pred)]
  
  out <- data.frame(
    u_bias = average.pred - average.truth,
    s_bias = (average.pred - average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  average.truth &&  average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - icate.truth)^2)) / sd(y),
    pehe.g1 = sqrt(mean((g1.pred - g1.truth)^2)) / sd(y),
    pehe.g2 = sqrt(mean((g2.pred - g2.truth)^2)) / sd(y), 
    model = .model, 
    tau = tau
  )
  
  return(out)
}


extract.stan4bart <- function(fit, .model){
  # patt
  samples.ppd.train <- stan4bart:::extract.stan4bartFit(fit, type = "ppd", sample = "train")
  samples.ppd.test <- stan4bart:::extract.stan4bartFit(fit, type = "ppd", sample = "test")
  samples.patt <- apply((samples.ppd.train[dat$treat ==1,] - samples.ppd.test[dat$treat ==1,]) * (2 * fit$frame[[fit$treatment]][dat$treat ==1] - 1), 2, mean)
  
  # icate
  samples.mu.train <- stan4bart:::extract.stan4bartFit(fit)
  samples.mu.test  <- stan4bart:::extract.stan4bartFit(fit, sample = "test")
  samples.icate <- (samples.mu.train[dat$treat ==1,]  - samples.mu.test[dat$treat ==1,]) * (2 * fit$frame[[fit$treatment]][dat$treat ==1] - 1)
  
  # results
  y.pred <- apply(samples.ppd.train, 1, mean)
  average.pred <- mean(samples.patt)
  interval <- c(average.pred - 1.96 * sd(samples.patt), average.pred + 1.96 * sd(samples.patt))
  icate.pred <- apply(samples.icate, 1, mean)
  
  g1 <- dat$g1[dat$treat ==1]
  g2 <- dat$g2[dat$treat ==1]
  
  g1.pred <- unlist(lapply(split(samples.icate, g1), mean))
  g1.pred <- g1.pred[!is.na(g1.pred)]
  
  g2.pred <- unlist(lapply(split(samples.icate, g2), mean))
  g2.pred <- g2.pred[!is.na(g2.pred)]
  
  out <- data.frame(
    u_bias = average.pred - average.truth,
    s_bias = (average.pred - average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  average.truth &&  average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - icate.truth)^2)) / sd(y),
    pehe.g1 = sqrt(mean((g1.pred - g1.truth)^2)) / sd(y),
    pehe.g2 = sqrt(mean((g2.pred - g2.truth)^2)) / sd(y), 
    model = .model, 
    tau = tau
  )
  
  return(out)
  
}



