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
  estimates <- summary(fit, target = 'sate')$estimates
  average.pred <- estimates$estimate
  interval <- c(estimates$ci.lower, estimates$ci.upper)
  samples.icate <- bartCause::extract(fit, 'icate')
  icate.pred <- apply(samples.icate, 2, mean)
  
  # get group level predictions 
  g1 <- dat$g1[dat$treat == 1]
  g2 <- dat$g2[dat$treat == 1]
  g1.pred <- unlist(lapply(split(icate.pred, g1), mean))
  g2.pred <- unlist(lapply(split(icate.pred, g2), mean))
  
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

interval.extract.bart <- function(fit, .model){
  samples.icate <- bartCause::extract(fit, 'icate')
  icate.pred <- apply(samples.icate, 2, mean)
  
  # get group level predictions 
  g1 <- dat$g1[dat$treat == 1]
  g2 <- dat$g2[dat$treat == 1]
  g1.pred <- unlist(lapply(split(icate.pred, g1), mean))
  g2.pred <- unlist(lapply(split(icate.pred, g2), mean))
  
  g1.splits <- split(as.data.frame(t(samples.icate)), g1)
  g1.splits <- lapply(names(g1.splits), function(i) apply(g1.splits[[i]], 2, mean))
  g1.intervals <- lapply(1:length(g1.splits), function(i)quantile(g1.splits[[i]], c(.025, .975)))
  g1.intervals <- bind_rows(g1.intervals)
  g1.intervals <- cbind(g1.intervals,g1.pred ,g1.truth)
  g1.intervals$u_bias <- g1.pred - g1.truth 
  g1.intervals$s_bias <- (g1.pred - g1.truth)/sd(y)
  g1.intervals$int.length <- g1.intervals$`97.5%` - g1.intervals$`2.5%`
  g1.intervals$model <- .model
  g1.intervals$group <- rownames(g1.intervals)
  g1.intervals$tau <- tau
  
  g2.splits <- split(as.data.frame(t(samples.icate)), g2)
  g2.splits <- lapply(names(g2.splits), function(i) apply(g2.splits[[i]], 2, mean))
  g2.intervals <- lapply(1:length(g2.splits), function(i)quantile(g2.splits[[i]], c(.025, .975)))
  g2.intervals <- bind_rows(g2.intervals)
  g2.intervals <- cbind(g2.intervals,g2.pred ,g2.truth)
  g2.intervals$u_bias <- g2.pred - g2.truth 
  g2.intervals$s_bias <- (g2.pred - g2.truth)/sd(y)
  g2.intervals$int.length <- g2.intervals$`97.5%` - g2.intervals$`2.5%`
  g2.intervals$model <- .model
  g2.intervals$group <- rownames(g2.intervals)
  g2.intervals$tau <- tau
  
  list(g1.intervals, g2.intervals)
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
  g1 <- dat$g1[dat$treat == 1]
  g1.pred <- unlist(lapply(split(icate.pred, g1), mean))
  
  g2 <- dat$g2[dat$treat == 1]
  g2.pred <- unlist(lapply(split(icate.pred, g2), mean))
  
  
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

interval.extract.stan <- function(fit, .model){
  samples.ppd.test <- extract(fit, type = "ppd", sample = "test")
  
  # Individual sample treatment effects
  samples.ite <- (dat$y - samples.ppd.test) * (2 * dat$treat - 1)
  samples.ite <- samples.ite[dat$treat == 1,]
  
  # Sample average treatment effect
  samples.sate <- apply(samples.ite, 2, mean)
  sate <- mean(samples.sate)
  sate.int <- quantile(samples.sate, c(.025, .975))
  
  
  # get group level predictions 
  g1 <- dat$g1[dat$treat == 1]

  g2 <- dat$g2[dat$treat == 1]

  # get group level intervals 
  g1.splits <- split(as.data.frame(samples.ite), g1)
  g1.splits <- lapply(names(g1.splits), function(i) apply(g1.splits[[i]], 2, mean))
  g1.pred <- unlist(lapply(g1.splits, mean))
  g1.intervals <- lapply(1:length(g1.splits), function(i)quantile(g1.splits[[i]], c(.025, .975)))
  g1.intervals <- bind_rows(g1.intervals)
  g1.intervals <- cbind(g1.intervals,g1.pred ,g1.truth)
  g1.intervals$u_bias <- g1.pred - g1.truth 
  g1.intervals$s_bias <- (g1.pred - g1.truth)/sd(y)
  g1.intervals$int.length <- g1.intervals$`97.5%` - g1.intervals$`2.5%`
  g1.intervals$model <- .model
  g1.intervals$group <- unique(g1)
  g1.intervals$tau <- tau
  
  
  g2.splits <- split(as.data.frame(samples.ite), g2)
  g2.splits <- lapply(names(g2.splits), function(i) apply(g2.splits[[i]], 2, mean))
  g2.pred <- unlist(lapply(g2.splits, mean))
  g2.intervals <- lapply(1:length(g2.splits), function(i)quantile(g2.splits[[i]], c(.025, .975)))
  g2.intervals <- bind_rows(g2.intervals)
  g2.intervals <- cbind(g2.intervals,g2.pred ,g2.truth)
  g2.intervals$u_bias <- g2.pred - g2.truth 
  g2.intervals$s_bias <- (g2.pred - g2.truth)/sd(y)
  g2.intervals$int.length <- g2.intervals$`97.5%` - g2.intervals$`2.5%`
  g2.intervals$model <- .model
  g2.intervals$group <- unique(g2)
  g2.intervals$tau <- tau
  
  list(g1.intervals, g2.intervals)
}







