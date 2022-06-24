loadData <- function(){
  dat <- readr::read_csv('data/wb_sample.csv', show_col_types = F)
  return(dat)
}

get_no_pool <- function(fit, .tau = tau){
  average.truth <- draw$average.truth
  icate.truth <- draw$icate.truth
  district.cate.truth <- draw$district.cate.truth$district.cate.truth
  school.cate.truth <- draw$school.cate.truth$school.cate.truth
  y <- draw$sim.dat$y
  
  average.pred <- fit$coefficients[2]
  se <- sqrt(diag(vcov(fit)))[2]
  interval <- c(average.pred - 1.96*se, average.pred + 1.96*se)
  icate.pred <- average.pred
  district.cate.pred <- average.pred
  school.cate.pred <- average.pred
  
  
  out <- data.frame(
    u_bias = average.pred - draw$average.truth,
    s_bias = (average.pred - draw$average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  draw$average.truth &&  draw$average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - draw$icate.truth)^2)) / sd(y),
    pehe.district = sqrt(mean((district.cate.pred - draw$district.cate.truth$district.cate.truth)^2)) / sd(y),
    pehe.school = sqrt(mean((school.cate.pred - draw$school.cate.truth$school.cate.truth)^2)) / sd(y), 
    model = 'no pooling no p.score', 
    tau = .tau
  )
  
  out <- list(out)

}


get_rbart <- function(fit, .tau = tau){
  average.truth <- draw$average.truth
  icate.truth <- draw$icate.truth
  district.cate.truth <- draw$district.cate.truth$district.cate.truth
  school.cate.truth <- draw$school.cate.truth$school.cate.truth
  y <- draw$sim.dat$y
  
  
  y.pred <- apply(bartCause::extract(fit, "mu.obs"), 2, mean)
  estimates <- summary(fit)$estimates
  average.pred <- estimates$estimate
  interval <- c(estimates$ci.lower, estimates$ci.upper)
  samples.icate <- bartCause::extract(fit, 'icate')
  icate.pred <- apply(samples.icate, 2, mean)
  
  
  district <- switch (fit$estimand,
    'ate' = draw$sim.dat$distric, 
    'att' = draw$sim.dat$district[draw$sim.dat$z ==1]
  )
  district.cate.pred <- unlist(lapply(split(samples.icate, district), mean))
  district.cate.pred <- district.cate.pred[!is.na(district.cate.pred)]
  
  school <- switch (fit$estimand,
                      'ate' = draw$sim.dat$schoolid, 
                      'att' = draw$sim.dat$schoolid[draw$sim.dat$z ==1]
  )
  school.cate.pred <- unlist(lapply(split(samples.icate, school), mean))
  school.cate.pred <- school.cate.pred[!is.na(school.cate.pred)]
  
  #if(fit$)
  out <- data.frame(
    u_bias = average.pred - draw$average.truth,
    s_bias = (average.pred - draw$average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  draw$average.truth &&  draw$average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - draw$icate.truth)^2)) / sd(y),
    pehe.district = sqrt(mean((district.cate.pred - draw$district.cate.truth$district.cate.truth)^2)) / sd(y),
    pehe.school = sqrt(mean((school.cate.pred - draw$school.cate.truth$school.cate.truth)^2)) / sd(y), 
    model = 'rbart no p.score', 
    tau = .tau
  )
  
  out <- list(out)
  
  return(out)
  
}


get_fixed_bart <- function(fit, .tau = tau){
  average.truth <- draw$average.truth
  icate.truth <- draw$icate.truth
  district.cate.truth <- draw$district.cate.truth$district.cate.truth
  school.cate.truth <- draw$school.cate.truth$school.cate.truth
  y <- draw$sim.dat$y
  
  
  y.pred <- apply(bartCause::extract(fit, "mu.obs"), 2, mean)
  estimates <- summary(fit)$estimates
  average.pred <- estimates$estimate
  interval <- c(estimates$ci.lower, estimates$ci.upper)
  samples.icate <- bartCause::extract(fit, 'icate')
  icate.pred <- apply(samples.icate, 2, mean)
  
  
  district <- switch (fit$estimand,
                      'ate' = draw$sim.dat$distric, 
                      'att' = draw$sim.dat$district[draw$sim.dat$z ==1]
  )
  district.cate.pred <- unlist(lapply(split(samples.icate, district), mean))
  district.cate.pred <- district.cate.pred[!is.na(district.cate.pred)]
  
  school <- switch (fit$estimand,
                    'ate' = draw$sim.dat$schoolid, 
                    'att' = draw$sim.dat$schoolid[draw$sim.dat$z ==1]
  )
  school.cate.pred <- unlist(lapply(split(samples.icate, school), mean))
  school.cate.pred <- school.cate.pred[!is.na(school.cate.pred)]
  
  #if(fit$)
  out <- data.frame(
    u_bias = average.pred - draw$average.truth,
    s_bias = (average.pred - draw$average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  draw$average.truth &&  draw$average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - draw$icate.truth)^2)) / sd(y),
    pehe.district = sqrt(mean((district.cate.pred - draw$district.cate.truth$district.cate.truth)^2)) / sd(y),
    pehe.school = sqrt(mean((school.cate.pred - draw$school.cate.truth$school.cate.truth)^2)) / sd(y), 
    model = 'fixed effect bart no p.score', 
    tau = .tau
  )
  
  out <- list(out)
  
  return(out)
  
}



get_vanilla_bart <- function(fit, .tau = tau){
  average.truth <- draw$average.truth
  icate.truth <- draw$icate.truth
  district.cate.truth <- draw$district.cate.truth$district.cate.truth
  school.cate.truth <- draw$school.cate.truth$school.cate.truth
  y <- draw$sim.dat$y
  
  
  y.pred <- apply(bartCause::extract(fit, "mu.obs"), 2, mean)
  estimates <- summary(fit)$estimates
  average.pred <- estimates$estimate
  interval <- c(estimates$ci.lower, estimates$ci.upper)
  samples.icate <- bartCause::extract(fit, 'icate')
  icate.pred <- apply(samples.icate, 2, mean)
  
  
  district <- switch (fit$estimand,
                      'ate' = draw$sim.dat$distric, 
                      'att' = draw$sim.dat$district[draw$sim.dat$z ==1]
  )
  district.cate.pred <- unlist(lapply(split(samples.icate, district), mean))
  district.cate.pred <- district.cate.pred[!is.na(district.cate.pred)]
  
  school <- switch (fit$estimand,
                    'ate' = draw$sim.dat$schoolid, 
                    'att' = draw$sim.dat$schoolid[draw$sim.dat$z ==1]
  )
  school.cate.pred <- unlist(lapply(split(samples.icate, school), mean))
  school.cate.pred <- school.cate.pred[!is.na(school.cate.pred)]
  
  #if(fit$)
  out <- data.frame(
    u_bias = average.pred - draw$average.truth,
    s_bias = (average.pred - draw$average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  draw$average.truth &&  draw$average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - draw$icate.truth)^2)) / sd(y),
    pehe.district = sqrt(mean((district.cate.pred - draw$district.cate.truth$district.cate.truth)^2)) / sd(y),
    pehe.school = sqrt(mean((school.cate.pred - draw$school.cate.truth$school.cate.truth)^2)) / sd(y), 
    model = 'vanilla bart no p.score', 
    tau = .tau
  )
  
  out <- list(out)
  
  return(out)
  
}


get_full_pool <- function(fit, .tau = tau){
  average.truth <- draw$average.truth
  icate.truth <- draw$icate.truth
  district.cate.truth <- draw$district.cate.truth$district.cate.truth
  school.cate.truth <- draw$school.cate.truth$school.cate.truth
  y <- draw$sim.dat$y
  
  average.pred <- fit$coefficients['z']
  se <- sqrt(diag(vcov(fit)))['z']
  interval <- c(average.pred - 1.96*se, average.pred + 1.96*se)
  icate.pred <- average.pred
  district.cate.pred <- average.pred
  school.cate.pred <- average.pred
  
  
  out <- data.frame(
    u_bias = average.pred - draw$average.truth,
    s_bias = (average.pred - draw$average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  draw$average.truth &&  draw$average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - draw$icate.truth)^2)) / sd(y),
    pehe.district = sqrt(mean((district.cate.pred - draw$district.cate.truth$district.cate.truth)^2)) / sd(y),
    pehe.school = sqrt(mean((school.cate.pred - draw$school.cate.truth$school.cate.truth)^2)) / sd(y), 
    model = 'full pooling no p.score', 
    tau = .tau
  )
  
  out <- list(out)
}

get_lme4 <- function(fit, .tau = tau){
  .singular <- isSingular(fit)
  y <- draw$sim.dat$y
  average.pred <- fit@beta[2]
  se <- sqrt(diag(vcov(fit)))['z']
  interval <- c(average.pred - 1.96*se, average.pred + 1.96*se)
  icate.pred <- average.pred
  district.cate.pred <- average.pred
  school.cate.pred <- average.pred
  
  out <- data.frame(
    u_bias = average.pred - draw$average.truth,
    s_bias = (average.pred - draw$average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  draw$average.truth &&  draw$average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - draw$icate.truth)^2)) / sd(y),
    pehe.district = sqrt(mean((district.cate.pred - draw$district.cate.truth$district.cate.truth)^2)) / sd(y),
    pehe.school = sqrt(mean((school.cate.pred - draw$school.cate.truth$school.cate.truth)^2)) / sd(y), 
    model = 'lme4 no p.score', 
    tau = .tau
  )
  
  out <- list(out, .singular)
  
  return(out)
  
}

get_stan4bart <- function(fit, .tau, type, estimand){
  # patt
  samples.ppd.train <- stan4bart:::extract.stan4bartFit(fit, type = "ppd", sample = "train")
  samples.ppd.test <- stan4bart:::extract.stan4bartFit(fit, type = "ppd", sample = "test")
  samples.patt <- apply((samples.ppd.train[sim.dat$z ==1,] - samples.ppd.test[sim.dat$z ==1,]) * (2 * fit$frame[[fit$treatment]][sim.dat$z ==1] - 1), 2, mean)
  # icate
  samples.mu.train <- stan4bart:::extract.stan4bartFit(fit)
  samples.mu.test  <- stan4bart:::extract.stan4bartFit(fit, sample = "test")
  samples.icate <- (samples.mu.train[sim.dat$z ==1,]  - samples.mu.test[sim.dat$z ==1,]) * (2 * fit$frame[[fit$treatment]][sim.dat$z ==1] - 1)
  
  # results
  y.pred <- apply(samples.ppd.train, 1, mean)
  average.pred <- mean(samples.patt)
  interval <- c(average.pred - 1.96 * sd(samples.patt), average.pred + 1.96 * sd(samples.patt))
  icate.pred <- apply(samples.icate, 1, mean)
  
  district <- switch (fit$estimand,
                      'ate' = draw$sim.dat$distric, 
                      'att' = draw$sim.dat$district[draw$sim.dat$z ==1]
  )
  
  school <- switch (estimand,
                    'ate' = draw$sim.dat$schoolid, 
                    'att' = draw$sim.dat$schoolid[draw$sim.dat$z ==1]
  )
  
  district.cate.pred <- unlist(lapply(split(samples.icate, district), mean))
  district.cate.pred <- district.cate.pred[!is.na(district.cate.pred)]
  
  school.cate.pred <- unlist(lapply(split(samples.icate, school), mean))
  school.cate.pred <- school.cate.pred[!is.na(school.cate.pred)]
  
  out <- data.frame(
    u_bias = average.pred - draw$average.truth,
    s_bias = (average.pred - draw$average.truth) /sd(y),
    ci_len =(interval[2] - interval[1]) / sd(y), 
    cover = as.double(interval[1] <=  draw$average.truth &&  draw$average.truth <= interval[2]),
    pehe = sqrt(mean((icate.pred - draw$icate.truth)^2)) / sd(y),
    pehe.district = sqrt(mean((district.cate.pred - draw$district.cate.truth$district.cate.truth)^2)) / sd(y),
    pehe.school = sqrt(mean((school.cate.pred - draw$school.cate.truth$school.cate.truth)^2)) / sd(y), 
    model = paste(type, 'stan4bart - no p.score'), 
    tau = .tau
  )
  
}


get_balance <- function(rawdata, treat,estimand="ATT"){
  if(missing(rawdata)) stop("rawdata is required")
  if(missing(treat)) stop("treatment vector (treat) is required")
  cat("Balance diagnostics assume that the estimand is the",estimand,"\n")
  #
  #raw.dat <- data.frame(rawdata, treat = treat)
  covnames <- colnames(rawdata)
  if (is.null(covnames)){
    cat("No covariate names provided.  Generic names will be generated.")
    covnames = paste("v",c(1:ncol(rawdata)),sep="")
  }
  K <- length(covnames)
  diff.means <- matrix(NA, K, 5)
  var.t <- numeric(K)
  var.c <- numeric(K)
  std.denom <- numeric(K)
  binary <- rep(1,K)
  
  for (i in 1:K) {
    # separate means by group
    diff.means[i, 1] <- mean(rawdata[treat==1, i])
    diff.means[i, 2] <- mean(rawdata[treat==0, i])
    # separate variances by group == only used as input to calculations below
    var.t[i] <- var(rawdata[(treat == 1), i])
    var.c[i] <- var(rawdata[(treat == 0), i])
    # denominator in standardized difference calculations
    if(estimand=="ATE"){std.denom[i] <- sqrt((var.t[i]+var.c[i])/2)}
    else{
      std.denom[i] <- ifelse(estimand=="ATT",sqrt(var.t[i]),sqrt(var.c[i]))
    }
    # difference in means
    diff.means[i, 3] <- diff.means[i, 1] - diff.means[i, 2]
    # standardized difference in means (sign intact)
    diff.means[i, 4] <- abs(diff.means[i, 3]/std.denom[i])
    if(length(unique(rawdata[,covnames[i]]))>2){
      binary[i] = 0
      diff.means[i, 5] <- sqrt(var.c[i]/var.t[i])
    }
  }
  
  dimnames(diff.means) <- list(covnames, c("treat", "control", "unstd.diff",
                                           "abs.std.diff", "ratio"))
  
  
  return(diff.means)
}

