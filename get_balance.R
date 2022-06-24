
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