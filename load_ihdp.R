load_ihdp <- function(rct = FALSE){
  dat <- read_csv('data/ihdp.csv')
  dat <- subset(dat, treat != 1 | momwhite != 0)
  dat$momage[dat$momage < 15] <- 15
  dat$momage[dat$momage > 40] <- 40
  dat$g1 <- as.factor(dat$momage)
  dat$g2 <- as.factor(dat$site.num)  
    
  covs.cont <- c("bw", "b.head", "preterm", "birth.o", "nnhealth")
  covs.cat  <- c("sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                   "cig", "first", "booze", "drugs", "work.dur", "prenatal")
  
  set.seed(2)
  if(isTRUE(rct)) ihdp$treat <- rbinom(nrow(dat), 1, .5)
    
  dat <- dat[, c(covs.cat, covs.cont, 'treat', 'g1', 'g2', 'momage')]
  
  
  return(dat)  
    
}
