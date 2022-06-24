load_ihdp <- function(){
  dat <- read_csv('data/ihdp.csv')
  dat <- subset(dat, treat != 1 | momwhite != 0)
  # truncate momage to group 1
    dat <- dat %>% 
    mutate(g1 = as.factor(if_else(momage < 15 , 15, momage)), 
           g1 = as.factor(if_else(momage > 40, 40 , momage)), 
           g2 = as.factor(site.num)) %>% 
    select(-momage)
    
    covs.cont <- c("bw", "b.head", "preterm", "birth.o", "nnhealth")
    covs.cat  <- c("sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                   "cig", "first", "booze", "drugs", "work.dur", "prenatal")
    
    dat <- dat[, c(covs.cat, covs.cont, 'treat', 'g1', 'g2')]
  
  
  return(dat)  
    
}
