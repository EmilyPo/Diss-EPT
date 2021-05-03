# summarize reinfections
# this function takes the simulation object and returns the median and IQR proportion of incident cases 
# that are reinfections over a given period of time
# two ways to measure: inf time to inf time, and previous diagnosis to inf time 
# stratified by sex 

reinf <- function(x, y=NULL, start=2000, end=NULL){
  
  nsims <- x$control$nsims
  
  if (is.null(start)) {start <- x$param$ept.start + 52*5}
  if(is.null(end)) {end <- x$control$nsteps }
  

  r12.f <- NULL
  r6.f <- NULL
  r3.f <- NULL
  r12.f.diag <- NULL
  r6.f.diag <- NULL
  r3.f.diag <- NULL
  
  r12.m <- NULL
  r6.m <- NULL
  r3.m <- NULL
  r12.m.diag <- NULL
  r6.m.diag <- NULL
  r3.m.diag <- NULL
  
  for (i in 1:nsims){
  # females
  r12 <- x$epi$reinf.12.female[[i]][start:end]
  r12.f <- cbind(r12.f, r12)
  
  r6 <- x$epi$reinf.6.female[[i]][start:end]
  r6.f <- cbind(r6.f, r6)
  
  r3 <- x$epi$reinf.3.female[[i]][start:end]
  r3.f <- cbind(r3.f, r3)
  
  r12.d <- x$epi$reinf.12.female.diag[[i]][start:end]
  r12.f.diag <- cbind(r12.f.diag, r12.d)
  r6.d <- x$epi$reinf.6.female.diag[[i]][start:end]
  r6.f.diag <- cbind(r6.f.diag, r6.d)
  r3.d <- x$epi$reinf.3.female.diag[[i]][start:end]
  r3.f.diag <- cbind(r3.f.diag, r3.d)
  

  # males
  r12 <- x$epi$reinf.12.male[[i]][start:end]
  r12.m <- cbind(r12.m, r12)
  r6 <- x$epi$reinf.6.male[[i]][start:end]
  r6.m <- cbind(r6.m, r6)
  r3 <- x$epi$reinf.3.male[[i]][start:end]
  r3.m <- cbind(r3.m, r3)
  
  r12.d <- x$epi$reinf.12.male.diag[[i]][start:end]
  r12.m.diag <- cbind(r12.m.diag, r12.d)
  r6.d <- x$epi$reinf.6.male.diag[[i]][start:end]
  r6.m.diag <- cbind(r6.m.diag, r6.d)
  r3.d <- x$epi$reinf.3.male.diag[[i]][start:end]
  r3.m.diag <- cbind(r3.m.diag, r3.d)

  }

  fems <- rbind(apply(r12.f, 2, mean, na.rm=T), 
                apply(r6.f, 2, mean, na.rm=T),
                apply(r3.f, 2, mean, na.rm=T),
                apply(r12.f.diag, 2, mean, na.rm=T),
                apply(r6.f.diag, 2, mean, na.rm=T),
                apply(r3.f.diag, 2, mean, na.rm=T))
  
  
  males <- rbind(apply(r12.m, 2, mean, na.rm=T), 
                 apply(r6.m, 2, mean, na.rm=T),
                 apply(r3.m, 2, mean, na.rm=T),
                 apply(r12.m.diag, 2, mean, na.rm=T),
                 apply(r6.m.diag, 2, mean, na.rm=T),
                 apply(r3.m.diag, 2, mean, na.rm=T))
  
  fems.mean <- apply(fems, 1, mean, na.rm=T)
  fems.se <- apply(fems, 1, sd, na.rm=T)/sqrt(ncol(fems))
  fems.low <- fems.mean - 2*fems.se 
  fems.high <- fems.mean + 2*fems.se 
  
  females <- data.frame(mean = fems.mean,
                         low = fems.low,
                         high = fems.high)
  
  males.mean <- apply(males, 1, mean, na.rm=T)
  males.se <- apply(males, 1, sd, na.rm=T)/sqrt(ncol(males))
  males.low <- males.mean - 2*males.se 
  males.high <- males.mean + 2*males.se 
  
  males <- data.frame(mean = males.mean,
                        low = males.low,
                        high = males.high)
  
  
  if (!is.null(y)) {
    x <- y
    r12.f2 <- NULL
    r6.f2 <- NULL
    r3.f2 <- NULL
    r12.f.diag2 <- NULL
    r6.f.diag2 <- NULL
    r3.f.diag2 <- NULL
    
    r12.m2 <- NULL
    r6.m2 <- NULL
    r3.m2 <- NULL
    r12.m.diag2 <- NULL
    r6.m.diag2 <- NULL
    r3.m.diag2 <- NULL
    
    for (i in 1:nsims){
      # females
      r12 <- x$epi$reinf.12.female[[i]][start:end]
      r12.f2 <- cbind(r12.f2, r12)
      
      r6 <- x$epi$reinf.6.female[[i]][start:end]
      r6.f2 <- cbind(r6.f2, r6)
      r3 <- x$epi$reinf.3.female[[i]][start:end]
      r3.f2 <- cbind(r3.f2, r3)
      
      r12.d <- x$epi$reinf.12.female.diag[[i]][start:end]
      r12.f.diag2 <- cbind(r12.f.diag2, r12.d)
      r6.d <- x$epi$reinf.6.female.diag[[i]][start:end]
      r6.f.diag2 <- cbind(r6.f.diag2, r6.d)
      r3.d <- x$epi$reinf.3.female.diag[[i]][start:end]
      r3.f.diag2 <- cbind(r3.f.diag2, r3.d)
      
      # males
      r12 <- x$epi$reinf.12.male[[i]][start:end]
      r12.m2 <- cbind(r12.m2, r12)
      r6 <- x$epi$reinf.6.male[[i]][start:end]
      r6.m2 <- cbind(r6.m2, r6)
      r3 <- x$epi$reinf.3.male[[i]][start:end]
      r3.m2 <- cbind(r3.m2, r3)
      
      r12.d <- x$epi$reinf.12.male.diag[[i]][start:end]
      r12.m.diag2 <- cbind(r12.m.diag2, r12.d)
      r6.d <- x$epi$reinf.6.male.diag[[i]][start:end]
      r6.m.diag2 <- cbind(r6.m.diag2, r6.d)
      r3.d <- x$epi$reinf.3.male.diag[[i]][start:end]
      r3.m.diag2 <- cbind(r3.m.diag2, r3.d)
      
      
      dat <- cbind(c(r12.f, r12.f2), c(r6.f, r6.f2), c(r3.f, r3.f2),
                   c(r12.f.diag, r12.f.diag2),  c(r6.f.diag, r6.f.diag2), c(r3.f.diag, r3.f.diag2), 
                   c(r12.m, r12.m2), c(r6.m, r6.m2), c(r3.m, r3.m2),
                   c(r12.m.diag, r12.m.diag2), c(r6.m.diag, r6.m.diag2), c(r3.m.diag, r3.m.diag2))
    }
  }

  return(list(females, males))

}

irate <- function(x, y=NULL, start=2000, end=NULL, denominator=100000){
  
  nsims <- x$control$nsims
  
  if (is.null(start)) {start <- x$param$ept.start + 52*5}
  if(is.null(end)) {end <- x$control$nsteps }
  
  
  ir1.f <- NULL
  ir2.f <- NULL
  ir3.f <- NULL
  ir4.f <- NULL
  ir5.f <- NULL
  ir6.f <- NULL

  ir1.m <- NULL
  ir2.m <- NULL
  ir3.m <- NULL
  ir4.m <- NULL
  ir5.m <- NULL
  ir6.m <- NULL
  
  for (i in 1:nsims){
    # females
    ir1.1 <- x$epi$i.female.1[[i]][start:end]
    ir1.f <- cbind(ir1.f, ir1.1)
    
    ir2.1 <- x$epi$i.female.2[[i]][start:end]
    ir2.f <- cbind(ir2.f, ir2.1)
    
    ir3.1 <- x$epi$i.female.3[[i]][start:end]
    ir3.f <- cbind(ir3.f, ir3.1)
    
    ir4.1 <- x$epi$i.female.4[[i]][start:end]
    ir4.f <- cbind(ir4.f, ir4.1)
    
    ir5.1 <- x$epi$i.female.5[[i]][start:end]
    ir5.f <- cbind(ir5.f, ir5.1)
    
    ir6.1 <- x$epi$i.female.6[[i]][start:end]
    ir6.f <- cbind(ir6.f, ir6.1)

    
    # males
    ir1.1 <- x$epi$i.male.1[[i]][start:end]
    ir1.m <- cbind(ir1.m, ir1.1)
    
    ir2.1 <- x$epi$i.male.2[[i]][start:end]
    ir2.m <- cbind(ir2.m, ir2.1)
    
    ir3.1 <- x$epi$i.male.3[[i]][start:end]
    ir3.m <- cbind(ir3.m, ir3.1)
    
    ir4.1 <- x$epi$i.male.4[[i]][start:end]
    ir4.m <- cbind(ir4.m, ir4.1)
    
    ir5.1 <- x$epi$i.male.5[[i]][start:end]
    ir5.m <- cbind(ir5.m, ir5.1)
    
    ir6.1 <- x$epi$i.male.6[[i]][start:end]
    ir6.m <- cbind(ir6.m, ir6.1)
  }
  
  
  fems <- rbind(apply(ir1.f, 2, mean, na.rm=T), 
                apply(ir2.f, 2, mean, na.rm=T),
                apply(ir3.f, 2, mean, na.rm=T),
                apply(ir4.f, 2, mean, na.rm=T),
                apply(ir5.f, 2, mean, na.rm=T),
                apply(ir6.f, 2, mean, na.rm=T))
  
  fems <- fems * 52* denominator
  
    males <- rbind(apply(ir1.m, 2, mean, na.rm=T), 
                  apply(ir2.m, 2, mean, na.rm=T),
                  apply(ir3.m, 2, mean, na.rm=T),
                  apply(ir4.m, 2, mean, na.rm=T),
                  apply(ir5.m, 2, mean, na.rm=T),
                  apply(ir6.m, 2, mean, na.rm=T))
    
    males <- males * 52* denominator
    
    fems.mean <- apply(fems, 1, mean, na.rm=T)
    fems.se <- apply(fems, 1, sd, na.rm=T)/sqrt(ncol(fems))
    fems.low <- fems.mean - 2*fems.se 
    fems.high <- fems.mean + 2*fems.se 
    
    females <- data.frame(mean = fems.mean,
                          low = fems.low,
                          high = fems.high)
    
    males.mean <- apply(males, 1, mean, na.rm=T)
    males.se <- apply(males, 1, sd, na.rm=T)/sqrt(ncol(males))
    males.low <- males.mean - 2*males.se 
    males.high <- males.mean + 2*males.se 
    
    males <- data.frame(mean = males.mean,
                        low = males.low,
                        high = males.high)
  
  if (!is.null(y)) {
    x <- y
    ir1.f2 <- NULL
    ir2.f2 <- NULL
    ir3.f2 <- NULL
    ir4.f2 <- NULL
    ir5.f2 <- NULL
    ir6.f2 <- NULL
    
    ir1.m2 <- NULL
    ir2.m2 <- NULL
    ir3.m2 <- NULL
    ir4.m2 <- NULL
    ir5.m2 <- NULL
    ir6.m2 <- NULL
    
    for (i in 1:nsims){
      # females
      ir1.1 <- x$epi$ir.uct.female.1519[[i]][start:end]
      ir1.f2 <- c(ir1.f2, ir1.1)
      
      ir2.1 <- x$epi$ir.uct.female.2024[[i]][start:end]
      ir2.f2 <- c(ir2.f2, ir2.1)
      
      ir3.1 <- x$epi$ir.uct.female.2529[[i]][start:end]
      ir3.f2 <- c(ir3.f2, ir3.1)
      
      ir4.1 <- x$epi$ir.uct.female.3034[[i]][start:end]
      ir4.f2 <- c(ir4.f2, ir4.1)
      
      ir5.1 <- x$epi$ir.uct.female.3539[[i]][start:end]
      ir5.f2 <- c(ir5.f2, ir5.1)
      
      ir6.1 <- x$epi$ir.uct.female.4044[[i]][start:end]
      ir6.f2 <- c(ir6.f2, ir6.1)
      
      # males
      ir1.1 <- x$epi$ir.uct.male.1519[[i]][start:end]
      ir1.m2 <- c(ir1.m2, ir1.1)
      
      ir2.1 <- x$epi$ir.uct.male.2024[[i]][start:end]
      ir2.m2 <- c(ir2.m2, ir2.1)
      
      ir3.1 <- x$epi$ir.uct.male.2529[[i]][start:end]
      ir3.m2 <- c(ir3.m2, ir3.1)
      
      ir4.1 <- x$epi$ir.uct.male.3034[[i]][start:end]
      ir4.m2 <- c(ir4.m2, ir4.1)
      
      ir5.1 <- x$epi$ir.uct.male.3539[[i]][start:end]
      ir5.m2 <- c(ir5.m2, ir5.1)
      
      ir6.1 <- x$epi$ir.uct.male.4044[[i]][start:end]
      ir6.m2 <- c(ir6.m2, ir6.1)
    }
      
    dat <- cbind(c(ir1.f, ir1.f2), c(ir2.f, ir2.f2), c(ir3.f, ir3.f2),
                 c(ir4.f, ir4.f2), c(ir5.f, ir5.f2), c(ir6.f, ir6.f2),
                 c(ir1.m, ir1.m2), c(ir2.m, ir2.m2), c(ir3.m, ir3.m2),
                 c(ir4.m, ir4.m2),c(ir5.m, ir5.m2), c(ir6.m, ir6.m2))
    } 
  
  return(list(females,males))
}

prev <- function(x, y=NULL, start=2000, end=NULL){
  
  nsims <- x$control$nsims
  
  if (is.null(start)) {start <- x$param$ept.start + 52*5}
  if(is.null(end)) {end <- x$control$nsteps }
  
  
  fem <- NULL
  male <- NULL
  
  for (i in 1:nsims){
    # females
    fem1 <- x$epi$i.prev.feml[[i]][start:end]
    fem <- cbind(fem, fem1)
    
    male1 <- x$epi$i.prev.male[[i]][start:end]
    male <- cbind(male, male1)

  }
  
  fs <- apply(fem, 2, mean, na.rm=T)
  fems.mean <- mean(fs)
  fems.se <- sd(fs)/sqrt(length(fs))
  fems.low <- fems.mean - 2*fems.se 
  fems.high <- fems.mean + 2*fems.se 
  
  ms <- apply(male, 2, mean, na.rm=T)
  males.mean <- mean(ms)
  males.se <- sd(ms)/sqrt(length(ms))
  males.low <- males.mean - 2*males.se 
  males.high <- males.mean + 2*males.se 
  
  females <- data.frame(mean = fems.mean,
                        low = fems.low,
                        high = fems.high)
  
  males <- data.frame(mean = males.mean,
                      low = males.low,
                      high = males.high)
  
  
  if (!is.null(y)) {
    x <- y
    
    fem2 <- NULL
    male2 <- NULL
    
    for (i in 1:nsims){
      # females
      fem2.1 <- x$epi$i.prev.feml[[i]][start:end]
      fem2 <- c(fem2, fem2.1)
      
      male2.1 <- x$epi$i.prev.male[[i]][start:end]
      male2 <- c(male2, male2.1)
    }
    
    dat <- cbind(c(fem, fem2), c(male, male2))
  } 

  
  return(list(females,males))
}
