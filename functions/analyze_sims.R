# summarize reinfections
# this function takes the simulation object and returns the median and IQR proportion of incident cases 
# that are reinfections over a given period of time
# two ways to measure: inf time to inf time, and previous diagnosis to inf time 
# stratified by sex 

reinf <- function(x, y=NULL, start=NULL, end=NULL){
  
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
  r12.f <- c(r12.f, r12)
  r6 <- x$epi$reinf.6.female[[i]][start:end]
  r6.f <- c(r6.f, r6)
  r3 <- x$epi$reinf.3.female[[i]][start:end]
  r3.f <- c(r3.f, r3)
  
  r12.d <- x$epi$reinf.12.female.diag[[i]][start:end]
  r12.f.diag <- c(r12.f.diag, r12.d)
  r6.d <- x$epi$reinf.6.female.diag[[i]][start:end]
  r6.f.diag <- c(r6.f.diag, r6.d)
  r3.d <- x$epi$reinf.3.female.diag[[i]][start:end]
  r3.f.diag <- c(r3.f.diag, r3.d)
  
  # males
  r12 <- x$epi$reinf.12.male[[i]][start:end]
  r12.m <- c(r12.m, r12)
  r6 <- x$epi$reinf.6.male[[i]][start:end]
  r6.m <- c(r6.m, r6)
  r3 <- x$epi$reinf.3.male[[i]][start:end]
  r3.m <- c(r3.m, r3)
  
  r12.d <- x$epi$reinf.12.male.diag[[i]][start:end]
  r12.m.diag <- c(r12.m.diag, r12.d)
  r6.d <- x$epi$reinf.6.male.diag[[i]][start:end]
  r6.m.diag <- c(r6.m.diag, r6.d)
  r3.d <- x$epi$reinf.3.male.diag[[i]][start:end]
  r3.m.diag <- c(r3.m.diag, r3.d)
  }

  dat <- cbind(r12.f, r6.f, r3.f, r12.f.diag, r6.f.diag, r3.f.diag,
              r12.m, r6.m, r3.m, r12.m.diag, r6.m.diag, r3.m.diag)
  
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
      r12.f2 <- c(r12.f2, r12)
      r6 <- x$epi$reinf.6.female[[i]][start:end]
      r6.f2 <- c(r6.f2, r6)
      r3 <- x$epi$reinf.3.female[[i]][start:end]
      r3.f2 <- c(r3.f2, r3)
      
      r12.d <- x$epi$reinf.12.female.diag[[i]][start:end]
      r12.f.diag2 <- c(r12.f.diag2, r12.d)
      r6.d <- x$epi$reinf.6.female.diag[[i]][start:end]
      r6.f.diag2 <- c(r6.f.diag2, r6.d)
      r3.d <- x$epi$reinf.3.female.diag[[i]][start:end]
      r3.f.diag2 <- c(r3.f.diag2, r3.d)
      
      # males
      r12 <- x$epi$reinf.12.male[[i]][start:end]
      r12.m2 <- c(r12.m2, r12)
      r6 <- x$epi$reinf.6.male[[i]][start:end]
      r6.m2 <- c(r6.m2, r6)
      r3 <- x$epi$reinf.3.male[[i]][start:end]
      r3.m2 <- c(r3.m2, r3)
      
      r12.d <- x$epi$reinf.12.male.diag[[i]][start:end]
      r12.m.diag2 <- c(r12.m.diag2, r12.d)
      r6.d <- x$epi$reinf.6.male.diag[[i]][start:end]
      r6.m.diag2 <- c(r6.m.diag2, r6.d)
      r3.d <- x$epi$reinf.3.male.diag[[i]][start:end]
      r3.m.diag2 <- c(r3.m.diag2, r3.d)
      
      dat <- cbind(c(r12.f, r12.f2), c(r6.f, r6.f2), c(r3.f, r3.f2),
                   c(r12.f.diag, r12.f.diag2),  c(r6.f.diag, r6.f.diag2), c(r3.f.diag, r3.f.diag2), 
                   c(r12.m, r12.m2), c(r6.m, r6.m2), c(r3.m, r3.m2),
                   c(r12.m.diag, r12.m.diag2), c(r6.m.diag, r6.m.diag2), c(r3.m.diag, r3.m.diag2))
    }
  }
  
  
  mean <- apply(dat, 2, mean, na.rm=T)
  sd <- apply(dat, 2, sd, na.rm=T)

  return(list(mean,sd))

}

irate <- function(x, y=NULL, start=NULL, end=NULL, denominator=100000){
  
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
    ir1.f <- c(ir1.f, ir1.1)
    
    ir2.1 <- x$epi$i.female.2[[i]][start:end]
    ir2.f <- c(ir2.f, ir2.1)
    
    ir3.1 <- x$epi$i.female.3[[i]][start:end]
    ir3.f <- c(ir3.f, ir3.1)
    
    ir4.1 <- x$epi$i.female.4[[i]][start:end]
    ir4.f <- c(ir4.f, ir4.1)
    
    ir5.1 <- x$epi$i.female.5[[i]][start:end]
    ir5.f <- c(ir5.f, ir5.1)
    
    ir6.1 <- x$epi$i.female.6[[i]][start:end]
    ir6.f <- c(ir6.f, ir6.1)
    
    # males
    ir1.1 <- x$epi$i.male.1[[i]][start:end]
    ir1.m <- c(ir1.m, ir1.1)
    
    ir2.1 <- x$epi$i.male.2[[i]][start:end]
    ir2.m <- c(ir2.m, ir2.1)
    
    ir3.1 <- x$epi$i.male.3[[i]][start:end]
    ir3.m <- c(ir3.m, ir3.1)
    
    ir4.1 <- x$epi$i.male.4[[i]][start:end]
    ir4.m <- c(ir4.m, ir4.1)
    
    ir5.1 <- x$epi$i.male.5[[i]][start:end]
    ir5.m <- c(ir5.m, ir5.1)
    
    ir6.1 <- x$epi$i.male.6[[i]][start:end]
    ir6.m <- c(ir6.m, ir6.1)
  }
  
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
    } else
      dat <- cbind(ir1.f, ir2.f, ir3.f, ir4.f, ir5.f, ir6.f,
                   ir1.m, ir2.m, ir3.m, ir4.m, ir5.m, ir6.m)
  
  dat <- dat*52*denominator
  
  mean <- apply(dat, 2, mean, na.rm=T)
  sd <- apply(dat, 2, sd, na.rm=T)
  
  return(list(mean,sd))
}
