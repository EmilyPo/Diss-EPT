# summarize reinfections
# this function takes the simulation object and returns the median and IQR proportion of incident cases 
# that are reinfections over a given period of time
# two ways to measure: inf time to inf time, and previous diagnosis to inf time 
# stratified by sex 

reinf <- function(x, start=NULL, end=NULL){
  
  if (is.null(start)) {start <- x$param$eptstart + 52*5}
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
  r12 <- x$epi$reinf.12.female[[i]][start:nsteps]
  r12.f <- c(r12.f, r12)
  r6 <- x$epi$reinf.6.female[[i]][start:nsteps]
  r6.f <- c(r6.f, r6)
  r3 <- x$epi$reinf.3.female[[i]][start:nsteps]
  r3.f <- c(r3.f, r3)
  
  r12.d <- x$epi$reinf.12.female.diag[[i]][start:nsteps]
  r12.f.diag <- c(r12.f.diag, r12.d)
  r6.d <- x$epi$reinf.6.female.diag[[i]][start:nsteps]
  r6.f.diag <- c(r6.f.diag, r6.d)
  r3.d <- x$epi$reinf.3.female.diag[[i]][start:nsteps]
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
  

  dat <- cbind(r12.f, r6.f, r3.f, r12.f.diag, r6.f.diag, r3.f.diag,
              r12.m, r6.m, r3.m, r12.m.diag, r6.m.diag, r3.m.diag)
  
  med <- apply(dat, 2, median)
  iqr <- apply(dat, 2, IQR)
  
  return(list(med,iqr))
  
  }
}

irate <- function(x, start=NULL, end=NULL){
  
  if (is.null(start)) {start <- x$param$eptstart + 52*5}
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
    r12 <- x$epi$reinf.12.female[[i]][start:nsteps]
    r12.f <- c(r12.f, r12)
    r6 <- x$epi$reinf.6.female[[i]][start:nsteps]
    r6.f <- c(r6.f, r6)
    r3 <- x$epi$reinf.3.female[[i]][start:nsteps]
    r3.f <- c(r3.f, r3)
    
    r12.d <- x$epi$reinf.12.female.diag[[i]][start:nsteps]
    r12.f.diag <- c(r12.f.diag, r12.d)
    r6.d <- x$epi$reinf.6.female.diag[[i]][start:nsteps]
    r6.f.diag <- c(r6.f.diag, r6.d)
    r3.d <- x$epi$reinf.3.female.diag[[i]][start:nsteps]
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
    
    
    dat <- cbind(r12.f, r6.f, r3.f, r12.f.diag, r6.f.diag, r3.f.diag,
                 r12.m, r6.m, r3.m, r12.m.diag, r6.m.diag, r3.m.diag)
    
    med <- apply(dat, 2, median)
    iqr <- apply(dat, 2, IQR)
    
    return(list(med,iqr))
    
  }
}