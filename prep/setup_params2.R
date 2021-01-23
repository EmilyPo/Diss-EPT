# params relating to the estimation of networks prior to simulation
# updated for chlamydia model 

setup_params <- function(ppopsize = 20000,
                         mcmcInterval = 1e5,
                         mcmcBurnin = 1e5, 
                         time.step = 7,
                         age.width = 45-15,
                         dissolution = ~offset(edges),
                         background_mortality = 2.52071e-05,
                         aging_out = 1/(age.width*(365/time.step)),
                         duration.marcoh = 847,
                         duration.other = 55,
                         shift.marcoh = 0.1511228,
                         shift.other = 0.14505
) {
  p <- list()
  p$ppopsize <- ppopsize
  p$mcmcInterval <- mcmcInterval
  p$mcmcBurnin <- mcmcBurnin
  p$time.step <- time.step
  p$age.width <- age.width
  p$dissolution <- dissolution
  p$background_mortality <- background_mortality
  p$aging_out <- aging_out
  p$duration.marcoh <- duration.marcoh
  p$duration.other <- duration.other
  p$shift.marcoh <- shift.marcoh
  p$shift.other <- shift.other
  p$mRate <- background_mortality + aging_out
  
  return(p)
}

