# calibration for open population 
# three networks - marriage/cohab, casual, inst
# calibration really only necessary for marcoh & casual (dur > 1)

library(here)
library(EpiModelHIV)
source(here("functions", "ddaf_functions.R"))
library(srvyr)
library(ggplot2)

n <- readRDS(here("fits", "netests.rds"))
n2 <- n

# 1: re-do EDA for log(D) instead of log(D-1)
n2[[1]]$coef.form[1] <- n[[1]]$coef.form.crude[1] - log(n[[1]]$coef.diss$duration)
n2[[2]]$coef.form[1] <- n[[2]]$coef.form.crude[1] - log(n[[2]]$coef.diss$duration)

# 2: adjust coefficients 
# mar
n2[[1]]$coef.form["nodefactor.agecat.1"] <- n[[1]]$coef.form["nodefactor.agecat.1"] + 1.5
n2[[1]]$coef.form["nodefactor.agecat.2"] <- n[[1]]$coef.form["nodefactor.agecat.2"] + 0.7
n2[[1]]$coef.form["nodefactor.agecat.3"] <- n[[1]]$coef.form["nodefactor.agecat.3"] 
#n2[[1]]$coef.form["nodefactor.age>30.TRUE"] <- n2[[1]]$coef.form["nodefactor.age>30.TRUE"] - 0.5

#casual
n2[[2]]$coef.form["nodefactor.ageM.15"]  <- n[[2]]$coef.form["nodefactor.ageM.15"] + 1.5
n2[[2]]$coef.form["nodefactor.ageM.18"]  <- n[[2]]$coef.form["nodefactor.ageM.18"] + 0.5

n2[[2]]$coef.form["nodefactor.ageF.15"]  <- n[[2]]$coef.form["nodefactor.ageF.15"] + 0.5

#sim <- netsim(n2,
#               param = param_het_ept(trackRels_sepNets = TRUE),
#               init = init_het_reldur(),
#               control = control_het_reldur(nsims = 1, nsteps=52*30))

#saveRDS(sim, here("fits", "calibration_sim_inprogress.rds"))

# degree plots 
plot(sim, y=c("meanDegMarcoh"), legend=T)
plot(sim, y=c("meanDegOther"), legend=T)

mdist <- get_dist(sim, network=1, ndynamic=2)
mdist[[1]];mdist[[2]]

odist <- get_dist(sim, network=2, ndynamic=2)
odist[[1]];odist[[2]]

# what is the mean age of extant ties in each network
# actually probably have to adjust marriages for if they started at beginning of sim
# those should ties should have previous duration 

r <- active_rels(sim)
m <- r[[1]]
c <- r[[2]]

mean(m$len)
mean(c$len)

# decrease underlying dissolution probability 
n2[[1]]$coef.diss$coef.adj <- log(n[[1]]$coef.diss$duration*5)
n2[[2]]$coef.diss$coef.adj <- log(n[[2]]$coef.diss$duration*1.03)

simDur <- netsim(n2,
              param = param_het_ept(trackRels_sepNets = TRUE),
              init = init_het_reldur(),
              control = control_het_reldur(nsims = 5, nsteps=52*30))

saveRDS(simDur, here("fits", "calibration_sim_duration.rds"))

plot(simDur, y=c("meanDegMarcoh"), legend=T)
plot(simDur, y=c("meanDegOther"), legend=T)

mdist <- get_dist(simDur, network=1, ndynamic=2)
mdist[[1]];mdist[[2]]

odist <- get_dist(simDur, network=2, ndynamic=2)
odist[[1]];odist[[2]]

# what is the mean age of extant ties in each network
# actually probably have to adjust marriages for if they started at beginning of sim
# those should ties should have previous duration 

r <- active_rels(simDur)
m <- r[[1]]
c <- r[[2]]

mean(m$len)
mean(c$len)
