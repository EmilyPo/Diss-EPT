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
n2[[3]]$coef.form[1] <- n[[3]]$coef.form.crude[1] - log(n[[3]]$coef.diss$duration)

# 2: adjust coefficients 
# mar
n2[[1]]$coef.form["nodefactor.age<30.TRUE"] <- n[[1]]$coef.form["nodefactor.age<30.TRUE"] + 1.8
n2[[1]]$coef.form["nodefactor.age>30.TRUE"] <- n[[1]]$coef.form["nodefactor.age>30.TRUE"] + .3

#cohab
n2[[2]]$coef.form["nodefactor.floor(age).15"] <- n2[[2]]$coef.form["nodefactor.floor(age).15"] -10 
n2[[2]]$coef.form["nodefactor.floor(age).16"] <- n2[[2]]$coef.form["nodefactor.floor(age).16"] -10
n2[[2]]$coef.form["nodefactor.floor(age).17"] <- n2[[2]]$coef.form["nodefactor.floor(age).17"] -10
n2[[2]]$coef.form["nodefactor.floor(age).21"] <- n2[[2]]$coef.form["nodefactor.floor(age).21"] + 0.5
n2[[2]]$coef.form["nodefactor.floor(age).22"] <- n2[[2]]$coef.form["nodefactor.floor(age).22"] + 1
n2[[2]]$coef.form["nodefactor.floor(age).35"] <- n2[[2]]$coef.form["nodefactor.floor(age).35"] - 5
n2[[2]]$coef.form["nodefactor.floor(age).36"] <- n2[[2]]$coef.form["nodefactor.floor(age).36"] - 5
n2[[2]]$coef.form["nodefactor.floor(age).37"] <- n2[[2]]$coef.form["nodefactor.floor(age).37"] - 5
n2[[2]]$coef.form["nodefactor.floor(age).38"] <- n2[[2]]$coef.form["nodefactor.floor(age).38"] - 5
n2[[2]]$coef.form["nodefactor.floor(age).39"] <- n2[[2]]$coef.form["nodefactor.floor(age).39"] - 5

#casual
n2[[3]]$coef.form["nodefactor.ageM.15"]  <- n[[3]]$coef.form["nodefactor.ageM.15"] + .6
n2[[3]]$coef.form["nodefactor.ageM.17"]  <- n[[3]]$coef.form["nodefactor.ageM.17"] + .6
n2[[3]]$coef.form["nodefactor.ageM.19"]  <- n[[3]]$coef.form["nodefactor.ageM.19"] 

n2[[3]]$coef.form["nodefactor.ageF.15"]  <- n[[3]]$coef.form["nodefactor.ageF.15"] + .4
n2[[3]]$coef.form["nodefactor.ageF.16"]  <- n[[3]]$coef.form["nodefactor.ageF.16"] + .3
n2[[3]]$coef.form["nodefactor.ageF.19"]  <- n[[3]]$coef.form["nodefactor.ageF.19"] - 1.5


sim2 <- netsim(n2,
              param = param_het_ept(trackRels_sepNets = TRUE),
              init = init_het_reldur(),
              control = control_het_reldur(nsims = 1, nsteps=52*60))

test <- netsim(n2,
               param = param_het_ept(trackRels_sepNets = FALSE),
               init = init_het_reldur(),
               control = control_het_reldur(nsims = 1, nsteps=2))

saveRDS(sim2, here("calibration", "calibration_sim_inprogress.rds"))

plot(sim2, y=c("meanDegMar", "meanDegCohab"), legend=T)
plot(sim2, y=c("meanDegOther"), legend=T)
plot(sim2, y=c("num.feml", "num.male"), legend=T)

# standardize ylim xlim

mdist <- get_dist(sim2,1)
mdist[[1]];mdist[[2]]

cdist <- get_dist(sim2,2)
cdist[[1]];cdist[[2]]

odist <- get_dist(sim2,3)
odist[[1]];odist[[2]]

# next steps
# also try making all olderpartnerC = 0
# calibrate duration



get_prob <- function(x){1/(1+exp(-x))}
get_logodds <- function(x){log(x/(1-x))}









# adjust nodecov<30 for m and fs
old_prob_f <- get_prob(n[[1]]$coef.form["nodecov.ageF<30"])
old_prob_m <- get_prob(n[[1]]$coef.form["nodecov.ageF<30"])


                       