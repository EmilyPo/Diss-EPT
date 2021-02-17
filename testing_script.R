library(here)
library(EpiModelHIV)
library(ddaf)

n <- readRDS(here("fits", "netests.rds"))

sim <- netsim(n,
              param = param_het_ept(trackRels_sepNets = FALSE),
              init = init_het_reldur(),
              control = control_het_reldur(nsims = 3, nsteps=1000))

saveRDS(sim, here("calibration", "og_sim.rds"))

plot(sim, y=c("meanDegMar", "meanDegCohab"), legend=T)
plot(sim, y=c("meanDegOther"), legend=T)
plot(sim, y=c("num.feml", "num.male"), legend=T)

library(srvyr)
library(ggplot2)

mar <- get_dist(sim, 1)
mar[[1]];mar[[2]]

cohab <- get_dist(sim, 2)
cohab[[1]];cohab[[2]]

casual<- get_dist(sim, 3)
casual[[1]];casual[[2]]