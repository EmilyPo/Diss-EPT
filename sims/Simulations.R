library(EpiModelHIV)
library(here)
library(EpiModelHPC)
library(parallel)
library(doParallel)
epistats <- readRDS(here("prep", "epistats.rds"))

#### Concurrency By Sex ####
nets <- readRDS(here("fits", "calibrated_netests.rds"))

sim.base.0 <- netsim(nets,
                   param = param_het_ept(epistats, ept.coverage = 0),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.0, here("sims", "may2", "simBase0.rds"))
rm(sim.base.0)

sim.base <- netsim(nets,
                   param = param_het_ept(epistats, ept.coverage = 0.25),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base, here("sims", "may2", "simBase.rds"))
rm(sim.base)

sim.base.ept <- netsim(nets,
                       param = param_het_ept(epistats, ept.coverage = 0.5),
                       init = init_het_ct(),
                       control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept, here("sims", "may2", "simBaseEPT.rds"))
rm(sim.base.ept)

sim.base.ept.high <- netsim(nets,
                            param = param_het_ept(epistats, ept.coverage = 0.75),
                            init = init_het_ct(),
                            control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept.high, here("sims", "may2", "simBaseEPThigh.rds"))
rm(sim.base.ept.high)

sim.base.ept.100 <- netsim(nets,
                            param = param_het_ept(epistats, ept.coverage = 1,
                                                  ept.provision.pers.ong = 1, ept.provision.main.ong = 1),
                            init = init_het_ct(),
                            control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept.100, here("sims", "may2", "simBaseEPT100.rds"))
rm(sim.base.ept.100)
rm(nets)

#### Average Concurrency ####
nets <- readRDS(here("fits", "calibrated_netests_avgconc.rds"))

sim.avg.0 <- netsim(nets,
                     param = param_het_ept(epistats, ept.coverage = 0),
                     init = init_het_ct(),
                     control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.0, here("sims", "may2", "simAvg0.rds"))
rm(sim.avg.0)

sim.avg <- netsim(nets,
                   param = param_het_ept(epistats, ept.coverage = 0.25),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg, here("sims", "may2", "simAvg.rds"))
rm(sim.avg)

sim.avg.ept <- netsim(nets,
                       param = param_het_ept(epistats, ept.coverage = 0.5),
                       init = init_het_ct(),
                       control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.ept, here("sims", "may2", "simAvgEPT.rds"))
rm(sim.avg.ept)

sim.avg.ept.high <- netsim(nets,
                            param = param_het_ept(epistats, ept.coverage = 0.75),
                            init = init_het_ct(),
                            control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.ept.high, here("sims", "may2", "simAvgEPThigh.rds"))
rm(sim.avg.ept.high)

sim.avg.ept.100 <- netsim(nets,
                           param = param_het_ept(epistats, ept.coverage = 1,
                                                 ept.provision.pers.ong = 1, ept.provision.main.ong = 1),
                           init = init_het_ct(),
                           control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.ept.100, here("sims", "may2", "simAvgEPT100.rds"))
rm(sim.avg.ept.100)
rm(nets)

#### No Concurrency ####
nets <- readRDS(here("fits", "calibrated_netests_noconc.rds"))

sim.nc.0 <- netsim(nets,
                    param = param_het_ept(epistats, ept.coverage = 0),
                    init = init_het_ct(),
                    control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.0, here("sims", "may2", "simNc0.rds"))
rm(sim.nc.0)

sim.nc <- netsim(nets,
                  param = param_het_ept(epistats, ept.coverage = 0.25),
                  init = init_het_ct(),
                  control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc, here("sims", "may2", "simNc.rds"))
rm(sim.nc)

sim.nc.ept <- netsim(nets,
                      param = param_het_ept(epistats, ept.coverage = 0.5),
                      init = init_het_ct(),
                      control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept, here("sims", "may2", "simNcEPT.rds"))
rm(sim.nc.ept)

sim.nc.ept.high <- netsim(nets,
                           param = param_het_ept(epistats, ept.coverage = 0.75),
                           init = init_het_ct(),
                           control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept.high, here("sims", "may2", "simNcEPThigh.rds"))
rm(sim.nc.ept.high)

sim.nc.ept.100 <- netsim(nets,
                          param = param_het_ept(epistats, ept.coverage = 1,
                                                ept.provision.pers.ong = 1, ept.provision.main.ong = 1),
                          init = init_het_ct(),
                          control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept.100, here("sims", "may2", "simNcEPT100.rds"))
rm(sim.nc.ept.100)
rm(nets)




## -------- NO SEX AFTER EPT FLAG -----------------
#### Concurrency By Sex ####
nets <- readRDS(here("fits", "calibrated_netests.rds"))

sim.base <- netsim(nets,
                   param = param_het_ept(epistats, ept.coverage = 0.25, noeptsex = TRUE),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base, here("sims", "may2",  "noeptsex", "simBase.rds"))
rm(sim.base)

nets <- readRDS(here("fits", "calibrated_netests.rds"))
sim.base100 <- netsim(nets,
                   param = param_het_ept(epistats, noeptsex = TRUE, ept.coverage = 1,
                                         ept.provision.pers.ong = 1, ept.provision.main.ong = 1),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base100, here("sims", "may2",  "noeptsex", "simBase100.rds"))
rm(sim.base)

rm(nets)

## -------- Add'l Runs -----------------
## All Previous Partners + All Current Partners 
## 100 ept coverage w/ same cascade of 