library(EpiModelHIV)
library(here)
library(EpiModelHPC)
library(parallel)
library(doParallel)
epistats <- readRDS(here("prep", "epistats.rds"))
# ----- w/ PT cascade ---------
#### Concurrency By Sex ####
nets <- readRDS(here("fits", "calibrated_netests.rds"))

sim.base.0 <- netsim(nets,
                   param = param_het_ept(epistats, ept.coverage = 0),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.0, here("sims", "may2", "noeptsex", "simBase0.rds"))
rm(sim.base.0)

sim.base <- netsim(nets,
                   param = param_het_ept(epistats, ept.coverage = 0.25),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base, here("sims", "may2", "noeptsex", "simBase2.rds"))
rm(sim.base)

sim.base.ept <- netsim(nets,
                       param = param_het_ept(epistats, ept.coverage = 0.5),
                       init = init_het_ct(),
                       control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept, here("sims", "may2", "noeptsex", "simBaseEPT.rds"))
rm(sim.base.ept)

sim.base.ept.high <- netsim(nets,
                            param = param_het_ept(epistats, ept.coverage = 0.75),
                            init = init_het_ct(),
                            control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept.high, here("sims", "may2", "noeptsex", "simBaseEPThigh.rds"))
rm(sim.base.ept.high)

sim.base.ept.highest <- netsim(nets,
                            param = param_het_ept(epistats, ept.coverage = 1),
                            init = init_het_ct(),
                            control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept.highest, here("sims", "may2", "noeptsex", "simBaseEPThighest.rds"))
rm(sim.base.ept.highest)

sim.base.ept.100.ong <- netsim(nets,
                            param = param_het_ept(epistats, ept.coverage = 1,
                                                  ept.provision.pers.ong = 1, ept.provision.main.ong = 1),
                            init = init_het_ct(),
                            control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept.100.ong, here("sims", "may2",  "noeptsex", "simBaseEPT100ong.rds"))
rm(sim.base.ept.100.ong)
rm(nets)

#### Average Concurrency ####
nets <- readRDS(here("fits", "calibrated_netests_avgconc.rds"))

sim.avg.0 <- netsim(nets,
                     param = param_het_ept(epistats, ept.coverage = 0),
                     init = init_het_ct(),
                     control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.0, here("sims", "may2",  "noeptsex",  "simAvg0.rds"))
rm(sim.avg.0)

sim.avg <- netsim(nets,
                   param = param_het_ept(epistats, ept.coverage = 0.25),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg, here("sims", "may2",  "noeptsex",  "simAvg.rds"))
rm(sim.avg)

sim.avg.ept <- netsim(nets,
                       param = param_het_ept(epistats, ept.coverage = 0.5),
                       init = init_het_ct(),
                       control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.ept, here("sims", "may2",  "noeptsex",  "simAvgEPT.rds"))
rm(sim.avg.ept)

sim.avg.ept.high <- netsim(nets,
                            param = param_het_ept(epistats, ept.coverage = 0.75),
                            init = init_het_ct(),
                            control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.ept.high, here("sims", "may2",  "noeptsex",  "simAvgEPThigh.rds"))
rm(sim.avg.ept.high)

sim.avg.ept.highest <- netsim(nets,
                           param = param_het_ept(epistats, ept.coverage = 1),
                           init = init_het_ct(),
                           control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.ept.highest, here("sims", "may2",  "noeptsex",  "simAvgEPThighest.rds"))
rm(sim.avg.ept.highest)


sim.avg.ept.100.ong <- netsim(nets,
                           param = param_het_ept(epistats, ept.coverage = 1,
                                                 ept.provision.pers.ong = 1, ept.provision.main.ong = 1),
                           init = init_het_ct(),
                           control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.ept.100.ong, here("sims", "may2",  "noeptsex",  "simAvgEPT100ong.rds"))
rm(sim.avg.ept.100.ong)
rm(nets)

#### No Concurrency ####
nets <- readRDS(here("fits", "calibrated_netests_noconc.rds"))

sim.nc.0 <- netsim(nets,
                    param = param_het_ept(epistats, ept.coverage = 0),
                    init = init_het_ct(),
                    control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.0, here("sims", "may2", "noeptsex", "simNc0.rds"))
rm(sim.nc.0)

sim.nc <- netsim(nets,
                  param = param_het_ept(epistats, ept.coverage = 0.25),
                  init = init_het_ct(),
                  control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc, here("sims", "may2", "noeptsex", "simNc.rds"))
rm(sim.nc)

sim.nc.ept <- netsim(nets,
                      param = param_het_ept(epistats, ept.coverage = 0.5),
                      init = init_het_ct(),
                      control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept, here("sims", "may2",  "noeptsex",  "simNcEPT.rds"))
rm(sim.nc.ept)

sim.nc.ept.high <- netsim(nets,
                           param = param_het_ept(epistats, ept.coverage = 0.75),
                           init = init_het_ct(),
                           control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept.high, here("sims", "may2",  "noeptsex",  "simNcEPThigh.rds"))
rm(sim.nc.ept.high)

sim.nc.ept.highest <- netsim(nets,
                          param = param_het_ept(epistats, ept.coverage = 1),
                          init = init_het_ct(),
                          control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept.highest, here("sims", "may2",  "noeptsex",  "simNcEPThighest.rds"))
rm(sim.nc.ept.highest)


sim.nc.ept.100ong <- netsim(nets,
                          param = param_het_ept(epistats, ept.coverage = 1,
                                                ept.provision.pers.ong = 1, ept.provision.main.ong = 1),
                          init = init_het_ct(),
                          control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept.100ong, here("sims", "may2",  "noeptsex",  "simNcEPT100ong.rds"))
rm(sim.nc.ept.100ong)
rm(nets)

## -------- NO PREVIOUS PARTERS, ALL CURRENT PARTNERS -----------------
#### Concurrency By Sex ####
nets <- readRDS(here("fits", "calibrated_netests.rds"))

sim.base.0 <- netsim(nets,
                     param = param_het_ept(epistats, ept.coverage = 0, 
                                           ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                           keep.window = 0),
                     init = init_het_ct(),
                     control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.0, here("sims", "may2", "current", "simBase02.rds"))
rm(sim.base.0)

sim.base <- netsim(nets,
                   param = param_het_ept(epistats, ept.coverage = 0.25,
                                         ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                         keep.window = 0),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base, here("sims", "may2", "current", "simBase2.rds"))
rm(sim.base)

sim.base.ept <- netsim(nets,
                       param = param_het_ept(epistats, ept.coverage = 0.5,
                                             ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                             keep.window = 0),
                       init = init_het_ct(),
                       control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept, here("sims", "may2", "current", "simBaseEPT2.rds"))
rm(sim.base.ept)

sim.base.ept.high <- netsim(nets,
                            param = param_het_ept(epistats, ept.coverage = 0.75,
                                                  ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                                  keep.window = 0),
                            init = init_het_ct(),
                            control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept.high, here("sims", "may2", "current", "simBaseEPThigh2.rds"))
rm(sim.base.ept.high)

sim.base.ept.100 <- netsim(nets,
                               param = param_het_ept(epistats, ept.coverage = 1,
                                                     ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                                     keep.window = 0),
                               init = init_het_ct(),
                               control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.base.ept.100, here("sims", "may2", "current", "simBaseEPT1002.rds"))
rm(sim.base.ept.100)
rm(nets)

#### Average Concurrency ####
nets <- readRDS(here("fits", "calibrated_netests_avgconc.rds"))

sim.avg.0 <- netsim(nets,
                    param = param_het_ept(epistats, ept.coverage = 0,
                                          ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                          keep.window = 0),
                    init = init_het_ct(),
                    control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.0, here("sims", "may2",  "current",  "simAvg0.rds"))
rm(sim.avg.0)

sim.avg <- netsim(nets,
                  param = param_het_ept(epistats, ept.coverage = 0.25,
                                        ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                        keep.window = 0),
                  init = init_het_ct(),
                  control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg, here("sims", "may2",  "current",  "simAvg.rds"))
rm(sim.avg)

sim.avg.ept <- netsim(nets,
                      param = param_het_ept(epistats, ept.coverage = 0.5,
                                            ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                            keep.window = 0),
                      init = init_het_ct(),
                      control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5,
                                               ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                               keep.window = 0))

saveRDS(sim.avg.ept, here("sims", "may2",  "current",  "simAvgEPT.rds"))
rm(sim.avg.ept)

sim.avg.ept.high <- netsim(nets,
                           param = param_het_ept(epistats, ept.coverage = 0.75,
                                                 ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                                 keep.window = 0),
                           init = init_het_ct(),
                           control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.ept.high, here("sims", "may2",  "current",  "simAvgEPThigh.rds"))
rm(sim.avg.ept.high)

sim.avg.ept.100 <- netsim(nets,
                              param = param_het_ept(epistats, ept.coverage = 1,
                                                    ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                                    keep.window = 0),
                              init = init_het_ct(),
                              control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.avg.ept.100, here("sims", "may2",  "current",  "simAvgEPT100.rds"))
rm(sim.avg.ept.100)
rm(nets)

#### No Concurrency ####
nets <- readRDS(here("fits", "calibrated_netests_noconc.rds"))

sim.nc.0 <- netsim(nets,
                   param = param_het_ept(epistats, ept.coverage = 0,
                                         ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                         keep.window = 0),
                   init = init_het_ct(),
                   control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.0, here("sims", "may2", "current", "simNc02.rds"))
rm(sim.nc.0)

sim.nc <- netsim(nets,
                 param = param_het_ept(epistats, ept.coverage = 0.25,
                                       ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                       keep.window = 0),
                 init = init_het_ct(),
                 control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc, here("sims", "may2", "current", "simNc2.rds"))
rm(sim.nc)

sim.nc.ept <- netsim(nets,
                     param = param_het_ept(epistats, ept.coverage = 0.5,
                                           ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                           keep.window = 0),
                     init = init_het_ct(),
                     control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept, here("sims", "may2",  "current",  "simNcEPT2.rds"))
rm(sim.nc.ept)

sim.nc.ept.high <- netsim(nets,
                          param = param_het_ept(epistats, ept.coverage = 0.75,
                                                ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                                keep.window = 0),
                          init = init_het_ct(),
                          control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept.high, here("sims", "may2",  "current",  "simNcEPThigh2.rds"))
rm(sim.nc.ept.high)

sim.nc.ept.100 <- netsim(nets,
                             param = param_het_ept(epistats, ept.coverage = 1,
                                                   ept.provision.pers.ong = 1, ept.provision.main.ong = 1,
                                                   keep.window = 0),
                             init = init_het_ct(),
                             control = control_het_ct(nsims = 5, nsteps=52*80, ncores=5))

saveRDS(sim.nc.ept.100, here("sims", "may2",  "current",  "simNcEPT1002.rds"))
rm(sim.nc.ept.100)



