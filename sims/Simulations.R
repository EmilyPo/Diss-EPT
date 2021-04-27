library(EpiModelHIV)
library(here)
epistats <- readRDS(here("prep", "epistats.rds"))

#### Concurrency By Sex ####
nets <- readRDS(here("fits", "calibrated_netests.rds"))

sim.base <- netsim(nets,
                   param = param_het_ept(epistats, keep.window = 52, ept.start = 1000, ept.coverage = 0.1,
                                         uct.tprob.mtf = 0.095*2, uct.tprob.ftm = 0.095*2.5,
                                         ct.test.agecat = c(0.363, 0.445, 0.358, 0.269, 0.161, 0.122)/3,
                                         ct.sympt.tx.int = 3, 
                                         condom.modifier = 0.2,
                                         acts.modifier = 4,
                                         female.modifier = 4,
                                         male.ntx.dur = 3*52, female.ntx.dur = 3*52,
                                         uct.sympt.prob.male = 0.1, uct.sympt.prob.female = 0.1
                   ),
                   
                   init = init_het_ct(prev.uct.male.young = 0.03*4, prev.uct.male.old = 0.0095*2,
                                      prev.uct.female.young = 0.0474*4, prev.uct.female.old = 0.0135*2,
                                      arrival.prev.male=0.03*2, arrival.prev.female = 0.0474*2),
                   
                   control = control_het_ct(nsims = 20, nsteps=52*50))


saveRDS(sim.base, here("sims", "fem.mod", "simBase.rds"))
rm(sim.base)

sim.base.ept <- netsim(nets,
                       param = param_het_ept(epistats, keep.window = 52, ept.start = 1000,
                                             uct.tprob.mtf = 0.095*2, uct.tprob.ftm = 0.095*2.5,
                                             ct.test.agecat = c(0.363, 0.445, 0.358, 0.269, 0.161, 0.122)/3,
                                             ct.sympt.tx.int = 3, 
                                             condom.modifier = 0.2,
                                             acts.modifier = 4,
                                             female.modifier = 4,
                                             male.ntx.dur = 3*52, female.ntx.dur = 3*52,
                                             uct.sympt.prob.male = 0.1, uct.sympt.prob.female = 0.1
                       ),
                       
                       init = init_het_ct(prev.uct.male.young = 0.03*4, prev.uct.male.old = 0.0095*2,
                                          prev.uct.female.young = 0.0474*4, prev.uct.female.old = 0.0135*2,
                                          arrival.prev.male=0.03*2, arrival.prev.female = 0.0474*2),
                       
                       control = control_het_ct(nsims = 20, nsteps=52*50))

saveRDS(sim.base.ept, here("sims", "fem.mod", "simBaseEPT.rds"))
rm(sim.base.ept)

rm(nets)

#### Average Concurrency ####
nets <- readRDS(here("fits", "calibrated_netests_avgconc.rds"))

sim.avg <- netsim(nets,
                  param = param_het_ept(epistats, keep.window = 52, ept.start = 1000, ept.coverage = 0.1,
                                        uct.tprob.mtf = 0.095*2, uct.tprob.ftm = 0.095*2.5,
                                        ct.test.agecat = c(0.363, 0.445, 0.358, 0.269, 0.161, 0.122)/3,
                                        ct.sympt.tx.int = 3, 
                                        condom.modifier = 0.2,
                                        acts.modifier = 4,
                                        female.modifier = 4,
                                        male.ntx.dur = 3*52, female.ntx.dur = 3*52,
                                        uct.sympt.prob.male = 0.1, uct.sympt.prob.female = 0.1
                  ),
                  
                  init = init_het_ct(prev.uct.male.young = 0.03*4, prev.uct.male.old = 0.0095*2,
                                     prev.uct.female.young = 0.0474*4, prev.uct.female.old = 0.0135*2,
                                     arrival.prev.male=0.03*2, arrival.prev.female = 0.0474*2),
                  
                  control = control_het_ct(nsims = 20, nsteps=52*50))


saveRDS(sim.avg, here("sims", "fem.mod", "simAvg.rds"))
rm(sim.avg)

sim.avg.ept <- netsim(nets,
                      param = param_het_ept(epistats, keep.window = 52, ept.start = 1000,
                                            uct.tprob.mtf = 0.095*2, uct.tprob.ftm = 0.095*2.5,
                                            ct.test.agecat = c(0.363, 0.445, 0.358, 0.269, 0.161, 0.122)/3,
                                            ct.sympt.tx.int = 3, 
                                            condom.modifier = 0.2,
                                            acts.modifier = 4,
                                            female.modifier = 4,
                                            male.ntx.dur = 3*52, female.ntx.dur = 3*52,
                                            uct.sympt.prob.male = 0.1, uct.sympt.prob.female = 0.1
                      ),
                      
                      init = init_het_ct(prev.uct.male.young = 0.03*4, prev.uct.male.old = 0.0095*2,
                                         prev.uct.female.young = 0.0474*4, prev.uct.female.old = 0.0135*2,
                                         arrival.prev.male=0.03*2, arrival.prev.female = 0.0474*2),
                      
                      control = control_het_ct(nsims = 20, nsteps=52*50))

saveRDS(sim.avg.ept, here("sims", "fem.mod", "simAvgEPT.rds"))
rm(sim.avg.ept)
rm(nets)
#### No Concurrency ####
nets <- readRDS(here("fits", "calibrated_netests_noconc.rds"))

#sim.nc <- netsim(nets,
#                  param = param_het_ept(epistats, keep.window = 52, ept.start = 1000, ept.coverage = 0.1,
#                                        uct.tprob.mtf = 0.095*2, uct.tprob.ftm = 0.095*2.5,
#                                        ct.test.agecat = c(0.363, 0.445, 0.358, 0.269, 0.161, 0.122)/3,
#                                        ct.sympt.tx.int = 3, 
#                                        condom.modifier = 0.2,
#                                        acts.modifier = 4,
#                                        female.modifier = 4,
#                                        male.ntx.dur = 3*52, female.ntx.dur = 3*52,
#                                        uct.sympt.prob.male = 0.1, uct.sympt.prob.female = 0.1
#                  ),
#                  
 #                 init = init_het_ct(prev.uct.male.young = 0.03*4, prev.uct.male.old = 0.0095*2,
#                                     prev.uct.female.young = 0.0474*4, prev.uct.female.old = 0.0135*2,
#                                     arrival.prev.male=0.03*2, arrival.prev.female = 0.0474*2),
#                  
#                  control = control_het_ct(nsims = 20, nsteps=52*50))


#saveRDS(sim.nc, here("sims", "fem.mod", "simNc.rds"))
#rm(sim.nc)

sim.nc.ept <- netsim(nets,
                      param = param_het_ept(epistats, keep.window = 52, ept.start = 1000,
                                            uct.tprob.mtf = 0.095*2, uct.tprob.ftm = 0.095*2.5,
                                            ct.test.agecat = c(0.363, 0.445, 0.358, 0.269, 0.161, 0.122)/3,
                                            ct.sympt.tx.int = 3, 
                                            condom.modifier = 0.2,
                                            acts.modifier = 4,
                                            female.modifier = 4,
                                            male.ntx.dur = 3*52, female.ntx.dur = 3*52,
                                            uct.sympt.prob.male = 0.1, uct.sympt.prob.female = 0.1
                      ),
                      
                      init = init_het_ct(prev.uct.male.young = 0.03*4, prev.uct.male.old = 0.0095*2,
                                         prev.uct.female.young = 0.0474*4, prev.uct.female.old = 0.0135*2,
                                         arrival.prev.male=0.03*2, arrival.prev.female = 0.0474*2),
                      
                      control = control_het_ct(nsims = 20, nsteps=52*50))

saveRDS(sim.nc.ept, here("sims", "fem.mod", "simNcEPT.rds"))
rm(sim.nc.ept)
rm(nets)