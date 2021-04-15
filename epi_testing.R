#epidemic_testing 

library(EpiModelHIV)
library(here)

nets <- readRDS(here("fits", "calibrated_netests.rds"))
epistats <- readRDS(here("prep", "epistats.rds"))

sim.base2 <- netsim(nets,
               param = param_het_ept(epistats, keep.window = 52, ept.start = 1000, ept.coverage = 0.1,
                                     uct.tprob.mtf = 0.095*2, uct.tprob.ftm = 0.095*2,
                                     ct.test.agecat = c(0.363, 0.445, 0.358, 0.269, 0.161, 0.122)/3,
                                     ct.sympt.tx.int = 3, 
                                     condom.modifier = 0.2,
                                     acts.modifier = 4,
                                     male.ntx.dur = 3*52, female.ntx.dur = 3*52,
                                     uct.sympt.prob.male = 0.1, uct.sympt.prob.female = 0.1
               ),
               
               init = init_het_ct(prev.uct.male.young = 0.03*4, prev.uct.male.old = 0.0095*2,
                                  prev.uct.female.young = 0.0474*4, prev.uct.female.old = 0.0135*2,
                                  arrival.prev.male=0.03*2, arrival.prev.female = 0.0474*2),
               
               control = control_het_ct(nsims = 10, nsteps=52*50))

saveRDS(sim.base2, here("sims", "simBase2.rds"))


plot(test, y= c("i.prev.1519", "i.prev.2024", "i.prev.2529", "i.prev.3034", "i.prev.3539", "i.prev.4044"), legend = T, main="Prev by Agecat, baseline transmissibility")
plot(test, y=c("i.prev.feml", "i.prev.male"), legend=T)
plot(test, y=c("inc.uct.male", "inc.uct.female"), legend = T)

plot(test, y= c("recov.uct.female.1519", "recov.uct.female.2024", "recov.uct.female.2529"), 
                        legend = T, main="Recovery")

plot(test, y= c("ir.uct.female.1519", "ir.uct.female.2024", "ir.uct.female.2529",
                         "ir.uct.female.3034", "ir.uct.female.3539", "ir.uct.female.4044"), 
                         legend = T, main="Incidence")





