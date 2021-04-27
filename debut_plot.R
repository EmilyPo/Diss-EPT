# 1. some evidence that younger women more susceptible to ct than older women
# 2. we don't have network terms for:
#  - concurrency by age (currently have concurrency by sex)

library(tidyverse)
# diagnosing 
x <- test

d <- NULL
for (i in 1:x$control$nsims) {
  debut <- x$attr[[i]]$effectiveDebut
  age <- floor(x$attr[[i]]$age)
  
  dat <- data.frame(age=as.factor(age), debut=debut)
  dat <- dat %>% group_by(age) %>% summarize(mean=mean(debut==1))
  
  d[[i]] <- dat
}


data <- cbind(d[[1]], d[[2]][,2], d[[3]][,2], d[[4]][,2], d[[5]][,2])

data$meanOverall <- apply(data[,2:6], 1, mean)


nsfg_debut <-  c(0.1454219, 0.2878825, 0.4257733, 0.5830055, 0.7019795, 0.7984001, 0.8471828, 0.8459053, 0.8804078, 0.9146240, 0.9414704, 0.9453254, 
                 0.9467009, 0.9639282, 0.9756667, 0.9820388, 0.9733068, 0.9727181, 0.9741881, 0.9809794, 0.9918823, 0.9832634, 0.9825848, 0.9836374,
                 0.9908677, 0.9896829, 0.9890290, 0.9871198, 0.9744849, 0.9877073)

plot(x=15:44, y=data$meanOverall, pch=8, main = "Proportion Sexually Active, Simulation vs NSFG",
     ylab="Proportion Sexually Active", xlab="Age", col="dodgerblue2")
points(x= 15:44, y=nsfg_debut, pch=19)
legend("bottomright", legend=c("Simulation", "NSFG"), col = c("dodgerblue2", "Black"), pch=c(8, 19))


mar <- cel[[1]][which(cel[[1]][,"ptype"]==1),]
females <- which(mar[,"s1"]==0)
mar[females, 1:7] <- mar[females, c(2,1,3,5,4,7,6)] 
meanAgeDiff <- mean(sqrt(mar$age1)-sqrt(mar$age2))

cas <- cel[[1]][which(cel[[1]][,"ptype"]==2),]
females <- which(cas[,"s1"]==0)
cas[females, 1:7] <- cas[females, c(2,1,3,5,4,7,6)] 
meanAgeDiff <- mean(sqrt(cas$age1)-sqrt(cas$age2))



x <- test.prevbump

agecat <- x$attr[[1]]$agecat 
deg.other <- x$attr[[1]]$deg.other 
deg.marcoh <- x$attr[[1]]$deg.marcoh
deg.all <- deg.marcoh + deg.other

d <- as.data.frame(cbind(agecat,deg.other, deg.marcoh, deg.all))
d$conc <- as.integer(deg.other>1)
d %>% group_by(agecat) %>% summarize(mean=mean(deg.all))
d %>% group_by(agecat) %>% summarize(mean=mean(conc))




# ok but what we really want if reinfection following treatment by sex 
