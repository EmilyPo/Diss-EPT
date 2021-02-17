# function writting in ddaf package but having trouble exporting
get_dist <- function(x, network, ndynamic){
  
  nsims <- x$control$nsims
  el <- get_mean_edgelist(x, network, byagesex=TRUE)
  pop <- get_agedist_ept(x)
  egodat <- egodata_meandegs_ept(ndynamic=ndynamic)
  
  
  ## grab correct target mean degs by age/sex of selected network
  if (ndynamic==2) {
    if (network==1){ #marcoh
      female <- egodat[[1]][,2]
      male <- egodat[[2]][,2]
    }
    
    if (network==2){ #casual
      female <- egodat[[3]][,2]
      male <- egodat[[4]][,2]
    }
  }
  
  if (ndynamic==3) {
    if (network==1){ # marriage
      female <- egodat[[1]][,2]
      male <- egodat[[2]][,2]
    }
    
    if (network==2){ # cohab
      female <- egodat[[3]][,2]
      male <- egodat[[4]][,2]
    }
    
    if (network==3){ #casual
      female <- egodat[[5]][,2]
      male <- egodat[[6]][,2]
    }
  }
  
  # edges counts by age for males and females
  females <- el[[1]]
  males <- el[[2]]
  
  # divide by mean pop counts
  females[,2] <- females[,2]/as.numeric(pop[[1]])
  males[,2] <- males[,2]/as.numeric(pop[[2]])
  
  # combine mean sim deg and mean egodata deg
  fdata <- as.data.frame(cbind(females, female))
  mdata <- as.data.frame(cbind(males, male))
  
  colnames(fdata) <- c("Ego Age", "Sim", "Egodata")
  colnames(mdata) <- c("Ego Age", "Sim", "Egodata")
  
  # plots
  fplot <- ggplot() +
    geom_bar(data = fdata, aes(x=`Ego Age`, y=Sim, fill="Sim"), stat = "identity",
             fill="#999999") +
    geom_smooth(fdata, mapping=aes(x=`Ego Age`, y=Egodata, color="blue")) +
    geom_point(fdata, mapping=aes(x=`Ego Age`, y=Egodata, color="blue"))  +
    xlab("Ego Age") +
    ylab("Mean Degree") +
    ggtitle("Mean Degree by Ego Age, Females") +
    scale_color_manual(name = "", labels = c("Egodata"), values = c("darkblue"))
  
  mplot <- ggplot() +
    geom_bar(data = mdata, aes(x=`Ego Age`, y=Sim, fill="Sim"), stat = "identity",
             fill="#999999") +
    geom_point(mdata, mapping=aes(x=`Ego Age`, y=Egodata, color="blue"))  +
    geom_smooth(mdata, mapping=aes(x=`Ego Age`, y=Egodata, color="blue"))  +
    xlab("Ego Age") +
    ylab("Mean Degree") +
    ggtitle("Mean Degree by Ego Age, Males") +
    scale_color_manual(name = "", labels = c("Egodata"), values = c("darkblue"))
  
  dat <- list(fplot, mplot, fdata, mdata)
  
  return(dat)
  
}
get_agedist_ept <- function(x){
  nsims <- x$control$nsims
  fages <- NULL
  mages <- NULL
  
  for (i in 1:nsims){
    f <- floor(x$attr[[i]]$ageF)
    fages <- c(fages, f)
    
    m <- floor(x$attr[[i]]$ageM)
    mages <- c(mages, m)
  }
  
  fages <- fages[fages<45 & fages>0]
  fmeans <- table(round(fages))/nsims
  mages <- mages[mages<45 & mages>0]
  mmeans <- table(round(mages))/nsims
  
  return(list(fmeans, mmeans))
}

get_mean_edgelist <- function(x, network, byagesex=TRUE){
  nsims <- x$control$nsims
  el <- x$el
  
  if (byagesex==TRUE){
    
    females <- NULL
    males <- NULL
    for (i in 1:nsims) {
      a <- floor(x$attr[[i]]$ageF[el[[i]][[network]]])
      b <- floor(x$attr[[i]]$ageM[el[[i]][[network]]])
      
      females <- c(females, a)
      males <- c(males,b)
      
      l <- list(table(females)/nsims, table(males)/nsims)
    }
    
    for (i in 1:2) {
      e <- l[[i]]
      d <- as.data.frame(e[-1])
      s <- setdiff(15:44, d[,1])
      
      if (length(s>0)){
        add <- data.frame(as.factor(s), rep(0,length(s)))
        colnames(add) <- colnames(d)
        d <- rbind(d, add)
        d[,1] <- as.numeric(as.character(d[,1]))
        d <- d[order(d[,1]),]
      }
      
      l[[i]] <- d
    }
    
    l[[1]]$females <- as.numeric(as.character(l[[1]]$females))
    l[[2]]$males <- as.numeric(as.character(l[[2]]$males))
    
    return(l)
    
  }
  else {
    el <- NULL
    for (i in 1:nsims) {
      el[[i]] <- x$el[[i]][[network]]
    }
    return(el)
  }
}
egodata_meandegs_ept <- function(path=NULL, ndynamic=NULL){
  
  if (ndynamic==3){
    
    if (is.null(path)) {
      fit <- readRDS("~/Dissertation/ExpeditedPartnerTherapy/old three networks/fits/fit.marriage.rds")
    } else {
      fit <- path
    }
    
    egos <- fit$egodata$egos
    egosvy <- as_survey(egos, weights="weight")
    
    # marriage
    mar_f <- egosvy %>% filter(male==0) %>% group_by(age) %>% summarise(mean=survey_mean(deg.mar))
    mar_m <- egosvy %>% filter(male==1) %>%group_by(age) %>% summarise(mean=survey_mean(deg.mar))
    
    # cohab
    cohab_f <- egosvy %>% filter(male==0) %>% group_by(age) %>% summarise(mean=survey_mean(deg.cohab))
    cohab_m <- egosvy %>% filter(male==1) %>%group_by(age) %>% summarise(mean=survey_mean(deg.cohab))
    
    # casual
    casual_f <- egosvy %>% filter(male==0) %>% group_by(age) %>% summarise(mean=survey_mean(deg.other))
    casual_m <- egosvy %>% filter(male==1) %>% group_by(age) %>% summarise(mean=survey_mean(deg.other))
    
    dat <- list(mar_f, mar_m, 
                cohab_f, cohab_m, 
                casual_f, casual_m)
  }
  
  if (ndynamic==2){
    if (is.null(path)) {
      fit <- readRDS("~/Dissertation/ExpeditedPartnerTherapy/fits/fit.marcoh.rds")
    } else {
      fit <- path
    }
    
    egos <- fit$egodata$egos
    egosvy <- as_survey(egos, weights="weight")
    
    # marriage
    mar_f <- egosvy %>% filter(male==0) %>% group_by(age) %>% summarise(mean=survey_mean(deg.marcoh))
    mar_m <- egosvy %>% filter(male==1) %>%group_by(age) %>% summarise(mean=survey_mean(deg.marcoh))
    
    # casual
    casual_f <- egosvy %>% filter(male==0) %>% group_by(age) %>% summarise(mean=survey_mean(deg.other))
    casual_m <- egosvy %>% filter(male==1) %>% group_by(age) %>% summarise(mean=survey_mean(deg.other))
    
    dat <- list(mar_f, mar_m, 
                casual_f, casual_m)
  }
  
  return(dat)
}
