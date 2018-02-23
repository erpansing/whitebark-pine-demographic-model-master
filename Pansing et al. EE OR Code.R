load("survival.Rda")

## _____________________________________________________________________________________________________
## Function to estimate proportion pilfered (95%CIs) and Odds ratios (95% CIs)
## _____________________________________________________________________________________________________
OR.frame.pilferage <-function(data, comparison, cache.characteristic){
  
  if(comparison == "all"){
    data <- data
    Comparison <- "All"
  }
  else if(comparison == "wc"){
    wc.area.subset <- subset(data, study.area =="wc")
    data <- wc.area.subset
    Comparison <- "WC"
  }
  else if(comparison == "tb"){
    tb.area.subset <- subset(data, study.area =="tb")
    data <- tb.area.subset
    Comparison <- "TB"
  }
  else if (comparison == "wc.sa"){
    wc.subset <- subset(data, study.area =="wc")
    wc.sa.subset <-subset(wc.subset, elevation.zone == "sa")
    data <- wc.sa.subset
    Comparison <- "WCSA"
  }
  else if (comparison == "wc.tl"){
    wc.subset <- subset(data, study.area =="wc")
    wc.tl.subset <-subset(wc.subset, elevation.zone == "tl")
    data <- wc.tl.subset
    Comparison <- "WCTL"
  }
  else if(comparison =="tb.sa"){
    wc.subset <- subset(data, study.area =="tb")
    wc.sa.subset <-subset(wc.subset, elevation.zone == "sa")
    data <- wc.sa.subset
    Comparison <- "TBSA"
  }
  else if (comparison == "tb.tl"){
    wc.subset <- subset(data, study.area =="tb")
    wc.tl.subset <-subset(wc.subset, elevation.zone == "tl")
    data <- wc.tl.subset
    Comparison <- "TBTL"
  }
  else{
    stop("invalid comparison")
  }
  
  
  if(cache.characteristic == "study.area"){
    characteristic <- data[,2]
    Characteristic <- "Study.Area"
    Ratio <- "TB/WC"
  } 
  else if(cache.characteristic == "elevation.zone"){
    characteristic <- data[,3]
    Characteristic <- "Elevation.Zone"
    Ratio <- "SA/TL"
  } 
  else if(cache.characteristic == "object"){
    characteristic <- data[,6]
    Characteristic <- "Object"
    Ratio <- c("N/R", "N/T","R/T")
  } 
  else if(cache.characteristic == "vegtype"){
    characteristic <- data[,9]
    Characteristic <- "Object"
    Ratio <- c("G/H","G/L","G/N","G/N","G/W","H/L","H/N","H/W","L/N","L/W",
               "N/W")
  } 
  else if(cache.characteristic == "cover"){
    characteristic <- data[,8]
    Characteristic <- "Object"
    Ratio <- c("1/2","1/3","1/4","2/3","2/4","3/4")
  }
  else if(cache.characteristic == "microsite"){
    characteristic <- data[,5]
    Characteristic <- "Microsite"
    Ratio <- c("LR/LT", "LR/O","LR/V", "LR/WR", "LR/WT",
               "LT/O","LT/V", "LT/WR", "LT/WT",
               "O/V", "O/WR",  "O/WT",
               "V/WR", "V/WT",
               "WR/WT")
  } 
  else if(cache.characteristic == "no.seeds.cached"){
    characteristic <- data[,4]
    Characteristic <- "No.Seeds.Cached"
    Ratio <- c("1/2", "1/3", "1/4", "1/5", "1/6", "1/7",
               "2/3", "2/4", "2/5", "2/6", "2/7",
               "3/4", "3/5", "3/6", "3/7",
               "4/5","4/6", "4/7",
               "5/6", "5/7",
               "6/7")
  } 
  else{
    stop("invalid cache characteristic")
  }
  
  n = as.matrix(tapply(data[,28], characteristic, length))
  k = as.matrix(tapply(data[,28], characteristic, sum, na.rm = T))
  pbar = k/n
  
  kn.frame <- cbind (n,k)
  print(kn.frame)
  
  
  
  conf.int <- apply(kn.frame, 1, function(x)  prop.test(x[2],x[1])$conf.int)
  conf.int <-t(conf.int)
  
  
  Proportion.Pilfered <- data.frame(LL = conf.int[,1], 
                                    Proportion.Pilfered = pbar,
                                    UL = conf.int[,2])
  print(Proportion.Pilfered)
  
  OR <- as.matrix(apply(combn(nrow(pbar), 2), 2, function(x) ((pbar[x[1],]/(1-pbar[x[1],])) / (pbar[x[2],]/(1-pbar[x[2],])))))
  
  # ## Calculate SE
  # 
  # germs <- tapply(survival[,28], survival$study.area, sum, na.rm = T)
  # no.germs <- tapply(survival[,28], survival$study.area, length)
  
  germs <- tapply(data[,28], characteristic, sum, na.rm = T)
  no.germs <- tapply(data[,28], characteristic, length)
  
  inverse.germs <- as.matrix(1/germs)
  inverse.no.germs <- as.matrix(1/no.germs)
  
  
  sum.germs <- as.matrix(apply(combn(nrow(inverse.germs), 2), 2, function(x) sum(inverse.germs[x[1],],inverse.germs[x[2],])))
  sum.no.germs <- as.matrix(apply(combn(nrow(inverse.no.germs), 2), 2, function(x) sum(inverse.no.germs[x[1],],inverse.no.germs[x[2],])))
  
  sum.total <- sum.germs + sum.no.germs
  SE <- sqrt(sum.total)
  
  alpha <- 0.05
  zalph <- qnorm(1-alpha/2)
  
  ln.or <- log(OR)
  or.ll <- exp(ln.or - zalph * SE)
  or.ul <- exp(ln.or + zalph * SE)
  
  OR.frame <- data.frame(Comparison = Comparison, Ratio = Ratio, Characteristic = Characteristic,
                         LL = or.ll,OR = OR, UL = or.ul)
  
  return(OR.frame)
  
  
}

## _____________________________________________________________________________________________________
## Function to estimate proportion germinated (95%CIs) and Odds ratios (95% CIs)
## _____________________________________________________________________________________________________

OR.frame <-function(data, comparison, cache.characteristic, germ.year){
  
  if(comparison == "all"){
    data <- data
    Comparison <- "All"
  }
  else if(comparison == "wc"){
    wc.area.subset <- subset(data, study.area =="wc")
    data <- wc.area.subset
    Comparison <- "WC"
  }
  else if(comparison == "tb"){
    tb.area.subset <- subset(data, study.area =="tb")
    data <- tb.area.subset
    Comparison <- "TB"
  }
  else if (comparison == "wc.sa"){
    wc.subset <- subset(data, study.area =="wc")
    wc.sa.subset <-subset(wc.subset, elevation.zone == "sa")
    data <- wc.sa.subset
    Comparison <- "WCSA"
  }
  else if (comparison == "wc.tl"){
    wc.subset <- subset(data, study.area =="wc")
    wc.tl.subset <-subset(wc.subset, elevation.zone == "tl")
    data <- wc.tl.subset
    Comparison <- "WCTL"
  }
  else if(comparison =="tb.sa"){
    wc.subset <- subset(data, study.area =="tb")
    wc.sa.subset <-subset(wc.subset, elevation.zone == "sa")
    data <- wc.sa.subset
    Comparison <- "TBSA"
  }
  else if (comparison == "tb.tl"){
    wc.subset <- subset(data, study.area =="tb")
    wc.tl.subset <-subset(wc.subset, elevation.zone == "tl")
    data <- wc.tl.subset
    Comparison <- "TBTL"
  }
  else{
    stop("invalid comparison")
  }
  if(cache.characteristic == "study.area"){
    characteristic <- data[,2]
    Characteristic <- "Study.Area"
    Ratio <- "TB/WC"
  } 
  else if(cache.characteristic == "elevation.zone"){
    characteristic <- data[,3]
    Characteristic <- "Elevation.Zone"
    Ratio <- "SA/TL"
  } 
  else if(cache.characteristic == "object"){
    characteristic <- data[,6]
    Characteristic <- "Object"
    Ratio <- c("N/R", "N/T","R/T")
  } 
  else if(cache.characteristic == "protection"){
    characteristic <- data[,7]
    Characteristic <- "Object"
    Ratio <- "P/U"
  } 
  else if(cache.characteristic == "cover"){
    characteristic <- data[,8]
    Characteristic <- "Object"
    Ratio <- c("1/2","1/3","1/4","2/3","2/4","3/4")
  }
  else if(cache.characteristic == "vegtype"){
    characteristic <- data[,9]
    Characteristic <- "Vegtype"
    Ratio <- c("G/H","G/L","G/N","G/W","H/L","H/N","H/W","L/N","L/W",
               "N/W")
  } 
  else if(cache.characteristic == "microsite"){
    characteristic <- data[,5]
    Characteristic <- "Microsite"
    Ratio <- c("LR/LT", "LR/O","LR/V", "LR/WR", "LR/WT",
               "LT/O","LT/V", "LT/WR", "LT/WT",
               "O/V", "O/WR",  "O/WT",
               "V/WR", "V/WT",
               "WR/WT")
  } 
  else if(cache.characteristic == "no.seeds.cached"){
    characteristic <- data[,4]
    Characteristic <- "No.Seeds.Cached"
    Ratio <- c("1/2", "1/3", "1/4", "1/5", "1/6", "1/7",
               "2/3", "2/4", "2/5", "2/6", "2/7",
               "3/4", "3/5", "3/6", "3/7",
               "4/5","4/6", "4/7",
               "5/6", "5/7",
               "6/7")
  } 
  else{
    stop("invalid cache characteristic")
  }
  
  
  if(germ.year == "germ.2013.y.n"){
    year <-data[,25]
    Year <- "2013"
  }
  else if(germ.year == "germ.2014.y.n"){
    year <-data[,26]
    Year <- "2014"
  }
  else if(germ.year == "germ.2013.2014"){
    year <-data[,27]
    Year <- "2013and2014"
  }  
  else if(germ.year == "survive.13.14"){
    year <- data[,20]
    Year <- "1stYearSurvival"
  }
  else if(germ.year == "living.2014.y.n"){
    year <- data[,23]
    Year <- "Living Seedlings 2014"
  }
  else{
    stop("invalid germ year")
  }
  
  
  n = as.matrix(tapply(year, characteristic, length))
  k = as.matrix(tapply(year, characteristic, sum, na.rm = T))
  pbar = k/n
  
  kn.frame <- cbind (n,k)
  print(kn.frame)
  
  
  
  conf.int <- apply(kn.frame, 1, function(x)  prop.test(x[2],x[1])$conf.int)
  conf.int <-t(conf.int)
  
  
  Proportion.Germinated <- data.frame(LL = conf.int[,1], 
                                      Proportion.Germinated = pbar,
                                      UL = conf.int[,2])
  print(Proportion.Germinated)
  
  OR <- as.matrix(apply(combn(nrow(pbar), 2), 2, function(x) ((pbar[x[1],]/(1-pbar[x[1],])) / (pbar[x[2],]/(1-pbar[x[2],])))))
  
  ## Calculate SE
  germs <- tapply(year, characteristic, sum, na.rm = T)
  no.germs <- tapply(year, characteristic, length)
  
  inverse.germs <- as.matrix(1/germs)
  inverse.no.germs <- as.matrix(1/no.germs)
  
  
  sum.germs <- as.matrix(apply(combn(nrow(inverse.germs), 2), 2, function(x) sum(inverse.germs[x[1],],inverse.germs[x[2],])))
  sum.no.germs <- as.matrix(apply(combn(nrow(inverse.no.germs), 2), 2, function(x) sum(inverse.no.germs[x[1],],inverse.no.germs[x[2],])))
  
  sum.total <- sum.germs + sum.no.germs
  SE <- sqrt(sum.total)
  
  alpha <- 0.05
  zalph <- qnorm(1-alpha/2)
  
  ln.or <- log(OR) 
  or.ll <- exp(ln.or - zalph * SE)
  or.ul <- exp(ln.or + zalph * SE)
  
  OR.frame <- data.frame(Year = Year, Comparison = Comparison, Ratio = Ratio, Characteristic = Characteristic,
                         LL = or.ll,OR = OR, UL = or.ul)
  
  return(OR.frame)
  
}
  
  
  ## _____________________________________________________________________________________________________
  ## Function to estimate proportion survived (95%CIs) and Odds ratios (95% CIs)
  ## _____________________________________________________________________________________________________
  
  OR.frame.survival <-function(data, comparison, cache.characteristic, germ.year){
    data <- subset(data, germ.2013.y.n==1)
    if(comparison == "all"){
      data <- data
      Comparison <- "All"
    }
    else if(comparison == "wc"){
      wc.area.subset <- subset(data, study.area =="wc")
      data <- wc.area.subset
      Comparison <- "WC"
    }
    else if(comparison == "tb"){
      tb.area.subset <- subset(data, study.area =="tb")
      data <- tb.area.subset
      Comparison <- "TB"
    }
    else if (comparison == "wc.sa"){
      wc.subset <- subset(data, study.area =="wc")
      wc.sa.subset <-subset(wc.subset, elevation.zone == "sa")
      data <- wc.sa.subset
      Comparison <- "WCSA"
    }
    else if (comparison == "wc.tl"){
      wc.subset <- subset(data, study.area =="wc")
      wc.tl.subset <-subset(wc.subset, elevation.zone == "tl")
      data <- wc.tl.subset
      Comparison <- "WCTL"
    }
    else if (comparison == "tb.sa"){
      wc.subset <- subset(data, study.area =="tb")
      wc.tl.subset <-subset(wc.subset, elevation.zone == "sa")
      data <- wc.tl.subset
      Comparison <- "TBSA"
    }
    else if (comparison == "tb.tl"){
      tb.subset <- subset(data, study.area =="tb")
      tb.tl.subset <-subset(tb.subset, elevation.zone == "tl")
      data <- tb.tl.subset
      Comparison <- "TBTL"
    }
    else{
      stop("invalid comparison")
    }
    
    
    if(cache.characteristic == "study.area"){
      characteristic <- data[,2]
      Characteristic <- "Study.Area"
      Ratio <- "TB/WC"
    } 
    else if(cache.characteristic == "elevation.zone"){
      characteristic <- data[,3]
      Characteristic <- "Elevation.Zone"
      Ratio <- "SA/TL"
    } 
    else if(cache.characteristic == "object"){
      characteristic <- data[,6]
      Characteristic <- "Object"
      Ratio <- c("N/R", "N/T","R/T")
    } 
    else if(cache.characteristic == "cover"){
      characteristic <- data[,8]
      Characteristic <- "Object"
      Ratio <- c("1/2","1/3","1/4","2/3","2/4","3/4")
    }
    else if(cache.characteristic == "microsite"){
      characteristic <- data[,5]
      Characteristic <- "Microsite"
      Ratio <- c("LR/LT", "LR/O","LR/V", "LR/WR", "LR/WT",
                 "LT/O","LT/V", "LT/WR", "LT/WT",
                 "O/V", "O/WR",  "O/WT",
                 "V/WR", "V/WT",
                 "WR/WT")
      
    } 
    else if(cache.characteristic == "protection"){
      characteristic <- data[,7]
      Characteristic <- "Object"
      Ratio <- "P/U"
    } 
    else if(cache.characteristic == "no.seeds.cached"){
      characteristic <- data[,4]
      Characteristic <- "No.Seeds.Cached"
      Ratio <- c("1/2", "1/3", "1/4", "1/5", "1/6", "1/7",
                 "2/3", "2/4", "2/5", "2/6", "2/7",
                 "3/4", "3/5", "3/6", "3/7",
                 "4/5","4/6", "4/7",
                 "5/6", "5/7",
                 "6/7")
    } 
    else{
      stop("invalid cache characteristic")
    }
    
    
    if(germ.year == "survive.13.14"){
      year <- data[,20]
      Year <- "1stYearSurvival"
    }
    
    
    else{
      stop("invalid germ year")
    }
    
    
    #   mean.survival = as.matrix(tapply(year, characteristic, mean, na.rm = T))
    #   k = as.matrix(tapply(year, characteristic, sum, na.rm = T))
    #   n = k/mean.survival
    
    
    n = as.matrix(tapply(year, characteristic, length))
    k = as.matrix(tapply(year, characteristic, sum, na.rm = T))
    pbar = k/n
    
    
    kn.frame <- cbind (n,k)
    print(kn.frame)
    
    
    conf.int <- apply(kn.frame, 1, function(x)  prop.test(x[2],x[1])$conf.int)
    conf.int <-t(conf.int)
    
    
    Proportion.Survived <- data.frame(LL = conf.int[,1], 
                                      Proportion.Survived = pbar,
                                      UL = conf.int[,2])
    print(Proportion.Survived)
    
    OR <- as.matrix(apply(combn(nrow(pbar), 2), 2, function(x) ((pbar[x[1],]/(1-pbar[x[1],])) / (pbar[x[2],]/(1-pbar[x[2],])))))
    
    ## Calculate SE
    germs <- tapply(year, characteristic, sum, na.rm = T)
    no.germs <- tapply(year, characteristic, length)
    
    inverse.germs <- as.matrix(1/germs)
    inverse.no.germs <- as.matrix(1/no.germs)
    
    
    sum.germs <- as.matrix(apply(combn(nrow(inverse.germs), 2), 2, function(x) sum(inverse.germs[x[1],],inverse.germs[x[2],])))
    sum.no.germs <- as.matrix(apply(combn(nrow(inverse.no.germs), 2), 2, function(x) sum(inverse.no.germs[x[1],],inverse.no.germs[x[2],])))
    
    sum.total <- sum.germs + sum.no.germs
    SE <- sqrt(sum.total)
    
    alpha <- 0.05
    zalph <- qnorm(1-alpha/2)
    
    ln.or <- log(OR)
    or.ll <- exp(ln.or - zalph * SE)
    or.ul <- exp(ln.or + zalph * SE)
    
    OR.frame <- data.frame(Year = Year, Comparison = Comparison, Ratio = Ratio, Characteristic = Characteristic,
                           LL = or.ll,OR = OR, UL = or.ul)
    
    return(OR.frame)
    
    
  }
  

  ##-------------------------------------------------------------------------------
  ##                        2013 Cohort
  ##-------------------------------------------------------------------------------
  
  ##-------------------------------------------------------------------------------
  ##                        Study Area
  ##-------------------------------------------------------------------------------
  
  load("2014_germ_survival.Rda")
  survival$pilf.total <- ifelse(survival$no.seeds.cached == survival$missing.2013, 1,0)
  
  germ.subset <- subset(survival, pilf.total==0)
  ##-----------------------------------
  ## 2013 Pilferage 
  ##-----------------------------------
  ##TB
  pilf.13.area <- OR.frame.pilferage(data = survival, comparison = "all",  cache.characteristic = "study.area")
  
  ##-----------------------------------
  ## 2013 Germination 
  ##-----------------------------------
  
  germ.13.area <-OR.frame(data = germ.subset, comparison = "all",  cache.characteristic = "study.area", germ.year = "germ.2013.y.n")
  
  ##-----------------------------------
  ## 2013 Survival 
  ##-----------------------------------
  
  survival.area <- OR.frame.survival(data = survival, comparison = "all",  cache.characteristic = "study.area", germ.year = "survive.13.14")
  
  
  ## Data Frame for Plotting
  
  study.area.OR.frame <- rbind(pilf.13.area, germ.13.area[,-1], survival.area[,-1]) %>%
    mutate(Stage = c("Pilferage","Germination","Survival")) %>%
    mutate(Stage = factor(Stage, levels = c("Pilferage","Germination", "Survival")))
  
  
  write.csv(study.area.OR.frame,"areaframe.csv")
  
  ##-------------------------------------------------------------------------------
  ##                        Elevation Zone
  ##-------------------------------------------------------------------------------
  
  ##-----------------------------------
  ## 2013 Pilferage 
  ##-----------------------------------
  ##TB
  tb.zone.pilf <- OR.frame.pilferage(data = survival, comparison = "tb",  cache.characteristic = "elevation.zone")
  
  ##WC
  wc.zone.pilf <- OR.frame.pilferage(data = survival, comparison = "wc",  cache.characteristic = "elevation.zone")
  
  ##-----------------------------------
  ## 2013 Germination 
  ##-----------------------------------
  ## TB
  tb.germ.13.zone <-OR.frame(data = germ.subset, comparison = "tb",  cache.characteristic = "elevation.zone", germ.year = "germ.2013.y.n")
  
  ## WC
  wc.germ.13.zone <-OR.frame(data = germ.subset, comparison = "wc",  cache.characteristic = "elevation.zone", germ.year = "germ.2013.y.n")
  
  ##-----------------------------------
  ## 2013 Survival 
  ##-----------------------------------
  ##TB
  tb.survival.zone <- OR.frame.survival(data = germ.subset, comparison = "tb",  cache.characteristic = "elevation.zone", germ.year = "survive.13.14")
  
  ##WC
  wc.survival.zone <- OR.frame.survival(data = survival, comparison = "wc",  cache.characteristic = "elevation.zone", germ.year = "survive.13.14")
  
  
  ## Data Frame for Plotting
  
  elevation.zone.OR.frame <- rbind(tb.zone.pilf, wc.zone.pilf, tb.germ.13.zone[,-1],
                                   wc.germ.13.zone[,-1], tb.survival.zone[,-1], wc.survival.zone[,-1])
  elevation.zone.OR.frame$Stage <- c(rep("Pilferage",2),rep("Germination",2),rep("Survival",2))
  
  
  LL<- 1/elevation.zone.OR.frame$UL
  UL<- 1/elevation.zone.OR.frame$LL
  OR<- 1/elevation.zone.OR.frame$OR
  
  inverse.elevation.zone.OR.frame <- data.frame("LL" = LL,"OR" = OR, "UL" = UL)
  elevation.zone.frame.labels<- elevation.zone.OR.frame[,-c(4:6)]
  
  elevation.zone.frame <- cbind(elevation.zone.frame.labels, inverse.elevation.zone.OR.frame) %>% 
    mutate(Ratio = c(rep("TL/SA",6))) %>%
    mutate(Stage = factor(Stage, levels = c("Pilferage","Germination", "Survival")))
  
  
  write.csv(elevation.zone.frame,"zoneframe.csv")
  
  ##-------------------------------------------------------------------------------
  ##                        Object
  ##-------------------------------------------------------------------------------
  
  ##-----------------------------------
  ## 2013 Pilferage 
  ##-----------------------------------
  
  ## TBSA
  tbsa.object.pilf <- OR.frame.pilferage(data = survival, comparison = "tb.sa",  cache.characteristic = "object")
  
  ## TBTL
  tbtl.object.pilf <- OR.frame.pilferage(data = survival, comparison = "tb.tl",  cache.characteristic = "object")
  
  ## WCSA
  wcsa.object.pilf <- OR.frame.pilferage(data = survival, comparison = "wc.sa",  cache.characteristic = "object")
  
  ## WCTL
  wctl.object.pilf <- OR.frame.pilferage(data = survival, comparison = "wc.tl",  cache.characteristic = "object")
  
  ##-----------------------------------
  ## 2013 Germination 
  ##-----------------------------------
  
  ## TBSA
  tbsa.germ.13.object <-OR.frame(data = germ.subset, comparison = "tb.sa",  cache.characteristic = "object", germ.year = "germ.2013.y.n")
  
  ## TBTL
  tbtl.germ.13.object <-OR.frame(data = germ.subset, comparison = "tb.tl",  cache.characteristic = "object", germ.year = "germ.2013.y.n")
  
  ## WCSA
  wcsa.germ.13.object <-OR.frame(data = germ.subset, comparison = "wc.sa",  cache.characteristic = "object", germ.year = "germ.2013.y.n")
  
  ## WCTL
  wctl.germ.13.object <-OR.frame(data = germ.subset, comparison = "wc.tl",  cache.characteristic = "object", germ.year = "germ.2013.y.n")
  
  ##-----------------------------------
  ## 2013 Survival 
  ##-----------------------------------
  
  ## TBSA
  tbsa.survival.object <- OR.frame.survival(data = survival, comparison = "tb.sa",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  ## TBTL
  tbtl.survival.object <- OR.frame.survival(data = survival, comparison = "tb.tl",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  ## WCSA
  wcsa.survival.object <- OR.frame.survival(data = survival, comparison = "wc.sa",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  ## WCTL
  wctl.survival.object <- OR.frame.survival(data = survival, comparison = "wc.tl",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  ##WC
  wc.survival.object <- OR.frame.survival(data = survival, comparison = "wc",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  
  ## Data Frame for Plotting
  
  object.OR.frame <- rbind(tbsa.object.pilf,tbtl.object.pilf,
                           wcsa.object.pilf,wctl.object.pilf,
                           tbsa.germ.13.object[,-1],tbtl.germ.13.object[,-1],
                           wcsa.germ.13.object[,-1],wctl.germ.13.object[,-1],
                           tbsa.survival.object[,-1],tbtl.survival.object[,-1],
                           wcsa.survival.object[,-1], wctl.survival.object[,-1])
  
  array1 <-array(data = NA,dim = c(36,3))
  LL <- NULL
  UL <- NULL
  OR <- NULL
  
  
  for(i in 1:nrow(object.OR.frame)){
    if(object.OR.frame[i,2] == "N/R"){
      LL = 1/object.OR.frame[i,6] 
      array1[i,1]<-LL
      UL = 1/object.OR.frame[i,4] 
      array1[i,3] <- UL
      OR = 1/object.OR.frame[i,5]
      array1[i,2] <- OR
    }
    else if(object.OR.frame[i,2] == "N/T"){
      LL = 1/object.OR.frame[i,6] 
      array1[i,1]<-LL
      UL = 1/object.OR.frame[i,4] 
      array1[i,3] <- UL
      OR = 1/object.OR.frame[i,5]
      array1[i,2] <- OR
    }
    else if(object.OR.frame[i,2] == "R/T"){
      LL = object.OR.frame[i,4] 
      array1[i,1]<-LL
      UL = object.OR.frame[i,6] 
      array1[i,3] <- UL
      OR = object.OR.frame[i,5]
      array1[i,2] <- OR
    }
  }
  
  object.OR.frame[,4:6]<-array1
  object.OR.frame$Ratio <-c(rep(c("R/N","T/N","R/T"),12))
  object.OR.frame$Stage <- c(rep("Pilferage",12),rep("Germination",12),rep("Survival",12))
  object.OR.frame$Stage <- factor(object.OR.frame$Stage,
                                  levels=c("Survival","Germination","Pilferage"))
  
  tb.object.frame <- subset(object.OR.frame, Comparison == "TBSA" |Comparison =="TBTL")
  tb.object.frame$Ratio <- factor(tb.object.frame$Ratio,
                                  levels=c("R/T","R/N","T/N"))
  tb.object.frame$Comparison <- factor(tb.object.frame$Comparison,
                                       levels=c("TBSA","TBTL"))
  
  wc.object.frame <- subset(object.OR.frame, Comparison == "WCSA" |Comparison =="WCTL")
  wc.object.frame$Ratio <- factor(wc.object.frame$Ratio,
                                  levels=c("R/T","R/N","T/N"))
  wc.object.frame$Comparison <- factor(wc.object.frame$Comparison,
                                       levels=c("WCSA","WCTL"))
  
  
  ##-------------------------------------------------------------------------------
  ##                        Cover
  ##-------------------------------------------------------------------------------
  
  ##-----------------------------------
  ## 2013 Pilferage 
  ##-----------------------------------
  
  cover.tbsa.pilf <- OR.frame.pilferage(data = survival, comparison = "tb.sa",
                     cache.characteristic = "cover")
  
  cover.tbtl.pilf <- OR.frame.pilferage(data = survival, comparison = "tb.tl",
                                   cache.characteristic = "cover")
  
  
  cover.wcsa.pilf <- OR.frame.pilferage(data = survival, comparison = "wc.sa",
                     cache.characteristic = "cover")
  
  cover.wctl.pilf <- OR.frame.pilferage(data = survival, comparison = "wc.tl",
                                   cache.characteristic = "cover")
  
  
  cover.OR.frame.pilf <- rbind(cover.tbtl.pilf, cover.tbsa.pilf, 
                          cover.wctl.pilf, cover.wcsa.pilf)
  
  ##-----------------------------------
  ## 2013 Germination 
  ##-----------------------------------
  
  cover.tbsa.germ <- OR.frame(data = germ.subset, 
                              comparison = "tb.sa",
                              cache.characteristic = "cover",
                              germ.year = "germ.2013.y.n")
  
  cover.tbtl.germ <- OR.frame(data = germ.subset, comparison = "tb.tl",
                                        cache.characteristic = "cover",
                              germ.year = "germ.2013.y.n")
  
  
  cover.wcsa.germ <- OR.frame(data = germ.subset, comparison = "wc.sa",
                                        cache.characteristic = "cover",
                              germ.year = "germ.2013.y.n")
  
  cover.wctl.germ <- OR.frame(data = germ.subset, comparison = "wc.tl",
                                        cache.characteristic = "cover",
                                        germ.year = "germ.2013.y.n")
  
  
  cover.OR.frame.germ <- rbind(cover.tbtl.germ, cover.tbsa.germ, 
                          cover.wctl.germ, cover.wcsa.germ)
  
  ##-----------------------------------
  ## 2013 Surviva 
  ##-----------------------------------
  
  cover.tbsa.survival <- OR.frame.survival(data = germ.subset, 
                              comparison = "tb.sa",
                              cache.characteristic = "cover",
                              germ.year = "germ.2013.y.n")
  
  cover.tbtl.survival <- OR.frame.survival(data = germ.subset, comparison = "tb.tl",
                              cache.characteristic = "cover",
                              germ.year = "germ.2013.y.n")
  
  
  cover.wcsa.survival <- OR.frame.survival(data = germ.subset, comparison = "wc.sa",
                              cache.characteristic = "cover",
                              germ.year = "germ.2013.y.n")
  
  cover.wctl.germ <- OR.frame.survival(data = survival, comparison = "wc.tl",
                              cache.characteristic = "cover",
                              germ.year = "germ.2013.y.n")
  
  
  cover.OR.frame.survival <- rbind(cover.tbtl.germ, cover.tbsa.germ, 
                               cover.wctl.germ, cover.wcsa.germ) 
  
  ##-----------------------------------
  ## 2013 Germination 
  ##-----------------------------------
  
  
  
  ##-----------------------------------
  ## 2012 - 2014 Survival 
  ##-----------------------------------
  
  ## Area
  survival.area.12.14 <- OR.frame(data = survival, comparison = "all",  cache.characteristic = "study.area", germ.year = "survive.13.14")
  
  
  ## TBSA
  tbsa.survival.object <- OR.frame(data = survival, comparison = "tb.sa",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  ## TBTL
  tbtl.survival.object <- OR.frame(data = survival, comparison = "tb.tl",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  ## WCSA
  wcsa.survival.object <- OR.frame(data = survival, comparison = "wc.sa",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  ## WCTL
  wctl.survival.object <- OR.frame(data = survival, comparison = "wc.tl",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  
  
  #-------------------------------------------------------------------------------
  ##                       2014  Study Area
  ##-------------------------------------------------------------------------------
  ##-------------------------------------------------------------------------------
  ##                        Study Area
  ##-------------------------------------------------------------------------------
  
  
  ##-----------------------------------
  ## 2014 Germination 
  ##-----------------------------------
  
  germ.14.area <-OR.frame(data = germ.subset, comparison = "all",  cache.characteristic = "study.area", germ.year = "germ.2014.y.n")
  
  ##-------------------------------------------------------------------------------
  ##                        Elevation Zone
  ##-------------------------------------------------------------------------------
  
  ##-----------------------------------
  ## 2014 Germination 
  ##-----------------------------------
  ## TB
  tb.germ.14.zone <-OR.frame(data = germ.subset, comparison = "tb",  cache.characteristic = "elevation.zone", germ.year = "germ.2014.y.n")
  
  ## WC
  wc.germ.14.zone <-OR.frame(data = germ.subset, comparison = "wc",  cache.characteristic = "elevation.zone", germ.year = "germ.2014.y.n")
  
  ##-------------------------------------------------------------------------------
  ##                        Object
  ##-------------------------------------------------------------------------------
  
  ##-----------------------------------
  ## 2014 Germination 
  ##-----------------------------------
  
  ## TBSA
  tbsa.germ.14.object <-OR.frame(data = germ.subset, comparison = "tb.sa",  cache.characteristic = "object", germ.year = "germ.2014.y.n")
  
  ## TBTL
  tbtl.germ.14.object <-OR.frame(data = germ.subset, comparison = "tb.tl",  cache.characteristic = "object", germ.year = "germ.2014.y.n")
  
  ## WCSA
  wcsa.germ.14.object <-OR.frame(data = germ.subset, comparison = "wc.sa",  cache.characteristic = "object", germ.year = "germ.2014.y.n")
  
  ## WCTL
  wctl.germ.14.object <-OR.frame(data = germ.subset, comparison = "wc.tl",  cache.characteristic = "object", germ.year = "germ.2014.y.n")
  ## WC
  wc.germ.14.object <-OR.frame(data = cohort.2014, comparison = "wc",  cache.characteristic = "object", germ.year = "germ.2014.y.n")
  
  
  
  ## Data Frame for Plotting
  
  cohort.2014.OR.frame <- rbind(germ.14.area,tb.germ.14.zone,wc.germ.14.zone,tbsa.germ.14.object,
                                tbtl.germ.14.object, wcsa.germ.14.object, wctl.germ.14.object)
  
  array1 <-array(0,0,dim = c(36,3))
  LL <- NULL
  UL <- NULL
  OR <- NULL
  
  
  for(i in 1:nrow(object.OR.frame)){
    if(object.OR.frame[i,2] == "N/R"){
      LL = 1/object.OR.frame[i,6] 
      array1[i,1]<-LL
      UL = 1/object.OR.frame[i,4] 
      array1[i,3] <- UL
      OR = 1/object.OR.frame[i,5]
      array1[i,2] <- OR
    }
    else if(object.OR.frame[i,2] == "N/T"){
      LL = 1/object.OR.frame[i,6] 
      array1[i,1]<-LL
      UL = 1/object.OR.frame[i,4] 
      array1[i,3] <- UL
      OR = 1/object.OR.frame[i,5]
      array1[i,2] <- OR
    }
    else if(object.OR.frame[i,2] == "R/T"){
      LL = object.OR.frame[i,4] 
      array1[i,1]<-LL
      UL = object.OR.frame[i,6] 
      array1[i,3] <- UL
      OR = object.OR.frame[i,5]
      array1[i,2] <- OR
    }
  }
  
  object.OR.frame[,4:6]<-array1
  object.OR.frame$Ratio <-c(rep(c("R/N","T/N","R/T"),12))
  object.OR.frame$Stage <- c(rep("Pilferage",12),rep("Germination",12),rep("Survival",12))
  object.OR.frame$Stage <- factor(object.OR.frame$Stage,
                                  levels=c("Survival","Germination","Pilferage"))
  
  tb.object.frame <- subset(object.OR.frame, Comparison == "TBSA" |Comparison =="TBTL")
  tb.object.frame$Ratio <- factor(tb.object.frame$Ratio,
                                  levels=c("R/T","R/N","T/N"))
  tb.object.frame$Comparison <- factor(tb.object.frame$Comparison,
                                       levels=c("TBSA","TBTL"))
  
  wc.object.frame <- subset(object.OR.frame, Comparison == "WCSA" |Comparison =="WCTL")
  wc.object.frame$Ratio <- factor(wc.object.frame$Ratio,
                                  levels=c("R/T","R/N","T/N"))
  wc.object.frame$Comparison <- factor(wc.object.frame$Comparison,
                                       levels=c("WCSA","WCTL"))
  
  
  
  #-------------------------------------------------------------------------------
  ##                          2013/14
  ##-------------------------------------------------------------------------------
  ##-------------------------------------------------------------------------------
  ##                        Study Area
  ##-------------------------------------------------------------------------------
  
  
  ##-----------------------------------
  ## 2013/14 Germination 
  ##-----------------------------------
  
  germ.13.14.area <-OR.frame(data = germ.subset, comparison = "all",  cache.characteristic = "study.area", germ.year = "germ.2013.2014")
  
  ##-------------------------------------------------------------------------------
  ##                        Elevation Zone
  ##-------------------------------------------------------------------------------
  
  ##-----------------------------------
  ## 2013/14 Germination 
  ##-----------------------------------
  ## TB
  tb.germ.13.14.zone <-OR.frame(data = germ.subset, comparison = "tb",  cache.characteristic = "elevation.zone", germ.year = "germ.2013.2014")
  
  ## WC
  wc.germ.13.14.zone <-OR.frame(data = germ.subset, comparison = "wc",  cache.characteristic = "elevation.zone", germ.year = "germ.2013.2014")
  
  ##-------------------------------------------------------------------------------
  ##                        Object
  ##-------------------------------------------------------------------------------
  
  ##-----------------------------------
  ## 2013/2014 Germination 
  ##-----------------------------------
  
  ## TBSA
  tbsa.germ.13.14.object <-OR.frame(data = germ.subset, comparison = "tb.sa",  cache.characteristic = "object", germ.year = "germ.2013.2014")
  
  ## TBTL
  tbtl.germ.13.14.object <-OR.frame(data = germ.subset, comparison = "tb.tl",  cache.characteristic = "object", germ.year = "germ.2013.2014")
  
  ## WCSA
  wcsa.germ.13.14.object <-OR.frame(data = germ.subset, comparison = "wc.sa",  cache.characteristic = "object", germ.year = "germ.2013.2014")
  
  ## WCTL
  wctl.germ.13.14.object <-OR.frame(data = germ.subset, comparison = "wc.tl",  cache.characteristic = "object", germ.year = "germ.2013.2014")
  ## WC
  wc.germ.14.object <-OR.frame(data = cohort.2014, comparison = "wc",  cache.characteristic = "object", germ.year = "germ.2014.y.n")
  
  
  
  ## Data Frame for Plotting
  
  cohort.2014.OR.frame <- rbind(germ.14.area,tb.germ.14.zone,wc.germ.14.zone,tbsa.germ.14.object,
                                tbtl.germ.14.object, wcsa.germ.14.object, wctl.germ.14.object)
  
  array1 <-array(0,0,dim = c(36,3))
  LL <- NULL
  UL <- NULL
  OR <- NULL
  
  
  for(i in 1:nrow(object.OR.frame)){
    if(object.OR.frame[i,2] == "N/R"){
      LL = 1/object.OR.frame[i,6] 
      array1[i,1]<-LL
      UL = 1/object.OR.frame[i,4] 
      array1[i,3] <- UL
      OR = 1/object.OR.frame[i,5]
      array1[i,2] <- OR
    }
    else if(object.OR.frame[i,2] == "N/T"){
      LL = 1/object.OR.frame[i,6] 
      array1[i,1]<-LL
      UL = 1/object.OR.frame[i,4] 
      array1[i,3] <- UL
      OR = 1/object.OR.frame[i,5]
      array1[i,2] <- OR
    }
    else if(object.OR.frame[i,2] == "R/T"){
      LL = object.OR.frame[i,4] 
      array1[i,1]<-LL
      UL = object.OR.frame[i,6] 
      array1[i,3] <- UL
      OR = object.OR.frame[i,5]
      array1[i,2] <- OR
    }
  }
  
  object.OR.frame[,4:6]<-array1
  object.OR.frame$Ratio <-c(rep(c("R/N","T/N","R/T"),12))
  object.OR.frame$Stage <- c(rep("Pilferage",12),rep("Germination",12),rep("Survival",12))
  object.OR.frame$Stage <- factor(object.OR.frame$Stage,
                                  levels=c("Survival","Germination","Pilferage"))
  
  tb.object.frame <- subset(object.OR.frame, Comparison == "TBSA" |Comparison =="TBTL")
  tb.object.frame$Ratio <- factor(tb.object.frame$Ratio,
                                  levels=c("R/T","R/N","T/N"))
  tb.object.frame$Comparison <- factor(tb.object.frame$Comparison,
                                       levels=c("TBSA","TBTL"))
  
  wc.object.frame <- subset(object.OR.frame, Comparison == "WCSA" |Comparison =="WCTL")
  wc.object.frame$Ratio <- factor(wc.object.frame$Ratio,
                                  levels=c("R/T","R/N","T/N"))
  wc.object.frame$Comparison <- factor(wc.object.frame$Comparison,
                                       levels=c("WCSA","WCTL"))
  
  
  
  
  ##-------------------------------------------------------------------------------
  ##                       Estimates for Nutcracker Dispersal Efficiency
  ##-------------------------------------------------------------------------------
  
  ## Number of seeds remaining after nutcrackers retrieve ~45% of seeds
  35000 - (35000 * 0.45)
  98000 - (98000 * 0.45)
  
  
  ## Low end = 19250 seeds remaining, high end = 44100
  
  
  ## Calculate number of 3-seed caches with 19250 seeds
  
  19250/3
  # ~ 6416.667 caches remaining after recovery
  
  ## Caches germinated after two years at each study area
  germ.13 <- OR.frame(data = survival, comparison = "all",  cache.characteristic = "study.area", germ.year = "living.2014.y.n")
  
  ## TB 
  .3798 * 19250
  #2437 caches germinated
  
  ## LL
  .3302 * 6417
  #2119
  
  ## UL
  .4319 * 6417
  # 2772
  
  ## WC
  0.2564 *6417
  #1645 caches germinated
  ## LL
  .2122 * 6417
  ##1362
  
  ## UL
  .3060 * 6417
  # 1964
  living.14.area <-OR.frame(data = survival, comparison = "all",  cache.characteristic = "study.area", germ.year = "living.2014.y.n")
  
  ## TB
  0.3579 * 6417
  #2297 caches living after two years
  ## LL
  0.3092 * 6417
  ## 1984
  ## UL
  0.4097 * 6417
  #2629
  
  ## WC
  0.071225 * 6417
  # 457 caches living after two years
  ## UL
  .0475 * 6417
  ## 305
  ## LL
  0.1047* 6417
  ## 672
  
  ## OR calculating as a sum of nutcracker and pilf
  
  # ~ 6416.667 caches remaining after nutcracker retreival
  
  6416.667-(6416.667*.40)
  #3850 caches remaining
  
  OR.frame(data = germ.subset, comparison = "all",  cache.characteristic = "study.area", germ.year = "germ.2013.2014")
  
  ##TB 
  .8479 * 3850
  # 3264
  
  OR.frame(data = survival, comparison = "all", cache.characteristic = "study.area", germ.year = "living.2014.y.n")
  
  #TBSA cumulative survival
  OR.frame.pilferage(data = survival, comparison = "tb.sa",  cache.characteristic = "object")
  OR.frame(data = survival, comparison = "tb.sa",  cache.characteristic = "object", germ.year = "germ.2013.y.n")
  OR.frame(data = survival, comparison = "tb.sa",  cache.characteristic = "object", germ.year = "survive.13.14")
  
  
  cumulative.survival.tbsa <- data.frame(Microsite = rep(c("Rock","Tree"), each = 3), 
                                         Estimate  = c(0.569, 0.483, 0.207,
                                                       .534,.279,.155),
                                         LL        = c(.433, .351, .116,
                                                       .400,.170,.078),
                                         UL        =c(.696, .617,.337,
                                                      .665, .411,.279),
                                         Stage     =rep(c("Pilferage", "Germination", "Survival"),  2),
                                         Year      =rep(c("2012", "2013", "2014"),  2)) %>%
    mutate(Stage = factor(Stage, levels = c("Pilferage", "Germination", "Survival")))
  

  
  
  library(ggplot2)
  
  ggplot(cumulative.survival.tbsa, aes(x     = Stage,
                          y     = Estimate, 
                          group = Microsite,
                          ymin  = LL, 
                          ymax  = UL))  + 
    geom_point(position = position_dodge(width = .3), 
               lwd      = 6, 
               col      = rep(c("black", "gray"),3)) + 
    geom_linerange(position = position_dodge(width = .3), 
                   lwd      = 1.5, 
                   col      = rep(c("black", "gray"),  3)) +
    labs(x = "", y = " Proportion survived") + 
    theme_bw() + 
    theme(panel.border = element_rect(fill     = NA, 
                                      color    = "black", 
                                      size     = 1.5, 
                                      linetype ="solid")) +
    scale_y_continuous(limits = c(0.0,1), 
                       breaks = c(0,.25,.5,.75,1)) +
    theme(plot.title = element_text(lineheight = 1.1, 
                                    size       = 25, 
                                    vjust      = 2)) +
    theme(axis.title.x = element_text(size  = 25,  
                                      vjust = 0)) + 
    theme(axis.text.x = element_text(size = 25)) + 
    theme(axis.title.y = element_text(size  = 25, 
                                      vjust = 2.75)) +
    theme(axis.text.y = element_text(size = 25)) + 
    theme(panel.grid.major = element_line(colour = "gray")) 
  
  
  OR.frame(data = survival, comparison = "tb.tl",  cache.characteristic = "object", germ.year = "living.2014.y.n")
  OR.frame(data = survival, comparison = "wc.sa",  cache.characteristic = "object", germ.year = "living.2014.y.n")
  OR.frame(data = survival, comparison = "wc.tl",  cache.characteristic = "object", germ.year = "living.2014.y.n")
  
  
  


