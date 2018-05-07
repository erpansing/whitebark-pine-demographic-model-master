##-----------------------------------------------------------##
##                                                           ##
##                 WBP Demographic Model                     ##
##                    Libby Pansing                          ##
##-----------------------------------------------------------##

# Note... having issues with masking of dplyr functions when loading
# plyr later in the script. IF you want to re-source the other
# scripts, you'll need to unload dplyr and plyr and reload plyr.
# Or restart R.
library(dplyr)
library(popbio)
library(tidyr)

## Import relevant parameter estimates calculated in other scripts.

## Required files:
## 1) GRIN parameter estimates.R
## 2) GRIN WBP Survival Estimates RMark Nest.R
## 3) 2017 YNP Data.Rda
## 4) survival.Rda

source("/Users/elizabethpansing/Box Sync/PhD/Code/WBP Demographic Model Master/WBP Model Active/GRIN parameter estimates.R") 

n <- c(62, 618, 79, 65, 91,  353, 30, 600, 50, 500, 120, 600)

## Write function for determining leaf area index (LAI) as a function
## of stage specific dbh and the number of trees on the landscape

##***********************************************************##
##                        Define LAI                         ##
##***********************************************************##

#mean dbh for each stage
d1 <- function() {   # dbh SEED1
  return(0)
}

d2 <- function() {   #dbh SEED2 
  return(0)
}

d3 <- function(){ #dbh CS
  return(0)
}

d4 <- function(){ #dbh SD1
  return(0)
}

d5 <- function(){ # dbh SAP
  return(6.69)
}

d6 <-function() { #dbh MA
  return(30.2)  # mean dbh taken from trees cored on HM in 2015. The average dbh of trees >12.5cm dbh (the sapling)
}

## Leaf area coefficients. Define the relationship between leaf area and diameter.
## Estimated via MLE assuming the general form y = ax^b

alpha1 <- function(){
  return(0.456)
}

alpha2 <- function(){
  return(0.0736)
}

alpha3 <- function(){
  return(2.070)
}


l <- function(){
  matrix(c(d1(),      # SEED1 do not contribute to LAI
    d2(),      # SEED2 do not contribute to LAI
    d3(),      # CS do not contribute to LAI
    alpha1(),  # SD1 don't have DBH but do contribute to LAI 
    alpha2() * d5() ^ alpha3(),
    alpha2() * d6() ^ alpha3(),
    d1(),      # SEED1 do not contribute to LAI
    d2(),      # SEED2 do not contribute to LAI
    d3(),      # CS do not contribute to LAI
    alpha1(),  # SD1 don't have DBH but do contribute to LAI 
    alpha2() * d5() ^ alpha3(),
    alpha2() * d6() ^ alpha3()), nrow = 6)
}


LAIb    <- function(){     # Background leaf area index. This is where competition can be incorporated...
  return(0)
}


LAI <- function(x) {       # LAI of the study area
  l <- l()
  c((t(matrix(l[,1], nrow = 6)) %*% x[1:6])/10000 + LAIb(), (t(matrix(l[,2], nrow = 6)) %*% x[7:12])/10000 + LAIb())
}

LAI(n)


## Now define the distributions from which survival and
## transition rates will be drawn to create stochastic demographic rates

##-----------------------------------------------------------##
##              Empirically derived survival                 ##
##              & transition distributions                   ##
##-----------------------------------------------------------##


##-----------------------------------------------------------##
##                          SEED1                            ##
##-----------------------------------------------------------##

s_SEED1 <- 0                         
                                     # Assume seeds either transition to SEED2, 
                                     # transition to CS (first year seedling)
                                     # or die


t1_SEED1    <- function(size = 1){       # survival probability of seeds (i.e., survive and transition to SEED2, but do not germinate)
  rbeta(n = size,                        # Drawn from a beta to give a 
        shape1 = SEED1_survive_alpha,    # probability of seeds transitioning
        shape2 = SEED1_survive_beta)     # to SEED2 stage
} 



t2_SEED1 <- function(size = 1){       # Germination probability of seeds
  rbeta(n = size,                     # Drawn from a beta to give a 
        shape1 = SEED1_germ_alpha,    # probability of seeds transitioning
        shape2 = SEED1_germ_beta)     # to SEED2 stage
}

##-----------------------------------------------------------##
##                          SEED2                            ##
##-----------------------------------------------------------##

s_SEED2 <- 0                         # Assume seeds either transition to 
                                     # CS (first year seedling) or die

t3_SEED2 <- function(size = 1){
  # Germination probability of seeds
  rbeta(n = size,                     # Drawn from a beta to give a 
        shape1 = SEED2_germ_alpha,    # probability of seeds transitioning
        shape2 = SEED2_germ_beta)     # to from SEED2 to CS stage (i.e., germination)
}                                     

##-----------------------------------------------------------##
##                            CS                             ##
##                     First Year Seedling                   ##
##-----------------------------------------------------------##

s_CS <- 0                            # Assume seeds transition to 
                                     # SD or die

t_CS     <- function(size = 1){          # Survival probability of first 
  rbeta(n = size,                    # year seedlings (cotyledon seedlings)
        shape1 = CS_survive_alpha,   # Drawn from a beta to give prob
        shape2 = CS_survive_beta)    # of transitioning to SD stage
}


##-----------------------------------------------------------##
##                          SD                               ##
##-----------------------------------------------------------##

s_SD     <- function(size = 1){         # survival probability of seedlings
  rbeta(n = size,                   # Drawn from a beta to give prob 
        shape1 = SD_survive_alpha,  # surviving any given year
        shape2 = SD_survive_beta)
}

##-----------------------------------------------------------##
##                           SAP                             ##
##-----------------------------------------------------------##

s_SAP     <- function(){         # Survival probability of saplings
  return(0.8)                    # Still looking for a distribution to 
}                                # use, so for now assuming constant

##-----------------------------------------------------------##
##                            MA                             ##
##-----------------------------------------------------------##

s_MA     <- function(){         # Survival rate of reproductively mature adults
  0.99                          # Assume limited death from senescence because 
}                               # of long lived nature of wbp (lifespan up to 1200 yrs)
                                # For now, assuming constant.

##-----------------------------------------------------------##
##                          DEFINE                           ##
##                         SURVIVAL                          ##
##                          VECTOR                           ##
##-----------------------------------------------------------##

survival_vector <- function(size = 1){   #survival vector
  c(
    s_SD(size = size),
    s_SAP(),
    s_MA())
}

survival_vector()

##-----------------------------------------------------------##
##                          DEFINE                           ##
##                       RESIDENCE TIME                      ##
##                          VECTOR                           ##
##-----------------------------------------------------------##

residence_SEED1 <- 1    # # years as seedling (SEED1)

residence_SEED2 <- 1    # # years as seedling (SEED2)

residence_CS    <- 1    # # years as seedling (CS)

residence_SD    <- 28   # # years as seedling (SD)

residence_SAP   <- 20   # # years as sapling (SAP)

residence_MA   <- Inf  # # years as reproductively mature
 
##-----------------------------------------------------------##
##                          DEFINE                           ##
##                       RESIDENCE TIME                      ##
##                          VECTOR                           ##
##-----------------------------------------------------------##


residence_vector <- 
  c(residence_SD,
    residence_SAP,
    residence_MA)


residence_vector


# ##-----------------------------------------------------------##
# ##                         FERTILITY                         ##
# ##-----------------------------------------------------------##
# 
# No_caches <- function(t, size = 1){    #Seed production in wbp is periodic
#   result <- NULL                   # with masting years every ~ 4 years
#                                    # so define cone production as a function
#   for(i in 1:size){                # of time described by cos with normally distributed error
#     value <- ((12.5*cos(1.5 * t) + 14 + rnorm(1, sd = 3.5)) * 45)/ 3
#     if( value >= 0){               # Max values and expected values from
#       result[i] <- value           # IGBST cone monitoring since 1980
#     } else if(value < 0){          # 
#       result[i] <- value - value   # # caches assumes 45 seeds/cone
#     }                              # and 3 seeds/cache. All available seeds cachek
#   }                                # Assumes 45% of caches created are left for regeneration.
#   return(result) 
# } 



##Here
##-----------------------------------------------------------##
##                         FERTILITY                         ##
##-----------------------------------------------------------##

No_seeds_per_cone <- 45

No_cones <- function(t, size = 1){ # Seed production in wbp is periodic
  result <- NULL                   # with masting years every ~ 4 years
                                   # so define cone production as a function
  for(i in 1:size){                # of time described by cos with normally distributed error
    value <- (12.5*cos(1.5 * t) + 14 + rnorm(1, sd = 3.5))
    if( value >= 0){               # Max values and expected values from
      result[i] <- value           # IGBST cone monitoring since 1980
    } else if(value < 0){          # 
      result[i] <- value - value   # # caches assumes 45 seeds/cone
    }                              # and 3 seeds/cache. All available seeds cachek
  }                                # Assumes 45% of caches created are left for regeneration.
  return(result) 
} 

##-----------------------------------------------------------##
##         Define variables assumed fixed & known            ##
##-----------------------------------------------------------##
Pfind  <- 0.55   # Proportion of seeds found by nutcrackers
Pcons  <- 0.3    # Proportion of seeds consumed by nutcracers (prior to caching?)
nBirds <- 3      # No. Clark's nutcrackers in the theoretical population
SpC    <- 3      # No. seeds per cache


e1 <- matrix(c(1,0,0,0,0,0,0,0,0,0,0,0))
e6 <- matrix(c(0,0,0,0,0,0,1,0,0,0,0,0))

No_caches_1 <- function(t, size = 1, x){
  (No_cones(t, 1) * No_seeds_per_cone * x[6] * (1-Pcons)* (1-Pfind)/3 * (1-dispersal12) +  
     No_cones(t, 1) * No_seeds_per_cone * x[12] * (1-Pcons)* (1-Pfind)/3 * dispersal21)* e1
}

No_caches_2 <- function(t, size = 1, x){
  (No_cones(t, 1) * No_seeds_per_cone * x[12] * (1-Pcons)* (1-Pfind)/3 * (1-dispersal21) +
     No_cones(t, 1) * No_seeds_per_cone * x[6] * (1-Pcons)* (1-Pfind)/3 * dispersal12)* e3
}

##-----------------------------------------------------------##
##                       Germination                         ##
##-----------------------------------------------------------##

##-----------------------------------------------------------##
##         Define variables dependent on time vars           ##
##-----------------------------------------------------------##

SpB    <- function(x){  # Number of seeds available to each bird
  x[1]/nBirds
}

## Define reduction factors. These variables reduce 
## 1) rALS decreases germination as light availability decreases
## 2) rCache increases caching propensity as seed availability increases


rALS   <- function(x){   
  1/(1 + exp(2*(LAI(x)-3))) 
}

# rCache1 <- function(x){
#   0.73/(1+ exp((31000-SpB(x[1:6]))/3000))  
# }
# 
# rCache2 <- function(x){
#   0.73/(1+ exp((31000-SpB(x[7:12]))/3000))  
# }

e3 <- matrix(c(0,0,1,0,0,0,0,0,0,0,0,0))
e9 <- matrix(c(0,0,0,0,0,0,0,0,1,0,0,0))

## Define germination probability
# r2 <- function(x){
#   as.vector(((1-Pfind)*(1-Pcons)/SpC) * rCache(x) * rALS(x)) * e3
# }

## Define germination probability
# r2 <- function(x){
#   as.vector(((1-Pfind)*(1-Pcons)/SpC) * rCache(x) * rALS(x)) * e3
# }


# No. seeds
germ1stpop1 <- function(t, size = 1, x){
  as.vector(x[1])*((1-Pcons)*rCache1(x)/3) * (1-Pfind) * as.vector(rALS(x)[1]) * rbeta(n = 1, shape1 = SEED1_germ_alpha, shape2= SEED1_germ_beta) * e3
}

germ1stpop2 <- function(t, size = 1, x){
  as.vector(x[7])*((1-Pcons)*rCache2(x)/3) * (1-Pfind) * as.vector(rALS(x)[2]) * rbeta(n = 1, shape1 = SEED1_germ_alpha, shape2= SEED1_germ_beta) * e9
}

germ2ndpop1 <- function(t, size = 1, x){
  as.vector(x[2] * rALS(x)) * rbeta(n = 1, shape1 = SEED2_germ_alpha, shape2 = SEED2_germ_beta) * e3 
}

germ2ndpop2 <- function(t, size = 1, x){
  as.vector(x[8] * rALS(x)) * rbeta(n = 1, shape1 = SEED2_germ_alpha, shape2 = SEED2_germ_beta) * e9 
}

##-----------------------------------------------------------##
##                  GET MATRIX ELEMENTS                      ##
##-----------------------------------------------------------##
si <- function(size = 1){                  # Gives probability of surviving and staying in the
  (1 - (1/residence_vector)) *         # same life stage for those life stages
    survival_vector(size = size)       # that have residence time > 1 (i.e., persist in the
}                                      # same life stage for > 1 year)

ti <- function(size = 1) {                 # Gives probability of surviving and transitioning
  (1/residence_vector) *               # to the next life stage for those life stages
    survival_vector(size = size)       # that have residence time > 1 (i.e., persist in the
}                                      # same life stage for > 1 year)

dispersal21 <- 0.05
dispersal12 <- 0.10


S <- function(t){
  matrix(c(0,            0,         0,        0,         0,            0,          0,            0,       0,         0,         0,  dispersal21*No_s, 
           t1_SEED1(1),  0,         0,        0,         0,            0,          0,            0,       0,         0,         0,            0,
           0,            0,         0,        0,         0,            0,          0,            0,       0,         0,         0,            0,
           0,            0,   t_CS(1),  si(1)[1],        0,            0,          0,            0,       0,         0,         0,            0,
           0,            0,         0,  ti(1)[1], si(1)[2],            0,          0,            0,       0,         0,         0,            0,
           0,            0,         0,        0,  ti(1)[2],     si(1)[3],          0,            0,       0,         0,         0,            0,
           #######################################################################################################################################
           0,            0,         0,        0,         0,  dispersal12,          0,            0,       0,         0,         0,            0, 
           0,            0,         0,        0,         0,            0, t1_SEED1(1),           0,       0,         0,         0,            0,  
           0,            0,         0,        0,         0,            0,          0,            0,       0,         0,         0,            0,  
           0,            0,         0,        0,         0,            0,          0,            0, t_CS(1),  si(1)[1],         0,            0,  
           0,            0,         0,        0,         0,            0,          0,            0,       0,  ti(1)[1],  si(1)[2],            0,  
           0,            0,         0,        0,         0,            0,          0,            0,       0,         0,  ti(1)[2],     si(1)[3]),  
         byrow = T, nrow = 6) 
}

t <- 1
S(t) 





##-----------------------------------------------------------##
##                                                           ## 
##              Fire return interval functions               ## 
##                                                           ## 
##-----------------------------------------------------------##
## Westerling et al. (2011) predict a decrease from historic fire 
## return intervals in the GYE from >120 years to <30 by the end 
## of the 21st century. The following functions describe that decrease
## in fire return interval. 

## We define current fire return intervals as ~230 years in 
## (Larson et al. 2009) in wbp forests. We set middle of the 
## century values following figure 3 in Westerling et al. 2011
## and end of century values as 30 years. 

## The following function uses those timeframes to estimate the shape
## of the declined in fire rates (i.e., lambda)

fire_return_decrease <- data.frame(Year = c(0, 40, 100), Interval = c(230, 75, 30)) 
fit <- lm(log(Interval)~Year, data = fire_return_decrease)
summary(fit)
predicted <- exp(fit$coefficients[1] + fit$coefficients[2]* fire_return_decrease$Year)
plot(Interval~Year, data = fire_return_decrease)
points(fire_return_decrease$Year, predicted, pch = 19)

## NEED TO: Try to put a distribution on the average rate
interval <- function(t){
  exp(fit$coefficients[1] + fit$coefficients[2]* t)
} 

## Function that determines whether fire occurs in current year
fire_current_year <- function(t, n = 1){
  c(rbinom(size = 1, n = 1, prob = 1/interval(t)),rbinom(size = 1, n = 1, prob = 1/interval(t))) 
}



##-----------------------------------------------------------##
##              Function that projects pop                   ##
##              sizes and incorporates fire                  ##
##-----------------------------------------------------------##

library(plyr)

n <- c(62, 580 + 38, 79, 65, 91,  353, 30, 600, 50, 500) # Arbitrary starting pop size vectors


##############################################################################################################################################            
############################################################################################################################################## 
############################################################################################################################################## 
############################################################################################################################################## 

project <- function(projection_time, n0, reps = 100, fire = T){       # stochastic projection function 
  # that tracks stage based pop sizes 
  # over time for reps number of iterations
  
  results <- stochastic_matrixes <-                          #Create null matrix that will hold stage
    array(0, dim = c(projection_time, length(n0) + 1, reps)) # based population sizes and year tracker
  
  for(j in 1:reps){       # Iterate through i years (projection_time) of population growth j times (iterations)
    # Incorporate FIRE
    if(j == 1){           
      fire_tracker   <- matrix(c(rep(1:reps, each = projection_time),rep(0, projection_time*reps*3)), 
                               nrow = projection_time * reps, ncol = 4, byrow = F)  
      LAI_tracker <- matrix(c(rep(1:reps, each = projection_time),rep(0, projection_time*reps*3)), 
                            nrow = projection_time * reps, ncol = 4, byrow = F)
    } else if(j != 1){
      fire_tracker <- fire_tracker
      LAI_tracker <- LAI_tracker
    }
    pops <- matrix(0, nrow = projection_time, ncol = length(n)) # Creates empty matrix to hold population sizes and LAI
    
    for(i in 1:projection_time){        
#--------------------------------------------------------------------------------------------------------------------------------------------      
      ## Set up initials for beginning of each iteration
      if (i == 1){
        n <- n0
        tSinceFire <- c(1,1)
      }else if(i != 1){
        n <- n
        tSinceFire <- tSinceFire +1
      }
      
      ## Update time counter for each time step
      t <-  i   # time counter
     
      ## Update LAI tracker
      LAI_tracker[j*projection_time - (projection_time)+i,2:4] <- c(i, LAI(n))
#--------------------------------------------------------------------------------------------------------------------------------------------      
      ##############################################################################################################################################            
      ##                                                                FIRE POSSIBLE                    
      ##############################################################################################################################################            
      if(fire == T){
        fire_current <- fire_current_year(t)  # Determine whether this iteration experiences a stand replacing fire
        fire_tracker[j*projection_time - (projection_time)+i,2:3] <- c(i, fire_current)
        
        ## 1) There's a fire. This kills the population. Assumes stand replacing burn that impacts entire population
        ##    And that no regeneration occurs the year of the fire
        ## 2) There's no fire and it's >1 year after fire. In this case, there are no modifications and the 
        ##    system proceeds as normal. If it's the year after fire, the only seed source is from the other subpopulation
        #--------------------------------------------------------------------------------------------------------------------------------------------              
        ## Fire can occur in different combinations
        ## 1) Fire in pop1 but not pop2: fire_current = 0,1
        ## 2) Fire in pop2 but not pop1: fire_current = 1,0
        ## 3) Fire in both populations: fire_current = 1,1
        
        ## 1) FIRE IN 1 BUT NOT 2
        
        if(fire_current[1] == T & fire_current[2] == F){                  
          
          tSinceFire <- C(0, tSinceFire[2]) 
          
          # Assuming stand replacing burn with no survival and no regeneration.
          # Most fires go out with first snow. e.g., Romme 1982
          mat <- S(t = t)
          out <- c(t(mat%*%n + germ1stpop1(t = t, size = 1, x = n) + 
                       germ1stpop2(t = t, size = 1, x = n) +
                       germ2ndpop1(t = t, size = 1, x = n) + 
                       germ2ndpop2(t = t, size = 1, x = n) +
                       No_seeds1(t = t, size = 1, x = n) +
                       No_seeds2(t = t, size = 1, x = n)))  # Defines the intermediate population size
          pops[i,] <- c(0,0,0,0,0,0, out[7:12])
          
          n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
          
        } 
        
        ## 2) FIRE IN 2 BUT NOT 1
        
        if(fire_current[1] == F & fire_current[2] == T){                  
          
          tSinceFire <- c(tSinceFire[1], 0)
          
          # Assuming stand replacing burn with no survival and no regeneration.
          # Most fires go out with first snow. e.g., Romme 1982
          mat <- S(t = t)
          out <- c(t(mat%*%n + 
                       germ1stpop1(t = t, size = 1, x = n) + 
                       germ1stpop2(t = t, size = 1, x = n) +
                       germ2ndpop1(t = t, size = 1, x = n) + 
                       germ2ndpop2(t = t, size = 1, x = n) +
                       No_seeds1(t = t, size = 1, x = n) +
                       No_seeds2(t = t, size = 1, x = n)))  # Defines the intermediate population size 
          
          pops[i,] <- c(out[1:6], 0,0,0,0,0,0)
          
          n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
      
        }
        
        ## 3) FIRE IN BOTH
        
        if(fire_current[1] == T & fire_current[2] == T){                  
          
          tSinceFire <- c(0, 0)
          
          # Assuming stand replacing burn with no survival and no regeneration.
          # Most fires go out with first snow. e.g., Romme 1982
          
          pops[i,] <- c(0,0,0,0,0,0, 0,0,0,0,0,0)
          
          n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
        
        }
        #--------------------------------------------------------------------------------------------------------------------------------------------              
        #                                                 No fire in current year in either subpopulation
        #--------------------------------------------------------------------------------------------------------------------------------------------              
        else if((fire_current[1] == F & fire_current[2] == F)){
          # tSinceFire <- tSinceFire + 1
          
          mat <- S(t = t)
          pops[i,]  <- c(t(mat%*%n + germ1st(t = t, size = 1, x = n) + germ2nd(t = t, size = 1, x = n) + No_seeds(t = t, size = 1, x = n)))  # Defines the intermediate population size 
          
          n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
         
        } 
        
      }
      ##############################################################################################################################################            
      ##                                                          FIRE = FALSE                    
      ############################################################################################################################################## 
      else if(fire == F){
        fire_tracker[j*projection_time - (projection_time)+i,2:3] <- c(i, 0)
        
        # tSinceFire <- tSinceFire + 1
        
        mat <- S(t = t)
        pops[i,] <- c(t(mat%*%n + germ1stpop1(t = t, size = 1, x = n) + 
                                  germ1stpop2(t = t, size = 1, x = n) +
                                  germ2ndpop1(t = t, size = 1, x = n) + 
                                  germ2ndpop2(t = t, size = 1, x = n) +
                                  No_seeds1(t = t, size = 1, x = n) +
                                  No_seeds2(t = t, size = 1, x = n)))  # Defines the intermediate population size
        
        n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1) 
      }
      ############################################################################################################################################## 
    }
    pops
    pops <- cbind(pops, rep(1:projection_time))  # Appends matrix to keep track of time during iteration
    
    
    results[, ,j] <- pops                        # Combines iterations into a j dimensional array
    
  }
  
  pop_sizes <- plyr::adply(results, 3)           # Changes array to dataframe so easier to manipulate later
  colnames(pop_sizes) <- c("Iteration", "SEED1_1", "SEED2_1", "CS_1", "SD_1", "SAP_1", "MA_1","SEED1_2", "SEED2_2", "CS_2", "SD_2", "SAP_2", "MA_2", "t")
  
  # fire_intervals <- cbind(iteration, intervals)  # Dataframe keeping track of fire return intervals
  
  results <- list(pop_sizes = pop_sizes, fire_tracker = fire_tracker, LAI_track = LAI_tracker) 
  
  return(results)
}
##############################################################################################################################################            
############################################################################################################################################## 
############################################################################################################################################## 
############################################################################################################################################## 





##-----------------------------------------------------------##
##                 Population projection                     ##
##-----------------------------------------------------------##

projection <- project(projection_time = 500, n0 = n, reps = 50) 

projection <- gather(projection$pop_sizes, Stage, Count, -Iteration, -t) %>%  
  filter(., !Stage == "SEED1") %>%          # Pop sizes in dataframe format
  filter(., !Stage == "SEED2") %>%          # and excluding seed numbers (most don't think of seeds)
  group_by(., Iteration, t) %>%             # as a part of the population, so presenting numbers as 
  summarise_at(., vars(Count), funs(sum)) %>% # number of living trees (i.e., post germination) is more
  ungroup(.)                                # intuitive


ggplot(data = projection, aes(x = t, y = Count, col = Iteration)) +  # plot pop sizes for all iterations.
  geom_line(lwd = 1) +
  theme(legend.position="none") +
  theme(axis.title.x=element_text( size=18, vjust=0)) +
  theme(axis.text.x=element_text(size=18))  +
  theme(axis.title.y=element_text( size=18, vjust=2.75, face = "bold")) +
  theme(axis.text.y=element_text(size = 18))


## Plot of projection iteration 1

projection1 <- projection %>% 
  filter(., Iteration == 1)

ggplot(data = projection1, aes(x = t, y = Count)) +  # plot pop sizes for all iterations.
  geom_line(lwd = 1) +
  theme(legend.position="none") +
  theme(axis.title.x=element_text( size=18, vjust=0)) +
  theme(axis.text.x=element_text(size=18))  +
  theme(axis.title.y=element_text( size=18, vjust=2.75, face = "bold")) +
  theme(axis.text.y=element_text(size = 18))

## Plot of projection iteration 1 for t = 30 years

projection1t30 <- projection1 %>% 
  filter(., t < 3)

ggplot(data = projection1t30, aes(x = t, y = Count)) +  # plot pop sizes for all iterations.
  geom_point() +
  geom_line(lwd = 1) +
  theme(legend.position="none") +
  theme(axis.title.x=element_text( size=18, vjust=0)) +
  theme(axis.text.x=element_text(size=18))  +
  theme(axis.title.y=element_text( size=18, vjust=2.75, face = "bold")) +
  theme(axis.text.y=element_text(size = 18))


hist(projection$Count[projection$t == 500], breaks = 20,
     main = "", xlab = "Population size at time = 500 years")

##-----------------------------------------------------------##
##              Get stochastic lambda and                    ##
##                    elasticities                           ##
##-----------------------------------------------------------##

## Create a list of 10,000 matrixes to use in stochastic lambda 
## and elasticity analyses

reps <- 10000
stochastic_matrixes <- array(0, dim = c(6,6,reps))

for(i in 1:reps){
  stochastic_matrixes[, , i] <- S(i)
}

A <- list(NULL)

for(i in 1:reps){
  mat <- stochastic_matrixes[,,i]
  A[[i]] <- mat
}

rm(stochastic_matrixes)

## Estimate elasticities
stoch.elast<- stoch.sens(A, tlimit = 500)
stoch.elast

## Estimate stochastic lambda
sgr <- stoch.growth.rate(A)
sgr_real <- exp(sgr$approx)





