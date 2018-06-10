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
library(ggplot2)

## Import relevant parameter estimates calculated in other scripts.

## Required files:
## 1) GRIN parameter estimates.R
## 2) GRIN WBP Survival Estimates RMark Nest.R
## 3) 2017 YNP Data.Rda
## 4) survival.Rda

source("/Users/elizabethpansing/Box Sync/PhD/Code/WBP Demographic Model Master/WBP Model Active/GRIN parameter estimates.R") 

n <- c(62, 618, 79, 65, 91,  353, 30, 600, 50, 500, 120, 600)
area <- 1e7

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
  matrix(c(#d1(),      # SEED1 do not contribute to LAI
    d2(),      # SEED2 do not contribute to LAI
    d3(),      # CS do not contribute to LAI
    alpha1(),  # SD1 don't have DBH but do contribute to LAI 
    alpha2() * d5() ^ alpha3(),
    alpha2() * d6() ^ alpha3(),
    #d1(),      # SEED1 do not contribute to LAI
    d2(),      # SEED2 do not contribute to LAI
    d3(),      # CS do not contribute to LAI
    alpha1(),  # SD1 don't have DBH but do contribute to LAI 
    alpha2() * d5() ^ alpha3(),
    alpha2() * d6() ^ alpha3()), nrow = 5)
}


LAIb    <- function(){     # Background leaf area index. This is where competition can be incorporated...
  return(0)
}


LAI <- function(x) {       # LAI of the study area
  l <- l()
  c((t(matrix(l[,1], nrow = 5)) %*% x[1:5])/area + LAIb(), (t(matrix(l[,2], nrow = 5)) %*% x[6:10])/area + LAIb())
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

s_MA     <- function(size = 1){ # Survival rate of reproductively mature adults
  rbeta(n = size,              # Assume limited death from senescence because 
        shape1 = MA_s_alpha,   # of long lived nature of wbp (lifespan up to 1200 yrs)
        shape2 = MA_s_beta)   # For now, assuming constant.
}

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


##-----------------------------------------------------------##
##                         FERTILITY                         ##
##-----------------------------------------------------------##

No_seeds_per_cone <- 45

rcones <- function(x){
  0.5/(1 + exp(5 *(LAI(x)-2.25)))
}

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

##-----------------------------------------------------------##
##               Define dispersal parameters                 ##
##-----------------------------------------------------------##

dispersal_25 <- function(mean_dispersal){   # Here we're assuming that the populations are separated by 25 km
  # x <- unlist(mean_dispersal[sample(1:nrow(mean_dispersal), size = 1),])
  pgamma(shape = dispersal_alpha, rate = dispersal_beta, 25, lower.tail = FALSE)

  # out <- pgamma(q = 25, shape = x$Alpha, rate = x$Beta, lower.tail = FALSE)
  # return(out)
}

e1 <- matrix(c(1,0,0,0,0,0,0,0,0,0,0,0))
e6 <- matrix(c(0,0,0,0,0,0,1,0,0,0,0,0))

## Do birds cache more in subpopulation 1 during years of fire in subpopulation 2? Or are they 
## lost to recruitment?

No_caches_1_nofire <- function(cones_1, cones_2, dispersal12, dispersal21, t, size = 1, x){
  (cones_1 * No_seeds_per_cone * x[5] * (1-Pcons)* (1-Pfind)/3 * (1-dispersal12) +  
    cones_2 * No_seeds_per_cone * x[10] * (1-Pcons)* (1-Pfind)/3 * dispersal21)
}

No_caches_1_fire1 <- 0 

No_caches_1_fire2 <- function(cones_1, t, size = 1, x){
     cones_1 * No_seeds_per_cone * x[5] * (1-Pcons)* (1-Pfind)/3
}

No_caches_2_nofire <- function(cones_1, cones_2,dispersal12, dispersal21, t, size = 1, x){
  (cones_1 * No_seeds_per_cone * x[10] * (1-Pcons)* (1-Pfind)/3 * (1-dispersal21) +
     cones_2 * No_seeds_per_cone * x[5] * (1-Pcons)* (1-Pfind)/3 * dispersal12)
}

No_caches_2_fire1 <- function(cones_2, t, size = 1, x){
  cones_2 * No_seeds_per_cone * x[10] * (1-Pcons)* (1-Pfind)/3 
}

No_caches_2_fire2 <- 0


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


rALS_germ   <- function(x){   
  1/(1 + exp(2*(LAI(x)-3))) 
}

rALS_sd     <- function(x){
  1/(1 + exp(2*(LAI(x)-4.5)))
}

# rCache1 <- function(x){
#   0.73/(1+ exp((31000-SpB(x[1:6]))/3000))  
# }
# 
# rCache2 <- function(x){
#   0.73/(1+ exp((31000-SpB(x[7:12]))/3000))  
# }

e2_1 <- matrix(c(0,1,0,0,0,0,0,0,0,0))
e2_2 <- matrix(c(0,1,0,0,0,0,0,0,0,0))
e7_1 <- matrix(c(0,0,0,0,0,0,1,0,0,0))
e7_2 <- matrix(c(0,0,0,0,0,0,1,0,0,0))

## Define germination probability
# r2 <- function(x){
#   as.vector(((1-Pfind)*(1-Pcons)/SpC) * rCache(x) * rALS(x)) * e3
# }

## Define germination probability
# r2 <- function(x){
#   as.vector(((1-Pfind)*(1-Pcons)/SpC) * rCache(x) * rALS(x)) * e3
# }


# No. germinated
germ1stpop1 <- function(caches1,t, size = 1, x){
  caches1* as.vector(rALS_germ(x)[1]) * rbeta(n = 1, shape1 = SEED1_germ_alpha, shape2= SEED1_germ_beta) * e2_1
}

germ1stpop2 <- function(caches2, t, size = 1, x){
  caches2 * as.vector(rALS_germ(x)[2]) * rbeta(n = 1, shape1 = SEED1_germ_alpha, shape2= SEED1_germ_beta) * e7_1
}


germ2ndpop1 <- function(t, size = 1, x){
  as.vector(x[1] * rALS_germ(x)[1]) * rbeta(n = 1, shape1 = SEED2_germ_alpha, shape2 = SEED2_germ_beta) * e2_2
}

germ2ndpop2 <- function(t, size = 1, x){
  as.vector(x[6] * rALS_germ(x)[2]) * rbeta(n = 1, shape1 = SEED2_germ_alpha, shape2 = SEED2_germ_beta) * e7_2
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

s1 <- si()
s2 <- si()
t1 <-ti()
t2 <-ti()

S <- function(t){
          # SEED2       CS                       SD           SAP           MA                SEED2_2       CS_2                     SD_2      SAP_2                  MA_2       
  matrix(c(
              0,         0,                       0,           0,  caches1*t1_SEED1(),             0,          0,                       0,         0,                   0,
              0,         0,                       0,           0,                   0,             0,          0,                       0,         0,                   0,
              0,    t_CS(),      s1[1]*rALS_sd(n)[1],          0,                   0,             0,          0,                       0,         0,                   0,
              0,         0,      t1[1]*rALS_sd(n)[1],      s1[2],                   0,             0,          0,                       0,         0,                   0,
              0,         0,                       0,       t1[2],               s1[3],             0,          0,                       0,         0,                   0,
           #############################################################################################################################################################################
           
              0,         0,                       0,           0,                   0,             0,          0,                       0,         0,   caches2*t1_SEED1(),  
              0,         0,                       0,           0,                   0,             0,          0,                       0,         0,   caches2*t2_SEED1(),
              0,         0,                       0,           0,                   0,             0,     t_CS(),     s2[1]*rALS_sd(n)[2],         0,                   0,  
              0,         0,                       0,           0,                   0,             0,          0,     t2[1]*rALS_sd(n)[2],     s2[2],                   0,  
              0,         0,                       0,           0,                   0,             0,          0,                       0,     t2[2],                s2[3]),  
         byrow = T, nrow = 10) 
}

t <- 1
S(t) 


S_fire_1 <- function(t){
  matrix(c(
                0,         0,        0,         0,            0,            0,       0,                       0,         0,                    0,
                0,         0,        0,         0,            0,            0,       0,                       0,         0,                    0,
                0,         0,        0,         0,            0,            0,       0,                       0,         0,                    0,
                0,         0,        0,         0,            0,            0,       0,                       0,         0,                    0,
                0,         0,        0,         0,            0,            0,       0,                       0,         0,                    0,
    ############################################################################################################################################
                0,         0,        0,         0,            0,            0,       0,                       0,         0,   caches2*t1_SEED1(),  
                0,         0,        0,         0,            0,            0,       0,                       0,         0,                    0,  
                0,         0,        0,         0,            0,            0, t_CS(1),     s2[1]*rALS_sd(n)[2],         0,                    0,  
                0,         0,        0,         0,            0,            0,       0,     t2[1]*rALS_sd(n)[2],     s2[2],                    0,  
                0,         0,        0,         0,            0,            0,       0,                       0,     t2[2],                 s2[3]),  
         byrow = T, nrow = 10) 
}

S_fire_1(t)


S_fire_2 <- function(t){
  matrix(c(
             0,         0,                       0,         0,   caches1*t1_SEED1(),                0,       0,         0,         0,            0,
             0,         0,                       0,         0,                    0,                0,       0,         0,         0,            0,
             0,   t_CS(1),     s1[1]*rALS_sd(n)[1],         0,                    0,                0,       0,         0,         0,            0,
             0,         0,     t1[1]*rALS_sd(n)[1],     s1[2],                    0,                0,       0,         0,         0,            0,
             0,         0,                       0,     t1[2],                s1[3],                0,       0,         0,         0,            0,
           ###################################################################################################################################################
      
             0,         0,                       0,         0,                    0,                0,         0,       0,         0,            0,
             0,         0,                       0,         0,                    0,                0,         0,       0,         0,            0,  
             0,         0,                       0,         0,                    0,                0,         0,       0,         0,            0,  
             0,         0,                       0,         0,                    0,                0,         0,       0,         0,            0,  
             0,         0,                       0,         0,                    0,                0,         0,       0,         0,            0),  
         byrow = T, nrow = 10) 
}

S_fire_2(t)

S_fire_both <- function(t){
  matrix(c(
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0,
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0,
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0,
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0,
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0,
        #############################################################################################################################
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0,
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0,  
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0,  
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0,  
            0,         0,        0,         0,            0,          0,            0,         0,        0,            0),  
         byrow = T, nrow = 10) 
}


S_fire_both(t)


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

fire_return_decrease <- data.frame(Year = c(0, 40, 100, 200), Interval = c(230, 75, 30, 30)) 
# fit <- lm(log(Interval)~Year, data = fire_return_decrease)
fit <- glm(Interval~Year, data = fire_return_decrease, family = inverse.gaussian)
summary(fit)
predicted_fire_return_decrease <- data.frame(Year = c(fire_return_decrease$Year, 500))
# predicted <- exp(fit$coefficients[1] + fit$coefficients[2]* fire_return_decrease$Year)
predicted <- predict.glm(object = fit, newdata = predicted_fire_return_decrease, type = "response")
# plot(Interval~Year, data = fire_return_decrease, xlim = c(0,550))
points(predicted_fire_return_decrease$Year, predicted, pch = 19)

## NEED TO: Try to put a distribution on the average rate
interval <- function(t){
  x <- data.frame(Year = t)
  predict.glm(object = fit, newdata = x, type = "response")
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

n <- c(800, 618, 79, 65, 91,  353, 30, 600, 50, 500, 120, 600)


##############################################################################################################################################            
############################################################################################################################################## 
############################################################################################################################################## 
############################################################################################################################################## 

project <- function(projection_time, n0, reps = 100, FRI_decrease = T, fire = T){       # stochastic projection function 
  # that tracks stage based pop sizes 
  # over time for reps number of iterations
  if(length(n0) != 10)
    stop("\nPopulation size vector must be of length 10")
  
  
  results <-                                                 #Create null matrix that will hold stage
    array(0, dim = c(projection_time, length(n0) + 1, reps)) # based population sizes and year tracker
  
  mats <-
    array(0, dim = c(length(n0), length(n0), projection_time, reps))
  
  for(j in 1:reps){       # Iterate through i years (projection_time) of population growth j times (iterations)
    # Incorporate FIRE
    if(j == 1){           
      fire_tracker <- matrix(c(rep(1:reps, each = projection_time),rep(0, projection_time*reps*3)), 
                             nrow = projection_time * reps, ncol = 4, byrow = F)  
      LAI_tracker  <- matrix(c(rep(1:reps, each = projection_time),rep(0, projection_time*reps*3)), 
                             nrow = projection_time * reps, ncol = 4, byrow = F)
      lambda       <- matrix(c(rep(1:reps, each = projection_time),rep(0, projection_time*reps*2)),
                             nrow = projection_time * reps, ncol = 3, byrow = F)
    } else if(j != 1){
      fire_tracker <- fire_tracker
      LAI_tracker  <- LAI_tracker
      lambda       <- lambda
    }
    pops <- matrix(0, nrow = projection_time, ncol = length(n)) # Creates empty matrix to hold population sizes and LAI
    
    for(i in 1:projection_time){        
      #--------------------------------------------------------------------------------------------------------------------------------------------      
      t <-  i    # time counter
      
      
      ## Update LAI tracker
      LAI_tracker[j*projection_time - (projection_time)+i,2:4] <- c(i, LAI(n))
      
      ## Update parameters drawn from distributions/samples that must remain constant during each year
      cones <- No_cones(t = t, size = 1) * rcones(n)
      
      
      dispersal12 <- dispersal_25(mean_dispersal = mean_dispersal)
      dispersal21 <- dispersal_25(mean_dispersal = mean_dispersal)
      
      s1 <- si(1)
      s2 <- si(1)
      
      t1 <- ti(1)
      t2 <- ti(1)
      
      ## Set up initials for beginning of each iteration
      if (i == 1){
        n          <- n0
        tSinceFire <- c(1,1)
        pops[i,]   <- n0
        n          <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
      }else if(i != 1){
        tSinceFire <- tSinceFire +1
        
        ## Update time counter for each time step
        
        #--------------------------------------------------------------------------------------------------------------------------------------------      
        ##############################################################################################################################################            
        ##                                                                FIRE POSSIBLE                    
        ##############################################################################################################################################            
        if(fire == T){
          if(FRI_decrease == T){
            fire_current <- fire_current_year(t)  # Determine whether this iteration experiences a stand replacing fire
          }else if(FRI_decrease == F){
            fire_current <- c(rbinom(size = 1, n = 1, prob = 1/230),rbinom(size = 1, n = 1, prob = 1/230))
          }
          fire_tracker[j*projection_time - (projection_time)+i,2:4] <- c(i, fire_current)
          
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
            
            tSinceFire <- c(0, tSinceFire[2]) 
            
            caches1 <- No_caches_1_fire1
            caches2 <- No_caches_2_fire1(cones[2], t = t, size = 1, x = n)
            
            # Assuming stand replacing burn with no survival and no regeneration.
            # Most fires go out with first snow. e.g., Romme 1982
            
            mat      <- S_fire_1(t = t)
            pops[i,] <- c(t(mat%*%n + 
                              # germ1stpop1(t = t, size = 1, x = n) + 
                              germ1stpop2(caches2 = caches2, t = t, size = 1, x = n) + 
                              # germ2ndpop1(t = t, size = 1, x = n) + 
                              germ2ndpop2(t = t, size = 1, x = n))) # +
                              # No_seeds1(t = t, size = 1, x = n) +
                              # No_caches_2_fire1(cones_2 = cones_2, t = t, size = 1, x = n)))  # Defines the intermediate population size
            
            n         <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
            
          } 
          
          ## 2) FIRE IN 2 BUT NOT 1
          
          if(fire_current[1] == F & fire_current[2] == T){                  
            
            tSinceFire <- c(tSinceFire[1], 0)
            
            
            caches1 <- No_caches_1_fire2(cones_1 = cones[1], t = t, size = 1, x = n)
            caches2 <- No_caches_2_fire2
            
            # Assuming stand replacing burn with no survival and no regeneration.
            # Most fires go out with first snow. e.g., Romme 1982
            mat      <- S_fire_2(t = t)
            pops[i,] <- c(t(mat %*% n + 
                              germ1stpop1(caches1 = caches1, t = t, size = 1, x = n) + 
                              # germ1stpop2(t = t, size = 1, x = n) +
                              germ2ndpop1(t = t, size = 1, x = n))) #+
                              # germ2ndpop2(t = t, size = 1, x = n) +
                              #No_caches_1_fire2(cones_1 = cones_1, t = t, size = 1, x = n))) #+
            # No_seeds2(t = t, size = 1, x = n)))  # Defines the intermediate population size
            
            n         <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
            
          }
          
          ## 3) FIRE IN BOTH
          
          if(fire_current[1] == T & fire_current[2] == T){                  
            
            cat(paste0("Extinction in iteration ", j, " year ",t,"\n"))
            mat      <- S_fire_both(t = t)
            pops[i,] <- c(t(mat %*% n))  # Defines the intermediate population size
            
            n         <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
            
          }
          #--------------------------------------------------------------------------------------------------------------------------------------------              
          #                                                 No fire in current year in either subpopulation
          #--------------------------------------------------------------------------------------------------------------------------------------------              
          else if(fire_current[1] == F & fire_current[2] == F){
            # tSinceFire <- tSinceFire + 1
            
            caches1 <- No_caches_1_nofire(cones_1 = cones[1], cones_2 = cones[2], dispersal12, dispersal21, t = t, size = 1, x = n)
            caches2 <- No_caches_2_nofire(cones_1 = cones[1], cones_2 = cones[2], dispersal12, dispersal21, t = t, size = 1, x = n)
            
            mat <- S(t = t)
            pops[i,]  <- c(t(mat%*%n + 
                               germ1stpop1(caches1 = caches1, t = t, size = 1, x = n) + 
                               germ1stpop2(caches2 = caches2, t = t, size = 1, x = n) + 
                               germ2ndpop1(t = t, size = 1, x = n) + 
                               germ2ndpop1(t = t, size = 1, x = n)))  # Defines the intermediate population size 
            
            n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
            
          } 
        }
        ##############################################################################################################################################            
        ##                                                          FIRE = FALSE                    
        ############################################################################################################################################## 
        else if(fire == F){
          fire_tracker[j*projection_time - (projection_time)+i,2:3] <- c(i, 0)
          
          # tSinceFire <- tSinceFire + 1
          
          caches1 <- No_caches_1_nofire(cones_1, cones_2, dispersal12, dispersal21, t = t, size = 1, x = n)
          caches2 <- No_caches_2_nofire(cones_1, cones_2, dispersal12, dispersal21, t = t, size = 1, x = n)
          
          mat <- S(t = t)
          pops[i,] <- c(t(mat%*%n + 
                            germ1stpop1(caches1 = caches1, t = t, size = 1, x = n) + 
                            germ1stpop2(caches2 = caches2, t = t, size = 1, x = n) + 
                            germ2ndpop1(t = t, size = 1, x = n) + 
                            germ2ndpop2(t = t, size = 1, x = n)))  # Defines the intermediate population size
          
          n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1) 
        }
        ############################################################################################################################################## 
      }
      if(i == 1){
        lambda[j*projection_time - (projection_time) + i, 2:3] <- c(i, NA)
      }else if(i != 1){
        lambda[j*projection_time - (projection_time) + i, 2:3] <- c(i, sum(pops[i,])/sum(pops[i-1,]))
      }
    }
    
    pops <- cbind(pops, rep(1:projection_time))  # Appends matrix to keep track of time during iteration
    
    
    results[, ,j] <- pops                        # Combines iterations into a j dimensional array
    
  }
  
  pop_sizes <- plyr::adply(results, 3)           # Changes array to dataframe so easier to manipulate later
  colnames(pop_sizes) <- c("Iteration", "SEED2_1", "CS_1", "SD_1", "SAP_1", "MA_1", "SEED2_2", "CS_2", "SD_2", "SAP_2", "MA_2", "t")
  
  # fire_intervals <- cbind(iteration, intervals)  # Dataframe keeping track of fire return intervals
  
  results <- list(pop_sizes = pop_sizes, fire_tracker = fire_tracker, LAI_track = LAI_tracker, lambda = lambda) 
  
  return(results)
}

##############################################################################################################################################            
############################################################################################################################################## 
############################################################################################################################################## 
############################################################################################################################################## 


##-----------------------------------------------------------##
##                 Population projection                     ##
##-----------------------------------------------------------##
# n <- c(50000, 20000, 50000, 40000, 75000,  80000, 30000, 5000, 60000, 50000, 70000, 60000)

# n <- c(20, 618, 79, 65, 91,  353, 30, 600, 50, 500, 120, 600)
n <- c(800, 300, 90, 100, 300, 700, 
       500, 600, 50, 500, 120, 600)

projection <- project(projection_time = 100, n0 = n, reps = 100, fire = TRUE, FRI_decrease = TRUE) 

pop_sizes <- gather(projection$pop_sizes, Stage, Count, -Iteration, -t) %>%  
  # filter(., !Stage == "SEED1_1") %>%          # Pop sizes in dataframe format
  # filter(., !Stage == "SEED1_2") %>% 
  filter(., !Stage == "SEED2_1") %>%          # and excluding seed numbers (most don't think of seeds)
  filter(., !Stage == "SEED2_2") %>%
  mutate(., Population = ifelse(substr(x = .$Stage, start = 3, stop = 4) == "_1", 1, 2)) %>% 
  # mutate(., Population = ifelse(grepl("1$", Stage), 1, 2))
  group_by(., Population, Iteration, t) %>%             # as a part of the population, so presenting numbers as 
  summarise_at(., vars(Count), funs(sum)) %>% # number of living trees (i.e., post germination) is more
  ungroup(.) %>%                                # intuitive
  mutate(., Density = Count/area)

population <- gather(projection$pop_sizes, Stage, Count, -Iteration, -t) %>% 
  mutate(., Population = ifelse(substr(x = .$Stage, start = 3, stop = 4) == "_1", 1, 2)) %>% 
  # mutate(., Population = ifelse(grepl("1$", Stage), 1, 2))
  group_by(., Population, Iteration, t) %>%             # as a part of the population, so presenting numbers as 
  summarise_at(., vars(Count), funs(sum)) %>% # number of living trees (i.e., post germination) is more
  ungroup(.)

lambdas <- as.data.frame(projection$lambda)
colnames(lambdas) <- c("Iteration", "Time", "Lambda")
median(lambdas$Lambda, na.rm = T)

Fire_years <- as.data.frame(projection$fire_tracker) %>% 
  dplyr::rename(., Iteration = V1, Year = V2, Pop1 = V3, Pop2 = V4) 

test <- Fire_years %>% 
  filter(., Pop1 == 1 & Pop2 == 1) %>% 
  group_by(Iteration) %>% 
  tally() %>% 
  ungroup() 

nrow(test)/max(pop_sizes$t)

pop_sizes %>% 
  filter(., Iteration %in% sample(1:max(as.numeric(as.character(pop_sizes$Iteration))), 25, replace = F)) %>% 
  ggplot(data = ., aes(x = t, y = Density, col = Iteration)) +  # plot pop sizes for all iterations.
  geom_line(lwd = 1) +
  theme(legend.position="none") +
  theme(axis.title.x=element_text( size=18, vjust=0)) +
  theme(axis.text.x=element_text(size=18))  +
  theme(axis.title.y=element_text( size=18, vjust=2.75, face = "bold")) +
  theme(axis.text.y=element_text(size = 18)) +
  facet_wrap(~Population)


pop_sizes %>% 
  filter(., Iteration %in% sample(seq(1,max(as.numeric(as.character(pop_sizes$Iteration))), 1), 4, replace = F)) %>%
  ggplot(data = ., aes( x = Density, fill = Iteration))+
  geom_histogram()+
  facet_wrap(~ Population + Iteration, nrow = 2)+
  theme(legend.position = "none")

pop_sizes %>% 
  filter(., Iteration %in% sample(seq(1,max(as.numeric(as.character(pop_sizes$Iteration))), 1), 5, replace = F)) %>%
  ggplot(data = ., aes(x = t, y = Density, col = Iteration))+
  geom_line() +
  facet_wrap(~Iteration + Population, nrow = 2, scale = "free")+
  theme(legend.position = "none")


pop_sizes %>% 
  group_by(., Population, Iteration) %>% 
  summarise_at(., vars(Density), funs(max)) %>% 
  ungroup(.) %>% 
  summarise_at(., vars(Density), funs(median))
  
##-----------------------------------------------------------##
##                 Population projection                     ##
##                  FRI constant                             ##
##-----------------------------------------------------------##
# n <- c(50000, 20000, 50000, 40000, 75000,  80000, 30000, 5000, 60000, 50000, 70000, 60000)

# n <- c(20, 618, 79, 65, 91,  353, 30, 600, 50, 500, 120, 600)
n <- c(800, 300, 90, 100, 300, 700, 
       500, 600, 50, 500, 120, 600)

projection1 <- project(projection_time = 100, n0 = n, reps = 1000, fire = TRUE, FRI_decrease = FALSE) 

pop_sizes1 <- gather(projection1$pop_sizes, Stage, Count, -Iteration, -t) %>%  
  filter(., !Stage == "SEED1_1") %>%          # Pop sizes in dataframe format
  filter(., !Stage == "SEED1_2") %>% 
  filter(., !Stage == "SEED2_1") %>%          # and excluding seed numbers (most don't think of seeds)
  filter(., !Stage == "SEED2_2") %>%
  mutate(., Population = ifelse(substr(x = .$Stage, start = 3, stop = 4) == "_1", 1, 2)) %>% 
  # mutate(., Population = ifelse(grepl("1$", Stage), 1, 2))
  group_by(., Population, Iteration, t) %>%             # as a part of the population, so presenting numbers as 
  summarise_at(., vars(Count), funs(sum)) %>% # number of living trees (i.e., post germination) is more
  ungroup(.) %>%                                # intuitive
  mutate(., Density = Count/area)


Fire_years1 <- as.data.frame(projection1$fire_tracker) %>% 
  dplyr::rename(., Iteration = V1, Year = V2, Pop1 = V3, Pop2 = V4) 

test1 <- Fire_years1 %>% 
  filter(., Pop1 == 1 & Pop2 == 1) %>% 
  group_by(Iteration) %>% 
  tally() %>% 
  ungroup() 

nrow(test1)/max(pop_sizes1$t)

pop_sizes1 %>% 
  filter(., Iteration %in% sample(1:max(as.numeric(as.character(pop_sizes1$Iteration))), 25, replace = F)) %>% 
  ggplot(data = ., aes(x = t, y = Density, col = Iteration)) +  # plot pop sizes for all iterations.
  geom_line(lwd = 1) +
  theme(legend.position="none") +
  theme(axis.title.x=element_text( size=18, vjust=0)) +
  theme(axis.text.x=element_text(size=18))  +
  theme(axis.title.y=element_text( size=18, vjust=2.75, face = "bold")) +
  theme(axis.text.y=element_text(size = 18)) +
  facet_wrap(~Population)


pop_sizes %>% 
  filter(., Iteration %in% sample(seq(1,max(as.numeric(as.character(pop_sizes$Iteration))), 1), 4, replace = F)) %>%
  ggplot(data = ., aes( x = Density, fill = Iteration))+
  geom_histogram()+
  facet_wrap(~ Population + Iteration, nrow = 2)+
  theme(legend.position = "none")

pop_sizes %>% 
  filter(., Iteration %in% sample(seq(1,max(as.numeric(as.character(pop_sizes$Iteration))), 1), 5, replace = F)) %>%
  ggplot(data = ., aes(x = t, y = Density, col = Iteration))+
  geom_line() +
  facet_wrap(~Iteration + Population, nrow = 2)+
  theme(legend.position = "none")


pop_sizes1 %>% 
  group_by(., Population, Iteration) %>% 
  summarise_at(., vars(Density), funs(max)) %>% 
  ungroup(.) %>% 
  summarise_at(., vars(Density), funs(median))


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
##                   LAI Diagnostics                         ##
##-----------------------------------------------------------##

LAI_data <- data.frame(DBH = c(0, 0, 0, 2.05, 12.5, 37), 
                       LA = c(0, 0, alpha1(), alpha2()*2.05^alpha3(), alpha2()*12.5^alpha3(), alpha2()*37^alpha3()))

ggplot(data = LAI_data, aes(x = DBH, y = LA))+
  geom_line()+
  geom_point()

lower <- function(x){
  quantile(x, prob = 0.025)
}

upper <- function(x){
  quantile(x, prob = 0.975)
}

LAI_values <- as.data.frame(projection[[3]]) %>% 
  dplyr::select(., Iteration = V1, Time = V2, LAI = V3) %>% 
  group_by(., Time) %>% 
  summarise_at(., vars(LAI), funs(mean, lower, upper)) %>% 
  dplyr::select(.,Time, LAI = mean,lower, upper)

test <- as.data.frame(projection[[3]]) %>% 
  dplyr::select(., Iteration = V1, Time = V2, LAI = V3) %>% 
  filter(., Time %in% seq(0,100, 5))

test %>% 
  # filter(., Iteration != 7) %>% 
  ggplot(data = ., aes(x = Time, y = LAI, col = as.factor(Iteration)))+
  geom_point() +
  facet_wrap(~Iteration, nrow = 2)+
  theme(legend.position = "none")


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





