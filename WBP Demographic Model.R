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
  c(d1(),      # SEED1 do not contribute to LAI
    d2(),      # SEED2 do not contribute to LAI
    d3(),      # CS do not contribute to LAI
    alpha1(),  # SD1 don't have DBH but do contribute to LAI 
    alpha2() * d5() ^ alpha3(),
    alpha2() * d6() ^ alpha3())
}


LAIb    <- function(){     # Background leaf area index. This is where competition can be incorporated...
  return(0)
}


LAI <- function(x) {       # LAI of the study area
  l <- l()
  return((t(l) %*% x)/10000 + LAIb())
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

s_SEED1 <- 0                         # Assume seeds either transition to SEED2, 
# transition to CS (first year seedling)
# or die



t1_SEED1    <- function(size = 1){        # survival probability of seeds
  rbeta(n = size,                     # Drawn from a beta to give a 
        shape1 = SEED1_survive_alpha, # probability of seeds transitioning
        shape2 = SEED1_survive_beta)  # to SEED2 stage
} 



t2_SEED1 <- function(size = 1){           # Germination probability of seeds
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


S <- function(t){  # Germination rates ([3,1] & [3,2]) and fecundity ([1,6]) are included
  # later as non-linear functions and added to the population vectors later
  matrix(c(               0,       0,        0,        0,        0,         0,
                          t1_SEED1(1),       0,        0,        0,        0,         0,
                          0,       0,        0,        0,        0,         0,
                          0,       0,  t_CS(1), si(1)[1],        0,         0,
                          0,       0,        0, ti(1)[1], si(1)[2],         0,
                          0,       0,        0,        0, ti(1)[2], si(1)[3]),
         byrow = T, nrow = 6, 
         dimnames = list(c("SEED1", "SEED2", "CS", "SD", "SAP", "MA"),
                         c("SEED1", "SEED2", "CS", "SD", "SAP", "MA"))) 
}
t <- 1
S(t) 

##-----------------------------------------------------------##
##                         FERTILITY                         ##
##-----------------------------------------------------------##

No_seeds_per_cone <- 45

No_cones <- function(t, size = 1){    #Seed production in wbp is periodic
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

e1 <- matrix(c(1,0,0,0,0,0))

No_seeds <- function(t, size = 1, x){
  No_cones(t,1) * No_seeds_per_cone  * x[6] * e1
}



##-----------------------------------------------------------##
##                       Germination                         ##
##-----------------------------------------------------------##
##-----------------------------------------------------------##
##         Define variables assumed fixed & known            ##
##-----------------------------------------------------------##
Pfind  <- 0.55 # Proportion of seeds found by nutcrackers
Pcons  <- 0.3 # Proportion of seeds consumed by nutcracers (prior to caching?)
nBirds <- 3   # No. Clark's nutcrackers in the theoretical population
SpC    <- 3.7 # No. seeds per cache

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

rCache <- function(x){
  0.73/(1+ exp((31000-SpB(x))/3000))  
}


e3 <- matrix(c(0,0,1,0,0,0))

## Define germination probability
# r2 <- function(x){
#   as.vector(((1-Pfind)*(1-Pcons)/SpC) * rCache(x) * rALS(x)) * e3
# }

germ1st <- function(t, size = 1, x){
  as.vector(No_seeds(t = t, size = 1, x = n)[1])*((1-Pcons)*rCache(x)/3.7) * (1-Pfind) * as.vector(rALS(x)) * rbeta(n = 1, shape1 = SEED1_germ_alpha, shape2= SEED1_germ_beta) * e3
}


germ2nd <- function(t, size = 1, x){
  as.vector(x[2]/3.7 * rALS(x)) * rbeta(n = 1, shape1 = SEED2_germ_alpha, shape2 = SEED2_germ_beta) * e3 
}

##-----------------------------------------------------------##
##              Function that projects pop                   ##
##              sizes and incorporates fire                  ##
##-----------------------------------------------------------##

library(plyr)

project <- function(projection_time, n0, reps = 100, fire = T){       # stochastic projection function 
  # that tracks stage based pop sizes 
  # over time for reps number of iterations
  
  results <- stochastic_matrixes <-                          #Create null matrix that will hold stage
    array(0, dim = c(projection_time, length(n0) + 1, reps)) # based population sizes
  
  for(j in 1:reps){       # Iterate through i years (projection_time) of population growth j times (iterations)
    # Incorporate FIRE
    if(j == 1){           
      intervals   <- NULL   # Create null vectors that will hold fire intervals for each iteration
      iteration   <- NULL   # Create null vectors that will show the iteration during which each fire occured
      LAI_tracker <- matrix(c(rep(1:reps, each = projection_time),rep(0, projection_time*reps*2)), 
                              nrow = projection_time * reps, ncol = 3, byrow = F)
    } else if(j != 1){
      intervals   <- intervals
      iteration   <- iteration
      LAI_tracker <- LAI_tracker
    }
    if(fire == TRUE){
    interval <- rnbinom(n = 4, size = 1, mu = 230) %>%  # Select fire years from a negative binomial distribution (assumes mean fire return interval = 230 yrs)
      cumsum(.)           # representing the waiting times btwn fires
    # cumsum gives the time = t years during which fires should occur in projection
    interval <- interval[-which(interval > projection_time)]  #trim to only include those within the projection time
    
    ## Creates dataframe result that tracks the fire years for each iteration.
    intervals <- append(intervals, interval, after = length(intervals))
    
    } else if(fire == FALSE){
      intervals <- NULL
      interval <- NULL
    }
    iteration <- append(iteration, rep(j, length(interval)), 
                        after = length(iteration))
    
    
    pops <- matrix(0, nrow = projection_time, ncol = length(n)) # Creates empty matrix to hold population sizes and LAI
    
    for(i in 1:projection_time){           # get population
        
      if (i == 1){
        n <- n0
      }else if(i != 1){
        n <- n
      }
      
      t <-  i                              # time counter
      tSinceFire <- ifelse(i == 1, 1, tSinceFire)
      fire_current <- ifelse(t %in% interval, TRUE, FALSE)  # Determine whether this iteration experiences a stand replacing fire
      
      # LAI_tracker_each_iteraton <- matrix(0, nrow = projection_time, ncol = 2)
      LAI_tracker[j*projection_time - (projection_time)+i,2:3] <- c(i, LAI(n))
      
      # Now, multiple possibilities
      # 1) There's a fire. This kills the population. Assumes stand replacing burn that impacts entire population
      #    And that no regeneration occurs the year of the fire
      # 2) It's a year after a fire, in which case we assume input from an outside population (i.e., system
      #    isn't entirely closed) and that the outside population is on a similar masting schedule as our 
      #    population
      # 3) There's no fire and it's >1 year after fire. In this case, there are no modifications and the 
      #    system proceeds as normal.
      
      if(fire_current == T){                  
        
        tSinceFire <- 0 
        
        pops[i,] <- c(0, 0, 0, 0, 0, 0)    # Assuming stand replacing burn with no survival and no regeneration. 
        n <- pops[i,]                      # Most fires go out with first snow. e.g., Romme 1982
        
      } else if(fire_current == F & tSinceFire == 1 & t != 1){
        
        tSinceFire <- tSinceFire + 1
        
        pops[i,] <- t(No_seeds(size = 1, t = t, x = n))
        n <- pops[i,]
        
        
      } else if((fire_current == F & tSinceFire != 1) |
                (fire_current == F & tSinceFire == 1 & t == 1)){
        tSinceFire <- tSinceFire + 1
        
        mat <- S(t = t)
         pops[i,]  <- t(mat%*%n + germ1st(t = t, size = 1, x = n) + germ2nd(t = t, size = 1, x = n) + No_seeds(t = t, size = 1, x = n))  # Defines the intermediate population size 

        n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
      
        
        # } else if(fire == F & tSinceFire == 1 & t == 2){
        #   tSinceFire <- tSinceFire + 1
        #   
        #   pops[i,] <- c(500000, 0, 0, 0, 0, 0)
        #   n <- as.matrix(pops[i,], nrow = length(pops[i,]), ncol = 1)
        }
          
    }
    
    pops <- cbind(pops, rep(1:projection_time))  # Appends matrix to keep track of time during iteration
    
    
    results[, ,j] <- pops                        # Combines iterations into a j dimensional array
    
  }
  
  pop_sizes <- plyr::adply(results, 3)           # Changes array to dataframe so easier to manipulate later
  colnames(pop_sizes) <- c("Iteration", "SEED1", "SEED2", "CS", "SD", "SAP", "MA", "t")
  
  fire_intervals <- cbind(iteration, intervals)  # Dataframe keeping track of fire return intervals
  
  results <- list(pop_sizes = pop_sizes, fire_intervals = fire_intervals, LAI_track = LAI_tracker)
  
  return(results)
}


##-----------------------------------------------------------##
##                 Population projection                     ##
##-----------------------------------------------------------##

n <- c(62, 580 + 38, 79, 65, 91,  353) # Arbitrary starting pop size vectors

projection <- project(projection_time = 500, n0 = n, reps = 10000, fire = TRUE) 

pops <- gather(projection$pop_sizes, Stage, Count, - Iteration, - t)

pop_sizes <- gather(projection$pop_sizes, Stage, Count, -Iteration, -t) %>%  
  filter(., !Stage == "SEED1") %>%          # Pop sizes in dataframe format
  filter(., !Stage == "SEED2") %>%          # and excluding seed numbers (most don't think of seeds)
  group_by(., Iteration, t) %>%             # as a part of the population, so presenting numbers as 
  summarise_at(., vars(Count), funs(sum)) %>% # number of living trees (i.e., post germination) is more
  ungroup(.)   %>%                              # intuitive
  mutate(., Density = Count/(10000))


fire_intervals <- as.data.frame(projection$fire_intervals)

ggplot(data = pop_sizes, aes(x = t, y = Count, col = Iteration)) +  # plot pop sizes for all iterations.
  geom_line(lwd = 1) +
  theme(axis.title.x=element_text( size=18, vjust=0)) +
  theme(axis.text.x=element_text(size=18))  +
  theme(axis.title.y=element_text( size=18, vjust=2.75, face = "bold")) +
  theme(axis.text.y=element_text(size = 18)) +
  labs(x = "Years", y = "no") + 
  theme(legend.position="none")+
  geom_hline(yintercept = 715)


## Plot of projection iteration 1

projection1 <- pop_sizes %>% 
  filter(., Iteration == 1)

ggplot(data = projection1, aes(x = t, y = Count)) +  # plot pop sizes for all iterations.
  geom_line(lwd = 1) +
  theme(legend.position="none") +
  theme(axis.title.x=element_text( size=18, vjust=0)) +
  theme(axis.text.x=element_text(size=18))  +
  theme(axis.title.y=element_text( size=18, vjust=2.75, face = "bold")) +
  theme(axis.text.y=element_text(size = 18))


(mean_density <- pop_sizes %>% 
  group_by(., Iteration) %>% 
  filter(., Count == max(Count)) %>% 
  ungroup(.) %>% 
  summarise(., mean(Density)))

(max_density <- pop_sizes %>% 
  group_by(., Iteration) %>% 
  filter(., Count == max(Count)) %>% 
  ungroup(.) %>% 
  summarise(., max(Density)))

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



LAI_values <- as.data.frame(projection[[3]]) %>% 
  dplyr::select(., Iteration = V1, Time = V2, LAI = V3) %>% 
  group_by(., Time) %>% 
  summarise_at(., vars(LAI), funs(mean, min, max)) %>% 
  dplyr::select(.,Time, LAI = mean, min, max)


ggplot(data = LAI_values, aes(x = Time, y = LAI))+
  geom_line() +
  geom_ribbon(data = LAI_values, aes(ymin = min, ymax = max), alpha = 0.5)


# ggplot(predictions, aes(x = time, y = prob)) +
#   # geom_point() +
#   geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5, fill = "#ef8a62") +
#   geom_line(size = 0.5) +
#   xlab("Year") +
#   ylab("Annual survival rate") +
#   scale_y_continuous(limits = c(0.75,1)) +
#   scale_x_continuous(limits = c(1990, 2017))


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



## Determining quasi extinction values

# when probability of observing a nutcracker dips below 10%
# solving for x gives the cone density index (ln(cones/ha))^2 when prob dips belows 10%

# ln(0.1/0.9) = 0.03883x - 1.5165

(cdi <- (log(0.18/(1-0.18)) + 1.5165)/0.03883)

(cdext <- exp(sqrt(cdi)))

odds <- exp(0.03*0 - 1.5165)
p <- odds/(1+odds)

# suggests that cone density that dips below 


quasiextinct <- stoch.quasi.ext(A, n0, Nx = 500, tmax = 100, maxruns = 10, nreps = 5000, 
                prob = NULL, sumweight = NULL, verbose=TRUE)

