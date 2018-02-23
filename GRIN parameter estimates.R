rm(list = ls())
load("survival.Rda")

##-----------------------------------------------------------##
##                                                           ##
##                 Define Functions                          ##
##             and script-specific vars                      ##
##                                                           ##
##-----------------------------------------------------------##


lower <- function(x){quantile(x, (0.025))} 

upper <- function(x){quantile(x, (0.975))}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

get_betas <- function(data, reps, size){
  output_vector <- NULL
  for(i in 1:reps){
    sample  <- sample(data, size = size, replace = T)
    output_vector[i] <- mean(sample)
  }
  
  mean  <- mean(output_vector)
  var  <- var(output_vector)
  lower <- lower(output_vector)
  upper <- upper(output_vector)
  
  betaparams <- estBetaParams(mu  = mean, 
                              var = var)
  alpha      <- betaparams$alpha
  beta       <- betaparams$beta
  
  
  result <- list(samp.dist  = output_vector, mean.samp.dist = mean,
                 var.samp.dist = var, lower.samp.dist = lower,
                 upper.samp.dist = upper, alpha = alpha, beta = beta)
  return(result)
}


get_gammas <- function(mu, var){
  alpha <- mu^2/var
  theta <- var/mu
  beta <- 1/theta
  
  result <- list(alpha = alpha, 
                 theta = theta,
                 beta = beta)
  return(result)
}




reps <- 10000
n <- nrow(survival)

##-----------------------------------------------------------##
##                                                           ##
##              First Year Seed Probabilities                ##
##                                                           ##
##-----------------------------------------------------------##

# 3 possible outcomes for a first year seed Pansing et al. 2017
# 1) Seed germinated 

SEED1_germ <- survival$germ.2013.y.n

# 2) Seed survives and becomes 2nd year seed
#    Bc cach is the sampling unit, germintaion in Year 1 means no survival

survival$intact.2013.y.n <- ifelse(SEED1_germ == 1, 0, survival$intact.2013.y.n)

SEED1_survive <- survival$intact.2013.y.n


# 3) Seed is pilfered, which we assume = death. 
#    OR
#    Seed dies (e.g., mold, inviable embryo, etc.)
# Here, we assume that all seeds that were not pilfered have the 
# ability to germinate in year 2. In other words, seed death in
# year 1 can ONLY result from cache pilferage

SEED1_die <- ifelse(SEED1_germ == 1 | SEED1_survive == 1, 0 ,1)

sum(mean(SEED1_germ), mean(SEED1_survive), mean(SEED1_die))

## ___________________________________________________________##
##                                                            ##
##     Create sampling distribution for SEED1 Germination     ##
##                                                            ##
## ___________________________________________________________##



SEED1_germ_params <- get_betas(data = SEED1_germ, reps = reps,
                               size = n)

SEED1_germ_alpha <- SEED1_germ_params$alpha
SEED1_germ_beta <- SEED1_germ_params$beta

# hist(SEED1_germ_params$samp.dist, 
#      xlab = "Proportion Germinated", main = "")
# abline(v = c(SEED1_germ_params$lower.samp.dist, 
#              SEED1_germ_params$upper.samp.dist), col = "blue", lwd = 1)


## Get sampling dist from beta parameters and plot

beta.dist <- rbeta(reps, shape1 = SEED1_germ_alpha, 
                   shape2 = SEED1_germ_beta)

# hist(beta.dist, xlab = "Proportion Germinated", main = "")
# abline(v = c(SEED1_germ_params$lower.samp.dist, 
#              SEED1_germ_params$upper.samp.dist), col = "blue", lwd = 1)

rm(SEED1_germ_params, beta.dist)

## ___________________________________________________________##
##                                                            ##
##     Create sampling distribution for SEED1 Survival        ##
##                                                            ##
## ___________________________________________________________##

SEED1_survival_params <- get_betas(data = SEED1_survive, size = n,
                                   reps = reps)

SEED1_survive_alpha      <- SEED1_survival_params$alpha
SEED1_survive_beta       <- SEED1_survival_params$beta


# hist(SEED1_survival_params$samp.dist, 
#      xlab = "Proportion Germinated", main = "")
# abline(v = c(SEED1_survival_params$lower.samp.dist, 
#              SEED1_survival_params$upper.samp.dist), col = "blue", lwd = 1)

## Get sampling dist from beta parameters and plot


beta.dist <- rbeta(reps, shape1 = SEED1_survival_params$alpha, 
                   shape2 = SEED1_survival_params$beta)

# hist(beta.dist, xlab = "Proportion Germinated", main = "")
# abline(v = c(SEED1_survival_params$lower.samp.dist, 
#              SEED1_survival_params$upper.samp.dist), col = "blue", lwd = 1)


rm(SEED1_germ, SEED1_survive, SEED1_die, SEED1_survival_params, beta.dist)

##-----------------------------------------------------------##
##                                                           ##
##              Second Year Seed Probabilities               ##
##                                                           ##
##-----------------------------------------------------------##

# Assuming 2 possibilities for 2nd year seeds:
# 1) Germination
# 2) Death (i.e., no germination). 
#
# We are assuming that if seeds do not germinate in year 1, they 
# will not. 

SEED2 <- survival %>%
  filter(., intact.2013.y.n == 1) 

SEED2$germ.2014.y.n[which(is.na(SEED2$germ.2014.y.n))] <- 0

SEED2_germ <- SEED2$germ.2014.y.n

## ___________________________________________________________##
##                                                            ##
##     Create sampling distribution for SEED2 Germination     ##
##                                                            ##
## ___________________________________________________________##

SEED_2_germ_params <- get_betas(data = SEED2_germ, reps = reps, size = n)

SEED2_germ_alpha      <- SEED_2_germ_params$alpha
SEED2_germ_beta       <- SEED_2_germ_params$beta


# hist(SEED_2_germ_params$samp.dist, 
#      xlab = "Proportion Germinated", main = "")
# abline(v = c(SEED_2_germ_params$lower.samp.dist, 
#              SEED_2_germ_params$upper.samp.dist), col = "blue", lwd = 1)

## Get sampling dist from beta parameters and plot


beta.dist <- rbeta(reps, shape1 = SEED_2_germ_params$alpha, 
                   shape2 = SEED_2_germ_params$beta)

# hist(beta.dist, xlab = "Proportion Germinated", main = "")
# abline(v = c(SEED_2_germ_params$lower.samp.dist, 
#              SEED_2_germ_params$upper.samp.dist), col = "blue", lwd = 1)


rm(SEED2, SEED2_germ, SEED_2_germ_params, beta.dist)

##-----------------------------------------------------------##
##                                                           ##
##                CS Survival Probabilities                  ##
##                                                           ##
##-----------------------------------------------------------##

CS <- survival %>%
  filter(., germ.2013.y.n == 1) %>%
  select(., tag, germ.2013.y.n, no.living.2013, living.2014.y.n, survive.13.14)

# Include those seedlings that germinated and died before being located.
CS$survive.13.14[CS$no.living.2013 == 0] <- 0

# Remove those caches that we couldn't find in 2014
CS <- CS %>%
  filter(., !is.na(survive.13.14))

# Define CS_survive vector

CS_survive <- CS$survive.13.14

# Redefine n bc # trials changed

n <- length(CS_survive)

CS_survive_params <- get_betas(data = CS_survive, size = n, reps = reps)

CS_survive_alpha      <- CS_survive_params$alpha
CS_survive_beta       <- CS_survive_params$beta


# hist(CS_survive_params$samp.dist, 
#      xlab = "Proportion Germinated", main = "")
# abline(v = c(CS_survive_params$lower.samp.dist, 
#              CS_survive_params$upper.samp.dist), col = "blue", lwd = 1)

## Get sampling dist from beta parameters and plot


beta.dist <- rbeta(reps, shape1 = CS_survive_params$alpha, 
                   shape2 = CS_survive_params$beta)

# hist(beta.dist, xlab = "Proportion Germinated", main = "")
# abline(v = c(CS_survive_params$lower.samp.dist, 
#              CS_survive_params$upper.samp.dist), col = "blue", lwd = 1)


rm(CS, CS_survive, CS_survive_params, beta.dist)

##-----------------------------------------------------------##
##                                                           ##
##             Seedling Survival Probabilities               ##
##                                                           ##
##-----------------------------------------------------------##



source("/Users/elizabethpansing/Box Sync/Yellowstone/88-Fires-Analysis/WBP Survival Estimates RMark Nest.R")


SD_beta <- estBetaParams(mu = SD_survival_mean, var = SD_survival_var)

SD_survive_alpha <- SD_beta$alpha
SD_survive_beta  <- SD_beta$beta

rm(SD_beta)


##-----------------------------------------------------------##
##                                                           ##
##              Sapling Survival Probabilities               ##
##                                                           ##
##-----------------------------------------------------------##



SAP_survival <- 0.8


rm(survival, n, reps, estBetaParams, get_betas, lower, upper)


##-----------------------------------------------------------##
##                                                           ##
##                        Fecundity                          ##
##                                                           ##
##-----------------------------------------------------------##

# Cone count estimates derived from IGBST cone count data from
# the GYE. https://www.usgs.gov/centers/norock/science/igbst-whitebark-pine-cone-production-annual-summaries?qt-science_center_objects=1#qt-science_center_objects



cones_1980_2016 <- c(25.69, 13.23, 16.98, 17.4, 6.43, 27.2, 1.37, 2.54, 2.40, 48.8,
                     1.54, 15.5, 15.38, 10.19, 2.03, 2.73, 25.05, 4.55, 8.4, 39.5,
                     5.7, 25.5, 2.4, 28.5, 6.9, 16.8, 34.4, 14.9, 8.6, 46.5, 
                     5.2, 19.8, 33, 5.2, 20.05, 15.89, 35.9)

mu.cones <- mean(cones_1980_2016)
var.cones <- var(cones_1980_2016)

cone_params <- get_gammas(mu = mu.cones, var = var.cones)

cone_alpha <- cone_params$alpha
cone_theta <- cone_params$theta

rm(cones_1980_2016, mu.cones, var.cones, cone_params)

##-----------------------------------------------------------##
##                                                           ##
##                      Fire Params                          ##
##                                                           ##
##-----------------------------------------------------------##

LarsonFireIntervals <- c(250, 350, 100)

mean_interval <- mean(LarsonFireIntervals)
var_interval <- var(LarsonFireIntervals)


fire_gamma <- get_gammas(mu = mean_interval, var = var_interval)

fire_alpha <- fire_gamma$alpha
fire_beta <- fire_gamma$beta

rm(LarsonFireIntervals, mean_interval, var_interval, fire_gamma, get_gammas)