rm(list = ls())
source(file = "GRIN parameter estimates.R")
library(popbio)
library(plyr)

##-----------------------------------------------------------##
##                                                           ##
##                 Define Functions                          ##
##             and script-specific vars                      ##
##                                                           ##
##-----------------------------------------------------------##


beta_parameters <- data.frame(alphas = c(SEED1_survive.alpha, SEED1_germ.alpha,
                                    SEED2_germ.alpha, CS_survive.alpha),
                         betas  = c(SEED1_survive.beta, SEED1_germ.beta,
                                    SEED2_germ.beta, CS_survive.beta))
gamma_parameters <- data.frame(alpha = cone_alpha, theta = cone_theta)

n_beta_params <- nrow(beta_parameters)
reps <- 10000

random_parameters_matrix <- matrix(0, nrow = reps, ncol = 6)
colnames(random_parameters_matrix) <- c("SEED1_survive", "SEED1_germ",
                                        "SEED2_germ", "CS_survive",
                                        "MA_survival", "Fecund")


for(i in 1:n_beta_params){
  random_parameters_matrix[,i] <- rbeta(reps, 
                                       shape1 = beta_parameters$alphas[i],
                                       shape2 = beta_parameters$betas[i])
}


random_parameters_matrix[,6] <- rgamma(n = reps, shape = cone_alpha, scale = cone_theta) * 45


random_parameters_matrix[ , 5] <- 1



random_parameters_matrix<- as.data.frame(random_parameters_matrix)


stochastic_matrixes <- array(0, dim = c(4,4,reps))

for(i in 1:reps){
  # create empty matrix
  mat <- matrix(data=0, nrow=4, ncol=4)
  
  mat[1,4] <- random_parameters_matrix$Fecund[i]
  mat[2,1] <- random_parameters_matrix$SEED1_survive[i]
  mat[3,1] <- random_parameters_matrix$SEED1_germ[i]
  mat[3,2] <- random_parameters_matrix$SEED2_germ[i]
  mat[4,3] <- random_parameters_matrix$CS_survive[i]
  mat[4,4] <- random_parameters_matrix$MA_survival[i]
  
  stochastic_matrixes[, , i] <- mat
}

lambdas <- NULL
el_SEED1_surivive <- NULL
el_SEED1_germ <- NULL
el_SEED2_germ <- NULL
el_CS_survive <- NULL
el_MA_survive <- NULL
el_Fecund <- NULL

for(i in 1:reps){
  
  eigenal <- eigen.analysis(stochastic_matrixes[,,i])
  
  lambdas[i] <- eigenal$lambda1
  el_SEED1_surivive[i] <- eigenal$elasticities[2,1]
  el_SEED1_germ[i] <- eigenal$elasticities[3,1]
  el_SEED2_germ[i] <- eigenal$elasticities[3,2]
  el_CS_survive[i] <- eigenal$elasticities[4,3]
  el_MA_survive[i] <- eigenal$elasticities[4,4]
  el_Fecund[i] <- eigenal$elasticities[1,4]
}

stochastic_matrixes <- alply(stochastic_matrixes, 3)

stochastic_elasticities <- stoch.sens(test_list)


