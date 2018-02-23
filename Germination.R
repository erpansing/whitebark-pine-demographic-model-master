load("survival.Rda")

library(dplyr)


survival <- survival %>%
  mutate(only.2014 = ifelse(.$germ.2013.y.n == 1, 0, .$germ.2014.y.n))

## _____________________________________________________________________________________________________
## Function to estimate proportion completely pilfered caches (# cached = # pilfered)
## _____________________________________________________________________________________________________
pilferage.frame <-function(data){
  
  characteristic <- data[,2]
  pilf <- data[,29]
  
  n = as.matrix(tapply(pilf[which(!is.na(pilf))], characteristic[which(!is.na(pilf))], function(x) sum(!is.na(x))))
  k = as.matrix(tapply(pilf, characteristic, sum, na.rm = T))
  
  kn.frame <- cbind (n,k)
  print(kn.frame)
  
  prop <- tapply(data[,29], characteristic, function(x) mean(x, na.rm = T))
  
  
  
  
  prop <- apply(kn.frame, 1, function(x)  prop.test(x[2],x[1]))
  
  
  Proportion.Pilfered <-  data.frame(LL         = c(prop$tb$conf.int[1],
                                                    prop$wc$conf.int[1]), 
                                      Prop.Pilf = c(prop$tb$estimate,
                                                    prop$wc$estimate),
                                      UL        = c(prop$tb$conf.int[2],
                                                    prop$wc$conf.int[2]))
  return(Proportion.Pilfered)
  
}


pilferage.frame(survival)[1,] 



## what if we make the assumption that pilferage is determined by the same
## processes across both study areas?


## _____________________________________________________________________________________________________
## Function to estimate proportion germinated (95%CIs) and Odds ratios (95% CIs)
## _____________________________________________________________________________________________________

germ.frame <-function(data, germ.year, subset){
  
  if(subset == "unpilfered"){
    data <- subset(data, data[,29] == 0) # select all unpilfered caches (pilf.total = 0)
                                         # completely pilfered caches = 1
  }
  else if(subset == "all"){           
    data <- data                         # select all caches created
  }
  # else if(subset == "only.2014"){
  #   data <- subset(data, data[, 30] == 1) # select all caches that ONLY germinated 
  #                                         # in 2014 (coded as 1s)
  # }
  else{
    stop("invalid subset")
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
  else if(germ.year == "only.2014"){
    year <- data[,30]
    Year <- "Germinated in only 2014"
  }
  else{
    stop("invalid germ year")
  }
  
  characteristic <- data[,2]
  

  n = as.matrix(tapply(year, characteristic, length))
  k = as.matrix(tapply(year, characteristic, sum, na.rm = T))
  
  kn.frame <- cbind (n,k)
  print(kn.frame)
  
  
  prop <- apply(kn.frame, 1, function(x)  prop.test(x[2],x[1]))
  
  
  Proportion.Germinated <- data.frame(LL        = c(prop$tb$conf.int[1],
                                                    prop$wc$conf.int[1]), 
                                      Prop.Germ = c(prop$tb$estimate,
                                                    prop$wc$estimate),
                                      UL        = c(prop$tb$conf.int[2],
                                                    prop$wc$conf.int[2]),
                                      n         = c(kn.frame[,1]))
  return(Proportion.Germinated)
}
  
germ.frame(survival, "germ.2013.y.n", subset = "all") 
germ.frame(survival, "germ.2013.y.n", subset = "unpilfered")

germ.frame(survival, "only.2014", subset = "all")
germ.frame(survival, "only.2014", subset = "unpilfered")

germ.frame(survival, "germ.2013.2014", subset = "all")
germ.frame(survival, "germ.2013.2014", subset = "unpilfered")


 
# how to incorporate the fact that probability of germination decreases
# as a function of time since caching?




## Assume that seeds WILL germinate after 2 years, otherwise they are dead
## Therefore, any seeds that weren't pilfered, or germinated after 2 years are dead

## Prob across 2 years....
## 1 = p.pilf + p.died + p.germ (2013 or 2014) 
p.pilf <- mean(survival$pilf.total)    ## Proportion pilfered in 2013

survival$germ.2013.2014[which(is.na(survival$germ.2013.2014))] <- 0
p.germ.both <- mean(survival$germ.2013.2014) # what germinated after 2 years
p.died <- 1 - p.pilf - p.germ.both  # what died in soil after 2 years

p.pilf+p.germ.both+p.died

p.died.2013 <- p.died/2  # assume half died in 2013 and half in 2014


## P.mort is pilfered plus those that were pilfered

p.mort <- p.pilf + p.died
p.germ <- germ.frame(survival, "germ.2013.y.n", "all")[1,2]
p.dorm <- germ.frame(survival, "only.2014", "all")[1,2]

p.mort + p.germ + p.dorm

n.seeds <- germ.frame(survival, "germ.2013.y.n", "all")[1,4]

(p.seed.mort   <- p.pilf + p.dead)
(mean.seed.mort <- p.seed.mort * n.seeds)
(var.seed.mort <- p.seed.mort * n.seeds * (1 - p.seed.mort))

## Method of Moments to get beta distribution parameters

# Binomial to beta
alpha <- function(n, mean, var){
  (n * mean - var)/(n*((var/mean)-mean-1)+mean)
}

test.alpha <- alpha(n = n.seeds, mean = mean.seed.mort, var = var.seed.mort)


beta <- beta <- function(alpha, mu){
  ((n-mean)*(n-(var/mean)))/ (n*((var/mean)-mean-1)+mean)
}

test.beta <- beta(n = n.seeds, mean= mean.seed.mort, var = var.seed.mort)

## Parameter estimates are negative, suggesting overdispersion...
## therefore try hypergeometric distribution


estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}



## MCMC Approach

lower <- function(x){quantile(x, (0.025))} 
upper <- function(x){quantile(x, (0.975))}

n              <- nrow(survival)
reps           <- 10000
germ.samp.dist <- NULL

for(i in 1:reps){
  sample            <- sample(survival$germ.2013.y.n, size = n, replace = T)
  germ.samp.dist[i] <- mean(sample)
}


mean.germ  <- mean(germ.samp.dist)
var.germ   <- var(germ.samp.dist)
lower.germ <- lower(germ.samp.dist)
upper.germ <- upper(germ.samp.dist)


hist(germ.samp.dist, xlab = "Proportion Germinated", main = "", 
     xlim = c(0.198,0.435))
abline(v = c(lower.germ, upper.germ), col = "blue", lwd = 1)


betaparams <- estBetaParams(mu = mean.germ, var = var.germ)
alpha      <- betaparams$alpha
beta       <- betaparams$beta

beta.dist <- rbeta(5000, shape1 = betaparams$alpha, 
                         shape2 = betaparams$beta)

hist(beta.dist, xlab = "Proportion Germinated", main = "", 
     xlim = c(0.198,0.435))
abline(v = c(lower.germ, upper.germ), col = "blue", lwd = 1)


## fitdistr

fitdistr(test1, densfun= "beta", start = list(shape1 = 200,shape2 = 500))

