# Conversion of normal to gamma using MOM


lambda <- function(mu){    #double check these.... might be using sd where it should be var
  mu1 <- mu
  mu2 <- mu^2 + mu*(1-mu)
  lambda <- mu1/(mu2 - mu1^2)
  return(lambda)
}

lambda(p.seed.mort)

alpha <- function(lambda, mu){
  lambda*mu
}

alpha(1.98913, p.seed.mort)


test <- rgamma(10000, shape = lambda(p.seed.mort), scale = alpha(lambda(p.seed.mort), p.seed.mort))

hist(test)