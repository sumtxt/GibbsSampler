rm(list=ls())
set.seed(12)
library(coda)


# ABC rejection sampler and ABC MCMC demonstration 
#  for simple linear regression with known variance
#
#	May 27, 2015
#
# see: 
# 	Approximate Bayesian Computation in Evolution and Ecology
#	Mark A. Beaumont in the Annu. Rev. Ecol. Evol. Syst. 2010. 41:379–406
#
# 	but also: Marjoram, Paul, John Molitor, Vincent Plagnol & Simon Tavare. 2003. 
#	“Markov chain Monte Carlo without likelihoods.” Proceedings of the National
#	 Academy of Sciences 100(26):15324–15328.
#




# Simulate data #
#############

N <- 100
X <- cbind(1, runif(N,-1,1))
truebeta <- c(-1,-0.45)
sigma <- 1
y <- X %*% truebeta + rnorm(N,0,sigma)

m1 <- summary(lm(y~X-1))


# Functions #
draw_pseudo_sample <- function(betadraw,X,N,sigma){
	ydraw <- X %*% betadraw + rnorm(N,0,sigma)
	return(ydraw)
	}

distance <- function(y,draw){
	return( sum( (y - draw)^2 ) )
	}

prior <- function(x) dnorm(x)



#  ABC rejection sampler #
###################

# see e.g. p. 381 (Beaumont)

sim <- 1000
epsilon <- 170

betastore <- matrix(NA, sim, 2)

for(i in 1:sim){
	repeat {
		betadraw <- rnorm(2,0,2) 
		ydraw <- draw_pseudo_sample(betadraw,X,N,sigma)
		if ( distance(y,ydraw) < epsilon ) break
		}
	cat(".")
	betastore[i,] <- betadraw
	}

apply(betastore, 2, function(x) c(mean(x),sd(x)))
coef(m1)

plot(as.mcmc(betastore))


# ABC Metropolis-Hastings sampler  #
############################

# see e.g. p. 387 (Beaumont)
# however, the ABC iterations are not counted towards the MCMC 
# iterations. 

mcmc <- 1000
burnin <- 500
epsilon <- 150

startvalue <- c(0,0) 

sample = array(NA, dim = c(mcmc+1,2))
sample[1,] <- startvalue

accept = 0
probab <- rep(NA,mcmc)

for (i in 1:mcmc){
	repeat {
		candidate <- rnorm(2,0,2)
		ydraw <- draw_pseudo_sample(candidate,X,N,sigma)
		if ( distance(y,ydraw) < epsilon ) break
		}
	ratio <- prior(candidate)/prior(sample[i,])
	probab[i] <- min(ratio, 1)
	if (runif(1) < probab[i]){
		sample[i+1,] <- candidate
		accept <- accept + 1
	}else{
		sample[i+1,] = sample[i,]
		}
	cat(".")
}

accept/mcmc
apply(sample, 2, function(x) c(mean(x),sd(x)))
coef(m1)

plot(as.mcmc(sample))




