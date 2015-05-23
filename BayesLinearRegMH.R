rm(list=ls())
library(MASS)
library(MCMCpack)
library(adaptMCMC)

# Linear Regression with Metropolis-Hasting algorithm,  
# multivariate normal random walk proposal
# 
# see:
# Chib/Greenberg. 1995. Understanding the Metropolis-Hastings Algorithm. The American Statistician. 
#
# MM, 30. June 2012
#
#


# Data 
######

alpha <- 2
beta  <- 1
sd <- 2
n <- 100

X <-  runif(n, -2,2)
y <-  alpha + beta * X + rnorm(n,0,sd=sd)


# Likelihood 
#######

lin.llik <- function(param,y,X){
	alphahat = param[1]
	betahat = param[2]
	sigma2 = abs(param[3])
	n <- length(y)
	pred <- y - (alphahat + betahat * X)
	logl <- -0.5*n*log(2*pi)-0.5*n*log(sigma2)-((t(pred)%*%pred)/(2*sigma2))
	return(logl)	
}


# Prior
#######

prior <- function(param){
	alpha = param[1]
	beta = param[2]
	sigma2 = param[3]
	p.alpha = dnorm(alpha, mean=0, sd = 4, log = T)
	p.beta = dnorm(beta, mean=0, sd = 4, log = T)
	p.sigma2 = dunif(sigma2, min=0, max=30, log = T)
	return(p.alpha+p.beta+p.sigma2)
}

# Proposal generator (symmetric, multivariate random walk)
##########

proposal <- function(param){
	candid <- param + rnorm(3,0,0.1)
	return(candid)
}

# RUN the random walk Metropolis-Hastings Algorithm 
#####

mcmc <- 50000
burnin <- 5000
# extraordinary bad stating values
startvalue <- c(4,0,10) 

sample = array(NA, dim = c(mcmc+1,3))
sample[1,] <- startvalue

accept = 0
probab <- rep(NA,mcmc)
for (i in 1:mcmc){
	candidate = proposal(sample[i,])
	# The probability to move, alpha(x,y) is defined as a ratio which includes the proposal density. Since the proposal density is symmetric, 
	# it can be ignored (see remark 4 on p. 329 Chib/Greenberg). For numerical reasons and the very fact that lin.llik / prior are defined on the log scale, 
	# the acceptance ratio is calculated on the log scale as well
	ratio <- exp(lin.llik(candidate,y,X)+ prior(candidate) - lin.llik(sample[i,],y,X)- prior(sample[i,]))
	probab[i] <- min(ratio, 1)
	if (runif(1) < probab[i]){
		sample[i+1,] <- candidate
		accept <- accept + 1
	}else{
		sample[i+1,] = sample[i,]
	}
}

posterior <- sample[-(1:burnin),]
posterior[,3] <- sqrt(posterior[,3])
accept/mcmc
colMeans(posterior)

par(mfrow=c(1,3))
plot(posterior[,1], type="l")
plot(posterior[,2], type="l")
plot(posterior[,3], type="l")

plot(posterior[1:1000,1],posterior[1:1000,2], type="l", main="Posterior (1:1000 draws)")


# CHECKs

# tailored MCMCpack Metropolis 
logposterior <- function(param, X, y) { lin.llik(param,y,X)+ prior(param) } 
post.samp <- MCMCmetrop1R(logposterior, theta.init=startvalue,
                          X=X, y=y,
                          thin=1, mcmc=mcmc, burnin=burnin,
                          tune=c(1.5, 1.5, 1.5),
                          verbose=0, logfun=TRUE)
posterior2 <- as.matrix(post.samp)
posterior2[,3] <- sqrt(posterior2[,3])
colMeans(posterior2)

# Robut Adaptive Rejection Metropolis 
post.samp <- adaptMCMC::MCMC(logposterior, n=mcmc, init=startvalue, adapt=FALSE, y=y, X=X)
posterior3 <- as.matrix(post.samp$samples)
posterior3 <- posterior3[-(1:burnin),]
posterior3[,3] <- sqrt(posterior3[,3])
colMeans(posterior3)




