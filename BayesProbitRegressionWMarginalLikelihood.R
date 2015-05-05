# R Code illustrating Marginal Likelihood Computation from Gibbs Sampler of a Probit Model as 
#  described in Chib (1995) and replicating Table 2 (row 2 to 8, col. 1/3 ) in the same paper
# 
# 
# References: 
# 	Chib, Siddhartha. 1995. “Marginal Likelihood from the Gibbs Output.” 
# 	Journal of the American Statistical Association 90(432):1313–1321.


rm(list=ls())

library(DPpackage)
data(nodal)
library(MASS)
library(msm)
library(mnormt)
library(coda)

probit.lik <- function(beta,y,X){
	zeta <- as.vector(X %*% beta)
	y.logic <- as.logical(y)
	n<-length(y)
	lgLik <- rep(NA,n)
	lgLik[y.logic] <- pnorm(zeta[y.logic], log.p = TRUE)
	lgLik[!y.logic] <- pnorm(zeta[!y.logic], lower.tail = FALSE, log.p=TRUE)
	logl<- sum(lgLik)
	return(logl)
	}

probit.bayes <- function(formula, data, burnin, mcmc, b0, B0, seed=10){

	set.seed(seed)
	X <- model.matrix(formula, model.frame(formula, data) )
	y <- data[, all.vars(formula)[1]]
	K <- ncol(X)
	n <- nrow(X)

	b0 <- rep(b0,K)
	B0 <- diag(K)*B0

	iter <- mcmc + burnin

	# Holding bins
	bhat <- bmu <- matrix(NA, iter, K)
	ystar <- matrix(NA, iter, n)
	lower <- rep(-Inf, n)
	lower[y==1] <- 0
	upper <- rep(Inf, n)
	upper[y==0] <- 0
	ystar[1,] <- rnorm(n)
	bhat[1, ] <- coef(glm(formula, data, family=binomial("probit")))

	# Compute constant 
	bsd <- solve(solve(B0) + t(X) %*% X) 

	# Sample 
	for(i in 2:iter){
		ystarmu <- X %*% bhat[(i-1), ]
		ystar[i, ] <- rtnorm(n, ystarmu, lower=lower, upper=upper)
		bmu[i, ] <- solve( solve(B0) + t(X) %*% X ) %*% (solve(B0) %*% b0 + t(X) %*% ystar[i, ])
		bhat[i,] <- mvrnorm(1, bmu[i,], bsd)
		}

	# Discard burnin 
	bmu <- bmu[(burnin+1):iter, ]
	bhat <- bhat[(burnin+1):iter, ]
	bmuhat <- apply(bhat, 2, mean)

	# > Calculate Marginal Likelihood using Chib 1995 #

	# 1. evaluate prior at posterior mean
	logprior <- log(dmnorm(bmuhat, mean=b0, varcov=B0 ))

	# 2. evaluate likelihood at posterior mean
	loglik <- probit.lik(bmuhat,y,X ) 

	# 3. evaluate posterior at posterior mean using its full conditional 
	logplugin <- rep(NA, nrow(bmu) )
	for(i in 1:nrow(bmu)){
		logplugin[i] <- dmnorm(bmuhat, mean=bmu[i,], varcov=bsd, log=FALSE )
		}

	# this is eq. 15 in the paper 
	chib95 <- logprior + loglik  - log(mean(logplugin))

	out <- list(bhat=as.mcmc(bhat), chib95=chib95, logprior=logprior, loglik=loglik)

	return(out)
	}

# Replicate Table 2 on p. 1318

mf2 <- ssln~age
mf3 <- ssln~log(acid)
mf4 <- ssln~xray
mf5 <- ssln~size
mf6 <- ssln~grade
mf7 <- ssln~log(acid) + size
mf8 <- ssln~log(acid) + xray + size
mf9 <- ssln~log(acid) + xray + size  + grade

b0 <- 0.75 
B0 <- 5^2 

m2_R <- probit.bayes(mf2, data=nodal, burnin=500, mcmc=5000, b0=b0, B0=B0 ) 
m3_R <- probit.bayes(mf3, data=nodal, burnin=500, mcmc=5000, b0=b0, B0=B0 ) 
m4_R <- probit.bayes(mf4, data=nodal, burnin=500, mcmc=5000, b0=b0, B0=B0 ) 
m5_R <- probit.bayes(mf5, data=nodal, burnin=500, mcmc=5000, b0=b0, B0=B0 ) 
m6_R <- probit.bayes(mf6, data=nodal, burnin=500, mcmc=5000, b0=b0, B0=B0 ) 
m7_R <- probit.bayes(mf7, data=nodal, burnin=500, mcmc=5000, b0=b0, B0=B0 ) 
m8_R <- probit.bayes(mf8, data=nodal, burnin=500, mcmc=5000, b0=b0, B0=B0 ) 
m9_R <- probit.bayes(mf9, data=nodal, burnin=500, mcmc=5000, b0=b0, B0=B0 ) 

# Results 

matrix(
	c(	m2_R$loglik, m2_R$chib95 , 
		m3_R$loglik, m3_R$chib95 , 
		m4_R$loglik, m4_R$chib95 , 
		m5_R$loglik, m5_R$chib95 , 
		m6_R$loglik, m6_R$chib95 ,
		m7_R$loglik, m7_R$chib95 ,
		m8_R$loglik, m8_R$chib95 ,
		m9_R$loglik, m9_R$chib95   ), 8, 3, byrow=TRUE)
