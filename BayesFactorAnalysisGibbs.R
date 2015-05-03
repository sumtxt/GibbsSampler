# Bayesian Factor Analysis in R #
#################################

#
# This script implements the Gibbs sampler 
# suggested in Greenberg 2009, p. 142 (Algorithm 8.6)
#
#
# n, N subjects / units 
# k, K indicators / items 
# d, D latent dimensions 

#
# Observed: 
# *********
# 
# X, indicator matrix, (N x K)
# 
# Unobserved:
# ***********
# 
# phi, latent factors, list of (N x D) matrices 
# lambda, factor loadings, list of (K x D) matrices 
# psi, precision matrix, list of (K x K) matrices 
#
#

rm(list=ls())
library(MCMCpack)
library(MASS)
library(pscl)
library(stringr)


# Get the data
data(swiss)
X <- swiss[,c("Agriculture", "Examination", "Education", "Catholic", "Infant.Mortality")]

# GET THE DEFAULT SOLUTION WITH MCMCPACK #
########################################## 

posterior <- MCMCfactanal(~Agriculture+Examination+Education+Catholic
                    +Infant.Mortality, factors=2,
                    lambda.constraints=list(Examination=list(1,"+"),
                       Examination=list(2,"-"), Education=c(2,0),
                       Infant.Mortality=c(1,0)), 
                    verbose=0, store.scores=TRUE, a0=1, b0=0.15,
                    data=swiss, burnin=5000, mcmc=50000, thin=20)

phi1_mcmc0 <-  posterior[,str_detect(colnames(posterior), "phi.*1")]
phi2_mcmc0 <-  posterior[,str_detect(colnames(posterior), "phi.*2")]

phi1_mu0  <- apply(phi1_mcmc0, 2, mean)
phi1_low0 <- apply(phi1_mcmc0, 2, quantile, probs=c(0.975), na.rm=TRUE)
phi1_up0  <- apply(phi1_mcmc0, 2, quantile, probs=c(0.025), na.rm=TRUE)

phi2_mu0  <- apply(phi2_mcmc0, 2, mean)
phi2_low0 <- apply(phi2_mcmc0, 2, quantile, probs=c(0.975), na.rm=TRUE)
phi2_up0  <- apply(phi2_mcmc0, 2, quantile, probs=c(0.025), na.rm=TRUE)



# REPLICATE #
##########

set.seed(12)

# Indices
K <- ncol(X)
N <- nrow(X)
D <- 2
M <- 5000

# Hparameters for priors 
L0j <- diag(D)
l0j <- matrix(0, D, 1 )
a0j <- 1
b0j <- 0.15

# Rescale data to zero mean / unit sd
for (i in 1:K) {
	X[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
	}
	
# Generate holding bins 
lambda <- list()
psi <- list()
phi <- list()

# Initialize with random starting values
lambda[[1]] <- matrix(rnorm(K*D), K, D )
psi[[1]] <- diag(K)
phi[[1]] <- matrix(NA,N,D)

# Make progress bar
pb <- txtProgressBar(min = 0, max = M, style = 3)

# Run MCMC 
for(m in 2:M){
	
	# -- Draw phi --
	phi[[m]] <- matrix(NA, N, D)
	for(i in 1:N){
		# Changes to Greenberg: 
		#  - for fi t(X[i,]) instead of X[i,]! otherwise not conformable 
		#  - I_n must be I_D
		Fi = solve( t(lambda[[m-1]]) %*% solve( psi[[m-1]] ) %*% lambda[[m-1]] + diag(D) )
		fi = Fi %*% ( t(lambda[[m-1]]) %*% solve( psi[[m-1]] ) %*% t(X[i, ]) ) 
		phi[[m]][i,] <- mvrnorm(1, mu=fi, Sigma=Fi)
		}
		
	# -- Draw lambda --
	lambda[[m]] <- matrix(NA, K, D)
	for(j in 1:K){
		# Greenberg has precision parameterization for the prior, but I use covariance.. 
		Lj = solve( (1/psi[[(m-1)]][j,j]) * t(phi[[m]]) %*% phi[[m]] + solve(L0j) )
		lj = Lj %*% (  (1/psi[[(m-1)]][j,j]) * t(phi[[m]]) %*% X[,j] + solve(L0j) %*% l0j) 		
		if (j != 2){
			lambda[[m]][j,] <- mvrnorm(1, mu=lj, Sigma=Lj)
			} else { 
			# Identification constraint for Examination: lambda[2,1] > 0 and lambda[2,2] < 0
			repeat { 
				lambda[[m]][j,] <- mvrnorm(1, mu=lj, Sigma=Lj)
			    if( lambda[[m]][j,1] > 0 &  lambda[[m]][j,2] < 0  ) break() 
				}
			}
		}
	# Identification constraint for Education / Infant Mortality: lambda[3,2] = 0 and lambda[5,1] = 0
	lambda[[m]][3,2] <- 0
	lambda[[m]][5,1] <- 0
	
	# -- Draw psi --
	psi[[m]] <- matrix(0, K, K)
	for(j in 1:K){
		aj <- a0j + N
		bj <- b0j + t(X[,j] - phi[[m]] %*% lambda[[m]][j,] ) %*% (X[,j] -  phi[[m]] %*% lambda[[m]][j,] ) 
		psi[[m]][j,j] <- rigamma(1, aj/2, bj/2)
		}
	setTxtProgressBar(pb, m)
	}
close(pb)
	
	
# Make array from lists	
lambda_mcmc <- array(unlist(lambda), dim=c(K,D,M) )
phi_mcmc <- array(unlist(phi), dim=c(N,D,M) )
psi_mcmc <- array(unlist(psi), dim=c(K,K,M) )

# # Traceplots
par(mfrow=c(1,3))
plot(lambda_mcmc[3,1,], type="l")
plot(phi_mcmc[3,1,], type="l")
plot(psi_mcmc[1,1,], type="l")

# Compare phi with phi from MCMCpack 
phi_mu <- apply(phi_mcmc, c(1,2), mean, na.rm=TRUE)
phi_up <- apply(phi_mcmc, c(1,2), quantile, probs=c(0.975), na.rm=TRUE)
phi_low <- apply(phi_mcmc, c(1,2), quantile, probs=c(0.025), na.rm=TRUE)

par(mfrow=c(1,2), pch=20, lwd=0.5)
plot(phi_mu[,1], phi1_mu0, ylim=c(-1.5,4), xlim=c(-1.5,4))
segments(phi_up[,1], phi1_mu0, phi_low[,1], phi1_mu0)
segments(phi_mu[,1], phi1_low0, phi_mu[,1], phi1_up0)
abline(a=0, b=1)
plot(phi_mu[,2], phi2_mu0, ylim=c(-1.5,2), xlim=c(-1.5,2))
segments(phi_up[,2], phi2_mu0, phi_low[,2], phi2_mu0)
segments(phi_mu[,2], phi2_low0, phi_mu[,2], phi2_up0)
abline(a=0, b=1)

# > same!
