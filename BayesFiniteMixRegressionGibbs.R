# Finite Mixture of a Linear Regression Model #
###############################################

#
# This script implements the Gibbs sampler 
# suggested in Frühwirth-Schnatter p. 253, Algorithm 8.1
# All coefficient are assumed to vary across mixtures. 
#


library(MCMCpack)
library(pscl)

set.seed(21)
rm(list=ls())

# SIMULATE SOME DATA 

KKK=3

beta1=c(-2, 2,  0.1)
beta2=c(2,  1,  0.4)
beta3=c(0.1,1.1,1.4)
betamatrix <- rbind(beta1,beta2, beta3)

nobs1=100
nobs2=150
nobs3=200

X1=cbind(1, runif(nobs1,-1,1),rbeta(nobs1,0.5, 1))
X2=cbind(1, runif(nobs2,-1,1),rbeta(nobs2,0.5, 1))
X3=cbind(1, runif(nobs3,-1,1),rbeta(nobs3,0.5, 1))

sigmavector <- c(1,0.4, 1.2)

y1=X1%*%beta1+rnorm(nobs1, 0, sqrt(sigmavector[1]))
y2=X2%*%beta2+rnorm(nobs2, 0, sqrt(sigmavector[2]))
y3=X3%*%beta3+rnorm(nobs3, 0, sqrt(sigmavector[3]))


y <- as.matrix(c(y1,y2,y3))
X <- rbind(X1,X2,X3)



# ESTIMATE WITH OWN CODE #


set.seed(12)

# Indices
K <- ncol(X)
N <- nrow(X)
D <- 3
M <- 5000

# Hparameters for priors 
b0 <- rep(0,K)
B0 <- diag(K) * 100
e0 <- 0.01
f0 <- 0.01
a0 <- 1

# Generate holding bins 
beta <- list()
sigma <- list()
eta <- list()

# Initialize with random starting values
sigma[[1]] <- matrix(runif(D, 0, 1), D, 1)
S <- sample(seq(1,D), N, replace=TRUE)

# Make progress bar
pb <- txtProgressBar(min = 0, max = M, style = 3)

# Run MCMC 
for(m in 2:M){
	
	# -- Draw: eta -- 
	a1 <- apply(matrix(seq(1,D),D,1), 1, function(x) sum(S==x) ) + a0
	eta[[m]] <- rdirichlet(1, alpha=a1)
	
	# -- Housekeeping --
	beta[[m]] <- matrix(NA, D, K)
	sigma[[m]] <- matrix(NA, D, 1)
	pS <- matrix(NA, N, D)
	
	for(i in 1:D){
		nk <- sum(S==i)
		Xk <- matrix(X[S==i, ], nk, K) 
		yk <- matrix(y[S==i], nk, 1)
		
		# -- Draw: beta --
		B1 <- solve( (1/sigma[[m-1]][i,]) * (t(Xk) %*% Xk) + solve(B0) )
		b1 <- B1 %*% ( (1/sigma[[m-1]][i,]) * (t(Xk) %*% yk ) + solve(B0) %*% b0 )
		beta[[m]][i,] <- mvrnorm(1, mu=b1, Sigma=B1)
		
		# -- Draw: gamma --
		e1 <- e0 + nk
		f1 <- f0 + t(yk - Xk %*% beta[[m]][i, ] ) %*% (yk - Xk %*% beta[[m]][i,] )
		sigma[[m]][i, ] <- rigamma(1, e1/2, f1/2)
		
		# -- Calculate p(S|.)
		pS[,i] <- eta[[m]][i] * dnorm(y, mean=(X%*%beta[[m]][i,]), sd=sqrt(sigma[[m]][i,]) )
		
		}
	
	# -- Draw: psi -- 
	S <- apply(pS, 1, function(x) sample(1:D, 1, prob=x/sum(x)) )
	
	setTxtProgressBar(pb, m)
	}
close(pb)
	
# Make array from lists	
eta_mcmc <- array(unlist(eta), dim=c(1,D,M) )
beta_mcmc <- array(unlist(beta), dim=c(D,K,M) )
sigma_mcmc <- array(unlist(sigma), dim=c(D,1,M) )


# # Traceplots
par(mfrow=c(3,2))
plot(beta_mcmc[1,1,], type="l")
plot(beta_mcmc[2,1,], type="l")
plot(beta_mcmc[1,2,], type="l")
plot(beta_mcmc[2,2,], type="l")
plot(beta_mcmc[1,3,], type="l")
plot(beta_mcmc[2,3,], type="l")

par(mfrow=c(1,2))
plot(sigma_mcmc[1,1,], type="l")
plot(sigma_mcmc[2,1,], type="l")

par(mfrow=c(1,2))
plot(eta_mcmc[1,1,], type="l", ylim=c(0,1))
abline(h=mean(eta_mcmc[1,1,]), col="red")
plot(eta_mcmc[1,2,], type="l", ylim=c(0,1))
abline(h=mean(eta_mcmc[1,2,]), col="red")

apply(beta_mcmc, c(1,2), quantile, probs=c(0.025,0.5,0.975))
apply(sigma_mcmc, c(1,2), quantile, probs=c(0.025,0.5,0.975))




# ESTIMATE WITH CANNED ROUTINE #

library(flexmix)
ex1 <- flexmix(y~X-1, k=KKK)
parameters(ex1)

# > same

