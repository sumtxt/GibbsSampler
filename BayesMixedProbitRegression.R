# Gibbs Samplers for Probit Regression 
# with Varying Intercept  
#############################

rm(list=ls())

library(MASS)
library(msm)
library(Rcpp)

library(rjags)
load.module("glm")

makeRunMean <- function(x) return(cumsum(x)/seq(along=x))

RMixedProbit <- function(K,sim,nburn,n,J,alpha_idx,B0,b0,X,y,nJ,e0,f0,betastart,omegastart){

	# HOLDING BINS
	bhat <- matrix(NA, sim, K)
	ystar <- matrix(NA, sim, n)
	lower <- rep(-Inf, n)
	lower[y==1] <- 0
	upper <- rep(Inf, n)
	upper[y==0] <- 0

	omega <- rep(NA,sim)
	alpha_draw <- matrix(NA, sim, J)

	# Start 
	ystar[1,] <- rnorm(n)
	bhat[1, ] <- betastart
	alpha <- rnorm(J)[alpha_idx]
	omega[1] <- omegastart

	# SAMPLE 
	for(i in 2:sim){
		bmu <- solve( solve(B0) + t(X) %*% X ) %*% (solve(B0) %*% b0 + t(X) %*% (ystar[(i-1),]-alpha) )
		bsd <- solve(solve(B0) + t(X) %*% X)
		bhat[i,] <- mvrnorm(1, bmu, bsd)
		ystarmu <- (X %*% bhat[i, ])+alpha
		ystar[i, ] <- rtnorm(n, ystarmu, lower=lower, upper=upper)
		alphavar <- 1/(omega[(i-1)] + nJ)
		ystarbar <- sapply(seq(1,J), function(z) mean(ystar[i,alpha_idx==z]-X[alpha_idx==z,] %*% bhat[i, ]) )
		mb <- alphavar * (ystarbar * nJ)
		alpha_draw[i,] <- rnorm(J,mb, sqrt(alphavar))
		alpha <- alpha_draw[i,alpha_idx]
		e1 <- e0 + (J/2)
		f1 <- f0 + (sum(alpha_draw[i,]^2)/2)
		omega[i] <- rgamma(1, e1, f1)
		cat(i, " - ")
		}

	bhat <- bhat[nburn:sim, ]
	omega <- omega[nburn:sim]

	out <- list(omega=omega,bhat=bhat)
	return(out)
	}


jagsMixedProbit <- "
model {
	for(n in 1:N){
		probit(p[n]) <- X[n,1]*beta[1] +
				     X[n,2]*beta[2] +
				     X[n,3]*beta[3] + alpha[alpha_idx[n]]		
		y[n] ~ dbern(p[n])
		}

	for(k in 1:K){
		beta[k] ~ dnorm(0,1/100)
		}
	
	for (j in 1:J) {
		alpha[j] ~ dnorm(0,omega)
		}
	omega ~ dgamma(K0,W0)

	}
"


# SIMULATE DATA #
K <- 3
n <- 500
X <- matrix(NA,n,K)
J <- 15
alpha_idx <- sample(seq(1,J,by=1), size=n, replace=TRUE)
tau <- 3
a <- rnorm(J,0,tau)
nJ <- as.vector(table(alpha_idx,useNA="no"))

tcoefs <- c(-0.5,-1,2)
X[,1] <- 1
X[,2] <- rnorm(n)
X[,3] <- runif(n,-1,1)
xB <- X %*% tcoefs + a[alpha_idx]
y <- rbinom(n, 1, pnorm(xB) )


# PRIOR 
b0 <- c(0,0,0)
B0 <- diag(K)*100

e0 <- 0.01
f0 <- 0.01

# ESTIMATE 
sim <- 10000
nburn <- 1000
nthin <- 1

betastart <- c(0,0,0)
omegastart <- 0.5

fromR <- RMixedProbit(K=K,sim=sim+nburn,nburn=nburn,n=n,J=J,alpha_idx=alpha_idx,
		B0=B0,b0=b0,X=X,y=y,nJ=nJ,e0=e0,f0=f0,betastart=betastart, 
		omegastart=omegastart)

forJags <- list(K0=e0,W0=f0,X=X,y=y,alpha_idx=alpha_idx,N=n,J=J,K=K)
jags <- jags.model(textConnection(jagsMixedProbit), data = forJags)
update(jags, nburn)
fromJags <- coda.samples(jags, variable.names=c("beta","omega"), 
	n.iter=10001,n.thin=nthin)

# RESULTS 
muR <- apply(fromR$bhat, 2, mean)
quR <- t(apply(fromR$bhat, 2, quantile, probs=c(0.025,0.975)))
mu_omegaR <- mean(fromR$omega)
qu_omegaR <- quantile(fromR$omega, probs=c(0.025,0.975))
( resR <- cbind(c(muR,mu_omegaR),rbind(quR,qu_omegaR)) )

summary(fromJags)
