# Gibbs Samplers for Probit Regression #
#############################

#
# This script implements the Gibbs sampler by
#  Albert/Chib 1993 using data augmentation
#  Imai and van Dyk (2005) using marginal augmentation 
#


library(MASS)
library(msm)

makeRunMean <- function(x) return(cumsum(x)/seq(along=x))

# SIMULATE DATA #

n <- 900
X <- matrix(NA,n,2)
tcoefs <- c(1,-3)
X[,1] <- rep(1,n)
X[,2] <- runif(n,-1,1)
xB <- X %*% tcoefs
y <- rbinom(n, 1, pnorm(xB))

sim <- 2000

# PRIOR 
b0 <- c(0,0)
B0 <- diag(2)*100



#########################
# SAMPLER by Albert/Chib 1993) #
#########################

# HOLDING BINS
bhat <- matrix(NA, sim, 2)
ystar <- matrix(NA, sim, n)
lower <- rep(-Inf, n)
lower[y==1] <- 0
upper <- rep(Inf, n)
upper[y==0] <- 0
ystar[1,] <- rnorm(n)
bhat[1, ] <- c(0,0)

# SAMPLE 
for(i in 2:sim){
	bmu <- solve( solve(B0) + t(X) %*% X ) %*% (solve(B0) %*% b0 + t(X) %*% ystar[(i-1), ])
	bsd <- solve(solve(B0) + t(X) %*% X)
	bhat[i,] <- mvrnorm(1, bmu, bsd)
	ystarmu <- X %*% bhat[i, ]
	ystar[i, ] <- rtnorm(n, ystarmu, lower=lower, upper=upper)
	# cat(i, " - ")
	}

# SUMMARIZE 

# Tracelot
bhat <- bhat[200:sim, ]
par(mfrow=c(1,2))
plot(bhat[, 1], type="l", col="grey", ylab=expression(beta[0]), xlab="Iteration")
lines(makeRunMean(bhat[,1]), col="blue")
abline(h=tcoefs[1], col="red")
plot(bhat[, 2], type="l", col="grey", ylab=expression(beta[1]), xlab="Iteration")
lines(makeRunMean(bhat[,2]), col="blue")
abline(h=tcoefs[2], col="red")




##########################
# SAMPLER by Imai/van Dyk (2005) #
##########################

bhat <- matrix(NA, sim, 2)
lower <- rep(-Inf, n)
lower[y==1] <- 0
upper <- rep(Inf, n)
upper[y==0] <- 0
bhat[1, ] <- c(0,0 )

v0 <- 1
alphasq0 <- 1

bsd <- solve(solve(B0) + t(X) %*% X)

alphasqA <- alphasqB <- rep(NA, sim)

for(i in 2:sim){
	#alphasq <- rigamma(1,v0/2,alphasq0/2)
	alphasq <- alphasq0/rchisq(1, v0)
	ystarmu <- (X %*% bhat[(i-1), ])*sqrt(alphasq)
	ystar <- rtnorm(n, mean=ystarmu, sd=sqrt(alphasq), lower=lower, upper=upper)
	
	bmu <- bsd %*% (t(X) %*% ystar)
	S <- sum((ystar- (X %*% bmu) )^2)
	v1 <- v0 + n
	alpha1sq <- as.numeric(alphasq0 + S + (t(bmu) %*% solve(B0) %*% bmu))
	
	alphasqA[i] <- rigamma(1,v1/2,alpha1sq/2)
	alphasq <- alphasqB[i] <- alpha1sq/rchisq(1, v1)
	betastar <- mvrnorm(1, mu=bmu, Sigma=(alphasq*bsd) )
	bhat[i,] <- betastar/sqrt(alphasq)

#	cat(i, " - ")
	}

bhat <- bhat[200:sim, ]
par(mfrow=c(1,2))
plot(bhat[, 1], type="l", col="grey", ylab=expression(beta[0]), xlab="Iteration")
lines(makeRunMean(bhat[,1]), col="blue")
abline(h=tcoefs[1], col="red")
plot(bhat[, 2], type="l", col="grey", ylab=expression(beta[1]), xlab="Iteration")
lines(makeRunMean(bhat[,2]), col="blue")
abline(h=tcoefs[2], col="red")




# ESIMATE WITH CANNED ROUTINE 

fit <- MCMCprobit(y~X[,2], burnin=200, mcmc=5000, thin=1, b0=b0, B0=B0)
fit <- as.matrix(fit)


plot(fit, type="n", ylab=expression(beta[1]), xlab=expression(beta[0]))
points(fit, pch=".")
for(i in 1:(sim-3200)){
	lines(fit[(i-1):i, ], lwd=0.4)
	}
