library(pscl)
library(coda)

set.seed(12)

# Simulate data #
#################

beta1=c(-2,2,0.1)
nobs=500

X=cbind(1, runif(nobs),rnorm(nobs))
e <- rnorm(nobs, 0, sqrt(1.8))
y=X%*%beta1+e


# Bayesian Linear Regression #
##############################

# Indices
K <- ncol(X)
N <- nrow(X)
M <- 5000

# Hyperparameters
b0 <- rep(0,K)
B0 <- diag(K) * 100
e0 <- 0.01
f0 <- 0.01

# Generate holding bins 
beta <- list()
sigmasq <- list()

# Starting values 
sigmasq[[1]] <- 0.1
beta[[1]] <- rnorm(K)

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = M, style = 3)

# :: Run MCMC ::
for(m in 2:M){
	
	# -- Draw: beta --
	B1 <- solve( (1/sigmasq[[m-1]]) * (t(X) %*% X) + solve(B0) )
	b1 <- B1 %*% ( (1/sigmasq[[m-1]]) * (t(X) %*% y ) + solve(B0) %*% b0 )
	beta[[m]] <- mvrnorm(1, mu=b1, Sigma=B1)
	
	# -- Draw: sigmasq --
	e1 <- e0 + N
	f1 <- f0 + t(y - X %*% beta[[m]] ) %*% (y - X %*% beta[[m]] )
	sigmasq[[m]] <- rigamma(1, e1/2, f1/2)
	
	setTxtProgressBar(pb, m)
	}
close(pb)
	
# Combine draws in coda mcmc-object
sigmasq_mcmc <- unlist(sigmasq)
beta_mcmc <- do.call(rbind, beta)
posterior <- cbind(beta_mcmc, sigmasq_mcmc)
colnames(posterior) <- c(paste(rep("beta"), seq(1,K), sep="_"), "sigmasq")
posterior <- as.mcmc(posterior)

# Traceplots
plot(posterior)

# Results 
summary(posterior)

# Compare to MCMCpack results 
summary(MCMCpack::MCMCregress(y~X-1,c0=e0, d0=f0))
