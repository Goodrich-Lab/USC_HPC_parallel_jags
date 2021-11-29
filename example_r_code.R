# script for R and RJAGS analysis
# Jesse Goodrich 111221

# Load Libraries 
library(rjags)
library(R2jags)
library(pbdMPI)

# Initiate connections
init()

#Set working directory 
setwd("/project/dconti_624/Users/jagoodri/example_jags_lin_mod")

# create simulated data -------------------------------------------------------


N <- 1000
alpha <- 2
beta <- 4
Y.sd <- 2
X <- rnorm(N, mean=0, sd=1) # generate standard normal independent variable
mu <- alpha + beta*X
Y <- rnorm(N, mean=mu, sd=Y.sd) # generate dependent variable

# Create function for running JAGS Model
jags_lin_mod <- function(N, Y, X){
  linear.model <- 
    "model {
  for(i in 1:N) {
    Y[i] ~ dnorm(mu[i], prec.sigma)
    mu[i] <- alpha + beta*X[i]
  }
  # prior on outcome variance
  prec.sigma <- 1/(sigma*sigma)
  sigma ~ dunif(0,3)
  
  alpha ~ dnorm(0, 1.0E-06) # prior on intercept
  beta ~ dnorm(0, 1.0E-06) # prior on slope

}"
  
  jdata <- list(N=N, Y=Y, X=X)
  var.s <- c("alpha", "beta", "sigma")
  model.fit <- jags.model(file=textConnection(linear.model),
                          data=jdata, 
                          n.chains=1, 
                          n.adapt=1000,
                          quiet=T)
  update(model.fit, n.iter=1000, progress.bar="none")
  model.fit <- coda.samples(model=model.fit, 
                            variable.names=var.s, 
                            n.iter=1000, 
                            thin=1, 
                            progress.bar="none")
  
  r <- summary(model.fit)
  
  colnames(r$statistics) <- c("Mean", "SD", "Naive_SE", "Time_series_SE")
  
  return(r$statistics)
}

# create function to run iterations of the function created above --------------
# This can also be used to iterate across different columns in a data frame
model_iterations <- function(i){
  # Run Model
  output <- jags_lin_mod(N = N, Y = Y, X = X)
  
  # Name Metabolite
  output$iteration = i
  #Return Output
  return(output)
}

# Set number of iterations -------------------------------------------------
n_iter <- 1280

# Run model ---------------------------------------------------------------
coefs <- pbdLapply(X = 1:n_iter, 
                   FUN = model_iterations, 
                   pbd.mode = "spmd")

# Save results  -------------------------------------------------------------
comm.write.csv(coefs,  file = "example_jags_lin_mod_code.csv")

# Write Success from each connection
message(paste("SUCCESS from rank", comm.rank()))

# End connections
finalize(mpi.finalize = TRUE)