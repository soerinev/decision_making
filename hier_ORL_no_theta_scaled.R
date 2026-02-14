install.packages("pacman")
pacman::p_load(R2jags, parallel)

set.seed(1983)

setwd('/work/exam/')

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#----------getting the data
# data from this paper: https://www.frontiersin.org/articles/10.3389/fpsyg.2014.00849/full
# available here: https://figshare.com/articles/dataset/IGT_raw_data_Ahn_et_al_2014_Frontiers_in_Psychology/1101324

#load control data
stationary4 <- read.table("data/stationary.txt",header=TRUE)

#----------prepare data for jags models - want trial x subject arrays for choice and gain & loss ----
# identify and count unique subject IDs
subIDs <- unique(stationary4$subjID)
nsubs <- length(subIDs)
ntrials_max <- 100

# all choices (x) and outcomes (X)
x_raw <- stationary4$deck
X_raw <- stationary4$gain + stationary4$loss #note the sign!

#--- assign choices and outcomes in trial x sub matrix

#different number of trials across subjects. We'll need to fix this by padding arrays of < 100
#this is just so we can make the array
#then we'll also need to record number of valid trials for each sub, 
#then run the JAGS model on only valid trials

# empty arrays to fill
ntrials_all <- array(0,c(nsubs))
x_all <- array(0,c(nsubs,ntrials_max))
X_all <- array(0,c(nsubs,ntrials_max))

for (s in 1:nsubs) {
  
  #record n trials for subject s
  ntrials_all[s] <- length(x_raw[stationary4$subjID==subIDs[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub <- x_raw[stationary4$subjID==subIDs[s]] 
  length(x_sub) <- ntrials_max
  
  X_sub <- X_raw[stationary4$subjID==subIDs[s]] 
  length(X_sub) <- ntrials_max
  
  # assign arrays
  x_all[s,] <- x_sub
  X_all[s,] <- X_sub
  
}

# Scaling the payoffs (cuz the learning parameter becomes less relevant for very large payoffs/losses)
X_all <- X_all/100

#----------testing our data curation by running JAGS on one subject

# Now we'll fit one subject just to make sure everything works

x <- x_all[1,]
X <- X_all[1,]

ntrials <- ntrials_all[1]

# set up jags and run jags model on one subject
data <- list("x","X","ntrials") 
params<-c("a_rew","a_pun","K","omega_f","omega_p","p")
samples <- jags.parallel(data, inits=NULL, params,
                model.file ="ORL_no_theta.txt",
                n.chains=4, n.iter=10000, n.burnin=2000, n.thin=1, n.cluster=4)

# let's look at the posteriors for the parameters
par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$a_rew))
plot(density(samples$BUGSoutput$sims.list$a_pun))
plot(density(samples$BUGSoutput$sims.list$K))
plot(density(samples$BUGSoutput$sims.list$omega_f))
plot(density(samples$BUGSoutput$sims.list$omega_p))

# Question: how would you expect the data to look on the basis of these posteriors?


###########################################################
#---------- run the hierarchical model on controls --------
###########################################################

x <- x_all
X <- X_all

ntrials <- ntrials_all

# set up jags and run jags model
data <- list("x","X","ntrials","nsubs") 
# NB! we're not tracking theta cuz we're not modelling it in order reduce complexity a bit (hence, we're just setting it to 1 in "hier_ORL.txt")
params<-c("mu_a_rew","mu_a_pun","mu_K","mu_omega_f","mu_omega_p") 

start_time = Sys.time()
samples <- jags.parallel(data, inits=NULL, params,
                         model.file ="hier_ORL_no_theta.txt",
                         n.chains=4, n.iter=10000, n.burnin=2000, n.thin=1, n.cluster=4)
end_time = Sys.time()
end_time - start_time

par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$mu_a_rew))
plot(density(samples$BUGSoutput$sims.list$mu_a_pun))
plot(density(samples$BUGSoutput$sims.list$mu_K))
plot(density(samples$BUGSoutput$sims.list$mu_omega_f))
plot(density(samples$BUGSoutput$sims.list$mu_omega_p))

# xlim scaled to Haines et al.
par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$mu_a_rew), xlim=c(0,0.4))
plot(density(samples$BUGSoutput$sims.list$mu_a_pun), xlim=c(0,0.08))
plot(density(samples$BUGSoutput$sims.list$mu_K), xlim=c(0,2))
plot(density(samples$BUGSoutput$sims.list$mu_omega_f), xlim=c(0,5))
plot(density(samples$BUGSoutput$sims.list$mu_omega_p), xlim=c(-4,2))


# WHAT NOW

# 1) access posterior samples
# The posterior samples are stored here:
samples$BUGSoutput$sims.list

# Each parameter is a vector of posterior samples
mu_a_rew <- samples$BUGSoutput$sims.list$mu_a_rew
mu_a_pun <- samples$BUGSoutput$sims.list$mu_a_pun
mu_K <- samples$BUGSoutput$sims.list$mu_K
mu_omega_f <- samples$BUGSoutput$sims.list$mu_omega_f
mu_omega_p <- samples$BUGSoutput$sims.list$mu_omega_p

# Check the length (should be 12000 samples)
length(mu_a_rew)



# 2) save data
# Option 1: Save as a dataframe
library(dplyr)

posterior_data <- data.frame(
  mu_a_rew = samples$BUGSoutput$sims.list$mu_a_rew,
  mu_a_pun = samples$BUGSoutput$sims.list$mu_a_pun,
  mu_K = samples$BUGSoutput$sims.list$mu_K,
  mu_omega_f = samples$BUGSoutput$sims.list$mu_omega_f,
  mu_omega_p = samples$BUGSoutput$sims.list$mu_omega_p
)

# Save as CSV
write.csv(posterior_data, "stationary5.csv", row.names = FALSE)

# Option 2: Save the entire samples object
saveRDS(samples, "stationary5.rds")

# Later you can load it back:
# samples <- readRDS("stationary_clean_samples.rds")


















