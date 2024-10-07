# Network Modeling - HS 2023
# C. Stadtfeld, A. Espinosa-Rada, A. Uzaheta
# Based on previous work from: K. Mepham, V. Amati
# Social Networks Lab
# Department of Humanities, Social and Political Sciences
# ETH Zurich
# 20 November 2023
#
# Assignment 2 - Task 2



# Task 2.1 ----------------------------------------------------------------
# The function "simulation" simulates the network evolution between 
# two time points. 
# Given the network at time t1, denoted by x1, the function simulates the 
# steps of the continuous-time Markov chain defined by a SAOM with outdegree,
# recip and transTrip statistics. Unconditional simulation is used.
# The function returns the network at time t2.
# The structure of the algorithm is described in the file
# _Simulating from SAOM.pdf_ available in
# the Lecture notes and additional material section on Moodle.

#' Simulate the network evolution between two time points
#'
#' @param n number of actors in the network
#' @param x1 network at time t1
#' @param lambda rate parameter
#' @param beta1 outdegree parameter
#' @param beta2 reciprocity parameter
#' @param beta3 transTrip parameter
#'
#' @return network at time t2
#'
#' @examples
#' netT1 <- matrix(c(
#'   0, 1, 0, 0, 0,
#'   0, 0, 0, 1, 0,
#'   0, 0, 0, 1, 1,
#'   1, 0, 1, 0, 0,
#'   0, 1, 1, 0, 1
#'   ), 
#'   nrow = 5, ncol =  5, byrow = TRUE)
#' netT2 <- simulation(5, netT1, 4, -2, 0.5, 0.05)
simulation <- function(n, x1, lambda, beta1, beta2, beta3) {
  t <- 0 # time
  x <- x1
  while (t < 1) {
    dt <- rexp(1, n * lambda)
    # --- MISSING ---
    i <- sample(1:n, 1)  # Select an actor at random
    #iterate through all the j's to calculate the evaluation function if the change was made by i,j
    vector_p_ij <-  vector("numeric", length = n) #an vector to store the change probability
    for (j in 1:n){
      if (j==i){change_stat <- c(0,0,0)}
      else{
        sign_change <- -2 * (x[i,j]) + 1 #map value 0 to 1, map 1 to -1, 
        num_two_path <- x %*% x
        num_two_path_i_j <- num_two_path[i,j] #number of two path from i to j
        num_common_pred <- sum(x[i,]*x[j,]) #number of predecessor of i and j
        num_direct_succ <- sum(x[,i]*x[,j]) #number of successor  of i and j
        change_stat <- c(sign_change, 
                         sign_change * x[j,i] , 
                         sign_change * (num_two_path_i_j + num_common_pred + num_direct_succ)) 
      }
      ev_fct = exp(sum(c(beta1,beta2,beta3)*change_stat)) #calculate the probability of every possible change
      vector_p_ij[j] <- ev_fct #store it in the vector of change
    }
  vector_p_ij_norm <- vector_p_ij/sum(vector_p_ij)#normalize the probability to have sum 1
  sample <- rmultinom(n=1, size=1, vector_p_ij_norm) #draw j from multinomial distribution using the calculate probability
  index_j <- which(sample == 1) #get the index of the actor to which we make the change
  if (i != index_j){
    x[i,j] <- 1- x[i,j] # if j is not i itself, we toggle the edge
  }
  else{
    x <- x # if j coincides with i, we make the network unchanged.
  }
  t <- t + dt
  }
  return(x)
}

#netT2 <- simulation(5, netT1, 4, -2, 0.5, 0.05)

# Task 2.2 ----------------------------------------------------------------
# Consider the two adjacency matrices in the files net1.csv and net2.csv.
# Estimate the parameters of the SAOM with outdegree, reciprocity and
# transitive triplets statistics using the function `siena07`.
# You can extract the estimated parameters from the `rate` and `theta`
#  components of the output object (e.g., res$rate and res$theta).

# ---MISSING---
library(RSiena)
library(sna)
setwd("C:/Users/zheng/Desktop/ETH Third Semester/Network Modeling/Assignment2")
net1 <- as.matrix(read.csv("net1.csv", header = FALSE))
net2 <- as.matrix(read.csv("net2.csv", header = FALSE))

nets <- sienaDependent( array(c(net1, net2), dim = c(22, 22, 2)))
mydata <- sienaDataCreate(nets)
mydata
print01Report(mydata, modelname = "nets")

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myeff

myAlgorithm <- sienaAlgorithmCreate(
  projname = "mynets",
  nsub = 4, n3 = 1000, seed = 2023
)
model2_2 <- siena07(
  myAlgorithm,
  data = mydata, effects = myeff,
  returnDeps = TRUE,
  useCluster = TRUE, nbrNodes = 3
)

model2_2

# Triad census
gof1.tc <- sienaGOF(model2_2, verbose = TRUE,
                    varName = "nets", TriadCensus)

plot(gof1.tc, center = TRUE, scale = TRUE)

source("printSiena.R")
source("siena07ToConverge.R")
printSiena(model2_2)


# Task 2.3 ----------------------------------------------------------------
# Conditioning on the first observation, generate 1,000 simulations of the 
# network evolution
# Compute the triad census counts for each simulated network.
# Save the results in an object, named `triadCensus`, in which rows are
# the index of a simulated network and columns are the type of triads.
# Column names should use the triad type name, e.g., "003", "012", "102", ... 

# ---MISSING---

num_simulations <- 1000
# Initialize matrix for storing triad census results
num_triad_types <- length(colnames(triad.census(net1))) 
triadCensus <- matrix(nrow = num_simulations, ncol = num_triad_types)
colnames(triadCensus) <- colnames(triad.census(net1))
# Run simulations
n_nodes <- 22
for (i in 1:num_simulations) {
  # Simulate network evolution
  print(i)
  netT2 <- simulation(n_nodes, net1, 4.097, -1.1067, 0.4750, 0.0763) # Use parameters of the previous model
  triad_census_result <- triad.census(netT2) # get the triad census of the simulation 
  triadCensus[i, ] <- triad_census_result # Store the result in the matrix
}



# Task 2.4 ----------------------------------------------------------------
## i. standardized the simulated network stats. ----
##   Name the resulting object as triadCensusStd

# ---MISSING---

column_mean <- as.vector(apply(triadCensus, MARGIN = 2, mean))
column_std <- as.vector(apply(triadCensus, MARGIN = 2, sd))

triadCensus_centered <- sweep(triadCensus, MARGIN = 2, STATS = column_mean, FUN = "-")
triadCensusStd <- sweep(triadCensus_centered, MARGIN = 2, STATS = column_std, FUN = "/")
triadCensusStd



## ii. variance-covariance matrix and its generalized inverse.         ----

# ---MISSING---
var_cov <- cov(triadCensusStd)
gen_inv <- MASS::ginv(var_cov)

## iii. standardized the observed values of the triad census counts    ----
##  in the second observation using values from i.

# ---MISSING---
triadCensus_net2 <- triad.census(net2)
triadCensus_net2 <- as.vector(triadCensus_net2)
triadCensus_net2_std <- (triadCensus_net2 - column_mean)/column_std


## iv. Monte-Carlo Mahalanobis distance computation                                ----
# Compute the Mahalanobis distance using the mhd function for 
# the auxiliar statistics of the simulated networks and the observed network.
# Remember to drop statistics with variance of 0 for the plot and
# Mahalanobis distance computation, report which statistics suffer this issue.

#' Compute the Mahalanobis distance
#'
#' @param auxStats numerical vector with the mean centered or standardized
#'   auxiliar statistics
#' @param invCov numerical matrix with the inverse of the variance-covariance
#'   matrix of the auxiliar statistics in the simulated networks
#'
#' @return numeric value with the Mahalanobis distance of auxiliar stats
#'
#' @examples
#' mhd(c(2, 4) - c(1.5, 2), solve(matrix(c(1, 0.8, 0.8, 1), ncol = 2)))
mhd <- function(auxStats, invCov) {
  t(auxStats) %*% invCov %*% auxStats
}

# ---MISSING---
#MHD of the simulations, return a vector of 1000 values
MHD_sim <- apply(triadCensusStd, MARGIN = 1, function(row) mhd(auxStats = row, invCov = gen_inv))
MHD_obs <- mhd(auxStats = triadCensus_net2_std, invCov = gen_inv)


## iii. Monte-Carlo p-value computation                                ----
# Compute the proportion of simulated networks where the distance is 
# equal or greater than the distance in the observed network.

# ---MISSING---
prop <- (sum(MHD_sim >= MHD_obs[1,1])) / length(MHD_sim)
prop

# violin plots ------------------------------------------------------------
# Fill out the missing part and run the code to obtain the violin plots

# install.packages(c("tidyverse", "ggplot2"))  # # run this line to install 
library(tidyverse)  # used: dplyr and tidyr
library(ggplot2)

# Given the array triadCensusStd, create a data frame from it in a long format, 
# do the same for the observed network statistics at time t2.
# Named the data frame "triadCensusDf" and "triadCensusObs".
# Drops statistics with variance of 0 for the plot.

triadCensusDf <- data.frame(triadCensusStd) |>
  select(where(~ var(.) > 0)) |>  # Drop statistics with zero variance
  pivot_longer(
    everything(),
    names_to = "triad", names_pattern = "^X(.+)$",
    values_to = "nnodes"
  )

# Compute the statistics of the observed network at time t2,
#  standardized using the stats from 2.4 literal i.
triadCensusObs <- triadCensus_net2_std |> 
  data.frame() |>
  pivot_longer(
    everything(),
    names_to = "triad", names_pattern = "^X(.+)$",
    values_to = "nnodes"
  ) |>
  filter(triad %in% unique(triadCensusDf$triad))

# The following code computes the 5% and the 95% quantiles
# of the triad counts by type
percTriad <- triadCensusDf |>
  group_by(triad) |>
  summarise(
    quant05 = quantile(nnodes, prob = 0.05),
    quant95 = quantile(nnodes, prob = 0.95)
  ) |>
  pivot_longer(
    starts_with("quant"),
    names_to = "quant", names_pattern = "quant(.+)",
    values_to = "nnodes"
  )


# The following code produces the violin plots
ggplot(triadCensusDf, aes(fct_inorder(triad), nnodes)) +
  geom_violin(trim = FALSE, scale = "width") +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_boxplot(width = 0.2, fill = "gray", outlier.shape = 4) +
  geom_point(data = triadCensusObs, col = "red", size = 2) +
  geom_line(
    data = triadCensusObs, aes(group = 1), col = "red", linewidth = 0.5
  ) +
  geom_line(
    data = percTriad, mapping = aes(group = quant),
    col = "gray", linetype = "dashed"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlab("triad type")

