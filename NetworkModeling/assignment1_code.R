#Assginment1
# setwd("C:/Users/zheng/Desktop/ETH Third Semester/Network Modeling/Assignment1")
library(sna)
library(network)
library(xtable)
library(ergm)

# Load the data
load("friend_net.Rda")

#task 1
#extract adjancency matrix
adj_matrix <- as.matrix(friend_net)
#extract sex covariate from the attribute 
gender <- get.node.attr(friend_net, "sex")
#set the x variable for first hypothesis
sameGender <- outer(gender,gender, "==") * 1

#perform QAP test
set.seed(1)
permutations <- 1000
nl1 <- netlogit(adj_matrix, sameGender, reps = permutations, nullhyp = "qapy")

nl1$names <- c("intercept","sameGender")
summary(nl1)

#task 1.2
gender <- attributes$sex
smoker <- attributes$smoke
activity <- attributes$activity
nNodes <- length(gender)

sender_gender <- matrix(gender, nNodes, nNodes, byrow= FALSE)
receiver_smoker <- matrix(smoker, nNodes, nNodes, byrow = TRUE)
same_activity <- outer(activity, activity, "==") * 1
nl2 <- netlogit(adj_matrix, list(sameGender, sender_gender, receiver_smoker, same_activity),
                reps = permutations, nullhyp = "qapspp")

nl2$names <- c("intercept","sameGender","sender_gender","receiver_smoker","same_activity")
summary(nl2)

# Task 1.5

#mark hockey as 1, other activity as 0
hockey <- activity
hockey[hockey == 2 | hockey == 3] <- 0
sender_hockey <- matrix(hockey, nNodes, nNodes, byrow = FALSE)
nl3 <- netlogit(adj_matrix, list(sameGender, sender_gender, receiver_smoker, same_activity, sender_hockey), reps = permutations, nullhyp = "qapspp")
nl3$names <- c("intercept","sameGender","sender_gender","receiver_smoker","same_activity","sender_hockey")
summary(nl3)


#task 2
num_edges <- sum(adj_matrix)

num_mutual <- sum(adj_matrix * t(adj_matrix))/2

num_homophily <- 0
# Iterate through all vertex pairs
for (i in 1:nNodes) {
  for (j in 1:nNodes) {
    # Check for directed tie and matching node attributes
    if (adj_matrix[i, j] == 1 && gender[i] == gender[j]) {
      num_homophily <- num_homophily + 1
    }
  }
}

obs_statistics <- c(num_edges,num_mutual,num_homophily)
obs_statistics

# MHstep ------------------------------------------------------------------
#' Simulate the next step of a network in Markov chain using Metropolis-Hasting
#' 
#' The function `MHstep` simulates the the Metropolis-Hastings step that defines
#' the Markov chain whose stationary distribution is the ERGM with 
#' edge, mutual and nodematch statistics
#'
#' @param net an object of class `matrix`. Adjacency matrix of the network.
#' @param nodeAttr a character or numeric vector. The node attribute. 
#' @param theta1 a numeric value. The value of the edge parameter of the ERGM.
#' @param theta2 a numeric value. The value of the mutual parameter of the ERGM.
#' @param theta3 a numeric value. The value of the nodematch parameter of the ERGM.
#'
#' @return next state of the Markov Chain
#'
#' @examples
#' MHstep(
#'   matrix(c(0, 1, 0, 0, 0, 0, 1, 1, 0), nrow = 3, ncol = 3),
#'   c("v", "g", "g"),
#'   -log(0.5), log(0.4), log(.8)
#' )
MHstep <- function(net, nodeAttr, theta1, theta2, theta3){
  
  # Number of vertices in the network
  nvertices <- nrow(net) 
  
  # Choose randomly two vertices, prevent loops {i,i} with replace = FALSE
  tie <- sample(1:nvertices, 2, replace = FALSE) 
  i <- tie[1]
  j <- tie[2]
  
  # Compute the change statistics
  
  #                --- MISSING---
  
  # Initialize counters for the statistics
  num_edges <- sum(net)
  num_mutual <- sum(net * t(net))/2 #we only consider every pair once
  num_homophily <- 0  #every pair should be cound twice
  # Iterate through all vertex pairs
  for (ii in 1:nvertices) {
    for (jj in 1:nvertices) {
      # Check for directed tie and matching node attributes
      if (net[ii, jj] == 1 && nodeAttr[ii] == nodeAttr[jj]) {
        num_homophily <- num_homophily + 1
      }
    }
  }
  #register the current states of the vector of statistics 
  current_stats <- c(num_edges, num_mutual, num_homophily)
  
  ##initialze the statiscis for the next possible state
  prop_net <- net
  prop_net[i, j] <- 1 - net[i, j] # if the tie was presented we remove. If the tie was lack, we add. 

  num_edges2 <- sum(prop_net)
  num_mutual2 <- sum(prop_net * t(prop_net))/2
  num_homophily2 <- 0
  # Iterate through all vertex pairs
  for (ii in 1:nvertices) {
    for (jj in 1:nvertices) {
      # Check for directed tie and matching node attributes
      if (prop_net[ii, jj] == 1 && nodeAttr[ii] == nodeAttr[jj]) {
        num_homophily2 <- num_homophily2 + 1
      }
    }
  }
  #the new statistics in the next possible state
  new_stats <- c(num_edges2, num_mutual2, num_homophily2)
  change_stats <- new_stats - current_stats # change of statistics
  #print(current_stats)
  #print(new_stats)
  
  
  # Compute the probability of the next state 
  # according to the Metropolis-Hastings algorithm
  
  #                --- MISSING---
  prob_trans = exp(sum(change_stats * c(theta1, theta2, theta3)))#probability of transition
  prob_trans = min(1, prob_trans)#by MH algorithm, we consider the minimum as transition probability
  
  # Select the next state: 
  
  #                --- MISSING---
  random_num <- runif(1)#introduce certain randomness to decide if we move to the next state
  if (random_num < prob_trans){ 
    net <- prop_net
  }
  #print(prob_trans)
  # Return the next state of the chain
  return(net)
}

# Markov Chain simulation -------------------------------------------------
#' The function MarkovChain simulates the networks from the ERGM with 
#' edge, mutual and nodematch statistics
#'
#' @param net an object of class `matrix`. Adjacency matrix of the network.
#' @param nodeAttr a character or numeric vector. The node attribute. 
#' @param theta1 a numeric value. The value of the edge parameter of the ERGM.
#' @param theta2 a numeric value. The value of the mutual parameter of the ERGM.
#' @param theta3 a numeric value. The value of the nodematch parameter of the ERGM.
#' @param burnin an integer value.
#'   Number of steps to reach the stationary distribution.
#' @param thinning an integer value. Number of steps between simulated networks.
#' @param nNet an integer value. Number of simulated networks to return as output.
#'
#' @return a named list:
#'   - netSim: an `array` with the adjancency matrices of the simulated networks.
#'   - statSim: a `matrix` with the value of the statistic defining the ERGM.
#'
#' @examples
#' MarkovChain(
#'   matrix(c(0, 1, 0, 0, 0, 0, 1, 1, 0), nrow = 3, ncol = 3),
#'   c("v", "g", "g"),
#'   -log(0.5), log(0.4), log(.8)
#' )
MarkovChain <- function(net, nodeAttr, theta1, theta2, theta3,
                        burnin = 10000, thinning = 1000, nNet = 1000){
  
  # Burnin phase: repeating the steps of the chain "burnin" times  
  nvertices <- nrow(net)
  burninStep <- 1 # counter for the number of burnin steps
  
  # Perform the burnin steps
  #                --- MISSING---
  while (burninStep <= burnin) {
    net <- MHstep(net, nodeAttr, theta1, theta2, theta3) # we perform MH step until we reach #burnin times 
    burninStep <- burninStep + 1
  }
  
  # After the burnin phase we draw the networks
  # The simulated networks and statistics are stored in the objects
  # netSim and statSim
  netSim <- array(0, dim = c(nvertices, nvertices, nNet))
  statSim <- matrix(0, nNet, 3)
  thinningSteps <- 0 # counter for the number of thinning steps
  netCounter <- 1 # counter for the number of simulated network
  
  #                --- MISSING---
  while (netCounter <= nNet){
    thinningSteps <- 0 # counter for the number of thinning steps
    while (thinningSteps <= thinning) {
      net <- MHstep(net, nodeAttr, theta1, theta2, theta3)
      thinningSteps <- thinningSteps + 1
    }
    netSim[,,netCounter] <- net
    num_edges <- sum(net)
    num_mutual <- sum(net * t(net)) / 2
    num_homophily <- 0
    # Iterate through all vertex pairs
    for (i in 1:nvertices) {
      for (j in 1:nvertices) {
        # Check for directed tie and matching node attributes
        if (net[i, j] == 1 && nodeAttr[i] == nodeAttr[j]) {
          num_homophily <- num_homophily + 1
        }
      }
    }
    statSim[netCounter,] <- c(num_edges, num_mutual, num_homophily)
    netCounter <- netCounter + 1
    #print(c(num_edges, num_mutual, num_homophily))
  }
  
  # Return the simulated networks and the statistics
  return(list(netSim = netSim, statSim = statSim))
}

zero_matrix <-array(0, dim = c(nNodes, nNodes))

plot_trace <- function(MC_res,obs_statistics,mc_param){
  par(xpd = NA,mfrow = c(3, 1))
  plot(MC_res$statSim[,1]-obs_statistics[1],type="l", main=paste("Trace of 1st statistic's diff",mc_param),xlab="sample",ylab="value")
  plot(MC_res$statSim[,2]-obs_statistics[2],type="l", main=paste("Trace of 2nd statistic's diff",mc_param),xlab="sample",ylab="value")
  plot(MC_res$statSim[,3]-obs_statistics[3],type="l", main=paste("Trace of 3nd statistic's diff",mc_param),xlab="sample",ylab="value")
}

MC_res1 <- MarkovChain(zero_matrix, gender, -2.76,  0.68,  1.21) # same
MC_res2 <-MarkovChain(zero_matrix, gender, -3.76,  -0.32,  0.21) # -1,-1,-1
MC_res3 <-MarkovChain(zero_matrix, gender, -3.76,  -0.32,  1.21) #-1,-1,+1
MC_res4 <-MarkovChain(zero_matrix, gender, -3.76,  1.68,  0.21)  # -1, +1, -1
MC_res5 <-MarkovChain(zero_matrix, gender, -3.76,  1.68,  2.21) # -1, +1, +1
MC_res6 <-MarkovChain(zero_matrix, gender, -1.76,  -0.32,  0.21)  # +1, -1, -1
MC_res7 <-MarkovChain(zero_matrix, gender, -1.76,  -0.32,  2.21)  # +1, -1, +1
MC_res8 <-MarkovChain(zero_matrix, gender, -1.26,  1.68,  0.21)  # +1, +1, -1
MC_res9 <-MarkovChain(zero_matrix, gender, -1.26,  1.68,  2.21)  # +1, +1, +1

plot_trace(MC_res1, obs_statistics, "-2.76,  0.68,  1.21")
plot_trace(MC_res2, obs_statistics, "-3.76,  -0.32,  0.21")
plot_trace(MC_res3, obs_statistics, "-3.76,  -0.32,  1.21" )
plot_trace(MC_res4, obs_statistics, "-3.76,  1.68,  0.21")
plot_trace(MC_res5, obs_statistics, "-3.76,  1.68,  2.21")
plot_trace(MC_res6, obs_statistics, "-1.76,  -0.32,  0.21")
plot_trace(MC_res7, obs_statistics, "-1.76,  -0.32,  2.21")
plot_trace(MC_res8, obs_statistics, "-1.26,  1.68,  0.21")
plot_trace(MC_res9, obs_statistics, "-1.26,  1.68,  2.21")


#task 3.1
set.seed(1)
#+gwesp(decay = 0.3, fixed = TRUE) 
model3.1 <- ergm(friend_net ~ edges + nodematch("sex"))
summary(model3.1)
theta1 <- model3.1$coef[1]
theta2 <- model3.1$coef[2]
# Computing the probability of a tie not reciprocating an existing tie is
oSameSex <- exp(theta1+theta2)
pSameSex <- oSameSex / (1 + oSameSex)
pSameSex


#task 3.2
set.seed(1)

model3.2 <- ergm(friend_net ~ edges + nodematch("sex") + mutual 
                 + gwesp(decay = 0.3, fixed = TRUE)
                 + gwodegree(decay = 0.3, fixed = TRUE) + gwidegree(decay = 0.3, fixed = TRUE) )

summary(model3.2)

# ERGM diagnostics and fit ------------------------------------------------
## Model convergence ------------------------------------------------------
mcmc.diagnostics(model3.2)



library(coda)
plot(model3.2$sample)


## Goodness of fit --------------------------------------------------------
model3.2gof <- gof(model3.2)
model3.2gof

par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))
plot(model3.2gof)

# ERGM interpretation -----------------------------------------------------
summary(model3.2)


#task 4.1
set.seed(1)
#set main.method to "Stochastic Approximation"
control.ergm(main.method = "Stochastic-Approximation")

#add smoke and activity attributes to friend_net
smoke = attributes$smoke
set.vertex.attribute(friend_net, "smoke", smoke)
activity = attributes$activity
set.vertex.attribute(friend_net, "activity", activity)
friend_net

#test without terms from 3.2
model4.1 <- ergm(friend_net ~ edges + nodeofactor("sex") + nodeifactor("smoke") + nodematch("activity"))
summary(model4.1)


#test with terms from 3.2
model4.1with <- ergm(friend_net ~ edges + nodeofactor("sex") + nodeifactor("smoke") + nodematch("activity")
                     + nodematch("sex") + mutual 
                     + gwesp(decay = 0.3, fixed = TRUE)
                     + gwodegree(decay = 0.3, fixed = TRUE) + gwidegree(decay = 0.3, fixed = TRUE))
summary(model4.1with)




#task 4.2
set.seed(1)

#add smoke attribute to friend_net
smoke = attributes$smoke
set.vertex.attribute(friend_net, "smoke", smoke)

model4.2 <- ergm(friend_net ~ edges + nodemix("smoke"))
summary(model4.2)

