library(RSiena)
library(sna)
library(here)
#install.packages(c("igraph", "ggraph", "ggplot2", "dplyr"))
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)





setwd("C:/Users/zheng/Desktop/ForthSemester/ND/my topic/Glasgow_data")
source("printSiena.R")
source("siena07ToConverge.R")

print(getwd())


load('Glasgow-demographic.RData')
load('Glasgow-friendship.RData')
load('Glasgow-geographic.RData')
load('Glasgow-lifestyle.RData')
load('Glasgow-selections.RData')
load('Glasgow-substances.RData')
load('Glasgow-various.RData')

friendship1_filtered <- friendship.1[selection129,selection129]
friendship2_filtered <- friendship.2[selection129,selection129]
friendship3_filtered <- friendship.3[selection129,selection129]

friendship1_filtered[friendship1_filtered == 2] <- 1
friendship2_filtered[friendship2_filtered == 2] <- 1
friendship3_filtered[friendship3_filtered == 2] <- 1



tobacco_filtered <- tobacco[selection129,]
alcohol_filtered <- alcohol[selection129,]
sex_filtered <-sex.F[selection129]
age_filtered <-age[selection129]
money_filtered <- money[selection129,]
romantic_filtered <- romantic[selection129,]
familysmoking_filtered <- familysmoking[selection129,]



###############################################
# Create graph objects for each wave
igraph_wave1 <- graph_from_adjacency_matrix(friendship1_filtered, mode = "directed", diag = FALSE)
igraph_wave2 <- graph_from_adjacency_matrix(friendship2_filtered, mode = "directed", diag = FALSE)

# Set seed for reproducibility and create layout for wave 1
set.seed(1)
mylayout_wave1 <- layout_with_fr(igraph_wave1)

# Assign the same layout to wave 2
mylayout_wave2 <- mylayout_wave1

# Visualization 1: Drinking behavior of wave 1
edges_wave1 <- geom_edge_link(alpha = .7, arrow = arrow(length = unit(1, 'mm'), type = "closed"), end_cap = circle(2, 'mm'))

nodes_wave1 <- geom_node_point(aes(x = mylayout_wave1[,1], y = mylayout_wave1[,2], 
                                   colour = as.numeric(alcohol_filtered[, "t1"]),
                                   shape = as.factor(sex_filtered),
                                   size = money_filtered[,1]), show.legend = TRUE)

plotlabs <- labs(colour = "Alcohol", shape = "Gender", size = "Money")
plottitle_wave1 <- ggtitle("Drinking behavior of wave 1")

ggraph(igraph_wave1, layout = mylayout_wave1) +
  edges_wave1 +
  nodes_wave1 +
  scale_color_gradient(low = "yellow", high = "red", na.value = "grey") +
  plotlabs +
  plottitle_wave1 +
  theme_graph()

# Visualization 2: Drinking behavior of wave 2
edges_wave2 <- geom_edge_link(alpha = .7, arrow = arrow(length = unit(1, 'mm'), type = "closed"), end_cap = circle(2, 'mm'))

nodes_wave2 <- geom_node_point(aes(x = mylayout_wave2[,1], y = mylayout_wave2[,2], 
                                   colour = as.numeric(alcohol_filtered[, "t2"]),
                                   shape = as.factor(sex_filtered),
                                   size = money_filtered[,2]), show.legend = TRUE)

plottitle_wave2 <- ggtitle("Drinking behavior of wave 2")

ggraph(igraph_wave2, layout = mylayout_wave2) +
  edges_wave2 +
  nodes_wave2 +
  scale_color_gradient(low = "yellow", high = "red", na.value = "grey") +
  plotlabs +
  plottitle_wave2 +
  theme_graph()

###############################################

friendship <- sienaDependent(array(c(friendship1_filtered,friendship2_filtered, friendship3_filtered ), dim = c(129,129,3)))

tobacco <-varCovar(tobacco_filtered)
alcohol <-varCovar(alcohol_filtered)
sex <- coCovar(sex_filtered)
age <- coCovar(age_filtered)
money <-varCovar(money_filtered)
romantic <- varCovar(romantic_filtered)
familysmoking <- varCovar(familysmoking_filtered)

mydata <- sienaDataCreate(friendship, tobacco, alcohol, sex, age, money, romantic, familysmoking)
mydata

print01Report(mydata, modelname = "glasgow")

myeff <- getEffects(mydata)

myeff <- includeEffects(myeff, egoX, altX, sameX, interaction1 = "sex")

myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "tobacco")

myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "alcohol")

myeff


myAlgorithm <- sienaAlgorithmCreate(projname = "glasgow_res", nsub = 4, n3 = 2000, seed = 1908)

model0 <- siena07(myAlgorithm,
                  data = mydata, effects = myeff, returnDeps = TRUE,
                  useCluster = TRUE, nbrNodes = 3
)

model0

model0 <- siena07(myAlgorithm,
                  data = mydata, effects = myeff, returnDeps = TRUE, prevAns = model0,
                  useCluster = TRUE, nbrNodes = 3
)

model0

siena07ToConvergence(myAlgorithm, dat = mydata, eff = myeff)

model0

printSiena(model0)


# Goodness of fit
# Indegree distribution
gof1.id <- sienaGOF(model0, verbose = FALSE,
                    varName = "friendship", IndegreeDistribution)
# Outdegree distribution
gof1.od <- sienaGOF(model0, verbose = FALSE,
                    varName = "friendship", OutdegreeDistribution)
# Triad census
gof1.tc <- sienaGOF(model0, verbose = FALSE,
                    varName = "friendship", TriadCensus)
# Geodesic distance
GeodesicDistribution <- function(i, data, sims, period, groupName,
                                 varName, levls = c(1:5, Inf), cumulative = TRUE, ...) {
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  require(sna)
  a <- sna::geodist(symmetrize(x))$gdist
  if (cumulative) {
    gdi <- sapply(levls, function(i) {
      sum(a <= i)
    })
  }
  else {
    gdi <- sapply(levls, function(i) {
      sum(a == i)
    })
  }
  names(gdi) <- as.character(levls)
  gdi
}
gof1.gd <- sienaGOF(model0, verbose = FALSE,
                    varName = "friendship", GeodesicDistribution)
plot(gof1.id)
plot(gof1.od)
plot(gof1.tc, center = TRUE, scale = TRUE)
plot(gof1.gd)

#####################################################
# Autocorrelation: Moran index
moran1 <- nacf(friendship1_filtered, tobacco_filtered[, 1], lag.max = 1, type = "moran",
               neighborhood.type = "out", mode = "digraph")
moran2 <- nacf(friendship2_filtered, tobacco_filtered[, 2], lag.max = 1, type = "moran",
               neighborhood.type = "out", mode = "digraph")
moran3 <- nacf(friendship3_filtered, tobacco_filtered[, 3], lag.max = 1, type = "moran",
               neighborhood.type = "out", mode = "digraph")

autocorr <- rbind(moran1,moran2,moran3)
autocorr[,2]

alcohol_filtered_imputed <- alcohol_filtered[, 1]
alcohol_filtered_imputed[is.na(alcohol_filtered_imputed)] <- mean(alcohol_filtered_imputed, na.rm = TRUE)
moran1 <- nacf(friendship1_filtered, alcohol_filtered_imputed, lag.max = 1, type = "moran",
               neighborhood.type = "out", mode = "digraph")

alcohol_filtered_imputed <- alcohol_filtered[, 2]
alcohol_filtered_imputed[is.na(alcohol_filtered_imputed)] <- mean(alcohol_filtered_imputed, na.rm = TRUE)
moran2 <- nacf(friendship2_filtered, alcohol_filtered_imputed, lag.max = 1, type = "moran",
               neighborhood.type = "out", mode = "digraph")

alcohol_filtered_imputed <- alcohol_filtered[, 3]
alcohol_filtered_imputed[is.na(alcohol_filtered_imputed)] <- mean(alcohol_filtered_imputed, na.rm = TRUE)
moran3 <- nacf(friendship2_filtered, alcohol_filtered_imputed, lag.max = 1, type = "moran",
               neighborhood.type = "out", mode = "digraph")

autocorr <- rbind(moran1,moran2,moran3)
autocorr[,2]

smokingbeh <- sienaDependent(tobacco, type = "behavior")
#drinkingbeh <- sienaDependent(alcohol, type = "behavior")

mydata <- sienaDataCreate(friendship, smokingbeh, sex, familysmoking)

print01Report(mydata, modelname = 'model1_vars')

myeff <- getEffects(mydata)
myeff

myeff <- includeEffects(myeff, sameX, interaction1 = "sex")


# behaviour effect explaining friendship evolution
myeff <- includeEffects(myeff, egoX, altX, simX,
                        name = "friendship", interaction1 = "smokingbeh")

#myeff <- includeEffects(myeff, egoX, altX, simX,
                        #name = "friendship", interaction1 = "drinkingbeh")

myeff <- includeEffects(myeff,  outdeg, indeg, avSim,
                        name = "smokingbeh", interaction1 = "friendship")

#myeff <- includeEffects(myeff, outdeg, indeg, avSim,
#                        name = "drinkingbeh", interaction1 = "friendship")

myeff <- includeEffects(myeff, effFrom,
                        name = "smokingbeh", interaction1 = "sex")

myeff


myAlgorithm <- sienaAlgorithmCreate(
  projname = "CoevKnecht",
  nsub = 4, n3 = 2000, seed = 1908
)

# Model estimation
model1 <- siena07(myAlgorithm,
                      data = mydata, effects = myeff,
                      returnDeps = TRUE, useCluster = TRUE, nbrNodes = 4, prevAns = TRUE
)
model1

siena07ToConvergence(myAlgorithm, dat = mydata, eff = myeff)

model1 <- siena07(myAlgorithm,
                  data = mydata, effects = myeff,
                  returnDeps = TRUE, useCluster = TRUE, nbrNodes = 4, prevAns = model.coev)
model1


printSienaCoev(model1)
