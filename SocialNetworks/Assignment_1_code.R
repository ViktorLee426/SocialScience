
library(sna)
library(igraph)
library(ggraph)

# Task 1: Describe and plot a social network ####

## a) load the affective network (2400 affective w1.csv) and the gender data (2400 sex.csv) ####

affective = read.csv("2400_affective_w1.csv")
affective=as.matrix(affective[,-1])
rownames(affective)=colnames(affective)


sex = read.csv("2400_sex.csv")
id = sex[,1]
sex=as.matrix(sex[,-1])
rownames(sex)=id



## b) recode the affective network of wave 1 into a friendship networks. (Hint: the value +2 stands #### 
#       for friendship ties – see the data description for details) 

# In order to recode the affective network into a friendship network we need to identify what are the friendship
# relationships (negative and possitive) and assign a weight to them representing the strength of their relationship
# Since the relationships vary from -2 to -1, we can assign different weights as follows:
# Let us fix person i and j, and denote the matrix of affection by A
#    - If 0<A[i,j]<A[j,i] (j is more affective towards i), then we assign an arrow from j to i.
#    - If A[i,j]<A[j,i]<0 (i hates more j), then we assign an arrow from i to j.
#    - If A[i,j]=A[j,i] (they think the same of each other) we assign an arrow both ways
#    - If A[i,j] = 0, then there is no realationship -> no arrow from i to j.
#    - 



# I will use a basic matrix for now, and change it later if its needed.
# Let us suppose that there is only friendship if there is a 2 in the matrix

friendship=affective
friendship[!friendship==2]=0
friendship[is.na(friendship)]=0
friendship[friendship==2]=1


# c) calculate basic network descriptives for this friendship network

# Network size (i.e., number of nodes):
ncol(friendship)
# Density: 
gden(friendship)
# average degree:
sum(friendship, na.rm = T)/dim(friendship)[1]   # how many ties each person has on average   
# reciprocity ratio:
grecip(friendship, measure = "dyadic.nonnull")
# gender composition in class
table(sex)
# count of same gender ties (i.e. both nodes have the same gender)
index_men=which(sex==1)
index_women=which(sex==2)
friendship_women=friendship[index_women,index_women]
friendship_men=friendship[index_men,index_men]

friendship_women
lower_women=friendship_women*lower.tri(friendship_women, diag = FALSE)
lower_women[upper.tri(friendship_women)]=NA
diag(lower_women)=NA
upper_women=friendship_women*upper.tri(friendship_women, diag = FALSE)
upper_women[lower.tri(friendship_women)]=NA
diag(upper_women)=NA
sum(lower_women==t(upper_women),na.rm=TRUE)


friendship_men
lower_men=friendship_men*lower.tri(friendship_men, diag = FALSE)
lower_men[upper.tri(friendship_men)]=NA
diag(lower_men)=NA
upper_men=friendship_men*upper.tri(friendship_men, diag = FALSE)
upper_men[lower.tri(friendship_men)]=NA
diag(upper_men)=NA
sum(lower_men==t(upper_men),na.rm=TRUE)



# Other meassure of your choice: Degree centrality
degree <- sna::degree(friendship, cmode = 'freeman') 
# freeman specifies the total of incoming and outgoing ties
max(degree) # using max(), we see the highest value
degree
degree == max(degree) 
which(degree == max(degree) )


# Briefly interpret the measures (where sensible)

#### TO BE DONE ####



# d) plot the friendship network.
#   .- The plot has to be informative
#   .- Color the nodes according to the gender of the person
#   .- the node size should be proportinal to the centrality meassure of your choice
friendship.igraph <- graph_from_adjacency_matrix(friendship,
                                             mode = "directed")
layout.graph <- create_layout(friendship.igraph, 
                              layout = 'fr')
edges <- geom_edge_link(alpha = .5)
edges <- geom_edge_link(alpha = .5, 
                        arrow = arrow(length = unit(1,'mm'), 
                                      type = "closed")) 
edges <- geom_edge_link2(alpha = .7,
                         arrow = arrow(length = unit(1, 'mm'), 
                                       type = "closed"),
                         end_cap = circle(2, 'mm')) # you may have to play with this a bit to get it right // add some space around the node in 
# order to better visualize the arrow-heads
degr.plot <- 
  ggraph(layout.graph)+ 
  edges+ 
  geom_node_point(aes(colour = sex,
                      size = degree))+
  theme_graph()+ 
  theme(legend.position = "none") + # remove the legends completely
  ggtitle("Degree centrality")
degr.plot


# e) Now also include the trust network (2400 trust w1.csv) and plot them both in one network plot (i.e., on top of each other).
# creating two random networks for this example
trust= read.csv("2400_trust_w1.csv")
trust=as.matrix(trust[,-1])
rownames(trust)=colnames(trust)

# set.seed(56)

net1 <- friendship
net2 <- trust
# a new matrix of the right size, prefilled with ties in net1
net3 <- net1
# where ties exist in the second network, we give a new value
net3[net2 == 1] <- 2
# where ties exist in both networks, we replace with a new, unique value
net3[net1 == 1 & net2 == 1] <- 3
# now you have 4 unique values in net3:
table(net3)

# 0 = no ties,
# 1 = tie in network 1,
# 2 = tie in network 2,
# 3 = tie in both networks

g=graph_from_adjacency_matrix(net3,mode="directed",weighted=TRUE)
layout.graph <- create_layout(g, 
                              layout = 'fr')
edges <- geom_edge_link2(alpha = .5,
                         arrow = arrow(length = unit(1, 'mm'), 
                                       type = "closed"),
                         end_cap = circle(2, 'mm'),aes(width = weight))
degr.plot <- 
  ggraph(layout.graph)+ 
  edges+ 
  geom_node_point(aes(colour = sex,
                      size = degree))+
  scale_edge_width(range = c(0.1, 1))+
  theme_graph()+ 
  theme(legend.position = "none") + # remove the legends completely
  ggtitle("Degree centrality")
degr.plot



# f) How large is overlap between the two networks? (i.e., how many ties between two people are present
#     in both networks)
sum(net1==net2,na.rm = TRUE)


# g) In a short paragraph (max. 250 words), describe what you see in the network plot and comment on the overlap between the two networks.




