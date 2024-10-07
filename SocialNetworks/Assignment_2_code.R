#  We get the matrix of the graphs from the High-Tech company: advice and friendship
library(mully)
library(RColorBrewer)
library(sna)
library(igraph)
library(ggraph)
library(networkdata)
library(ggforce)
net1 <-as.matrix(ht_advice, matrix.type = c("adjacency"))
net2 <-as.matrix(ht_friends, matrix.type = c("adjacency"))

# a new matrix of the right size, prefilled with ties in net1
net3 <- net1
# where ties exist in the second network, we give a new value
net3[net2 == 1] <- 2
# where ties exist in both networks, we replace with a new, unique value
net3[net1 == 1 & net2 == 1] <- 3
# So we have that the least weight will correspond to the interactions where only
# advice is being asked, but there is no friendship. In the case of weight 2, 
# the person asking for advice considers the other person as a friend, but does not ask for advice from him.

g=graph_from_adjacency_matrix(net3,mode="directed",weighted=TRUE)

# We set the attributes from the original graphs
g <- set_vertex_attr(g, "age", value = vertex_attr(ht_advice, "age"))
g <- set_vertex_attr(g, "level", value = vertex_attr(ht_advice, "level"))
g <- set_vertex_attr(g, "dept", value = vertex_attr(ht_advice, "dept"))


set.seed(1)
layout.graph <- create_layout(g, layout = 'stress')
# we color the edges and make the width proportional to the weight in order to make better
edges <- geom_edge_link2(alpha = .45,
                         arrow = arrow(length = unit(1.5, 'mm'),
                                       type = "closed"),
                         end_cap = circle(2, 'mm'),aes(width = weight,color=factor(weight)))

#nodes are colored according to gender and size proportional to the degree
nodes <- geom_node_point(aes(shape=factor(level),size = degree(g,mode="in"),color = age)) 
plotlabs <- labs(shape = "Hierarchy level ",size = "In-degree",color ="Age")
# we draw the plot
degr.plot <-
  ggraph(layout.graph)+
  labs(edge_color = "Friendship") +
  edges+
  nodes+
  geom_node_text(aes(label = dept,size =
                       15),nudge_x=-0.05, nudge_y=0.05)+
  plotlabs+
  scale_edge_color_manual(values =c("yellow",'green',"blue"),label = c("Only asks for advice","Consider as friend but not ask for advice","Asks for advice and friendship")) +
  scale_color_gradient( low = "blue",high = "red")+
  scale_edge_width(guide="none",range = c(0.5, 0.8))+
  theme_graph()+theme(legend.key.size = unit(1, 'cm'),legend.title = element_text(size=15),  legend.text = element_text(size=13))+
  scale_shape_discrete(label = c("CEO","Vice President","Manager"))+
  ggtitle("High-tech Managers (Advice)") 
degr.plot


####2a  walktrap algorithm:
library("concaveman")
ht.wt <- cluster_walktrap(ht_advice)

length(ht.wt) #How many communities were identified by the algorithm?
sizes(ht.wt) #What are their sizes?
member.wt <- membership(ht.wt) #Who belongs to which community?
modularity(ht.wt)

set.seed(3)
g <- ggraph(ht_advice, layout = 'fr') +
  geom_edge_arc(color="gray89",strength=0.15,arrow=arrow(length=unit(2.5,'mm')),end_cap=circle(3, 'mm')) + 
  labs(shape="Hierarchy level",color ="Age", size="Indegree") +
  geom_node_point(aes(shape=factor(level) ,size = igraph::degree(ht_advice,mode="in"),color=age)) +
  theme_graph() + 
  #scale_color_gradient2( low = "blue",high = "red")+
  geom_node_text(aes(label = dept,size = 15),nudge_x=-0.05, nudge_y=0.05)+
  scale_shape_discrete(label = c("CEO","Vice President","Manager"))+
  ggtitle("High-tech Managers (Advice)") 
g +geom_mark_hull(aes(x,y,fill = factor(member.wt)), 
                  colour = NA,show.legend = FALSE)+scale_fill_brewer(palette = "Dark2")


########2a edge-betweenness algorithm:

ht.wt <- cluster_edge_betweenness(ht_advice)

length(ht.wt) #How many communities were identified by the algorithm?
sizes(ht.wt) #What are their sizes?
member.wt <- membership(ht.wt) #Who belongs to which community?
modularity(ht.wt)

nb.cols <- 21
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

set.seed(3)
g <- ggraph(ht_advice, layout = 'fr') +
  geom_edge_arc(color="gray89",strength=0.15,arrow=arrow(length=unit(2.5,'mm')),end_cap=circle(3, 'mm')) + 
  labs(shape="Hierarchy level",color ="Age", size="Indegree") +
  geom_node_point(aes(shape=factor(level) ,size = igraph::degree(ht_advice,mode="in"),color=age)) +
  theme_graph() + 
  #scale_color_gradient2( low = "blue",high = "red")+
  geom_node_text(aes(label = dept,size = 15),nudge_x=-0.05, nudge_y=0.05)+
  scale_shape_discrete(label = c("CEO","Vice President","Manager"))+
  ggtitle("High-tech Managers (Advice)") 
g +geom_mark_hull(aes(x,y,fill = factor(member.wt)), 
                  colour = NA,show.legend = FALSE)+scale_fill_manual(values = mycolors)




###2b test1
rm(list = ls()) #clear the environment

#rename graph for convenience
ht <- ht_advice

#get params for Erdos-Renyi simulation & print 
nodes <- length(ht) # Nodes on the original graph
edges <- gsize(ht) # Edges on the original graph
nodes
edges


# First we compute the proportion of same-department edges on the original graph:
same_department_edges_original=sum(V(ht)$dept[ends(ht, E(ht))[,1]] == V(ht)$dept[ends(ht, E(ht))[,2]])/edges

# Secondly, we generate 10000 Erdos-Renyi graphs and compute the same quantity of interest
num_iterations <- 10000
same_department_edges_sim=numeric(num_iterations)

set.seed(1)

#main loop 
for (i in 1:num_iterations) {
  #create erdos-renyi stylised graph with 21 nodes and 190 edges in G(N, M) mode, directed and without loops
  e <- erdos.renyi.game(
    nodes,
    edges,
    "gnm",
    directed = TRUE,
    loops = FALSE
  )
  # Fix the original attributes of the nodes:
  e <- set_vertex_attr(e, "dept", value = vertex_attr(ht, "dept"))
  # And compute the proportion of same-department edges on the simulated graphs:
  same_department_edges_sim[i] = sum(V(e)$dept[ends(e, E(e))[,1]] == V(e)$dept[ends(e, E(e))[,2]])/edges
}

# Create a histogram and plot the results, add in dashed line the position of the original value of
# the quantity of interest 
df = data.frame(same_department_edges=same_department_edges_sim)
g = ggplot(df,aes(x=same_department_edges)) + geom_histogram(binwidth = 0.01)
g + geom_vline(xintercept=same_department_edges_original, linetype="dashed", color = "#a91b0d",size=2)+
  annotate("text", x=0.25, y=1600, label="Same department edges original", angle=90,color = "#a91b0d")


# Compute p-value
cdf <- ecdf(sort(same_department_edges_sim)) 
cdf_value <- cdf(same_department_edges_original)
(p_value <- 1 - cdf_value)
# 0.1824


#####2b test2

rm(list = ls()) #clear the environment

#rename graph for convenience
ht <- ht_advice

#get params for Erdos-Renyi simulation & print 
nodes <- length(ht) # Nodes on the original graph
edges <- gsize(ht) # Edges on the original graph
nodes
edges


# First we compute the proportion of same-level edges on the original graph:
same_level_edges_original=sum(V(ht)$level[ends(ht, E(ht))[,1]] == V(ht)$level[ends(ht, E(ht))[,2]])/edges

# Secondly, we generate 10000 Erdos-Renyi graphs and compute the same quantity of interest
num_iterations <- 10000
same_level_edges_sim=numeric(num_iterations)


set.seed(1)

#main loop 
for (i in 1:num_iterations) {
  #create erdos-renyi stylised graph with 21 nodes and 190 edges in G(N, M) mode, directed and without loops
  e <- erdos.renyi.game(
    nodes,
    edges,
    "gnm",
    directed = TRUE,
    loops = FALSE
  )
  e <- set_vertex_attr(e, "level", value = vertex_attr(ht, "level"))
  
  same_level_edges_sim[i] = sum(V(e)$level[ends(e, E(e))[,1]] == V(e)$level[ends(e, E(e))[,2]])/edges
}

# Create histogram and compute p-value
df = data.frame(same_level_edges=same_level_edges_sim)
g = ggplot(df,aes(x=same_level_edges)) + geom_histogram(binwidth = 0.01)
g + geom_vline(xintercept=same_level_edges_original, linetype="dashed", color = "#a91b0d",size=2)+
  annotate("text", x=0.54, y=1400, label="Same level edges original", angle=90,color = "#a91b0d")


# Compute p-value
cdf <- ecdf(sort(same_level_edges_sim))
cdf_value <- cdf(same_level_edges_original)
(cdf_value)   # Probability on the left side -> p-value
# 0.029

