library(igraph)
library(evolqg)
library(EvomodR)
library(ggplot2)
library(reshape2)
library(bipartite)
library(readr)
library(matrixcalc)


burn = ReadRawPop(folder = "/home/diogro/projects/evomod/c-gsl/output/burn_in", "burnin")
#stab = ReadRawPop(folder = "/home/diogro/projects/evomod/c-gsl/output/CorrelatedBurin_Results/burn_in/", "stab")
divsel = ReadRawPop(folder = "/home/diogro/projects/evomod/c-gsl/output/burn_in/DivSel-0.004", "DivSel-0.004")

divsel$Graph = alply(divsel$B, 3, graph.incidence)
burn$Graph = alply(burn$B, 3, graph.incidence)

divsel$B[,,2] %*% t(divsel$B[,,2])
t(divsel$B[,,2]) %*% divsel$B[,,2]


plotAdj <- function (bg) {
  pr=bipartite.projection(bg)
  l <-layout.circle(pr$proj2) 
  gg.com <- fastgreedy.community(pr$proj2)
  V(pr$proj2)$color <- gg.com$membership + 1
  plot(pr$proj2, 
       edge.width=(E(pr$proj2)$weight/200)^2,
       edge.color="black", 
       vertex.label=V(pr$proj2)$name,
       vertex.color = V(pr$proj2)$color,
       layout = l)
  return(pr)
}

bg = graph.incidence(burn$B[,,3])
#sg = graph.incidence(stab$B[,,3])
dg = graph.incidence(divsel$B[,,3])
b = plotAdj(bg)
#s = plotAdj(sg)
d = plotAdj(dg)

hist(apply(burn$B, 1:2, mean))
hist(apply(divsel$B, 1:2, mean))

get.adjacency(b$proj2,sparse=FALSE,attr="weight")
#get.adjacency(s$proj2,sparse=FALSE,attr="weight")
get.adjacency(d$proj2,sparse=FALSE,attr="weight")


mod = cbind(rep(c(1, 0), each = 5),
            rep(c(0, 1), each = 5))
TestModularity(get.adjacency(b$proj2, sparse = FALSE, attr="weight"), mod)
TestModularity(get.adjacency(d$proj2, sparse = FALSE, attr="weight"), mod)

degree = data_frame(index = 1:513, 
                    burn = c(degree.distribution(bg), rep(0, 15)), 
                    divsel = degree.distribution(dg))

networklevel(burn$B[,,3], index = 'H2')

degreedistr(burn$B[,,1])
degreedistr(divsel$B[,,10])

x = data_frame(modular = E(d$proj2)$weight, 
               burn    = E(b$proj2)$weight)
ggplot(melt(x), aes(value, color = variable)) + geom_histogram(binwidth = 10, aes(fill = variable)) + theme_classic()

num_graphs = 1000
shared_edges = matrix(0, num_graphs, num_graphs)
for (i in 1:num_graphs){
   for (j in i:num_graphs){
     shared_edges[i, j] <- sum(E(divsel$Graph[[i]]) %in% E(divsel$Graph[[j]]))
     shared_edges[j, i] <- shared_edges[i, j]/ecount(divsel$Graph[[i]])
   }
}
shared_edges_burn = matrix(0, num_graphs, num_graphs)
for (i in 1:num_graphs){
  for (j in i:num_graphs){
    shared_edges_burn[i, j] <- sum(E(burn$Graph[[i]]) %in% E(burn$Graph[[j]]))
    shared_edges_burn[j, i] <- shared_edges_burn[i, j]/ecount(burn$Graph[[i]])
  }
}
x = data_frame(burn    = shared_edges_burn[upper.tri(shared_edges_burn)], 
               modular = shared_edges     [upper.tri(shared_edges)])
ggplot(melt(x), aes(variable, value, color = variable)) + geom_violin() + theme_classic()
