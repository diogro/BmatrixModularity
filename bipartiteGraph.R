library(igraph)
library(evolqg)
library(EvomodR)


burn = ReadRawPop(folder = "/home/diogro/projects/evomod/c-gsl/output/burn_in", "burnin")
stab = ReadRawPop(folder = "/home/diogro/projects/evomod/c-gsl/output/CorrelatedBurin_Results/burn_in/", "stab")
divsel = ReadRawPop(folder = "/home/diogro/projects/evomod/c-gsl/output/burn_in/DivSel-0.004", "DivSel-0.004")


divsel$B[,,2] %*% t(divsel$B[,,2])
t(divsel$B[,,2]) %*% divsel$B[,,2]


plotAdj <- function (bg) {
  pr=bipartite.projection(bg)
  l <-layout.reingold.tilford(pr$proj2) 
  l[1:5,1] = c(0, -1, -2, -1, 0)
  l[1:5,2] = 1:5
  l[6:10,1] = c(1, 2, 3, 2, 1)
  l[6:10,2] = 1:5  
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
sg = graph.incidence(stab$B[,,3])
dg = graph.incidence(divsel$B[,,3])
b = plotAdj(bg)
s = plotAdj(sg)
d = plotAdj(dg)


hist(apply(burn$B, 1:2, mean))
hist(apply(divsel$B, 1:2, mean))

get.adjacency(b$proj2,sparse=FALSE,attr="weight")
get.adjacency(s$proj2,sparse=FALSE,attr="weight")
get.adjacency(d$proj2,sparse=FALSE,attr="weight")


mod = cbind(rep(c(1, 0), each = 5),
            rep(c(0, 1), each = 5))
# TestModularity(b, mod)
# TestModularity(d, mod)

da = graph.adjacency(d$proj2, mode= "undirected", weighted = TRUE)
gg.com <- fastgreedy.community(d$proj2)
V(d)$color <- gg.com$membership + 1
plot(da)

