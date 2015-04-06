library(igraph)
library(evolqg)
library(EvomodR)


burn = ReadRawPop(folder = "/home/diogro/projects/evomod/c-gsl/output/burn_in", "burnin")
divsel = ReadRawPop(folder = "/home/diogro/projects/evomod/c-gsl/output/burn_in/DivSel-0.004", "DivSel-0.004")


plotAdj <- function (bg) {
  pr=bipartite.projection(bg)
  l <-layout.reingold.tilford(pr$proj2) 
  l[1:5,1] = 0
  l[1:5,2] = 1:5
  l[6:10,1] = 1
  l[6:10,2] = 1:5  
  plot(pr$proj2, 
       edge.width=E(pr$proj2)$weight^2/30000,
       edge.color="black", 
       vertex.label=V(pr$proj2)$name,
       layout = l)
  return(get.adjacency(pr$proj2,sparse=FALSE,attr="weight"))
}

bg = graph.incidence(burn$B[,,1])
dg = graph.incidence(divsel$B[,,1])
b = plotAdj(bg)
d = plotAdj(dg)
mod = cbind(rep(c(1, 0), each = 5),
            rep(c(0, 1), each = 5))
TestModularity(b, mod)
TestModularity(d, mod)