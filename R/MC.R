#' Computes the maximal cliques of  decomposable graph
#'
#' @param M the lower triangular matrix defining the undirected decomposable graph.
#'
#' @return MCli  list of maximal cliques
#' @export

MC <- function(M)
  {
  M = as.matrix(M)
MCli <- list()
MCli[[1]] = vector(mode="integer",length=0L)

#maximal cliques - add on the new cliques from each round

n = nrow(M)
if(nrow(M)==1)
  {
MCli[[1]] <- 1
return(MCli)
}

SU = sum(M)
if(SU == 0.5*n*(n-1))
{
  MCli[[1]] <- seq(1,n)
  return(MCli)
}

#the if statements deals with M itself being an MC. We have now dealt with
#(a) the whole CC is a single vertex and (b) the whole CC is complete.

counterMC <- 0

g <- function(x)
{
  A = B[[x]]
  #vertex x together with neighbours
  I = M[A,A]
  #adjacency matrix for subgraph
  m <- length(A)
  h=1*(isTRUE(sum(I)== 0.5*m*(m-1)))
  #if so, A is complete
  return(h)
}
#g decides if the vertex is simplicial at this round.
#if so, vertex plus neighbours is a max clique unless
#it is (part of a) separator (complete set in final round)


repeat
{
B <- nhd(M)$C
  #list of (vertex plus neighbours) those for
  #simplicial vertices are maximal cliques

  g1=sapply(1:nrow(M),g)
  MCL = which(g1==1)

#g1[x] == 1 if and only if B[[x]] is a MC
  #which implies that x is a simplicial vertex
#hence MCL lists the simplicial vertices at this round

hC <- function(x){B[[x]]}
Cnew = unique(lapply(MCL,hC))
mark <- vector(mode="integer",length=length(Cnew))
for(i in 1:length(Cnew))
{
f <- function(j){1*(all(Cnew[[i]] %in% MCli[[j]]))}
f1 <- sapply(1:(length(MCli)),f)
mark[i] = 1 - max(f1)
#so mark[i]=1 if and only if all f1 are 0 which
#means that Cnew[i] is a new MC
}
set = which(mark == 1)

if(length(set) >= 1)
{
for(j in 1:length(set))
{
MCli[[counterMC+j]] = Cnew[[set[j]]]
}
MCli <- unique(MCli)
}
counterMC <- length(MCli)


M[MCL,] = 0
M[,MCL]=0

#having found the MCs at this round, remove the simplicial nodes


if(sum(M)==0)
{
return(MCli)
#and we terminate if, having removed
#the simplicial nodes, the graph is empty}
}
}
}
