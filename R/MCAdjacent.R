#' Adjacency Matrix for Maximal Cliques
#'
#' This function computes the max clique
#' adjacency matrix, whereby
#' two MCs are adjacent if  there is a junction tree in
#' which they are adjacent
#'
#' @param JTr a junction tree computed by Kruskal's algorithm
#' @param W the weighted clique graph (edge weights give number of vertices in
#' intersections)
#' @return Adj - lower triangular matrix 1 if the corresponding
#' MCs are adjacent in some junction tree and 0 otherwise
#'
#'
#'@export
#'
#'
MCAdjacent <- function(JTr,W){
  #JTr is lower triangular giving edges of junction tree
  #W is lower triangular, for j > i W[j,i] number of vertices
  #in intersection C[j],C[i]

  #first task is to get a directed version of the J Tree
  #any will do, so take 1 as the source vertex.

  Wsym = W + t(W)
  d = ncol(JTr)
  DT = matrix(0,d,d)
  #initialisation DT[i,j] means from i to j
  JTOT = JTr + t(JTr)
  DT[1,] = JTOT[1,]
  JTOT[,1]=0
  JTOT[1,]=0
  A = which(DT[1,]==1)

  repeat
  {
    B=c()

    DT[A,]= JTOT[A,]
    JTOT[A,]=0
    JTOT[,A]=0
    for(j in A)
    {
    B = c(B,which(DT[j,]==1))
    }
  if(sum(JTOT)==0)
  {
    break
  }
    A = B
  }

#DT is a directed version of the junction tree, all paths pointing away from 1
  #each vertex has at most one parent.

#Now find paths back to the root
  A = vector(mode="list",length=d)
  f = function(j){A[[j]] = j}
  A = lapply(1:d,f)
  for(k in 2:d)
  {
    a = k
    repeat
    {
    b = which(DT[,a]==1)
    if(length(b)==0)
    {
      break
    }
    A[[k]] <- c(A[[k]],b)
    a=b
    }
  }

  Adj <- matrix(0,d,d)

 for(i in 2:d)
 {
   for(j in 1:(i-1))
   {
     S = intersect(A[[i]],A[[j]])
     if(i %in% S)
     {
       path = c(setdiff(A[[j]],A[[i]]),i)
     }
     if(j %in% S)
     {
       path = c(setdiff(A[[i]],A[[j]]),j)
     }
     if(length(intersect(S,c(i,j)))==0)
     {
       Bi = setdiff(A[[i]],A[[j]])
       Bj = rev(setdiff(A[[j]],A[[i]]))
       node = S[1]
       path = c(Bi,node,Bj)
     }
       len = length(path)-1
       f <- function(k){Wsym[path[k],path[k+1]]}
       f1 = min(sapply(1:len,f))
       if(Wsym[i,j] >= f1)
       {
         Adj[i,j]=1
       Adj[j,i]=1
       }


   }
 }
Adj[upper.tri(Adj,diag=TRUE)]=0
   return(Adj)
}
