#' Edges which can be added to a decomposable graph to give a decomposable graph
#'
#'For an input M, lower triangular matrix of 1's and 0's, where 1's are
#'undirected edges of a DECOMPOSABLE graph, addedgemat(M) returns lower
#'triangular matrix where entry [y,x] is 1 if and only if there is no edge
#'and adding the edge between y and x returns a decomposable graph.
#'
#' @param M lower triangular matrix of 1's and 0's, where 1's denote presence
#' of edges in a decomposable graph
#'
#'
#' @return
#' Out - Lower triangular matrix where 1 denotes edge not present, but can be
#' added to give a decomposable graph.
#'
#' @export



addedgemat <- function(M){

  #recall that an edge x-y can be added if and only if
  #x and y belong to MCs which are adjacent in some junction tree.
  #M is an n by n lower triangular adjacency matrix
  #for an undirected decomposable graph M[x,y]= 1 if there is an
  #edge x-y and x > y
  #m is the maximal clique size permitted
  #the output is a matrix nrow(M) by ncol(M) with 1's where an edge can
  #be added to give a decomposable graph with max clique size <= m and 0 where it can't
  #(either the graph already has an edge there, or else the new edge creates a graph
  #that is not decomposble).

  #simp(M)$MC gives the maximal cliques of the graph: S the maximal separators
  #(i.e. no repetition; also if S1 is a subset of S2, only S2 recorded)

  C <- MC(M)
  #this gives the collection of maximal cliques of the graph defined by M

  W = WeightedCliqueGraph(C)
  JTr = kruskal(C)
  J <- MCAdjacent(JTr,W)
  #this gives the adjacency matrix for the maximal cliques

  AUE <- matrix(0,nrow(M),nrow(M))
  Nsym = matrix(0,nrow(M),nrow(M))

  for(i in 1:(nrow(J)-1))
      {
        for(j in (i+1):nrow(J))
        {
          if(J[j,i]==1)
          {
            A = setdiff(C[[i]],C[[j]])
        B = setdiff(C[[j]],C[[i]])
        Nsym[A,B]=1
        Nsym[B,A]=1
          }
        }
  }
   for(i in 2:nrow(M))
   {
     for(j in 1:(i-1))
     {
       if(Nsym[i,j]==1 & M[i,j]==0)
       {
         AUE[i,j]=1
       }
     }
   }


  return(AUE)
}
