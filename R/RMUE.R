 #'Removes undirected edge
#'
#' @param CG - an essential graph object
#' @param x = (x[1],x[2]) two vertices in the same CC
#' @param v - the layer
#' @param w - the CC
#'
#' @return: Out
#' \itemize{
#' \item Out$CG is the new essential graph
#' \item Out$Ratio is the ratio between reverse and forward proposals
#' }
#'@export

RMUE <- function(CG,x,v,w)
  {
  #removing an undirected edge.
  #If the removal leaves two
  #connected components, then we
  #have two different CCs and we direct to RMUESplitCC
  #
  #if the edge is 'removable' in the sense
  #that we have a connected decomposable
  #CC, we remove the edge and choose among the
  #different possibilities for immoralities out of
  #the created vee structures.
  #
  #If removing the edge means that the
  #graph is no longer triangulated, then some
  #vee structures must be converted to immoralities
  #to return an essential graph.
  #
  #consider all maximal cliques containing
  #the edge, then choose at random either one of them or
  #none where vee structures remain undirected; for
  #all the others resolve by making the vee structures
  #into immoralities.



  #We'll put the output graph into NewCG, initially this is CG.

  NewCG <- CG

  #Firstly - check if removal of the undirected edge leads to splitting the CC
  #We only need to check if (x[1],x[2]) is itself one of the maximal
  #cliques (i.e. the MC is an edge).



  Cmat <-CG$CCued[[w]]
  y = match(x,CG$CC[[w]])

  H <- MC(Cmat)
  f <- function(j){1*(setequal(y,H[[j]]))}
  f1 = sapply(seq(1,length(H)),f)
  k = which(f1==1)
  if(length(k)==0){
    mark <-0
    }
  if(length(k) >= 1)
    {
    mark <- 1
  Cel <- k[1]
  }

  #so if mark == 0, edge removal in and of itself
  #does not split CC.
  #if mark == 1 then it does and we label the clique
  #from decomposition given by simp as Cel.



#now deal with the case where edge removal leads directly
  #to two CCs - compute the CCs and their undirected edges
  #then dump the current parent sets and choose new parent
  #sets for each of the two CCs, both i2ip1 and DEex

  if(mark == 1)
    {
    #i.e. the edge to be removed is an MC and hence CC gets split
     Out <- RMUESplitCC(CG,x,v,w)
    return(Out)
  }


    if(mark == 0)
    {
    Out = RMUESameCC(CG,x,v,w)
    return(Out)
    }
  }
