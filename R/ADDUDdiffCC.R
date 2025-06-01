#' Add undirected edge x[1],x[2] different CCs, same parents
#'
#'this is called from ADDSamParDiffCC when the
#'decision is made to add an undirected edge <x[1],x[2]>
#'x[1] and x[2] in different CCs which have the same parents.
#'
#' @param CG - an essential graph object
#' @param x = (x[1],x[2]) (two vertex labels)
#' @param v   is a positive integer giving the layer (same layer
#' if same parents)
#' @param w = (w[1],w[2]) (CC labels for CCs containing x[1] resp. x[2])
#'
#'
#' @return
#' OUT is the new chain graph after the edge has been added (and legal mergers
#' carried out)
#'
#'@export




ADDUDdiffCC <- function(CG,x,v,w)
{
  #this routine is called from ADDSamParDiffCC if the add undirected edge
  #option is chosen - we only compute the new CG.
  #adds undirected edge, hence v[1]=v[2] (=v), but different CCs (same layer).
  #This move requires TOTPAR[[w[1]]] = TOTPAR[[w[2]]] (equal parent sets).
  #We merge the two CCs
  #We then have to do possible legal mergers, just in case the additional
  #undirected edge creates a complete set so that the Studeny conditions are violated.

  NewCG <- CG


  #firstly, deal with adding extra edge. When we merge
  #the two CCs, the new CC will take the label min(w[1],w[2])

  o = sort(w)



  NewCG$CC[[o[1]]] <-sort(c(CG$CC[[o[1]]],CG$CC[[o[2]]]))

  #this merges the CCs
  #now deal with the undirected edges.

  l <- length(NewCG$CC[[o[1]]])
  NewCG$CCued[[o[1]]] <- matrix(0,l,l)
  NewCG$CCued[[o[1]]][match(CG$CC[[o[1]]],NewCG$CC[[o[1]]]),match(CG$CC[[o[1]]],NewCG$CC[[o[1]]])] <- CG$CCued[[o[1]]]
  NewCG$CCued[[o[1]]][match(CG$CC[[o[2]]],NewCG$CC[[o[1]]]),match(CG$CC[[o[2]]],NewCG$CC[[o[1]]])] <-
    CG$CCued[[o[2]]]

  #so far we have put in the existing edges; now the new undirected edge between x = (x[1],x[2])

  y <- match(x,NewCG$CC[[o[1]]])
  z<-sort(y,decreasing=TRUE)

  NewCG$CCued[[o[1]]][z[1],z[2]] <- 1

  #the new CC and its undirected edges have now been sorted out.
  #choose undirected edges for the new CC - valid subset from the union of the current subsets

  NewCG$DEi2ip1[[o[1]]] = CG$DEi2ip1[[o[1]]]
  NewCG$DEex[[o[1]]] = CG$DEi2ip1[[o[1]]]
  NewCG$TOTPAR[[o[1]]] = CG$TOTPAR[[o[1]]]


  #now remove the redundant CC, labelled o[2]

  a = seq(1,length(NewCG$CC))
  aprime = a[-o[2]]
  NewCG$CC <- NewCG$CC[-o[2]]
  NewCG$CCued <- NewCG$CCued[-o[2]]
  NewCG$DEi2ip1 <- NewCG$DEi2ip1[-o[2]]
  NewCG$DEex <- NewCG$DEex[-o[2]]
  NewCG$TOTPAR <- NewCG$TOTPAR[-o[2]]
  NewCG$ordlay[[v]] <- NewCG$ordlay[[v]][-match(o[2],NewCG$ordlay[[v]])]

  NewCG <- NEWLAY(NewCG)
  NewCG <-LEGMERGE(NewCG)


  return(NewCG)
}
