#' Indicates if directed edge is possible without creating partially directed cycles
#'
#' When x = (x[1],x[2]) is chosen, for x[1] in CC(w1) in layer v1 and
#' x[2] in CC(w2) in layer v2, this algorithm tells us if a directed
#'  edge from x[1] to x[2] is possible without creating a directed cycle
#'  in the CC DAG.
#'
#' @param CG essential graph object
#' @param v1, v2 (layers of x[1] resp x[2])
#' @param w1, w2 (CCs containing x[1] resp x[2])
#'
#' @return EDGEYES which is 1 if directed edge x[1] to x[2] is possible
#' 0 otherwise
#'@export


edgeposs <- function(CG,v1,v2,w1,w2)
{
  #v1 and v2 respective layers of x[1] and x[2]
  #w1 and w2 respective CCs of x[1] and x[2]
  if(v1<=v2)
  {
    EDGEYES = 1
  }
if(v1 >= v2+1)
{
  #check if move is possible. In other words, adding
  #the edge does not create a cycle in the DAG of CCs

  CClistt = vector(mode="integer", length = 0L)
  #list of CCs which are ancestors of w1

  vertlist = vector(mode="integer",length = 0L)
  #list of ancestor vertices of w1

  CClistt <- sort(c(CClistt,w1))
  CClistcurrent <- CClistt
  vertlist <- sort(c(vertlist,CG$TOTPAR[[w1]]))

  sett = 1:length(CG$CC)

  repeat
  {
    f <- function(j){1*(length(intersect(vertlist,CG$CC[[j]])) >=1)}
    f1 = sapply(sett,f)
    #locate CCs which contain ancestor vertices

    f2 = which(f1==1)

    CClisttold <- CClistt
    CClistt <- sort(unique(c(CClisttold,f2)))
    CClistcurrent = setdiff(CClistt,CClisttold)
    if(length(CClistcurrent)==0)
    {
      break
    }
    newvertlist = unlist(CG$TOTPAR[f2])
    vertlist <- sort(unique(c(vertlist,newvertlist)))
  }

  EDGEYES =1 - 1*(w2 %in% CClistt)
}
return(EDGEYES)
}

