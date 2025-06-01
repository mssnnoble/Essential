#' Computes probability of extra parents
#'
#' Computes probabilities of 'extra edge' configuration (i.e. directed edges
#' from 2 or more layers back) for each CC.
#'
#' Routine is called from GraphPR
#'
#' @param CG - essential graph object
#' @param p - 'extra edge' probability
#'
#' @return L - vector of length equal to number of CCs
#' where L[k] gives the probability for configuration of DEex[[k]] (number of
#' parents for CC[[k]] that are 2 or more layers back)
#'@export

DEexPr <-function(CG,p){
  DEextra = CG$DEex
  prVIL = CG$prVIL
  ordlay = CG$ordlay
  CCs = CG$CC
  CCued = CG$CCued
  DEi2ip1 = CG$DEi2ip1
  L <- rep(1,length(CCs))
  #L[i] gives the probability of the set of 'extra edges' pointing into CCs[i]
  #for layers 1 and 2, there are no extra parents, probability is 1, therefore
  #no need to change it.
if(length(ordlay)>=3)
{
  for(j in 3:length(ordlay))
    {
    for(k in 1:length(ordlay[[j]]))
      {
      if(length(DEextra[[ordlay[[j]][k]]])>=1)
        {
        m = length(DEextra[[ordlay[[j]][k]]])
        n = length(prVIL[[j-2]])
        L[ordlay[[j]][k]] = (p^m)*((1-p)^(n-m))
        }
      else
        {
          #deal with the case of no extra edges. This can be obtained in two ways.
          #either no extra edges generated, or else a structure that violates Studeny
          #created
          if(checkcompletek(DEi2ip1[[ordlay[[j]][k]]],CCs,ordlay,CCued,j-1)[1]==0 )
            {
        n = length(prVIL[[j-2]])
        L[ordlay[[j]][k]] = (1-p)^n
        }
        else
          {
            e = checkcompletek(DEi2ip1[[ordlay[[j]][k]]],CCs,ordlay,CCued,j-1)
        n = length(prVIL[[j-2]])
        t = length(DEi2ip1[[ordlay[[e[2]]][e[3]]]]) + length(DEextra[[ordlay[[e[2]]][e[3]]]])
        L[ordlay[[j]][k]] = (1-p)^n + (p^t)*((1-p)^(n-t))
        }
        }
      }
  }
  return(L)
}
  if(length(ordlay)<=2)
  {
    return(L)
  }
}
