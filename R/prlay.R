#'Computes the probability of the layering, given the CCs
#'
#'@param CG an essential graph object
#' @param alpha - Chinese Restaurant concentration parameter for layers
 #'
#' @return P1 probability of the layering
#'
#'@export


prlay <- function(CG,alpha){
  ordlay = CG$ordlay
  CCs = CG$CC

  if(is.list(ordlay)==FALSE)
  {
    a = list()
    a[[1]] = ordlay
    ordlay = a
  }

  if(is.list(CCs)==FALSE)
  {
    a = list()
  a[[1]] = CCs
  CCs = a
  }



  sl <- sapply(ordlay,length)
CClen = sapply(CCs,length)
  #deal with cases. First case is whether or not it could have arisen from
  #merging two layers, the other case, it couldn't.
  #
  #It can only arise from merging of two layers if
  #all CCs have either 1 or 2 vertices
  #ground layer contains 2 CCs, all other layers contain one CC

A = 1*(sl[1]==2 & length(ordlay) == length(CCs)-1 & max(CClen) <=2)

  if(A==1)
  {
    #these three conditions define and characterise the situation

    s2 = length(ordlay)
    r = rep(1,s2+1)
    #first term in P1 - assignment to layers (2,1,1....). Then choice for ground
    #is determined, others assigned at random
    #second term in P1 - assignment to layers (1,1,1,....) each layer singleton.
    #probability of subset of 2 chosen randomly for merger from s2+1
    #s2 - 1 then arranged randomly
    A1 = ((1/(factorial(s2 -1)))*(prs(sl,alpha)))
    B1 = ((prs(r,alpha))*(2/(s2*(s2+1)))*(1/(factorial(s2-1))))
    P1 = A1+B1
    return(P1)
  }
if(A ==0)
  {
    k = 0
  for(i in 1:length(ordlay))
    {
    if(length(ordlay[[i]]) >=2)
    {
      k <- k+1
      }
    if(length(ordlay[[i]]) <=1)
      {
      if(length(CCs[[ordlay[[i]][1]]] >=3))
      {
        k<-k+1
      }
    }
    }
  P1 <- (prs(sl,alpha))*(1/k)*(1/(factorial(length(ordlay)-1)))
  return(P1)
  }
}
