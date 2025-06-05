#'Computes probability of undirected edge configuration of the CCs
#'
#'@param CG an essential graph object
#'CG$CCs listing of the CCs
#'CG$CCued listing of lower triangular matrices indicating undirected edges
#'@param p extra edge probability parameter
#'
#'@return PC list of the edge configuration probabilities for the CCs
#'@export

prCC <-function(CG,p){
  #  CCs a listing of the CCs
  #CCued lower triangular matrices, one for each CC,
  #OUTPUT: given the CCs, this is the probability of the
  #undirected edge configuration

  CCs = CG$CC
  CCued = CG$CCued
  ordlay = CG$ordlay

  a <- as.integer(sapply(CCued,sum))
  #a gives number of edges in each CC
  b <- as.integer(sapply(CCs,length))
  #b gives number of vertices in each CC
  e <- a-b+1

  #e lists number of extra edges over and above spanning tree

  f <- rep(0,length(CCs))
  if(length(ordlay[[1]]) == 1 & length(ordlay)>=2){f[ordlay[[1]][1]] = 1}
  #f indicates if CC is required to be non-complete

  g <- function(j){dbinom(e[j], 0.5*(b[j]-1)*(b[j]-2)-f[j],p)}
  pr = sapply(1:length(CCs),g)

  PC <- rep(1,length(CCs))
  for(i in 1:length(CCs))
    {
    if(e[i]<=19 & b[i] <=8)
      {
          PC[i] = pr[i]*(1/decomp[e[i]+1,b[i]])
    }
    else{
      fnoug = (b[i])^(b[i]-2)
      diiv = (fnoug - 1)/((0.5*(b[i]-1)*(b[i]-2)))
      invp = fnoug - diiv*e[i]
      PC[i] = pr[i]/(invp)
    }
  }
  return(PC)
}
