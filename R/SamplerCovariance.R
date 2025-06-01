#' Metropolis  Sampler: input is covariance
#'
#' @param covar is a d by d covariance matrix with diagonal elements all 1.
#' @param n is the number of observations (for which covariance is likelihood raised
#' to n/2)
#' @param p is the edge probability parameter
#' @param alpha is the layering CCs parameter for Chinese Restaurant Process
#' @param theta is the assign vertices to CCs parameter for CRP
#' @param N is the number of MCMC iterations
#'
#' @return
#' OUT contains a list of N, each with three components
#' OUT$CG is the essential graph (CG = chain graph)
#' OUT$Post is the posterior score for the graph
#' OUT$Prior is the prior graph probability
#' OUT$Lik is the likelihood score
#'
#'@export

SamplerCovariance <- function(covar,n,N,theta,alpha,p)
{


  Covarmat = covar
  d = ncol(Covarmat)

  CGinit = essgraph(theta,alpha,p,d)

  #we use a random initialisation


  MCOut <- list()
  MCOut[[1]] <- list()
  MCOut[[1]]$OldCG = CGinit
  MCOut[[1]]$NewCG = CGinit



  for(i in 2:(N+1))
  {
    Out = list()
    CGPRPN <- MCMCmove(MCOut[[i-1]]$NewCG)
    CGOLD = CGPRPN$OldCG
    CGProp = CGPRPN$NewCG
    LikOld = LikScore(CGOLD,n,Covarmat)
    LikNew = LikScore(CGProp,n,Covarmat)
    R3=LikNew - LikOld
    RatioProp = CGPRPN$Ratio
    ScoreNew = prod(unlist(GraphPR(CGProp,theta,alpha,p)))
    ScoreOld = prod(unlist(GraphPR(CGOLD,theta,alpha,p)))
    RatioGraphs = (ScoreNew)/(ScoreOld)
    A = min(1,RatioProp*RatioGraphs*(exp(R3)))
    if(is.nan(A)==TRUE)
    {A=1}
    print(A)
    accept = runif(1,0,1)

    Out$OldCG = CGOLD
    Out$RatioProp = RatioProp
    Out$RatioGraphs = RatioGraphs


    if(accept <= A)
    {
      Out$move = 1L
      Out$NewCG = CGProp
      Out$logLikelihood = LikNew
      Out$GraphScore = log(ScoreNew)
      Out$WholeScore = (LikNew) + (log(ScoreNew))
    }
    if(accept > A)
    {
      Out$move=0L
      Out$NewCG = CGOLD
      Out$logLikelihood = LikOld
      Out$GraphScore = log(ScoreOld)
      Out$WholeScore = (LikOld) + (log(ScoreOld))
    }
    MCOut[[i]] = Out
    print(Out$move)
    print("i=")
    print(i)
    save(MCOut, file="savedgraphs.Rda")
    }
  MCOut <- MCOut[2:(N+1)]


  return(MCOut)

}
