#' Metropolis  sampler: input is data matrix
#'
#' @param data is an n by d data matrix n instantiations of d-variate vector,
#'      which we centre and scale. We wrote our own (simple) centre/scale routine
#'      because the built-in routine seemed to give the wrong answer.
#' @param N is the number of MCMC iterations
#' @param theta is the assign vertices to CCs parameter for CRP
#' @param alpha is the layering CCs parameter for Chinese Restaurant Process
#' @param p is the edge probability parameter
#'
#' @return
#' OUT a list of length N, where OUT[[j]] has following components:
#' \itemize{
#' \item OUT[[j]]$OldCG current graph (before proposal)
#' \item OUT[[j]]$NewCG new graph
#' \item OUT[[j]]$logLikelihood log of the likelihood score (missing 2 pi)
#' \item OUT[[j]]$GraphScore log of probability of new graph
#' \item OUT[[j]]$WholeScore sum of GraphScore and logLikelihood
#' \item OUT[[j]]$move 1 if proposal accepted, 0 otherwise
#'}
#'@export

SamplerData <- function(data,N,theta,alpha,p)
{
#just in case we need a generalised inverse of the covariance matrix

  datamat = centrescale(data)

  #we insist on centring and scaling the data

  Covarmat = as.matrix(cov(datamat))
  d = ncol(datamat)
  n = nrow(datamat)

  CGinit = essgraph(theta,alpha,p,d)

  #we use a random initialisation


  MCOut <- list()

  MCOut[[1]] <- list()
  MCOut[[1]]$OldCG = CGinit
  MCOut[[1]]$NewCG = CGinit
  MCOut[[1]]$theta = theta
  MCOut[[1]]$alpha = alpha
  MCOut[[1]]$p = p



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
      Out$theta = theta
      Out$alpha = alpha
      Out$p = p
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
  }
  MCOut <- MCOut[2:(N+1)]



  return(MCOut)

}
