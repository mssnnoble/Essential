#' Metropolis sampler for essential graphs from prior distribution
#'
#' This algorithm simply runs a Metropolis Hastings sampler through essential
#' graphs using 'prior' for 'posterior' and thus may be used for sampling
#' according to the given 'prior'.
#'
#' @param N is the number of MCMC iterations
#' @param d is the number of vertices
#' @param alpha is the layering CCs parameter for Chinese Restaurant Process
#' @param theta is the assign vertices to CCs parameter for CRP
#' @param p is the edge probability for the prior
#'
#' @return
#' OUT contains a list of N, each with the following components
#' \itemize{
#' \item OUT[[j]]$OldCG is the input CG from the previous move
#' \item OUT[[j]]$NewCG is the important one - the new CG.
#'   which is the proposal if the proposal is accepted
#'   or the OldCG if the proposal is rejected
#' \item OUT[[j]]$RatioProp is   P(NewCG -> OldCG)/P(OldCG -> NewCG)
#'   which is the ratio of probabilities of proposing old from new to new
#'   from old
#' \item OUT[[j]]$RatioGraphs is the ratio of the prior probabilities of new to old graphs
#' \item OUT[[j]]$Score is the logaritm of the (prior) probability of NewCG
#' }
#'
#'@export

Sampler <- function(N,d,theta,alpha,p)
{
CGinit = essgraph(theta,alpha,p,d)

  #we use a random initialisation


MCOut <- list()
MCOut[[1]] <- list()
MCOut[[1]]$CG = CGinit
MCOut[[1]]$NewCG = CGinit



   for(i in 2:(N+1))
  {
     Out = list()
    CGPRPN <- MCMCmove(MCOut[[i-1]]$NewCG)
    CGOLD = CGPRPN$OldCG
    CGProp = CGPRPN$NewCG
    RatioProp = CGPRPN$Ratio
    ScorePrev = log(prod(unlist(GraphPR(CGOLD,theta,alpha,p))))
    ScoreNew = log(prod(unlist(GraphPR(CGProp,theta,alpha,p))))
    RatioGraphs = (prod(unlist(GraphPR(CGProp,theta,alpha,p))))/(prod(unlist(GraphPR(CGOLD,theta,alpha,p))))

    A = min(1,RatioProp*RatioGraphs)
    if(is.nan(A)==TRUE)
    {A=1}
    print(A)
    accept = runif(1,0,1)

    Out$OldCG = CGOLD
    Out$RatioProp = RatioProp
    Out$RatioGraphs = RatioGraphs
    Out$PropCG = CGProp

    if(accept <= A)
   {
    Out$move = 1L
    Out$NewCG = CGProp
    Out$Score = ScoreNew
     }
  if(accept > A)
  {
    Out$move=0L
    Out$NewCG = CGOLD
    Out$Score = ScorePrev
  }
    MCOut[[i]] = Out
print(Out$move)
print("i=")
print(i)
   }
MCOut <- MCOut[2:(N+1)]

 return(MCOut)

}
