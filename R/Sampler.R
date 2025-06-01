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
#' OUT contains a list of N, each with three components
#' OUT$OldCG is the input CG from the previous move
#' OUT$NewCG is the important one - the new CG.
#'   This is the proposal if the proposal is accepted
#'   or the OldCG if the proposal is rejected
#' OUT$RatioProp is the ratio P(NewCG -> OldCG)/P(OldCG -> NewCG)
#'   numerator probability of proposing OldCG from new
#'   denominator probability of proposing NewCG from old
#'OUT$RatioGraphs is the ratio of the prior probabilities of new to old graphs
#'OUT$Score is the prior probability of NewCG
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
print(MCOut$CG)



   for(i in 2:N)
  {
     Out = list()
    CGPRPN <- MCMCmove(MCOut[[i-1]]$NewCG)
    CGOLD = CGPRPN$OldCG
    CGProp = CGPRPN$NewCG
    RatioProp = CGPRPN$Ratio
    ScorePrev = prod(unlist(GraphPR(CGOLD,theta,alpha,p)))
    ScoreNew = prod(unlist(GraphPR(CGProp,theta,alpha,p)))
    RatioGraphs = (prod(unlist(GraphPR(CGProp,theta,alpha,p))))/(prod(unlist(GraphPR(CGOLD,theta,alpha,p))))
    #R3 = CGprop$LikRat
    A = min(1,RatioProp*RatioGraphs)
    print(A)
    accept = runif(1,0,1)

    Out$OldCG = CGOLD
    Out$RatioProp = RatioProp
    Out$RatioGraphs = RatioGraphs
    Out$PropCG = CGProp
    Out$x = x

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


 return(MCOut)

}
