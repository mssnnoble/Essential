#' Adding edge in a CC if edge is addable
#'
#'This is called from MCMCmove if x[1],x[2] are in the same
#'CC, there is no edge and the edge is addable.
#'All possibilities are considered of an undirected edge, also all
#'possibilities for
#'x[1] -> CC2 <- CC1 are also considered where CC1 - CC2 edges exist and
#'are directed to form immoralities which vanish by legal merger when
#'x[1] -> CC2 is removed
#'
#' @param CG an essential graph
#' @param x = (x[1],x[2]) two vertices
#' @param v the layer in which x[1] and x[2] appear
#' @param w the label of the CC in which x[1] and x[2] appear
#'
#'
#' @return
#' Out where
#'Out$CG gives new CG (after addition of edge)
#'Out$Ratio gives ratio of forward to reverse proposals
#'
#'@export


ADDEdgeSameCCAddable <- function(CG,x,v,w)
{
  # CG - chain graph object
  #x = x[1],x[2] - between x[1] and x[2] to be added
  #v[1]=v[2]=v (same layer)
  #w[1] = w[2] = w (same CC) and (x[1],x[2]) undirected edge is addable


  f <- function(j){length(CG$CC[[j]])}
  f1 <- sapply(1:length(CG$CC),f)
  d = sum(f1)

  Out = list()
  NewCG = CG

  lCC2 = length(CG$CC[[w]])

  x2 = match(x[2],CG$CC[[w]])
  x1 = match(x[1],CG$CC[[w]])

  #this gets x[1] and x[2] into CC coordinates

   Cmat <- CG$CCued[[w]]
  CCold <- CG$CC[[w]]

  nh = nhd(Cmat)$L

  y = match(x,CG$CC[[w]])
  z = sort(y,decreasing = TRUE)

#This algorithm is called if the undirected edge is
#addable.

#Firstly, find the directed possibilities:

DEPOSSIB = CCWholeSplitADD(x1,x2,Cmat)

count = DEPOSSIB$count
tott = count + 1
whichpick = sample(tott,1)

if(whichpick == tott)
{
 # in this case we simply add the undirected edge

NewCG$CCued[[w]][z[1],z[2]] <- 1

  NewCG <- LEGMERGE(NewCG)
  Out$CG = NewCG


  #now we need to compute the forward and backward proposal probabilities.

   forwardpr = (2/(d*(d-1)))*(1/tott)



  #now compute reverse probability when undirected edge is added; we suppose
  #NewCG has been
  #created by adding an undirected edge. We need the number of ways of
  #getting new graphs when creating immoralities (that are zapped by
  #addition of undirected x[1] - x[2])

  f <- function(j){1*(all(x %in% NewCG$CC[[j]]))}
  f1 = sapply(1:length(NewCG$CC),f)
  k = which(f1==1)
  MCl = MC(NewCG$CCued[[k]])
  y = match(x,NewCG$CC[[k]])
  g <- function(j){1*(all(y %in% MCl[[j]]))}
  g1=sapply(1:length(MCl),g)
  mcl = which(g1==1)

  totsets = list()
  totsets[[1]] = vector(mode="integer",length=0L)
  for(j in mcl)
  {
    PossNonIm <- setdiff(MCl[[j]],y)
    aaa = length(PossNonIm)
    all_subsets <- lapply(1:aaa,function(k) combn(aaa,k,simplify=FALSE))
    all_subsets <- unlist(all_subsets,recursive=FALSE)
    all_subsets = as.list(all_subsets)
    f <- function(k){PossNonIm[all_subsets[[k]]]}
    new_subsets = lapply(1:length(all_subsets),f)
    totsets = unique(append(totsets,new_subsets))
  }
count = length(totsets)
reversepr = (2/d*(d-1))*(1/count)

ratio = reversepr/forwardpr

  Out$CG <- NewCG
  Out$Ratio <- ratio

  return(Out)
}
if(whichpick <= (tott-1))
{
  PointOut = DEPOSSIB$CC1[[whichpick]]
  PointIn = DEPOSSIB$CC2[[whichpick]]
  CC2 = PointIn

  MaCl = MC(Cmat)
  NewCCPar <- Compelled(Cmat,MaCl,PointOut,PointIn)
  NNewCC = NewCCPar$NNewCC
  Par = NewCCPar$Par

  #We now get the total parent sets of the new CCs. PA denotes the parents of the new CCs.
  #Each vertex in a CC has the same parents

  nnewCC = length(NNewCC)
  sss <- 1:nnewCC



  PA <- vector(mode="list", length=nnewCC)
  PA <-lapply(PA,as.integer)

  f <- function(j){sort(c(CG$CC[[w]][Par[[j]]],CG$TOTPAR[[w]]))}
  PA <- lapply(seq_along(sss),f)

  old = length(CG$CC)

  NewCG$CC[[w]] <- CG$CC[[w]][NNewCC[[1]]]
  NewCG$CCued[[w]] <- as.matrix(Cmat[NNewCC[[1]],NNewCC[[1]],drop = FALSE])
  NewCG$TOTPAR[[w]] <- PA[[1]]

  for(j in 1:(nnewCC - 1))
  {
    NewCG$CC[[old+j]] =  CG$CC[[w]][NNewCC[[j+1]]]
  }

  for(j in 1:(nnewCC-1))
  {
    NewCG$CCued[[old+j]] <- as.matrix(Cmat[NNewCC[[j+1]],NNewCC[[j+1]],drop=FALSE])
  }

  for(j in 1:(nnewCC-1))
  {
    NewCG$TOTPAR[[old+j]] <-PA[[j+1]]
  }


  #note here that the TOTPAR are correct - we get the correct i2ip1 and DEex
  #after the new layering

  NewCG = NEWLAY(NewCG)
  Out$CG = NewCG

  #Now we have to compute the ratio of reverse to forward proposal.


  CANI <- function(z){CCWholeSplitADD(x1,z,Cmat)$count}
  N = sapply(CC2,CANI)
  sumthis <- sum(1/(N+1))
  forwardpr =(1/(d*(d-1)))*sumthis

  #For reverse, we remove the directed edge from whole of S

  Spr = CG$CC[[w]][CC2]
  f <- function(j){1*(all(Spr %in% NewCG$CC[[j]]))}
  f1=sapply(1:length(NewCG$CC),f)
  k=which(f1==1)

  frev <- function(z){CCSplit(z,NewCG$CCued[[k]])$count}
  Nrev = sapply(1:length(NewCG$CC[[k]]),frev)


  reversepr = (2/(d*(d-1)))*sum(1/Nrev)
  Ratio = reversepr/forwardpr
  Out$Ratio= Ratio
  return(Out)
}
}
