#' Adding an edge, same CC, if edge is not addable
#'
#'This is called from MCMCmove if x[1],x[2] are in the same
#'CC, there is no edge and the edge not addable.
#'All possibilities are considered of
#'x[1] -> CC2 <- CC1 where CC1 and CC2 are subsets of CC(w) where
#'x[2] is in CC2, CC2 connected, CC1 complete set and all CC1 vertices
#'connected to CC2 vertices, no directed cycles formed when compelled edges
#' are directed)
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
#'export


ADDEdgeSameCCNotAdd <- function(CG,x,v,w)
{
  #x[1] and x[2] are in the same CC, w, same layer v
  #but no edge x[1] - x[2]. First check if edge is addable.
  #If addable, then undirected edge is one of the possibilities
  #
  #Compute all possibilities for adding directed edge, x[1] to part of CC[[w]]
  #containing x[2]
  #and create new immoralities.
  #

  f <- function(j){length(CG$CC[[j]])}
  f1 <- sapply(1:length(CG$CC),f)
  d = sum(f1)
  #total number of vertices

  Out <- list()


  NewCG = CG

  lCC2 = length(CG$CC[[w]])

  x2 = match(x[2],CG$CC[[w]])
  x1 = match(x[1],CG$CC[[w]])

  #get x[1] and x[2] into CC coordinates

  Cmat = CG$CCued[[w]]

  ANI = CCWholeSplit(x1,x2,Cmat)
  count = ANI$count

  if(count == 0)
  {
    Out = list()
    Out$CG = CG
    Out$Ratio = 0
    return(Out)
  }


  whichpick = sample(count,1)

  CC1 = ANI$CC1[[whichpick]]
  CC2 = ANI$CC2[[whichpick]]
  CC3 = ANI$CC3[[whichpick]]

  CC11 = intersect(CC1,nhd(Cmat)$L[[x2]])

  PointOut <- CC11
  PointIn <- CC2



  #now CC[[w]] is split into new CCs and we run the routine to do so
  #with initialisation PointIn = CC2 and PointOut = x1 union CC1.
  #Then CC2 will be a CC
  #where i2ip1 are y and x[1], etc .....

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


  CANI <- function(z){CCWholeSplit(x1,z,Cmat)$count}
  N = sapply(CC2,CANI)
  sumthis <- sum(1/N)
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
