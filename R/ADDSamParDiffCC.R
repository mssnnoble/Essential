#' Add edge, x[1], x[2] from different CCs, but have same parents.
#'
#'deals with adding edge (directed or undirected) when x[1] and x[2] are from
#'different CCs with the same parent sets.
#'For an ordered pair (x[1],x[2]) consider all possibilities of x[1] - x[2]
#'or x[1] -> CC2 <- CC1
#'(if undirected edge then *only* to vertex x[2])
#' @param CG - an essential graph object
#' @param x = (x[1],x[2]) (two vertex labels)
#' @param v = (v[1],v[2]) (layers of x[1] resp x[2]) if layers
#' different; v is a positive integer if both are in same layer
#' @param w = (w[1],w[2]) (CC labels for CCs containing x[1] resp. x[2])
#'
#'
#' @return
#' OUT contains a list where
#' \itemize{
#' \item OUT$CG is the new essential graph (after x[1] -> CC2 <- CC1 added
#' CC1 and CC2 suitable subsets of CC(w[2]) or undirected edge x[1]-x[2]
#' has been added)
#'  \item OUT$Ratio is ratio of forward proposal probability divided by
#' reverse proposal probability
#'}
#'
#'@export


ADDSamParDiffCC <- function(CG,x,v,w)
{
#x[1] and x[2] are in different CCs, w[1] resp w[2]
  #they are in the same layer v
  #
  #sent to this routine if TOTPAR[w[1]] = TOTPAR[w[2]]
  #
  #EITHER add undirected edge x[1] - x[2]
  #and merge the CCs
  #
  #OR adds a directed edge, x[1] to part of CC[[w[2]]]
  #and create new immoralities; a neighbour y of x[2] is
  #chosen in w[2] and a connected component of N(y)\y containing x[2]
  #the only new immoralities are y > z < x[1] for each z in
  #the connected component, which becomes a CC. Other compelled
  #edges are directed.
  #
  #This enables a reverse move for the situation where a directed
  #edge is removed from the whole of a CC, thus removing
  #immoralities created in the forward move; legal mergers return the original.

f <- function(j){length(CG$CC[[j]])}
f1 <- sapply(1:length(CG$CC),f)
d=sum(f1)

NewCG = CG
Cmat = CG$CCued[[w[2]]]

x2 = match(x[2],CG$CC[[w[2]]])
if(length(CG$CC[[w[2]]])==1)
{
  count=0
}

if(length(CG$CC[[w[2]]]) >=2)
{
SPLITIT = CCSplit(x2,Cmat)
count = SPLITIT$count

CC1 = SPLITIT$CC1

f<- function(j){1*(length(CC1[[j]])>=1)}
f1 <- sapply(1:count,f)
g1 = which(f1==1)
CC1 = SPLITIT$CC1[g1]
CC2 = SPLITIT$CC2[g1]
CC3 = SPLITIT$CC3[g1]
count = sum(f1)
}
#recall: if parent sets are the same, we need CC1 non empty to get immoralities
#to ensure new directed edges don't vanish with legal merge

whichpick = sample((count + 1),1 )


if(whichpick == count+1)
{
  mark = "u"
   NewCG = ADDUDdiffCC(CG,x,v,w)
}
if(whichpick <= count)
{
  mark = "d"
  CC1 = CC1[[whichpick]]
  CC2 = CC2[[whichpick]]
  CC3 = CC3[[whichpick]]

  if(length(CC3)==0)
  {
    ln = length(CG$CC)
     NewCG$CC[[w[2]]] = CG$CC[[w[2]]][CC1]
    NewCG$CC[[ln+1]] = CG$CC[[w[2]]][CC2]
    NewCG$CCued[[w[2]]] = as.matrix(CG$CCued[[w[2]]][CC1,CC1])
    NewCG$CCued[[ln+1]] = as.matrix(CG$CCued[[w[2]]][CC2,CC2])
    NewCG$TOTPAR[[ln+1]] = sort(unique(c(CG$TOTPAR[[w[2]]],CG$CC[[w[2]]][CC1],x[1])))
    NewCG = LEGMERGE(NewCG)
  }

if(length(CC3)>=1)
{
  PointOut <- CC1
  PointIn <- CC2

  MC = MC(Cmat)
  NEWW = Compelled(Cmat,MC,PointOut,PointIn)
  NNewCC = NEWW$NNewCC
  Par = NEWW$Par

  lccb = length(NNewCC)

  NewCG$CC = NewCG$CC[-w[2]]
  NewCG$CCued = NewCG$CCued[-w[2]]
  NewCG$DEi2ip1 = NewCG$DEi2ip1[-w[2]]
  NewCG$DEex = NewCG$DEex[-w[2]]
  NewCG$TOTPAR = NewCG$TOTPAR[-w[2]]

  f <- function(j){1*(length(intersect(NNewCC[[j]],CC2))>=1)}
  f1 <- sapply(1:lccb,f)

  C2Cl = which(f1==1)
  C3Cl = which(f1==0)

  a = length(NewCG$CC)
  for(j in 1:lccb)
  {
    NewCG$CC[[a+j]] = vector(mode="integer",length=0L)
    NewCG$CCued[[a+j]] = as.matrix(matrix(0,1,1))
    NewCG$TOTPAR[[a+j]]=vector(mode="integer",length=0L)
  }
  for(j in C2Cl)
  {
    NewCG$CC[[a+j]] = CG$CC[[w[2]]][NNewCC[[j]]]
    NewCG$CCued[[a+j]] = as.matrix(Cmat[NNewCC[[j]],NNewCC[[j]]])
    NewCG$TOTPAR[[a+j]] = sort(c(x[1],CG$TOTPAR[[w[2]]],CG$CC[[w[2]]][Par[[j]]]))

  }
  for(j in C3Cl)
  {
    NewCG$CC[[a+j]] = CG$CC[[w[2]]][NNewCC[[j]]]
    NewCG$CCued[[a+j]] = as.matrix(Cmat[NNewCC[[j]],NNewCC[[j]]])
    NewCG$TOTPAR[[a+j]] = sort(c(CG$TOTPAR[[w[2]]],CG$CC[[w[2]]][Par[[j]]]))

  }
  #note: TOTPARs have to be correct; NEWLAY doesn't use DEex or DEi2ip1
  NewCG = NEWLAY(NewCG)
  NewCG = LEGMERGE(NewCG)
}
}



#now compute reverse and forward probabilities

if(mark == "u")
{
  x1 =  match(x[1],CG$CC[[w[1]]])
  Cmat1 = CG$CCued[[w[1]]]
  SPOP = CCSplit(x1,Cmat1)
  f <- function(j){1*(length(SPOP$CC1[[j]] >=1))}
  coo = SPOP$count
  f1 <- sapply(1:coo,f)
  cooo = 1+sum(f1)
  forwardpr = (1/(d*(d-1)))*((1/cooo)+(1/(1 + count)))


   #now for backward

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

}

if(mark=="d")
{
  #can obtain the same graph as a possibility for any x2 in CC2
  tot = 0
  for(z in CC2)
  {
    hh = CCSplit(z,Cmat)$CC1
    f <- function(a){1*(length(hh[[a]])>=1)}
    lspln = length(hh)
    f1 = sapply(1:lspln,f)
    suuu = 1 + sum(f1)
    tot = tot + (1/suuu)
  }
forwardpr = (1/(d*(d-1)))*(1/tot)

#now for backward probability. Remove directed edge x[1] -> CC
#where CC is the CC in NewCG containing x[2]
#each vertex in this CC must be considered

f <- function(j){1*(x[2] %in% NewCG$CC[[j]])}
lll = length(NewCG$CC)
f1 <- sapply(1:lll,f)
g1 = which(f1 == 1)[1]
N <- vector(mode="integer",length = length(NewCG$CC[[g1]]))
cc <- function(z){CCSplit(z,NewCG$CCued[[g1]])$count}
N = sapply(1:length(NewCG$CC[[g1]]),cc)

reversepr = (1/(d*(d-1)))*(sum(1/N))
  }



ratio = reversepr/forwardpr

Out <- list()
Out$CG <- NewCG
Out$Ratio <- ratio
return(Out)

}


