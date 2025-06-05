#' Add directed edge, different CCs and different parent sets
#'
#'deals with adding directed edge when x[1] and x[2] are from CCs with
#'different parent sets.
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
#' \item OUT$CG is the new essential graph (after x[1] -> S added
#' where S is a chosen subset of CC(w[2]))
#' \item OUT$ratio is ratio of forward proposal probability divided by
#' reverse proposal probability
#'}
#'
#'@export


ADEDiffPar <- function(CG,x,v,w)
  {
  #x[1] in CC[w[1]], x[2] in CC[w[2]]. Sent here if
  #1) w[1] != w[2] and
  #2) TOTPAR[[w[1]]] != TOTPAR[[w[2]]]
  #Firstly, check whether x[1] -> w[2] and x[2] -> w[1]
  #are both possible. If not, re-order x so that
  #x[1] -> CC(w[2]) is possible (and the other isn't)
  #If both are possible, we do not re-order (and this is taken into
  #account when computing proposal probability)
  #add a directed edge, x[1] -> S2
  #We have all splits of CC(w[2]) into disjoint sets S1,S2,S3
  #where S1 is complete, S2 connected, all connections between
  #S1 and S2 vertices exist.
  #create immoralities x[2] -> S2 <- S1 which are removed when
  #reverse move is applied.


f<-function(j){length(CG$CC[[j]])}
f1=sapply(1:length(CG$CC),f)
d=sum(f1)

#number of vertices

NewCG = CG

Out = list()


#edgeposs = 1 means that x[1] -> CC[w[2]] can be added without
#creating a partially cycle.
#edgepossrev = 1 means that x[2] -> CC[w[1]] can be added without
#creating a cycle.
#NewCG will be output CG

if(length(v)==1)
{
  #we had v[1] == v[2] so this was resolved as v = v[1]
  edgepossib = 1
  edgepossibrev = 1
}

if(length(v)==2)
{
if(v[1] >= v[2]+1)
{
  edgepossibrev = 1
  edgepossib = edgeposs(CG,v[1],v[2],w[1],w[2])
  }

if(v[2] >= v[1]+1)
{
  edgepossib = 1
  edgepossibrev = edgeposs(CG,v[2],v[1],w[2],w[1])
}

if(v[1] >= v[2]+1)
{
if(edgepossib == 0)
{
  v = c(v[2],v[1])
  w = c(w[2],w[1])
  x = c(x[2],x[1])
}
}
}


Cmat = CG$CCued[[w[2]]]
x2 = match(x[2],CG$CC[[w[2]]])

SPLITIT = CCSplit(x2,Cmat)
count = SPLITIT$count

whichpick = sample(count,1)

CC1 = SPLITIT$CC1[[whichpick]]
CC2 = SPLITIT$CC2[[whichpick]]
CC3 = SPLITIT$CC3[[whichpick]]

if(length(CC1)+length(CC3)==0)
{
  #in this case add DE to whole CC
  NewCG$DEi2ip1[[w[2]]]= sort(c(NewCG$DEi2ip1[[w[2]]],x[1]))
  NewCG$TOTPAR[[w[2]]]= sort(c(NewCG$TOTPAR[[w[2]]],x[1]))
  NewCG = LEGMERGE(NewCG)
  Out$CG = NewCG
}

if(length(CC1)>=1 & length(CC3)==0)
{
  ln = length(CG$CC)
  NewCG$CC[[w[2]]] = CG$CC[[w[2]]][CC1]
  NewCG$CC[[ln+1]] = CG$CC[[w[2]]][CC2]
  NewCG$CCued[[w[2]]] = as.matrix(CG$CCued[[w[2]]][CC1,CC1])
  NewCG$CCued[[ln+1]] = as.matrix(CG$CCued[[w[2]]][CC2,CC2])
  NewCG$TOTPAR[[ln+1]] = sort(unique(c(CG$TOTPAR[[w[2]]],CG$CC[[w[2]]][CC1],x[1])))
  NewCG = LEGMERGE(NewCG)
  Out$CG = NewCG
}

if(length(CC1)>=1 & length(CC3) >= 1)
{
PointOut <- CC1
PointIn <- CC2

MC = MC(Cmat)
NEWW = Compelled(Cmat,MC,PointOut,PointIn)
NNewCC = NEWW$NNewCC
Par = NEWW$Par

NewCG$CC = NewCG$CC[-w[2]]
NewCG$CCued = NewCG$CCued[-w[2]]
NewCG$DEi2ip1 = NewCG$DEi2ip1[-w[2]]
NewCG$DEex = NewCG$DEex[-w[2]]
NewCG$TOTPAR = NewCG$TOTPAR[-w[2]]

lncc = length(NNewCC)
f <- function(j){1*(length(intersect(NNewCC[[j]],CC2))>=1)}
f1 <- sapply(1:lncc,f)

C2Cl = which(f1==1)
C3Cl = which(f1==0)

a = length(NewCG$CC)
for(j in 1:lncc)
{
  NewCG$CC[[a+j]] = vector(mode="integer",length=0L)
  NewCG$CCued[[a+j]]= as.matrix(matrix(0,1,1))
  NewCG$TOTPAR[[a+j]] = vector(mode="integer",length=0L)
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
}

if(length(CC1)==0 & length(CC3) >= 1)
{
  MCl = MC(Cmat)
  POINTIN = CC2
  POINTOUT = vector(mode="integer",length=0L)
  Dir = Compelled(Cmat,MCl,POINTOUT,POINTIN)
  ln = (length(Dir$NNewCC))
  f <- function(j){1*(setequal(CC2,Dir$NNewCC[[j]]))}
  f1 <- sapply(1:ln,f)
  g1 = which(f1==1)
  NewCG$CC[[w[2]]] = CG$CC[[w[2]]][CC2]
  NewCG$CCued[[w[2]]] = as.matrix(CG$CCued[[w[2]]][CC2,CC2])
  NewCG$TOTPAR[[w[2]]] = sort(c(x[1],CG$TOTPAR[[w[2]]]))
  a = length(CG$CC)
  NNewCC = Dir$NNewCC[-g1]
  Par = Dir$Par[-g1]
  for(j in 1:(ln-1))
  {
    NewCG$CC[[j+a]] = CG$CC[[w[2]]][NNewCC[[j]]]
    NewCG$CCued[[j+a]] = as.matrix(Cmat[NNewCC[[j]],NNewCC[[j]],drop=FALSE])
    NewCG$TOTPAR[[j+a]] = sort(c(CG$TOTPAR[[w[2]]],CG$CC[[w[2]]][Par[[j]]]))
    NewCG$DEi2ip1[[j+a]] = CG$CC[[w[2]]]
  }
  NewCG = NEWLAY(NewCG)
}



#this is the new graph. We also need the proposal and the reverse-proposal
#probabilities.

N = vector(mode="integer",length=length(CC2))
f <- function(z){CCSplit(z,CG$CCued[[w[2]]])$count}
N = sapply(CC2,f)

prprel = (1/(d*(d-1)))*(sum(1/N))
if(edgepossib + edgepossibrev == 2)
{
  forwardpr = prprel
}
if(edgepossib + edgepossibrev == 1)
{
  forwardpr = 2*prprel
}



#Now we consider reverse move. Find k such that NewCG$CC[[k]] = CC2
#probability of the move is the probability of (a) selecting (x[1],z)
#(or (z,x[1])) for a z in NewCC[[k]] and (b) conditioned on this selection
#electing to remove edge from the whole CC.
#recall that we consider partial deletions where the CC can be split in two
#(which get merged when the edge is added)

if(1*(is.list(NewCG$CC)==TRUE))
{
f <- function(j){1*(x[2] %in% NewCG$CC[[j]])}
ln = length(NewCG$CC)
f1 = sapply(1:ln,f)
g1 = which(f1==TRUE)[1]
}
if(1*(is.list(NewCG$CC)==FALSE))
{
  g1==1
}
CNewMat = NewCG$CCued[[g1]]
NC2 = length(NewCG$CC[[g1]])
g <- function(j){CCSplit(j,CNewMat)$count}
g1 = sapply(1:NC2,g)

reversepr = (2/(d*(d-1)))*sum(1/g1)



  ratio = reversepr/forwardpr


Out$CG = NewCG
Out$Ratio = ratio
return(Out)
}
