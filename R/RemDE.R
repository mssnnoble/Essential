#'Removes directed edge
#'
#'deals with removing a directed edge from x[1] to S, subset CC(w[2]).
#'We find all splits of CC(w[2]) into disjoint subsets CC2, CC1, CC3 such
#'that their union is CC(w[2]), x[2] in CC2, CC2 connected,
#'CC1 complete, all connections between CC1 and CC2 present. Then we direct
#'CC1 -> CC2
#'
#'Adding edge x[1] to CC2 reverses the move.
#'
#' @param CG - an essential graph object
#' @param x = (x[1],x[2]) (two vertex labels)
#' @param v = (v[1],v[2]) (layers of x[1] resp x[2]) if layers
#' different; v is a positive integer if both are in same layer
#' @param w = (w[1],w[2]) (CC labels for CCs containing x[1] resp. x[2])
#'
#'
#' @return
#' OUT contains a list where
#' OUT$CG is the new essential graph (after x[1] -> S removed
#' S suitable subset of CC(w[2]))
#' OUT$ratio is ratio of forward proposal probability divided by
#' reverse proposal probability
#'
#'@export



RemDE <- function(CG,x,v,w)
  {

  Out <- list()

  f <- function(j){length(CG$CC[[j]])}
  f1 <- sapply(1:length(CG$CC),f)
  d=sum(f1)


  NewCG <- CG





  if(v[1] > v[2])
  {
    #we re-order so that we add directed edge *from* x[1] *to*
    #some part of w[2]

      x <-c(x[2],x[1])
    v <- c(v[2],v[1])
    w <- c(w[2],w[1])
  }

  y=match(x[2],CG$CC[[w[2]]])

  Cmat = CG$CCued[[w[2]]]
  ALLCC = 1:nrow(Cmat)


      #first task is to find splits

    if(length(CG$CC[[w[2]]])==1)
  {
    CC1 = list()
    CC1[[1]] = vector(mode = "integer", length = 0L)
    CC2 = list()
    CC2[[1]] = ALLCC
    CC3 = list()
    CC3[[1]]= vector(mode = "integer", length = 0L)
    count = 1
    POSSIB = count
  }

  if(length(CG$CC[[w[2]]]) >=2)
{
  SPLITIT = CCSplit(y,Cmat)
  CC1 = SPLITIT$CC1
CC2 = SPLITIT$CC2
CC3 = SPLITIT$CC3
count = SPLITIT$count
}

#now select one of the splits at random

a = sample(count,1)
POSSIB = count
CC1 = CC1[[a]]
CC2 = CC2[[a]]
CC3 = CC3[[a]]

S = CG$CC[[w[2]]][CC2]

      if(length(CC1)+length(CC3)==0)
      {
        #this means remove whole x[1] to CC edge.

         NewCG$CC[[w[2]]] = CG$CC[[w[2]]]
        NewCG$CCued[[w[2]]] = CG$CCued[[w[2]]]
         NewCG$TOTPAR[[w[2]]] <-as.integer(NewCG$TOTPAR[[w[2]]][-match(x[1],NewCG$TOTPAR[[w[2]]])])
        NewCG <- NEWLAY(NewCG)
        NewCG <- LEGMERGE(NewCG)
      }

      if(length(CC1)>=1 & length(CC3)==0)
      {
        KK2 = CC2
        KK1 = CC1


       NewCG$CC[[w[2]]] = CG$CC[[w[2]]][KK1]

       ln = length(CG$CC)

       NewCG$CC[[ln+1]] = CG$CC[[w[2]]][KK2]
       NewCG$CCued[[w[2]]] = as.matrix(CG$CCued[[w[2]]][KK1,KK1])
       NewCG$CCued[[ln+1]] = as.matrix(CG$CCued[[w[2]]][KK2,KK2])
       NewCG$DEi2ip1[[ln+1]] = CG$CC[[w[2]]][KK1]
       oldpar = CG$TOTPAR[[w[2]]][-match(x[1],CG$TOTPAR[[w[2]]])]
       NewCG$TOTPAR[[ln+1]] = sort(c(CG$CC[[w[2]]][KK1],oldpar))
       NewCG$DEex[[ln+1]] = setdiff(NewCG$TOTPAR[[ln+1]],NewCG$DEi2ip1[[ln+1]])
       NewCG = NEWLAY(NewCG)
       NewCG = LEGMERGE(NewCG)
      }

if(length(CC1)==0 & length(CC3)>=1)
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
NewCG$TOTPAR[[w[2]]] = CG$TOTPAR[[w[2]]][-match(x[1],CG$TOTPAR[[w[2]]])]
a = length(CG$CC)
NNewCC = Dir$NNewCC[-g1]
Par = Dir$Par[-g1]

for(j in 1:(ln-1))
{
  NewCG$CC[[a+j]] = vector(mode="integer", length = 0L)
  NewCG$CCued[[a+j]] = as.matrix(matrix(0,1,1))
  NewCG$TOTPAR[[a+j]] = vector(mode="integer",length=0L)
}

for(j in 1:(ln-1))
{
 NewCG$CC[[j+a]] = CG$CC[[w[2]]][NNewCC[[j]]]
  NewCG$CCued[[j+a]] = as.matrix(Cmat[NNewCC[[j]],NNewCC[[j]],drop=FALSE])
  NewCG$TOTPAR[[j+a]] = sort(c(CG$TOTPAR[[w[2]]],CG$CC[[w[2]]][Par[[j]]]))
 }
NewCG = NEWLAY(NewCG)
NewCG = LEGMERGE(NewCG)



}

if(length(CC1)>=1 & length(CC3)>=1)
{
  MCl = MC(Cmat)

  POINTIN = CC2
  POINTOUT = CC1
  Dir = Compelled(Cmat,MCl,POINTOUT,POINTIN)
  NNewCC = Dir$NNewCC
  Par = Dir$Par

  #now decide which of the new CCs are in CC2 and which are in CC3.
  #By construction, all in one or the other

  NewCG$CC = NewCG$CC[-w[2]]
  NewCG$CCued = NewCG$CCued[-w[2]]
  NewCG$DEi2ip1 = NewCG$DEi2ip1[-w[2]]
  NewCG$DEex = NewCG$DEex[-w[2]]
  NewCG$TOTPAR = NewCG$TOTPAR[-w[2]]

  a1 = length(NewCG$CC)

  lccb = length(NNewCC)

  for(j in 1:lccb)
  {
    NewCG$CC[[a1+j]] = vector(mode="integer",length=0L)
    NewCG$CCued[[a1+j]] = as.matrix(matrix(0,1,1))
    NewCG$TOTPAR[[a1+j]] = vector(mode = "integer", length = 0L)
  }

  f <- function(j){1*(length(intersect(NNewCC[[j]],CC2))>=1)}
  f1 <- sapply(1:lccb,f)

  C2Cl = which(f1==1)
  C3Cl = which(f1==0)


    for(j in C2Cl)
  {
    NewCG$CC[[a1+j]] = CG$CC[[w[2]]][NNewCC[[j]]]
    NewCG$CCued[[a1+j]] = as.matrix(Cmat[NNewCC[[j]],NNewCC[[j]]])
    NewCG$TOTPAR[[a1+j]] = sort(c(CG$TOTPAR[[w[2]]][-match(x[1],CG$TOTPAR[[w[2]]])],CG$CC[[w[2]]][Par[[j]]]))
    }
  for(j in C3Cl)
  {
    NewCG$CC[[a1+j]] = CG$CC[[w[2]]][NNewCC[[j]]]
    NewCG$CCued[[a1+j]] = as.matrix(Cmat[NNewCC[[j]],NNewCC[[j]]])
    NewCG$TOTPAR[[a1+j]] = sort(c(CG$TOTPAR[[w[2]]],CG$CC[[w[2]]][Par[[j]]]))
       }
  #note: TOTPARs have to be correct; NEWLAY doesn't use DEex or DEi2ip1
  NewCG = NEWLAY(NewCG)
  NewCG = LEGMERGE(NewCG)

}


Out$CG = NewCG

#now compute forward and reverse proposal probabilities.
#forward
#assume rearranged so that existing arrow is from x_1 to
#CC[w_2]. Move is remove part (or all) of this arrow.
#
#for forward move, we can obtain this by choosing
#any (x1,x2) pair where x2 in S

N = vector(mode="integer",length=length(S))
g <- function(z){CCSplit(z,Cmat)$count}
N <- sapply(CC2,g)
Ninv = 1/N

forwardpr = (2/(d*(d-1)))*sum(Ninv)

#Now consider the reverse probability, that of proposing CG from NewCG.
#Let us locate the NewCG$CC's which contain x[2] and x[1]

#general idea: we get it by choosing (x[1],z) for any z in the
#CC containing x[2].
#if both ways are possible, we don't multiply by 2, if
#only one way is possible, we multiply by 2.
#if both ways are possible, we have to additionally
#check whether or not parents are different (i.e
#include possibility of undirected edge chosen)
#
#also, after directed edge removal and legal merge, they may both be
#in same CC.

f <- function(j){1*(x[2] %in% NewCG$CC[[j]])}
f1 <- sapply(1:length(NewCG$CC),f)
k2 = which(f1==1)
S = NewCG$CC[[k2]]
Cmat = NewCG$CCued[[k2]]

g <- function(j){1*(x[1] %in% NewCG$CC[[j]])}
g1 <- sapply(1:length(NewCG$CC),g)
k1 = which(g1==1)

if(k1==k2)
{
  z1 = match(x[1],NewCG$CC[[k1]])
  #x[1],x[2] both in NewCG$CC[[k1]]
  set = match(CG$CC[[w[2]]][CC2],NewCG$CC[[k2]])
  f <- function(z){CCWholeSplit(z1,z,Cmat)$count}
  f1 <- sapply(set,f)
  reversepr = (1/(d*(d-1)))*(sum(1/f1))

}
if(k1 != k2)
{
 la2 <- function(j){(1*(k2 %in% NewCG$ordlay[[j]]))}
la2li <- sapply(1:length(NewCG$ordlay),la2)
v2 = which(la2li == 1)

la1 <- function(j){(1*(k1 %in% NewCG$ordlay[[j]]))}
la1li <- sapply(1:length(NewCG$ordlay),la1)
v1 = which(la1li == 1)


lowlay = 1*(v1 < v2)

if(v1 == v2)
{
  epx1tox2 = 1
  epx2tox1 = 1
   if(1*(setequal(NewCG$TOTPAR[[k2]],NewCG$TOTPAR[[k1]]))==0)
  {
    ueposs = 0
  }
   if(1*(setequal(NewCG$TOTPAR[[k2]],NewCG$TOTPAR[[k1]]))==1)
     {
     ueposs = 1
   }
}
if(v1 < v2)
{
  epx1tox2 = 1
  epx2tox1 = edgeposs(NewCG,v2,v1,k2,k1)
  ueposs = 0
}
if(v1 > v2)
{
  epx2tox1 = 1
epx1tox2 = edgeposs(NewCG,v1,v2,k1,k2)
ueposs = 0
}

f <- function(z){CCSplit(z,Cmat)$count}
f1 <- sapply(1:nrow(Cmat),f)


if(ueposs == 1)
{
  ff1 = 1+f1
  reversepr = (1/(d*(d-1)))*(sum(1/(ff1)))
}

if(ueposs ==0 & epx1tox2 + epx2tox1 ==2)
{
  reversepr = (1/(d*(d-1)))*(sum(1/f1))
}

if(ueposs ==0 & epx1tox2+epx2tox1==1)
{
  reversepr = (2/(d*(d-1)))*(sum(1/f1))
  }
}


Ratio = reversepr/forwardpr
Out$CG = NewCG
    Out$Ratio <- Ratio

    return(Out)
}
