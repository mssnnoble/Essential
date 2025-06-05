#' Finds possible splits of a CC for (partial) removal of directed edge
#'
#'This is called from RemDE, ADEDiffPar, ADDSamParDiffCC for the purpose of
#'finding possible connected CC2 (add DE for add moves, remove DE for RemDE)
#'and a CC1 which is a complete set and all CC2 vertices connected to each
#'CC1 vertex. (CC3 is simply the remaining vertices not in CC1 or CC2) These
#'edges become directed CC1 -> CC2 in the add/rem algorithms.
#'This includes empty set for CC1 - such possibilities have to be removed
#'for the ADDSamParDiffCC (otherwise directed edges not compelled)
#'
#'
#' @param y a vertex which has to be in CC2
#' @param Cmat lower triangular, 1's and 0's, 1's denote undirected
#' edges within the chain component.
#'
#'
#'
#' @return
#' Out where
#' \itemize{
#' \item Out$CC2 lists all possible connected subsets containing y
#' \item Out$CC1
#' \item Out$CC3
#' \item where (Out$CC1[[j]], Out$CC2[[j]], Out$CC3[[j]]) is jth
#'possible split.
#' \item Out$count gives number of ways.
#'}
#'
#'@export



CCSplit <- function(y,Cmat)
{
   lCC2 = nrow(Cmat)
   ALLCC = seq(1,lCC2)

b = lCC2
  if(b == 1)
  {
    #we have a singleton CC
    CC2 = list()
    CC2[[1]]=ALLCC
    count=1
    Out = list()
    Out$count = count
    CC1 = list()
    CC1[[1]] = vector(mode = "integer", length = 0L)
    CC3 = list()
    CC3[[1]] = vector(mode = "integer", length = 0L)

    Out$CC1 = CC1
    Out$CC2 = CC2
    Out$CC3 = CC3

    return(Out)
  }
  if(b >= 2)
  {
    set = list()
    set[[1]] = seq(1,3)
    all_splits = unname(as.matrix(expand.grid(rep(set,lCC2))))
    #we have all possible splits into three disjoint subsets

    ap = which(all_splits[,y]==2)
    all_splits = all_splits[ap,]
    #we have removed those for which y is not in CC2

    leng = nrow(all_splits)
    a = rep(1,leng)

#now we insist that B is connected

    for(i in 1:leng)
    {
B = which(all_splits[i,]==2)

MTRXB = as.matrix(Cmat[B,B,drop=FALSE])
lE = nrow(MTRXB)
Ctrav <- list()
Ctrav[[1]] <- MTRXB + t(MTRXB) + diag(lE)
if(lE >=2)
{
  for(k in 2:lE){
    Ctrav[[k]] <- Ctrav[[k-1]] %*% (MTRXB+t(MTRXB))
  }
}
Econn <- Reduce('+',Ctrav)
Econn <- sign(Econn)

#Econn[i,k]==1 if and only if i and k are in the same
#connected components
if(sum(Econn) != lE*lE)
{
  a[i]=0
}
}


    ap = which(a==1)
    all_splits = as.matrix(all_splits[ap,,drop=FALSE])

  #now B is connected
leng = nrow(all_splits)
a = rep(1,leng)

#next, we insist that the set A is complete
#and we insist that all elements of A are
#neighbours of each element of B

for(i in 1:leng)
{
  A = which(all_splits[i,]==1)
  B = which(all_splits[i,]==2)
if(length(A)>=1)
{
  lA = length(A)
  if(sum(Cmat[A,A])!=0.5*lA*(lA-1))
  {
    a[i]=0
  }
  if(a[i]==1)
  {
f <- function(z){length(setdiff(A,nhd(Cmat)$L[[z]]))}
f1 = sapply(B,f)
if(sum(f1)>=1)
{
  a[i]=0
}
  }
}
}
ap = which(a==1)
all_splits = as.matrix(all_splits[ap,,drop=FALSE])


#we now have all possible splits into three subsets

count = nrow(all_splits)

#count gives the total number of legal splits

CC1 <- vector(mode="list",length=count)
CC1 <-lapply(CC1,as.integer)
CC2 <- vector(mode="list",length=count)
CC2 <-lapply(CC2,as.integer)
CC3 <- vector(mode="list",length=count)
CC3 <-lapply(CC3,as.integer)

fcc1 <- function(j){which(all_splits[j,]==1)}
CC1 <- lapply(1:count,fcc1)
fcc2 <- function(j){which(all_splits[j,]==2)}
CC2 <- lapply(1:count,fcc2)
fcc3 <- function(j){which(all_splits[j,]==3)}
CC3 <- lapply(1:count,fcc3)

Out = list()
Out$count = count
Out$CC1 = CC1
Out$CC2 = CC2
Out$CC3 = CC3
 return(Out)
}
}
