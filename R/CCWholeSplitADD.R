#' Finds possible splits of a CC for adding within-CC directed edge
#'
#'This is called from ADDEdgeSameCCNotAdd for the purpose of finding possible
#'splits of a CC into CC2, CC1, CC3 where CC2 is connected and contains x[2],
#'CC1 is connected and contains a complete set CC11 where all connections
#'between CC11 and
#'CC2 exist, CC3 remaining vertices.
#'
#'
#'@param x1 a vertex which has to be in CC1
#' @param x2 a vertex which has to be in CC2
#' @param Cmat lower triangular, 1's and 0's, 1's denote undirected
#' edges within the chain component.
#'
#'
#'
#' @return
#' Out where
#'Out$CC2 lists all possibilities for CC2
#'Out$CC1 lists all possiblilities for CC1
#'Out$CC3 list of sets of remaining vertices
#'where Out$CC1[[j]], Out$CC2[[j]], Out$CC3[[j]] is jth
#'possible split.
#'Out$count gives number of ways.
#'


CCWholeSplitADD <- function(x1,x2,Cmat)
{
oot = intersect(unlist(nhd(Cmat)$L[x1]),unlist(nhd(Cmat)$L[x2]))
withinmc = sort(c(x2,oot))
mbcc1 = sort(c(x1,oot))

#mbcc1 contains vertices that MC containing x1 on junction tree must have
#withinmc contains vertices that adjacent junction tree containing x2 must
#have (junction trees not unique - may be several candidates)

#oot is the separator between the two

MaCl = MC(Cmat)

h <- function(j){1*(all(mbcc1 %in% MaCl[[j]]))}
h1 <- sapply(1:length(MaCl),h)
h2 = which(h1==1)


f <- function(j){1*(all(withinmc %in% MaCl[[j]]))}
f1 <- sapply(1:length(MaCl),f)
immmcs = which(f1==1)
#lists candidates for MC containing x2, adjacent to MC containing x2

g <- function(j){1*(length(setdiff(MaCl[[j]],withinmc))>=1)}
g1 <- sapply(immmcs,g)
ssuumm = sum(g1)

if(ssuumm==0)
{
  count = 0

  CC1 <- list()
  CC1[[1]] = vector(mode = "integer", length = 0L)
  CC2 <- list()
  CC2[[1]] = vector(mode="integer",length=0L)
  CC3 <- list()
  CC3[[1]] = vector(mode = "integer",length = 0L)

  Out = list()
  Out$count = 0
  Out$CC1 = CC1
  Out$CC2 = CC2
  Out$CC3 = CC3
  return(Out)
}


g2 = which(g1==1)
goodcls = immmcs[g2]

#goodcls gives the labels of the MCs that can be used to generate
#immoralites (so that x[1] -> CC2 when added will be compelled)
#To get a directed edge, we need an MC containing withinmc vertices
#which contains additional vertices to make the immorality.


lCC2 = nrow(Cmat)

set = list()
set[[1]] = seq(1,3)
all_splits = unname(as.matrix(expand.grid(rep(set,lCC2))))

leng = nrow(all_splits)
f <- function(j){1*(all(all_splits[j,mbcc1]==1))*(all_splits[j,x2]==2)}
f1 <- sapply(1:leng,f)
ap = which(f1==1)
all_splits = as.matrix(all_splits[ap,,drop=FALSE])

#we have all possible splits into three disjoint subsets where CC1
#contains x1 and CC2 contains x2.
#we want CC1 the PointOut vertices. These include x1 (these denoted by 1).
#CC2 the PointIn, connected, containing x2, directed edge x1 -> CC2 added
#CC3 all the other vertices

leng = nrow(all_splits)
a = seq(1,leng)

f <- function(j)
{
an = 1
CC1 = which(all_splits[j,]==1)
CC2 = which(all_splits[j,]==2)
diffCC1 = setdiff(CC1,mbcc1)
inanMC = sort(c(diffCC1,CC2))

if(length(diffCC1)==0)
{
  an=0
}
if(an==1)
{
h <- function(k){1*(all(inanMC %in% MaCl[[k]]))}
h1 = sapply(goodcls,h)
if(sum(h1) == 0)
{
  an=0
}
}
return(an)
}

a = sapply(1:leng,f)
ap = which(a==1)
all_splits = as.matrix(all_splits[ap,,drop=FALSE])
count = nrow(all_splits)

CC1 <- vector(mode="list",length=count)
CC1 <-lapply(CC1,as.integer)
CC2 <- vector(mode="list",length=count)
CC2 <-lapply(CC2,as.integer)
CC3 <- vector(mode="list",length=count)
CC3 <-lapply(CC2,as.integer)

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
