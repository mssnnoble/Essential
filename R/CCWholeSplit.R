#' Finds possible splits of a CC for adding within-CC directed edge
#'
#'This is called from ADDEdgeSameCCNotAdd for the purpose of finding possible
#'splits of a CC into CC2, CC1, CC3 where CC2 is connected and contains x[2],
#'CC1 is   a complete set  where all connections
#'between CC1 and
#'CC2 exist, CC3 remaining vertices.
#'
#'CCWholeSplit deals with <x[1],x[2]> not addable, it sends to CCWholeSplitADD
#'in the case where <x[1],x[2]> is addable.
#'
#'
#'@param x1 a vertex from which new directed edge proceeds.
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
#'@export


CCWholeSplit <- function(x1,x2,Cmat)
{
  xpar = sort(c(x1,x2),decreasing = TRUE)
  litmus = addedgemat(Cmat)[xpar[1],xpar[2]]

  if(litmus == 0)
  {
    #litmus = 0 means undirected edge not addable - x1 and x2 are not in
    #adjacent MCs on any junction tree
  lCC2 = nrow(Cmat)
   MaCl = MC(Cmat)
   JTr = kruskal(MaCl)
   xx2 = c(x2)
   path = MCPath(x1,xx2,MaCl,JTr)
   lep = length(path)

#get path from CC containing x1 to CC containing x2
#the MC path algorithm ensures that neither x1 nor x2 appear
#in intervening MCs on path.

i2 = path[2]
#first that does not contain x1
    i3 = path[3]
   JTrmod = JTr
   JTrmod[i2,i3] = 0
   JTrmod[i3,i2]=0

   lE = nrow(JTrmod)
   Ctrav <- list()
   Ctrav[[1]] <- JTrmod + t(JTrmod) + diag(lE)
   if(lE >=2)
   {
     for(k in 2:lE){
       Ctrav[[k]] <- Ctrav[[k-1]] %*% (JTrmod+t(JTrmod))
     }
   }
   Econn <- Reduce('+',Ctrav)
   Econn <- sign(Econn)
   concomp = unique(Econn)
   part1 = which(concomp[1,]==1)
   part2 = which(concomp[2,]==1)
comp1 = sort(unique(unlist(MaCl[part1])))
comp2 = sort(unique(unlist(MaCl[part2])))
f = 1*(x1 %in% comp1)
g = 1*(x1 %in% comp2)
if(f==1)
{
  Acand = comp1
}
if(g==1)
{
  Acand = comp2
}

#Acand contains x1

# task is to find all possible decompositions of
#1:length(CC) into CC1* and CC2 where CC1* = CC1 union CC12, x1 in CC12,
#x2 in CC2 and
  #1) there is an edge from each CC1 to each CC2 and
  #2) no edges from CC12 to any CC2 or CC3
  #3) CC1 is a complete set
  #4) nhd(CC11) does not contain x1


A1 = sort(intersect(nhd(Cmat)$L[[x2]],Acand))

  #neighbours of x2 not including x2. All the i2ip1 parents for
  #CC2 must be contained in this. Each non-empty subset of POSSCC11
  #gives a candidate CC11


    set = list()
    set[[1]] = seq(1,3)
    all_splits = unname(as.matrix(expand.grid(rep(set,lCC2))))
    #we have all possible splits into three disjoint subsets

    leng = nrow(all_splits)
    f <- function(j){1*(all_splits[j,x2]==2)*(all_splits[j,x1]==1)}
    f1 <- sapply(1:leng,f)
    ap = which(f1==1)
    all_splits = as.matrix(all_splits[ap,,drop=FALSE])
    #we have removed those for which x2 is not in CC2 and x1 is not in CC1
    #hence both CC1 and CC2 are non-empty



    f <- function(i){1*(all(Acand %in% which(all_splits[i,]==1)))}
f1 = sapply(1:nrow(all_splits),f)
ap = which(f1==1)
all_splits = as.matrix(all_splits[ap,,drop=FALSE])

leng = nrow(all_splits)
a = rep(1,leng)
   for(i in 1:leng)
{
  A = which(all_splits[i,]==1)
B = which(all_splits[i,]==2)
C = which(all_splits[i,]==3)

f <- function(z){intersect(unlist(nhd(Cmat)$L[z]),A)}
L <- lapply(B,f)
Lp <- unique(L)

if(length(Lp)>=2)
{
  a[i]=0
}
if(a[i]==1)
{
  A1 = Lp[[1]]
  b = length(A1)
  if(b==0)
  {
    a[i]=0
  }
  if(b >=1)
  {
  if(sum(as.matrix(Cmat[A1,A1,drop=FALSE]))!=0.5*b*(b-1))
  {
    a[i]=0
  }
}
   }
}
ap = which(a==1)
all_splits = as.matrix(all_splits[ap,,drop=FALSE])

leng = nrow(all_splits)
a = rep(1,leng)
for(i in 1:leng)
{
  A = which(all_splits[i,]==1)
  B = which(all_splits[i,]==2)
  C = which(all_splits[i,]==3)

  A1 = sort(intersect(nhd(Cmat)$L[[x2]],A))

  path = MCPath(x1,B,MaCl,JTr)
  lep = length(path)
  f <- function(j){1*(x1 %in% MaCl[[path[j]]])}
  f1 <- sapply(1:lep,f)


  sett = sort(unique(c(A1,B)))
  if(1*(all(MaCl[[path[lep]]] %in% sett))==0)
  {
    a[i]=0
  }
  if(a[i]==1)
  {
  i2 = path[lep]
  #first that does not contain x1
  i3 = path[lep-1]
  JTrmod = JTr
  JTrmod[i2,i3] = 0
  JTrmod[i3,i2]=0

  lE = nrow(JTrmod)
  Ctrav <- list()
  Ctrav[[1]] <- JTrmod + t(JTrmod) + diag(lE)
  if(lE >=2)
  {
    for(k in 2:lE){
      Ctrav[[k]] <- Ctrav[[k-1]] %*% (JTrmod+t(JTrmod))
    }
  }
  Econn <- Reduce('+',Ctrav)
  Econn <- sign(Econn)
  concomp = unique(Econn)
  part1 = which(concomp[1,]==1)
  part2 = which(concomp[2,]==1)
  comp1 = sort(unique(unlist(MaCl[part1])))
  comp2 = sort(unique(unlist(MaCl[part2])))
  f = 1*(x1 %in% comp1)
  g = 1*(x1 %in% comp2)
  if(f==1)
  {
    Acand = comp1
  }
  if(g==1)
  {
    Acand = comp2
  }
  if(1*(setequal(Acand,A))==0)
  {
    a[i]=0
  }

  }



  if(a[i]==1)
      {
        f <- function(j){length(intersect(B,MaCl[[path[j]]]))}
        f1 <- sapply(1:(lep - 1),f)
        g1 = sum(f1)
        if(g1 >=1)
        {
          a[i]=0
        }
      }
    #intersects with CC2
      if(a[i]==1)
      {
      bigset = sort(c(A1,B))
    SettDifff = length(setdiff(MaCl[[i2]],bigset))
    if(SettDifff >= 1)
    {
      a[i]=0
    }
      }

  }



ap = which(a==1)
all_splits = as.matrix(all_splits[ap,,drop=FALSE])


#we have not yet checked CC2 connected

leng = nrow(all_splits)
a = rep(1,leng)

for(i in 1:leng)
{
  CC2 = which(all_splits[i,]==2)
  CMX = as.matrix(Cmat[CC2,CC2,drop=FALSE])
  lE = nrow(CMX)
  Ctrav <- list()
  Ctrav[[1]] <- CMX + t(CMX) + diag(lE)
  if(lE >=2)
  {
    for(k in 2:lE){
      Ctrav[[k]] <- Ctrav[[k-1]] %*% (CMX+t(CMX))
    }
  }
  Econn <- Reduce('+',Ctrav)
  Econn <- sign(Econn)
  if(sum(Econn) != lE*lE)
  {
    a[i]=0
  }
}
ap = which(a==1)
all_splits = as.matrix(all_splits[ap,,drop=FALSE])

count = nrow(all_splits)

CC1 <- vector(mode="list",length=count)
CC1 <-lapply(CC1,as.integer)
CC2 <- vector(mode="list",length=count)
CC2 <-lapply(CC2,as.integer)
CC3 <- vector(mode="list",length=count)
CC3 <-lapply(CC2,as.integer)

fcc1 <- function(j){intersect(which(all_splits[j,]==1),nhd(Cmat)$L[[x2]])}
CC1 <- lapply(1:count,fcc1)
fcc2 <- function(j){which(all_splits[j,]==2)}
CC2 <- lapply(1:count,fcc2)
fcc3 <- function(j)
  {
  sort(c(which(all_splits[j,]==3),setdiff(which(all_splits[j,]==1),nhd(Cmat)$L[[x2]])))
  }
CC3 <- lapply(1:count,fcc3)

Out = list()
Out$count = count
Out$CC1 = CC1
Out$CC2 = CC2
Out$CC3 = CC3

return(Out)
  }
if(litmus == 1)
{
 Out = CCWholeSplitADD(x1,x2,Cmat)
}
}
