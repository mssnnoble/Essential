#' Directs compelled edges
#'
#'After generating new immoralities by adding a directed edge and directing
#'some undirected edges, this forces direction of compelled edges (previously
#'undirected) to ensure no other immoralities and no directed or partially
#'directed cycles
#'
#' @param Cmat this is the lower triangular directed edge matrix of a CC
#' @param MC This is the list of maximal cliques corresponding to Cmat
#' @param PointOut This is list of vertices from which arrows are pointing out to
#' form the new immoralities
#' @param PointIn This is vertices into which arrows are pointing
#'
#'
#' @return
#' \itemize{
#' \item Out$NewCC this gives a list of the new CCs generated when compelled
#' edges are appropriately directed
#' \item  Out$Par This gives a list of the additional parents (due to directing
#' compelled edges)
#' }
#'
#' @export


Compelled <- function(Cmat,MC,PointOut,PointIn)
{

S <- PointIn
len = ncol(Cmat)

f <- function(j){1*(length(intersect(PointIn,MC[[j]]))>=1)}
f1 <- sapply(1:length(MC),f)
INTSCL = which(f1==1)

JTree <- kruskal(MC)

Par1 <- vector(mode="list",length=len)
Par1 <-lapply(Par1,as.integer)

#Par1 is the *within CC* parents of the vertices of the CC.

MCliq = seq(1,length(MC))

neighbours = nhd(Cmat)$L
Cmat[PointOut,PointIn] <- 0
Cmat[PointIn,PointOut] <- 0
POINT = sort(c(PointIn,PointOut))
for(j in S)
{
  Par1[[j]] =  PointOut
  for(k in setdiff(neighbours[[j]],POINT))
  {
    npark = intersect(neighbours[[k]],POINT)
    Par1[[k]] <- unique(sort(c(Par1[[k]],npark)))
    PointIn <- sort(c(PointIn,k))
    Cmat[npark,k] = 0
    Cmat[k,npark] = 0
    POINT <- sort(c(POINT,k))
  }
}

unused = seq(1,length(MC))
used = vector(mode="integer",length=0L)

k= INTSCL[1]

jt = JTree + t(JTree)
#symmetrise junction tree matrix to make it easier to see
#UNUSED MCs which are reachable in one step from USED MCs.

repeat
{
  if(length(unused)==0)
  {
    break
  }
  if(length(unused) <= length(MC)-1)
  {
    f <- function(j){intersect(which(jt[used[j],]==1),unused)}
    f2 <- lapply(1:length(used),f)
    g2 = unique(unlist(f2))
    if(length(intersect(g2,INTSCL))>=1)
    {
      k = intersect(g2,INTSCL)[1]
    }
    if(length(intersect(g2,INTSCL))==0)
    {
      k = g2[1]
    }
    #k gives label of an unused MC, which is adjacent to
    #a used MC on some junction tree. We're about to use it
    #so add it to used and remove from unused.
  }
  used <-sort(unique(c(used,k)))
  unused <- unused[-match(k,unused)]
  if(length(unused)<=(length(MC)-2))
  {
    jt[used,k]=0
    jt[k,used]=0
  }
  Point = unique(c(PointIn,PointOut))
  subset = as.integer(intersect(MC[[k]],Point))
  subset2 = as.integer(intersect(MC[[k]],PointIn))
  if(length(subset2)>=1)
  {
    PINew = setdiff(MC[[k]],subset)

    Cmat[subset,PINew]= 0
    Cmat[PINew,subset]= 0
    PointIn = sort(c(PointIn,PINew))

    f <- function(j){unique(sort(union(Par1[[j]],subset)))}
    Par1[PINew] = lapply(PINew,f)
  }
}

#at this stage, we have the within-old-CC parents (of the new CCs)
#and Cmat gives the undirected edges.
#We now compute the connected components

lE <- nrow(Cmat)
Cmat <- Cmat + diag(lE)
Ctrav <- list()
Ctrav[[1]] <- Cmat + t(Cmat)
if(lE >=2)
{
  for(j in 2:lE){
    Ctrav[[j]] <- Ctrav[[j-1]] %*% (Cmat+t(Cmat))
  }
}
Econnect <- Reduce('+',Ctrav)
Econnect <- sign(Econnect)

#Econnect[i,j]==1 if and only if i and j are in the same
#connected components

lenn = nrow(Econnect)
f <- function(j){which(Econnect[j,]==1)}
NNewCC <- unique(lapply(1:lenn,f))

#Par1 is the list of parents *for each vertex* we now construct a list
#for each CC

ln = length(NNewCC)
Par = vector(mode = "list", length = ln)
Par = lapply(Par,as.integer)
f <- function(j){Par1[[NNewCC[[j]][1]]]}
Par = lapply(1:ln,f)

Out <- list()
Out$NNewCC <- NNewCC
Out$Par <- Par
return(Out)
}
