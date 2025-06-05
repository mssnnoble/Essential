#'Removes undirected edge where removal does not split the CC into two CCs
#'
#' @param CG - an essential graph object
#' @param x = (x[1],x[2]) two vertices in the same CC
#' @param v - the layer
#' @param w - the CC
#'
#' @return: Out where
#' \itemize{
#' \item Out$CG is the new essential graph
#' \item Out$Ratio is the ratio between reverse and forward proposals
#' }
#'@export

RMUESameCC <- function(CG,x,v,w)
{
  #Here removal of an undirected edge does not split the CC.
  #If edge is removable, all vee structures introduced come from same
  #MC (in graph before removal)
  #If Nvee is the number of vee structure vertices, there are 2^Nvee
  #possibilities for deciding on vee structures.
  #If edge is not removable, it appears in more than one MC.
  # vee structures are introduced, some of which must be resolved as
  #immoralities. One of these MC (at most) can be chosen for non-immoralities
  #if there are MCvee in that MC, we have 2^MCvee possible ways for choosing
  #non-immoralities, all the rest are immoralities.

f <- function(j){length(CG$CC[[j]])}
f1 <- sapply(1:length(CG$CC),f)
d = sum(f1)


  Out <- list()


  NewCG <- CG
  CCold = seq(1,length(CG$CC[[w]]))

  Par <- vector(mode="list",length = length(CCold))
  Par <- lapply(Par,as.integer)

  #Par[[j]] will be the within-old-CC parents of vertex j
  #(so that all vertices within a new CC have same parents)

  y <- sort(match(x,CG$CC[[w]]))
  z <- sort(y,decreasing = TRUE)

  #CCold is CC co-ordinates for GC$CC[[w]]
  #y = (y[1],y[2]) is x in these co-ordinates


  Cmat <- CG$CCued[[w]]

  MC <- MC(Cmat)
  #get maximal cliques of Cmat

  Cmat[z[1],z[2]] <- 0


  JTree <- kruskal(MC)
  SEQ <- seq(1,length(MC))

  #get junction tree and its adjacency matrix
  #this is necessary for dealing with compelled
  #edges following Meek's rules (no flags, no partially
  #directed cycles)



  f<-function(j){1*(all(y %in% MC[[j]]))}
  f1 <- sapply(1:length(MC),f)
  MCvee <- which(f1==1)

  #MCvee gives 1 for those MCs which contain both x1 and x2
#exactly one value in MCvee if (and only if) edge is removable

  totsets = list()
  totsets[[1]] = vector(mode="integer",length=0L)
  for(j in MCvee)
  {
    PossNonIm <- setdiff(MC[[j]],y)
    aaa = length(PossNonIm)
    all_subsets <- lapply(1:aaa,function(k) combn(aaa,k,simplify=FALSE))
    all_subsets <- unlist(all_subsets,recursive=FALSE)
    all_subsets = as.list(all_subsets)
    f <- function(j){PossNonIm[all_subsets[[j]]]}
    new_subsets = lapply(1:length(all_subsets),f)
    totsets = unique(append(totsets,new_subsets))
  }

#the totsets[[1]] entry adds the case where all vee structures are immoralities.

  ln = length(totsets)
  POSSIB = ln

  jay = sample(ln,1)
  SS0 = totsets[[jay]] # these are the vertices that are non-immoralities
  #for vee-structures x1 - z - x2, z in SS0

  #We now find SS1


SETTT = sort(unique(unlist(MC[MCvee])))


 SS1 = setdiff(SETTT,sort(c(y,SS0))) #immoralities


  PointIn <- SS1
  PointOut <- c(SS0,y)

  NEWW = Compelled(Cmat,MC,PointOut,PointIn)
  NNewCC = NEWW$NNewCC
  Par = NEWW$Par

nnewCC = length(NNewCC)
  sss <- 1:nnewCC
  PA <- vector(mode="list", length=nnewCC)
  PA <-lapply(PA,as.integer)

  f <- function(j){sort(c(CG$CC[[w]][Par[[j]]],CG$TOTPAR[[w]]))}
  PA <- lapply(seq_along(sss),f)




NewCG$CC[[w]] <- CG$CC[[w]][NNewCC[[1]]]
NewCG$CCued[[w]] <- as.matrix(Cmat[NNewCC[[1]],NNewCC[[1]]])
NewCG$TOTPAR[[w]] <- PA[[1]]
NewCG$DEi2ip1[[w]] <- PA[[1]]
NewCG$DEex[[w]] <- PA[[1]]

if(nnewCC >=2)
{
old <- length(CG$CC)
  for(j in 1:(nnewCC-1))
  {
    NewCG$CC[[old+j]] <- CG$CC[[w]][NNewCC[[j+1]]]
  }

  for(j in 1:(nnewCC-1))
  {
    NewCG$CCued[[old+j]] <- as.matrix(Cmat[NNewCC[[j+1]],NNewCC[[j+1]]])
  }



  for(j in 1:(nnewCC-1))
  {
    NewCG$TOTPAR[[old+j]] <-PA[[j+1]]
    NewCG$DEi2ip1[[old+j]] <-PA[[j+1]]
    NewCG$DEex[[old+j]] <- PA[[j+1]]
  }
}

  #we now have the total parent sets for the updated CC list.
  #now to get them into the correct layerings,
  #the i2ip1 and DEex are wrong, but have been added so that the
  #object is of the correct format to apply NEWLAY (which only deals
  #with TOTPAR)

NewCG <- NEWLAY(NewCG)

#now the forward and reverse probabilities


forwardpr = (2/(d*(d-1)))*(1/POSSIB)

#now compute the reverse proposal, which is adding undirected
#edge to NewCG, to get CG.
#
#note different rules apply if x1 and x2 are in the same CC or
#different CCs.
#
#if they are in the same CC, then the edge is addable - so there
#is exactly one possible move - add undirected edge. This happens if
#either (x1,x2) or (x2,x1) are chosen.
#
#if x1 and x2 are in different CCs in NewCG, these CCs have same
#parents; TOTPAR for both is CG$TOTPAR[[w]]

#firstly, find the CCs which contain x[1] resp. x[2]
f <- function(j){1*(x[1] %in% NewCG$CC[[j]])}
f1 <- sapply(1:length(NewCG$CC),f)
k1 = which(f1==1)
#k1 the index of CC which contains x[1]

g <- function(j){1*(x[2] %in% NewCG$CC[[j]])}
g1 <- sapply(1:length(NewCG$CC),g)
k2 = which(g1==1)
#k2 the index of CC which contains x[2]

if(k1==k2)
{
  reversepr = 2/(d*(d-1))
}

if(k1!=k2)
{
  y1 = match(x[1],NewCG$CC[[k1]])
  y2 = match(x[2],NewCG$CC[[k2]])
  Cmat1 = NewCG$CCued[[k1]]
  Cmat2 = NewCG$CCued[[k2]]
  ANI1 = CCSplit(y1,Cmat1)
  count1 = ANI1$count
  ANI2 = CCSplit(y2,Cmat2)
  count2 = ANI2$count

  f <- function(j){1*(length(ANI1$CC1[[j]])>=1)}
  f1 <- sapply(1:count1,f)
  N1 = 1 + sum(f1)

  #note: with same parents, need non-trivial CC1 to do add directed edge move

  g <- function(j){1*(length(ANI2$CC1[[j]])>=1)}
  g1 <- sapply(1:count2,g)
  N2 = 1+sum(g1)

reversepr = (1/(d*(d-1)))*((1/N1)+(1/N2))
}

ratio = reversepr/forwardpr
 Out$CG = NewCG
  Out$Ratio = ratio
  return(Out)
}
