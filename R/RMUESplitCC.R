#'Removes undirected edge  where removing splits the CC into two CCs.
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

RMUESplitCC <- function(CG,x,v,w)
{


  Out <- list()

  f <- function(j){length(CG$CC[[j]])}
  f1 <-sapply(1:length(CG$CC),f)
  d=sum(f1)

  NewCG <- CG
  #NewCG will be the output CG




#now establish the ordlay, the CCs and the CCueds after the split.

l = length(CG$CC[[w]])
CClc <- seq(1,l)

#CC in local co-ordinates

l2 = length(CG$CC)

NewCG$ordlay[[v]] <- sort(c(l2+1,CG$ordlay[[v ]]))

#CC is split therefore layer v has an extra CC, which we label l+1
#We now compute the membership of the two CCs and also both the CCueds

y = sort(match(x,CG$CC[[w]]),decreasing=TRUE)

Cmat <- CG$CCued[[w]]
Cmat[y[1],y[2]] = 0

#edge removed. Now compute the connected components

dim = nrow(Cmat)
E1 = Cmat + t(Cmat)+diag(dim)

E = list()
E[[1]] = E1
for(j in 2:dim)
{
  E[[j]] = E[[j-1]]%*%E1
}
F = Reduce('+',E)
F = sign(F)

  #computes neighbourhoods after edge removal

  Con1 <- which(F[1,]==1)
  #neighbours of vertex labelled 1 (in CC co-ordinates)

  Con2 <- sort(setdiff(CClc,Con1))

  #here we know a priori that there are two connected
  #components.

  NewCG$CC[[w]] <- CG$CC[[w]][Con1]
  NewCG$CC[[l2+1]] <- CG$CC[[w]][Con2]

  #We now have the two new CCs

  len1 = length(NewCG$CC[[w]])
  len2 = length(NewCG$CC[[l2+1]])

  NewCG$DEi2ip1[[l2+1]] = CG$DEi2ip1[[w]]
  NewCG$DEex[[l2+1]] = CG$DEex[[w]]
  NewCG$TOTPAR[[l2+1]] = CG$TOTPAR[[w]]

  NewCG$CCued[[w]] <- as.matrix(Cmat[Con1,Con1])
  NewCG$CCued[[l2+1]] <- as.matrix(Cmat[Con2, Con2])

  #we now have the two new CCueds

  forwardpr = 2/(d*(d-1))

  #reversepr is probability of adding an undirected edge.

  g1 = 1*(x[2] %in% NewCG$CC[[w]])
  g2 = 1*(x[2] %in% NewCG$CC[[l2+1]])

  wst2 = w*g1 + (l2+1)*g2
  #wst2 is label of CC containing x[2]

  h1 = 1*(x[1] %in% NewCG$CC[[w]])
  h2 = 1*(x[1] %in% NewCG$CC[[l2+1]])
  wst1 = w*h1 + (l2+1)*h2
  #wst1 is label of CC containing x[1]
  #these should have different labels.

  Cmat1 = NewCG$CCued[[wst1]]
  x1 = match(x[1],NewCG$CC[[wst1]])
  Cmat2 = NewCG$CCued[[wst2]]
  x2 = match(x[2],NewCG$CC[[wst2]])



  ANI1 = CCSplit(x1,Cmat1)
  count1 = ANI1$count
  ANI2 = CCSplit(x2,Cmat2)
  count2 = ANI2$count

  f <- function(j){1*(length(ANI1$CC1[[j]])>=1)}
  f1 <- sapply(1:count1,f)

  g <- function(j){1*(length(ANI2$CC1[[j]])>=1)}
  g1 <- sapply(1:count2,g)

N1 = 1+sum(f1)
N2 = 1+sum(g1)


   reversepr = (1/(d*(d-1)))*((1/N1)+(1/N2))

  ratio = reversepr/forwardpr


    Out$CG <- NewCG
  Out$Ratio <- ratio
   return(Out)
}
