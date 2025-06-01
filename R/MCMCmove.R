#' Randomly selects two different vertices and computes a proposal.
#'
#' @param CG an essential graph object
#'
#' @return Out  a list where
#' Out$OldCG = CG (the original - input - essential graph)
#' Out$NewCG <- NewCG (the proposed graph)
#' Out$Ratio <- Ratio (ratio of reverse proposal to forward proposal)
#'
#'
#'@export


MCMCmove <- function(CG){


  d =  sum(unlist(lapply(CG$CC,length)))

  #so that d is the number of vertices

  SCC = 0
  TPSAME = 0
  x <- sample(1:d,2,replace=FALSE)
  print(x)

   #We add an edge if there is no x[1],x[2] edge (choose between
  #possibilities that do not involve changing parent sets of the CCs).
  #We delete if there is an edge.

  #We first compute v[1] and v[2] the layers of x[1] and x[2] resp.
  #and also w[1] and w[2] the CCs of x[1] and x[2] resp.
  #
  #If x[1] -> w[2], delete edge from x[1] to x[2] when v[1] != v[2] means
  #delete the whole vertex to CC directed edge
  #(so, when computing the proposal probability, same proposal for any
  #x[2] in w[2])

v <- vector(mode="integer",length = 2)
w <- vector(mode="integer",length=2)


#v for layers, w for CCs

#get layers
 f1 <-function(b){1*(x[1] %in% CG$VIL[[b]])}
 f2 <- function(b){1*(x[2] %in% CG$VIL[[b]])}
 g1 = sapply(1:length(CG$VIL),f1)
 g2 = sapply(1:length(CG$VIL),f2)

v[1] <-which(g1==1)
v[2] <- which(g2==1)

#get CCs
h1 <- function(b){1*(x[1] %in% CG$CC[[b]])}
h2 <- function(b){1*(x[2] %in% CG$CC[[b]])}
h1a <- sapply(1:length(CG$CC),h1)
h2a <- sapply(1:length(CG$CC),h2)
w[1] <- which(h1a==1)
w[2] <- which(h2a==1)

print(v)
print(w)

if(v[1]!=v[2])
{
  SL = 0
}
if(v[1]==v[2])
{
  v = v[1]
  SL=1
#SL = 1 denotes same layer
if(w[1]!=w[2])
{
SCC = 0
#SCC = 1 denotes same CC SCC = 0 denotes different CCs
    TPSAME = 1*(setequal(CG$TOTPAR[[w[1]]],CG$TOTPAR[[w[2]]]))
}

  #TPSAME means that the parent sets of w[1] and w[2] are the same.
if(w[1] == w[2])
  {
  w=w[1]
  SCC = 1
  }
}

#now to decide on addition or removal

if(SCC==1)
  {
y <- match(x,CG$CC[[w]])
z = sort(y,decreasing = TRUE)
if(CG$CCued[[w]][z[1],z[2]]==1)
{
   #hence there is an undirected edge
  print("RMUE")
  NEW = RMUE(CG,x,v,w)
  NewCG = NEW$CG
  Ratio = NEW$Ratio
}
if(CG$CCued[[w]][z[1],z[2]]==0)
{
  #hence no edge
  if(addedgemat(CG$CCued[[w]])[z[1],z[2]] ==1)
  {
    #so edge addable
    NEW = ADDEdgeSameCCAddable(CG,x,v,w)
    NewCG = NEW$CG
    Ratio = NEW$Ratio
  }
  if(addedgemat(CG$CCued[[w]])[z[1],z[2]]==0)
  {
    #edge not addable
    NEWCGTHING = ADDEdgeSameCCNotAdd(CG,x,v,w)
    NewCG = NEWCGTHING$CG
    Ratio = NEWCGTHING$Ratio
  }
}
}
if(SCC==0 & TPSAME==1)
{
 #so same layer, different CCs, same parent sets
  NEW = ADDSamParDiffCC(CG,x,v,w)
  NewCG = NEW$CG
  Ratio = NEW$Ratio
}
if(SCC==0 & TPSAME ==0)
{
 addremde = 1*(x[1] %in% CG$TOTPAR[[w[2]]]) + 1*(x[2] %in% CG$TOTPAR[[w[1]]])
 if(addremde == 1)
 {
   #this means that there is a directed edge which is to be removed
   NEW = RemDE(CG,x,v,w)
   NewCG = NEW$CG
   Ratio = NEW$Ratio
 }
  if(addremde == 0)
  {
    #there is no edge; we have to add a directed edge
    NEW = ADEDiffPar(CG,x,v,w)
    NewCG = NEW$CG
    Ratio = NEW$Ratio
  }
}

Out = list()
Out$OldCG = CG
  Out$NewCG <- NewCG
  Out$Ratio <- Ratio
  return(Out)
  }

