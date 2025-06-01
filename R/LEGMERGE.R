#' Legal Merger: performs all possible legal mergers of a quasi-essential graph
#'
#' Start with a quasi-essential graph, perform all possible legal mergers
#' between CCs to obtain the corresponding essential graph
#'
#' @param CG a quasi-essential graph object
#'
#'
#' @return NG - the corresponding essential graph.
#' @export



LEGMERGE <- function(CG){

  NG = NEWLAY(CG)

  #j marks layers, i marks the CC within the layer.
  #take it a layer at a time (from layer 2 onwards)
  #a CC from layer j may be merged into a CC from layer
  #j-1.
  #During the i cycle, we first mark all CCs that can be
  #merged into the lower layer, merge them, then re-compute
  #the layering.
  #Since this has a 'knock-on' effect, we therefore repeat
  #the same layer (since a merger then re-computation can)
  #knock a layer j+1 CC down to layer j) until the layer
  #is stable, then move on to the next layer.

  j <- 2
  repeat
  {
  repeat
    {
      reliis <- vector(mode="integer",length=0L)

      #this denotes the indices of the CCs that should be
      #removed from layer j.

      gtot <- 0
      if(length(NG$ordlay) >=j)
      {
      for(i in 1:length(NG$ordlay[[j]]))
      {
        g <- function(k){1*(all(NG$DEi2ip1[[NG$ordlay[[j]][i]]] %in% NG$CC[[NG$ordlay[[j-1]][k]]]))}
        g1  = sapply(1:length(NG$ordlay[[j-1]]),g)

        #for layer j, CC i, g1[k]==1 if all i2ip1 parents are in CC k of layer j-1.
        #hence at most one of the g1's will be 1 - they satisfy one of the Studeny
        #rules for merging.

        if(sum(g1)==1)
        {
          #i.e all i2ip1 parents are in the same CC

          k = which(g1==1)

          #we use k to denote the label of this CC in layer j-1
          A = NG$DEi2ip1[[NG$ordlay[[j]][i]]]
          B <- match(A,NG$CC[[NG$ordlay[[j-1]][k]]])
          counter = NG$ordlay[[j-1]][k]
          M = as.matrix(NG$CCued[[counter]][B,B])
          ind = 0.5*length(A)*(length(A)-1)
          if(sum(M)==ind)
          {
            if(setequal(NG$DEex[[NG$ordlay[[j]][i]]],NG$TOTPAR[[NG$ordlay[[j-1]][k]]]) == TRUE)
            {
              w2 <- NG$ordlay[[j]][i]
              #index of CC in layer j that has to be removed from layer j
              gtot <- gtot + 1
              #total number to be removed from layer j
              reliis <- c(reliis,w2)
              #listing of the CCs that have to be removed from layer j
              w1 <- NG$ordlay[[j-1]][k]
              #this is the CC in the lower layer (j-1) with which w2 is merged
              newCC <- sort(c(NG$CC[[w1]],NG$CC[[w2]]))
              l1 = length(newCC)
              newCCued <- matrix(0,l1,l1)
              V1 = match(NG$CC[[w1]],newCC)
              V2 = match(NG$CC[[w2]],newCC)
              edge = match(NG$DEi2ip1[[w2]],newCC)
              newCCued[V1,V1] <- NG$CCued[[w1]]
              newCCued[V2,V2] <- NG$CCued[[w2]]
              newCCued[edge,V2] <- 1
              newCCued <- newCCued + t(newCCued)
              newCCued[upper.tri(newCCued,diag=TRUE)] <- 0
              #we have the merged CC and its undirected edges

              NG$CC[[w1]] <- newCC
              NG$CCued[[w1]] <- as.matrix(newCCued)
              NG$CC[[w2]] <- vector(mode="integer",length=0L)
              NG$CCued[[w2]] <- matrix(0,1,1)
              NG$TOTPAR[[w2]] <- vector(mode="integer",length=0L)
              NG$DEi2ip1[[w2]] <- vector(mode="integer",length=0L)
              NG$DEex[[w2]] <- vector(mode="integer",length=0L)
            }
          }
        }
      }
      }
      if(length(reliis)>=1)
      {
        NG$CC <- NG$CC[-reliis]
        NG$ordlay[[j]] <-NG$ordlay[[j]][-match(reliis,NG$ordlay[[j]])]
        NG$CCued <- NG$CCued[-reliis]
        NG$DEi2ip1 <- NG$DEi2ip1[-reliis]
        NG$DEex <- NG$DEex[-reliis]
        NG$TOTPAR <-NG$TOTPAR[-reliis]
        NG <- NEWLAY(NG)
      }
      if(gtot == 0)
      {
        break
      }
    }
    j <- j+1
    if(length(NG$ordlay)<=j-1)
    {
      break
    }

  }
  return(NG)
}
