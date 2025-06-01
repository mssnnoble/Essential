#' Computes layering of a graph when CCs and their parents are established
#'
#' Input is output in the style of an essential graph object,
#' where the  CCs, CCueds, TOTPAR are correct and  are the CCs, CCueds, TOTPARs
#' of an essential graph. This routine computes
#' the layering, the i2ip1 parents, the DEex parents, the VIL and prVIL to give
#' a full essential graph object
#'
#' @param CG - essential graph object where CCs, CCueds, TOTPARs are correct
#'
#'
#' @return L - the correct essential graph object
#' @export



NEWLAY <- function(CG){


  L <- CG

  #L will become the output

  whole <- seq(1,length(L$CC))
  LEN <- length(whole)
  f <- function(j){length(L$TOTPAR[[j]])}
  f1 <- sapply(whole,f)
  A <- which(f1==0)
  L$DEi2ip1 = vector(mode="list",length = LEN)
  L$DEi2ip1 = lapply(L$DEi2ip1,as.integer)
  L$DEex = vector(mode="list",length=0L)
  L$DEex = lapply(L$DEex,as.integer)

  #this gets the CCs with no parents - hence these
  #are in Layer 1.

  for(j in A)
  {
    L$DEi2ip1[[j]] <- vector(mode="integer",length=0L)
    L$DEex[[j]] <- vector(mode="integer",length=0L)
    L$TOTPAR[[j]] <- vector(mode="integer",length=0L)
  }
  L$ordlay[[1]] <- A
  #A gives labels of CCs in Layer 1

  L$VIL[[1]] <- vector(mode = "integer", length=0L)
   L$VIL[[1]] <- sort(unlist(L$CC[A]))
  L$prVIL[[1]] <- L$VIL[[1]]

  whole <- whole[-match(A,whole)]
  LEN <- length(whole)

  if(LEN == 0)
  {
    LVI <- vector(mode="list",length=1)
    LVI <-lapply(LVI,as.integer)
    LVI[[1]] = L$VIL[[1]]
    L$VIL = LVI
    L$prVIL = LVI
    LOR <- vector(mode="list",length=1)
    LOR <-lapply(LOR,as.integer)
    LOR[[1]] = A
    L$ordlay = LOR
    return(L)
  }

  if(LEN >=1)
  {
  counter <- 1L
  repeat{
    counter <- counter + 1L
    f <- function(j){1*(all(L$TOTPAR[[j]] %in% L$prVIL[[counter-1]]))}
    f1 = sapply(whole,f)
    A <- whole[which(f1==1)]
    hip1 <- function(j){intersect(L$TOTPAR[[j]],L$VIL[[counter-1]])}
    hdex <-function(j){setdiff(L$TOTPAR[[j]],L$VIL[[counter-1]])}
    for(j in A)
    {
      L$DEi2ip1[[j]] <- as.integer(hip1(j))
      L$DEex[[j]] <- as.integer(hdex(j))
    }
    L$ordlay[[counter]] <- A
    L$VIL[[counter]] = vector(mode="integer",length = 0L)
    L$VIL[[counter]] = sort(unlist(L$CC[A]))
    L$prVIL[[counter]]=sort(c(L$prVIL[[counter-1]],L$VIL[[counter]]))
    whole = whole[-match(A,whole)]
    LEN <- length(whole)
    if(LEN==0)
    {
      break
    }
  }
set = 1:counter
L$ordlay <- L$ordlay[set]
L$VIL <- L$VIL[set]
L$prVIL <- L$prVIL[set]


  return(L)
}
}
