#' Generates extra parents
#'
#' This routine generates the 'extra' directed edges - those from two or more
#' layers back for each CC - ensuring that the i2ip1 edges remain directed
#' (i.e. not adding edges that will violate Studeny conditions)
#'
#' If, for a CC, a configuration is chosen that violates Studeny conditions
#' then we instead keep 'extra edges' empty.
#'
#' Routine is called from essgraph (generating essential graph)
#'
#' @param i2ip1 - directed edges from previous layer (already chosen)
#' @param CC - list where CC[[k]] denotes vertices in the kth CC
#' @param CCued - list of matrices CCued[[k]] denotes undirected edges
#' for CC[[k]]
#' @param prVIL where prVIL[[k]] denotes the collection of all vertices in
#' layers 1, ... ,k
#' @param ordlay where ordlay [[j]] gives the labels of the CCs in layer j
#' @param p - the 'extra edge' probability
#'
#' @return M - list of vectors,
#' where M[[k]] gives the 2 or more layers back parents for CC[[k]]
#'
#'@export

DEextra <-function(i2ip1,CC,CCued,prVIL,ordlay,p){

  M <- vector(mode="list",length=length(CC))
  M <-lapply(M,as.integer)

  if(length(ordlay) <= 2)
  {
    return(M)
    }
  else{
    for(i in 3:length(ordlay))
      {
      a = rbinom(length(ordlay[[i]]),length(prVIL[[i-2]]),p)
      for(j in 1:length(ordlay[[i]]))
        {
        set = 1:length(prVIL[[i-2]])
        z1 = sort(sample(set,a[j],replace = FALSE))
        z = prVIL[[i-2]][z1]
        f <- checkvalid(z,i,j,i2ip1,M,CC,CCued,ordlay)
        if(f==1)
        {
          M[[ordlay[[i]][j]]] <- z
          }
      }
    }
    return(M)
  }
}
