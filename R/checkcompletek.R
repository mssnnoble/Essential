#' Checks if a set  in layer k is a complete set.
#'
#' @param x - set to be checked
#' @param CC - chain components
#' @param ordlay - layering of the CCs
#' @param CCued - collection of undirected edge matrices for the CCs
#' @param k - layer.
#'
#' @return
#' \itemize{
#' \item a = c(1,k,i) if x is complete in layer k, CC labelled i
#' \item a = c(0,0,0) if x is not complete
#' }
#'@export

checkcompletek <- function(x,CC,ordlay,CCued,k){
  #check if x is a complete set in layer k.
  #return 0 if not complete
  # return 1 and also CC of which it is a subset if complete
  a = rep(0,3)
  for(i in 1:length(ordlay[[k]])){
    if(setequal(intersect(x,CC[[ordlay[[k]][i]]]),x) == TRUE)
      {
      z = match(x,CC[[ordlay[[k]][i]]])
      t <- length(z)
      y <-CCued[[ordlay[[k]][i]]][z,z]
      if(sum(y) == 0.5*t*(t-1))
      {
        a = c(1,k,i)
        }
      else{
        a=c(0,0,0)
        }
      }
    else{
      a=c(0,0,0)
      }
  }
  return(a)
  }
