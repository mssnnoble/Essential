#' checks if a set of vertices is complete
#'
#' @param matr lower triangular matrix of 0s and 1s indicating undirected
#' edges of a CC
#' @param set a set of vertices within the CC (in CC coordinates)
#'
#' @return ret = 1 if set is a complete set and 0 if it is not a complete set
#'
#'@export

complete = function(matr,set){

  m = matr[set,set]
  l = length(set)
  ret = 0
  if(sum(matr) == 0.5*l*(l-1))
  {
    ret = 1
    }
  return(ret)
  }
