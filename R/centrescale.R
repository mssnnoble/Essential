#' Centres and scales a data matrix
#'
#' @param new: n by p data set each row is an instantiation of a p-variate
#' observtion.
#'
#' @return out - n by p matrix where out[i,j] obtained from new[i,j] by
#'subtracting column average and dividing through by
#' column standard deviation
#'
#'@export

centrescale <- function(new)
{
  new = as.matrix(new)
  me = colSums(new)/nrow(new)
  H = matrix(0,nrow(new),ncol(new))
  for(i in 1:nrow(H))
  {
    H[i,] = new[i,]-me
  }
  G = H^2
  b = sqrt(colSums(G)/(nrow(G)-1))

  out = matrix(0,nrow(new),ncol(new))
  for(i in 1:nrow(out))
  {
    out[i,]=H[i,]/b
  }
  return(out)
}
