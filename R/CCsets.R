#' List of customers at each table for Chinese Restaurant Process
#'
#'This algorithm generates the list of customers at each table
#'of a Chinese Restaurant Process. This is used to generate the chain
#'components. It is also used to generate the layering for the chain
#'components
#'
#' @param CRP this is the list of table numbers for each vertex (customer)
#'
#'
#' @return
#' L - for  table i, L[[i]] gives the customers at the table.
#' @export



CCsets <- function(CRP){
  #we take the output of CRP (i.e. a list of table numbers)
  #and use this to generate a list of elements (vertices or CCs)
  #for each table (CC or layer)
  L <- vector(mode="list",length = max(CRP))
  L <- lapply(L,as.integer)
  f <- function(j){which(CRP == j)}
  L <- lapply(1:max(CRP),f)
  return(L)
}
