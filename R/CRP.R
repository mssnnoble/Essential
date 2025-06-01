#' Chinese Restaurant Process
#'
#'assigns customers to tables following a Chinese Restaurant Process
#'
#' @param conc weight of orange ball; probability of sitting at an empty
#' table for customer j is conc /(conc + (j-1))
#' @param size total number of customers
#'
#'
#' @return x = vector, length = size. x[j] gives table assignment for
#' customer j.
#'
#'We use this twice: firstly for CCs, secondly
#' for layers.
#'
#' Function is used to generate CCs and to assign CCs to raw layers
#' @export


CRP <- function(conc,size){

  x<-rep(1,size)
  if (size > 1){
    numTables <-1
    for (i in 2:size){
      y <- runif(1)
      if (y <= conc/(conc + i - 1)){
        numTables <- numTables + 1
        x[i] <- numTables
      }
      else {
        numbs = seq.int(1,i-1)
        index <- sample(x[numbs],1)
        x[i] = x[index]
      }
    }
  }
  return(x)
}
