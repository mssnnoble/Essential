#' Generates a Pruefer code for each CC
#'
#'takes the CCs and generates a pruefer code for each
#' returning a single entry of 0 if CC has one or two nodes.

#'
#' @param CC - list of CCs (we only use length of each CC)

#'
#'
#' @return pruf - list of pruefer codes, one for each CC, pruefer code is 0
#' if CC has length 1 or 2.
#'@export

createpruf <- function(CC){
  #

  pruf <- vector(mode="list",length=length(CC))
  pruf <-lapply(pruf,as.integer)

  for(i in 1:length(CC)){if(length(CC[[i]]) <= 2){pruf[[i]] = vector(mode="integer",length=0L)}
    else{pruf[[i]]=sample(1:length(CC[[i]]),length(CC[[i]])-2,replace = TRUE)}}
  return(pruf)
}
