#' Checks if a candidate set of extra parents is valid
#'
#'
#' @param z is the set we want to check (returned from DEextra)
#' @param i is the layer
#' @param j number within layer
#' @param i2ip1 set of previous layer parents (obtained for all CCs)
#' @param M candidate 'extra' parents
#' @param CC list of all CCs
#' @param CCued list of all adjacency matrices for the CCs
#' @param ordlay layering of the CCs
#'
#' @return a which is 1 if valid, 0 if not valid
#' @export


checkvalid <- function(z,i,j,i2ip1,M,CC,CCued,ordlay){



  #Recall the conditions for admissible parents: need
  #1. p_1(\tau) is contained in a particular \sigma
  #2. p_1(\tau) \cap \sigma is complete
  #3. p_2(\tau) = p_1(\sigma) union p_2(\sigma)

  a=1
  c=0

  #c will be the index of the CC in the lower layer

  d=0
  e=0

  #z candidate set for p_2(\tau)

  f <- function(j){sort(c(i2ip1[[j]],M[[j]]))}
  TOTPAR = lapply(1:length(CC),f)

  uf <- function(b){
    1*(isTRUE(setequal(z,TOTPAR[[ordlay[[i-1]][b]]])))
    }

  #clique labelled CC[[ordlay[[i-1]][b]]] is candidate for \sigma we get union of p_1 and p_2
  #of \sigma and check if they are equal (condition 3.)

  ug <- sapply(1:length(ordlay[[i-1]]),uf)
  c = which(ug==1)

  if(length(c)==0)
    {
    a <- 1
  return(a)
    }

  #if one of the b's satisfies we pull it out and call it c, otherwise the set z is admissible

  else{
    if(length(c)==1)
  {
      label = c[1]
    e = 1
  }

      #e = 1 means that 3. is satisfied and also 1. is satisfied

    if(e ==1)
      {
      k = match(i2ip1[[ordlay[[i]][j]]],CC[[ordlay[[i-1]][label]]])
      y <- CCued[[ordlay[[i-1]][label]]]
      if(complete(y,k)==1)
      {
        a = 0
        }
      }
  }
  return(a)
}
