#' Adjusts Layering to ensure suitable ground layer exists
#'
#'This algorithm adjusts the layering to ensure that there exists a layer
#'that is suitable for a ground layer (i.e. so that it is possible to generate
#'immoralities between layer 1 and layer 2 vertices)
#'
#' @param z collection of CCs
#' @param q their layers
#'
#'
#' @return a - the new layering
#' @export

#we adjust the layers to ensure there is a suitable ground layer
#i.e. ensure that at least one layer has possibility of not being complete


CCsetsADJ <- function(z,q){ #z the collection of CCs, q their layers
  #here we deal with the exceptional case where none of the layers are suitable
  #for a ground layer. We pick two layers and merge them. This layer will be used as
  #the ground layer in ordersets.
  a <- q
   if(max(sapply(z,length)) < 3 & max(sapply(a,length)) < 2)
  { #i.e. CCs have max 2, layers are all 1
    a[[length(a)-1]] <- c(a[[length(a)]],a[[length(a) - 1]])
    a[[length(a)]] <- NULL
  }
  return(a) #a is the new layering
}
