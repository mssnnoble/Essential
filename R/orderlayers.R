#' orders layers at random ensuring a suitable ground layer is chosen
#'
#' The function orderlayers orders the layers at random, ensuring that a suitable
#' ground layer is chosen (i.e one that can admit immoralities for layer 2 CCs).
#' This means taking the specially created layer if this was necessary,
#' otherwise choosing a ground layer at random from among
#' suitable layers (i.e. not a layer which is a single CC which is either a
#' singleton or double). Other layers selected at random from the unselected layers
#'
#' @param z - this is the list of CCs, output of CCsets
#' @param q - adjusted, unordered layering of the CCs (adjusted to ensure
#' there is a suitable layer for ground layer)
#'
#'
#' @return b - the ordered layers, b[[1]] the ground layer. b[[j]] lists labels
#' of CCs in layer j.
#' @export



orderlayers = function(z,q){

   m=0
  ones <- rep(1,length(q))
  ones2 <- rep(1,length(q))
  b <- q
  #
  #eventually ones will give a 1 if the layer is suitable for ground layer and
  #0 otherwise
  #
  if(length(b)==1){
    return(b)
  }
  else{
    for (i in 1:length(q)){
      if (length(q[[i]]) == 1 & length(z[[q[[i]][1]]]) <=2 )
        {
        ones[i]=0}
      }
    #
    #length(q[i])==1 means that (unordered) layer i contains exactly one cc
    #length(z[[q[[i]][1]]]) gives number of vertices in the only CC in this layer
    #ones[i]=0 means that this layer won't be considered for ground layer
    #

    possible = which(ones==1)

    samp <- sample(1:length(possible),1)
    b[[1]]<- q[[possible[samp]]]
    #q[[possible[samp]]] is the ground layer
    rrem <- q[-possible[samp]]
    c = sample(seq.int(1,(length(q)-1)),(length(q)-1))
    for(i in 2:length(q)){
      b[[i]] = rrem[[c[i-1]]]}
    return(b)
     }
}
