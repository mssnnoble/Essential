#'Gives a path on the junction between x1 and CC2
#'
#'@param x1 - a vertex
#'@param CC2 - a set of vertices
#' @param C collection of maximal cliques
#' @param J junction tree
#'
#' @return L - path from MC containing x1 to MC which intersects with CC2
#'
#'
#'@export


MCPath <- function(x1,CC2,C,J)
  {

  #J the junction tree (encoded as a matrix)
  #  C is the collection of maximal cliques

 Linit <- function(x1,CC2,C,J)
  {
  L <- vector(mode="list",length=length(C))
  L <-lapply(L,as.integer)


  #L[[j]] path from an MC containing y[1] to MC C[j]

  f1 <- function(j){1*(x1 %in% C[[j]])}
  f2 <- function(j){1*(length(intersect(CC2, C[[j]])) >=1)}

  g1 <- sapply(1:length(C),f1)
  #MCs containing x1
  g2 <- sapply(1:length(C),f2)
#MCs that intersect with CC2
  y1MC <- which(g1==1)

  #MCs containing x1

  y2MC <- which(g2 == 1)
  #MCs which intersect with CC2

  SEQ <- seq(1,length(C))
  L[[y1MC[1]]] <- y1MC[1]
  USED <- y1MC[1]

    #MCs that have already been dealt with
  UNUSED <- SEQ[-match(y1MC[1],SEQ)]
  #MCs where path hasn't been computed

  #initialisation - choose an MC containing x1

  SYMM <- J + t(J)

  N <- which(SYMM[y1MC[1],]==1)
  #MCs which are neighbours in Junction Tree of y1MC[1]
  Ll = length(N)
  #Ll number of MCs that are neighbours.

  for(i in 1:Ll)
  {
    L[[N[i]]]= c(L[[y1MC[1]]],N[i])
  }
  USED <- c(USED,N)
  UNUSED <- UNUSED[-match(N,UNUSED)]
  SYMM[L[[y1MC[1]]],N] <- 0
  SYMM[N,L[[y1MC[1]]]] <- 0
  #remove edges that have been transversed

  repeat
  {
    f = function(j){1*(SYMM[USED[j],]==1)}
    f1 = lapply(1:length(USED),f)
    #f1[[j]][i] == 1 if corresponding MC is one step away from USED[j]
    #and 0 otherwise.
    #(row of 0's if j is a leaf)
    f = function(j){which(f1[[j]]==1)}
    ONESTEP = sapply(1:length(f1),f)

        #ONESTEP[[i]] lists MCs one step away from USED[i]

    g = function(j){intersect(ONESTEP[[j]],y2MC)}
    g1 = lapply(1:length(f1),g)

    #g1[[j]] gives MCs containing y2 which are one step away from USED[j]

    h = function(j){sum(g1[[j]])}
    sss = sapply(1:length(g1),h)
    summit = sum(sss)
    if(summit >=1)
    {
      a = which(sss >=1)[1]
      #so g1[[a]][1] is one step away from a used MC
      os = g1[[a]][1]
      b = intersect(which(SYMM[g1[[a]][1],]==1), USED)[1]
      L[[os]] = c(L[[b]],os)
      return(L[[os]])
    }
    else
    {
      for(i in 1:length(f1))
      {
        if(length(ONESTEP[[i]])>=1)
        {
          N <- sort(unique(c(N,ONESTEP[[i]])))
          for(j in 1:length(ONESTEP[[i]]))
          {
            a = intersect(which(SYMM[ONESTEP[[i]][j],]==1),USED[i])
            L[[ONESTEP[[i]][j]]] = c(L[[USED[i]]],ONESTEP[[i]][j])
          }
        }
      }
    }
    SYMM[USED,N] = 0
    SYMM[N,USED]=0
    USED <- sort(unique(c(USED,N)))
    UNUSED <- SEQ[-match(N,SEQ)]
  }
  if(length(UNUSED) == 0)
  {
    return(L)
    }
}

L = Linit(x1,CC2,C,J)
f = function(j){1*(x1 %in% C[[L[j]]])}
f1 = sapply(1:length(C),f)
lowind = max(which(f1==1))

g = function(j){1*(length(intersect(CC2,C[[L[j]]]))>=1)}
g1 = sapply(1:length(C),g)
highind = min(which(g1==1))

path = L[lowind:highind]

return(path)

}
