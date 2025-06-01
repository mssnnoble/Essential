#' Likelihood score for covariance factorised along the graph
#'
#' @param CG a chain graph object
#' @param n number of instantiations in an n by d data set
#' @param Covarmat - d by d empirical covariance matrix for a data set
#'
#' @return LIK log of likelihood score factorising along chain graph,
#' assuming normality and omitting the contribution from (2 pi) term
#'@export

LikScore <- function(CG,n,Covarmat)
{
d = ncol(Covarmat)



L <- vector(mode = "numeric", length = length(CG$CC))

# we'll compute the scores for each CC (chain component), conditioned on its
# parents. Resulting score will be the product of these.
# CC - find perfect order, C_1,....,C_m with seps S_1,...,S_{m-1} find
# cov of C_i \ S_{i-1}|S_{i-1} cup p(CC) (covariance of C_i without S_{i-1},
#  conditioned on separator S_{i-1} plus CC parents

for(j in 1:length(CG$CC))
  {
    C = MC(CG$CCued[[j]])
    J = kruskal(C)
    PO = perford(J)
    SET = CG$TOTPAR[[j]]

if(length(SET) >=1)
{
     a = det(as.matrix(Covarmat[SET,SET,drop=FALSE]))
    }

if(length(SET)==0)
{
  a = -1
  #this is just a marker to indicate that the set is empty
  }

    #note that R coerces things - a 1 by 1 matrix is no
    #longer a matrix (and has to be treated differently)

    K = vector(mode="numeric",length=nrow(J))

   for(k in 1:nrow(J))
  {
     vert = CG$CC[[j]][C[[PO[1,k]]]]
#vertices of MC labelled PO[1,k]

         if(PO[2,k]==0)
           {
           if(length(SET)==0)
           {
            K[k] = det(as.matrix(Covarmat[vert,vert,drop=FALSE]))
           }
     if(length(SET) >=1 & a==0)
     {
       INV = MASS::ginv(Covarmat[SET,SET,drop=FALSE])
     }

     if(length(SET)>=1 & a != 0)
     {
       INV = solve(Covarmat[SET,SET,drop=FALSE])
     }
      if(length(SET) >= 1)
      {
           Bmat = (Covarmat[vert,SET,drop=FALSE])%*%INV%*%(Covarmat[SET,vert,drop=FALSE])
      }
      if(length(SET)==0)
      {
        Bmat = matrix(0,length(vert),length(vert))
      }
           K[k] = det(Covarmat[vert,vert,drop=FALSE] - Bmat)
    }
    if(PO[2,k]!=0)
    {
      ns =  intersect(CG$CC[[j]][C[[PO[1,k]]]],CG$CC[[j]][C[[PO[2,k]]]])
      set1 = sort(c(ns,SET))
      set2 = sort(setdiff(CG$CC[[j]][C[[PO[1,k]]]],CG$CC[[j]][C[[PO[2,k]]]]))
       if(length(set1)>=1)
      {
        reqi = det(as.matrix(Covarmat[set1,set1,drop=FALSE]))
      }
      if(reqi==0)
      {
        A = Covarmat[set2,set1,drop=FALSE] %*% (MASS::ginv(Covarmat[set1,set1,drop=FALSE])) %*%
          Covarmat[set1,set2,drop=FALSE]
      }
      if(reqi != 0)
      {
        A = Covarmat[set2,set1,drop=FALSE] %*% (solve(Covarmat[set1,set1,drop=FALSE])) %*%
          Covarmat[set1,set2,drop=FALSE]
      }



      K[k] = det(Covarmat[set2,set2,drop=FALSE] - A)
    }
   }
    L[j] <- prod(K)
  }
  LIK <- -(n/2)*(log(prod(L)))
  return(LIK)
}
