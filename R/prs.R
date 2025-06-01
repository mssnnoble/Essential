#' Computes probability of table assignments according to Chinese Restaurant Process
#'
#' @param s: vector of lengths of the tables
#' @param theta: weight of 'orange ball' for Polya urn
#'
#' @return pr: probability of the table assignments according to
#' Ewens sampling formula with parameter theta.
#' @export


prs <- function(s,theta){
  n <- sum(s)
  a = rep(1,n)
  b = rep(1,n)
  f <- function(i){sum(s==i)}
  a <- sapply(1:n,f)
  g <- function(i){i*a[i]}
  b <- sapply(1:n,g)
  h <- function(i){(i/(theta + i - 1))}
  c <- sapply(1:n,h)
  dfn <- function(i){(theta^(a[i]))/((i^(a[i]))*(factorial(a[i])))}
  d <- sapply(1:n,dfn)



  pr = (prod(c))*(prod(d))
  return(pr)
}
