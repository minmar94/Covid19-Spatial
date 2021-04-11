# Packages ----------------------------------------------------------------

library(MASS)
require(tidyverse)
require(magrittr)


# Functions - Richards -----------------------------------------------------

dRich <- function(ti, pars)
{
  r <- pars[1]
  h <- pars[2]
  p <- pars[3]
  s <- pars[4]
  
  out <- exp(log(r)+log( (1+exp(h*(p-ti)))^(-s) -  (1+exp(h*(p-ti+1)))^(-s)  ))
  
  return(out)
}

dRichLin <- function(ti, pars)
{
  r <- pars[1]
  h <- pars[2]
  p <- pars[3]
  s <- pars[4]
  
  out <- exp(log(r)+log(s)+log(h)+h*(p-ti)-(s+1)*log(1+exp(h*(p-ti))))
  
  return(out)
}

dRichLinDiff <- function(ti, pars)
{
  h <- pars[1]
  p <- pars[2]
  s <- pars[3]
  
  out <- exp(-h-(s+1)*(log(1+exp(h*(p-ti)))-log(1+exp(h*(p-ti+1)))))
  
  return(out)
}
