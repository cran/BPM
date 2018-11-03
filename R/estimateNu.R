#' Estimate noise intensity (nv) for non-DMCs,
#' using maximum likelihood estmiation.
#'
#' @param  z       A matrix. Observated mixed turmor samples.
#' @param  phi     mode of beta-values of each row in pure nomral samples y.
#' @param  maxit   A postive integer. The iteration number used in maximum likelihood.
#' @param  beginP  A number, where the method start to search from for root.
#'
#' return  estimated nv (noise intensity)
#'
#'

estimateNu<-function(z, phi, maxit = 50,beginP=20) {
  m = nrow(z) # row number
  n= ncol(z) #  column number # cliff
  c1 <- 0
  for (i in 1 : m) {
    c1 <- c1 + sum(phi[i] * log(z[i,]/(1-z[i,])),na.rm = TRUE)
  }
  c1 <- c1 + sum(log(1-z),na.rm=TRUE)
  c1 <- c1 / n

  ### solve c1=logLikFun(phim,nu) using bisection method
  nu <- beginP #start point
  nu_high <- nu
  nu_low <- nu
  for (tt in 1 : maxit)
    { temp = .logLikFun(phi, nu)

    if (temp == c1)
      {break }
    else if (temp < c1)
      { if (nu_high == nu)
        { nu_low <- nu
          nu <- nu * 2
          nu_high <- nu
        }
        else
        { nu_low <- nu
          nu <- (nu + nu_high) / 2
        }
      }
    else
      { if (nu_low == nu)
        { nu_high <- nu
          nu <- nu / 2
          nu_low <- nu
        }
        else
        { nu_high <- nu
          nu <- (nu + nu_low) / 2
        }
      }
   }
  nu
}


### ---- logLikFun  (3.13)/or formula (6) in supplementary S2.
.logLikFun <-function (phi, nu) {
  K <- length(phi)
  ll <- sum(phi * digamma(phi*nu+1),na.rm=TRUE)
  ll <- ll + sum((1-phi) * digamma((1-phi)*nu+1),na.rm=TRUE)
  ll <- ll - K * digamma(nu + 2)
  ll
}

