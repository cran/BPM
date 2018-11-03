#' Sampling xi and alpha (tumor purity)
#' @param  y        A matrix, observed pure normal samples
#' @param  z        A matrix, observed mixed tumor samples
#' @param mstates   A matrix, hyper/hypo of dataset
#' @param xprior    A matrix, prior knowledge about purity
#' @param maxit     A number, maximum iteraction
#' @param burnin    A number, "burn-in" sample
#' @param xpar      Logistic, default is FALSE
#' @param n_ab0     initial value of n_ab
#' @param alp0      initial value of alpha
#' @param xbar0     initial value of xbar
#' @param trace     Logisitc, check the values in code, default is FALSE
#' @param verbose   Logistic, output the message,default is FALSE
#' @return x_bar    x_mode,
#' x_last   x2
#' x_sample x_sample
#' xpar     xprior2,
#' nab      n_ab2,
#' alp      alp2
#'@export
#'
#'
fullSampler <-function (y, z, mstates, xprior = NULL,
                         maxit = 1000, burnin = maxit,
                         xpar = FALSE, n_ab0 = NULL,
                         alp0 = NULL, xbar0 = NULL,
                         trace = FALSE, verbose = FALSE) {
  #---- sample  xi alpha from Top-K DMCs
  # y normal; z mixed/cancer; mstats hyper/hypo;
  m <- nrow(y) #loci
  n <- ncol(y) #sample Number
  if (is.null(xprior)) {
    #xprior <- matrix(c(10,1,1,10),2)
    xprior <- matrix(c(0.5,0.5,0.5,0.5),2)

  }

  hyper <- (mstates == 1)
  hypo <- (mstates == 2)

  # Initialization
  if (!is.null(n_ab0)) {
    n_ab <- n_ab0
  }else{
    n_ab <- 50
  }
  if (!is.null(alp0)) {
    alp <- alp0
  }else{
    alp <- stats::rbeta(n, 5, 5)
  }

  x_bar <- rep(NA, m)
  m1 <- rowMeans(y, na.rm = T)
  v1 <- apply(y, 1, stats::var, na.rm = T)
  y_bar <- m1 + (2*m1 - 1)*v1 / (m1*(1-m1) - 3*v1)
  #print(v1)
  #print(y_bar)
  y_bar[y_bar >= 1 | y_bar <= 0] <- m1[y_bar >= 1 | y_bar <= 0]

  if (is.null(xbar0)) {
    x_bar[hyper] <- .sampleXbar(y_bar[hyper], z[hyper,], alp, n_ab, xprior[1,])
    x_bar[hypo] <- .sampleXbar(y_bar[hypo], z[hypo,], alp, n_ab, xprior[2,])
  } else {
    x_bar <- xbar0
  }

  if (xpar) {
    xprior[1,] <- .sampleXprior(x_bar[hyper])
    xprior[2,] <- .sampleXprior(x_bar[hypo])
  }
  # Burn-in
  track <- list()
  if (trace) {
    if (xpar) {
      track$xprior <- array(0, dim = c(2,2,burnin))
    }
    track$nab <- rep(0, burnin)
  }


  for (it in 1 : burnin) {
    if (it %% 10 == 1 & verbose) {
      cat("Burn-in at iteration: ", it, "\n")
      if (xpar) {
        print(xprior)
      }
    }

    n_ab <- .sampleNab(x_bar, y_bar, z, alp, n_ab)
    for (k in 1 : n) {
      alp[k] <- .sampleAlp(x_bar, y_bar, z[,k], alp[k], n_ab)
    }

    x_bar[hyper] <- .sampleXbar(y_bar[hyper], z[hyper,], alp, n_ab, xprior[1,], x_bar[hyper])
    x_bar[hypo] <- .sampleXbar(y_bar[hypo], z[hypo,], alp, n_ab, xprior[2,], x_bar[hypo])

    if (xpar) {
      xprior[1,] <- .sampleXprior(x_bar[hyper], xprior[1,])
      xprior[2,] <- .sampleXprior(x_bar[hypo], xprior[2,])
    }
    if (trace) {
      if (xpar) {
        track$xprior[,,it] <- xprior
      }
      track$nab[it] <- n_ab
    }
  }

  it <- 1
  n_ab2 <- rep(0, maxit)
  alp2 <- matrix(0, maxit, n)
  x2 <- x_bar
  xprior2 <- array(0, dim = c(2,2,maxit))
  x_mean <- matrix(0, length(x_bar), maxit/10)
  x_mean[,1] <- x2
  nrec <- 100
  x_sample <- matrix(0, maxit/10, nrec)
  x_sample[1,] <- c(x_bar[1:nrec])

  n_ab2[it] <- .sampleNab(x_bar, y_bar, z, alp, n_ab)
  for (k in 1 : n) {
    alp2[it,k] <- .sampleAlp(x_bar, y_bar, z[,k], alp[k], n_ab2[it])
  }
  x2[hyper] <- .sampleXbar(y_bar[hyper], z[hyper,], alp2[it,], n_ab2[it], xprior[1,], x2[hyper])
  x2[hypo] <- .sampleXbar(y_bar[hypo], z[hypo,], alp2[it,], n_ab2[it], xprior[2,], x2[hypo])
  if (xpar) {
    xprior2[1,,it] <- .sampleXprior(x_bar[hyper], xprior[1,])
    xprior2[2,,it] <- .sampleXprior(x_bar[hypo], xprior[2,])
  }

  for (it in 2 : maxit) {
    if (it %% 10 == 1) {
      if (verbose)
      cat("Sampling at iteration: ", it, "\n")
      if (xpar) {
        print(xprior2[,,it-1])
      }
      x_mean[,(it+9)/10] <- x2
      x_sample[(it+9)/10,] <- c(x2[1:nrec])
    }

    n_ab2[it] <- .sampleNab(x2, y_bar, z, alp2[it-1,], n_ab2[it-1])
    for (k in 1 : n) {
      alp2[it,k] <- .sampleAlp(x2, y_bar, z[,k], alp2[it-1,k], n_ab2[it])
    }

    if (xpar) {
      x2[hyper] <- .sampleXbar(y_bar[hyper], z[hyper,], alp2[it,], n_ab2[it], xprior2[1,,it-1], x2[hyper])
      x2[hypo] <- .sampleXbar(y_bar[hypo], z[hypo,], alp2[it,], n_ab2[it], xprior2[2,,it-1], x2[hypo])
      xprior2[1,,it] <- .sampleXprior(x2[hyper], xprior2[1,,it-1])
      xprior2[2,,it] <- .sampleXprior(x2[hypo], xprior2[2,,it-1])
    } else {
      x2[hyper] <- .sampleXbar(y_bar[hyper], z[hyper,], alp2[it,], n_ab2[it], xprior[1,], x2[hyper])
      x2[hypo] <- .sampleXbar(y_bar[hypo], z[hypo,], alp2[it,], n_ab2[it], xprior[2,], x2[hypo])
    }
  }

  x_mode <- rep(NA, m)
  m1 <- rowMeans(x_mean, na.rm = T)
  v1 <- apply(x_mean, 1, stats::var, na.rm = T)
  x_mode <- m1 + (2*m1 - 1)*v1 / (m1*(1-m1) - 3*v1)
  x_mode[x_mode >= 1 | x_mode <= 0] <- m1[x_mode >= 1 | x_mode <= 0]

  if (trace) {
    return(list(x_bar = x_mode, x_last = x2, x_sample = x_sample, xpar = xprior2, nab = n_ab2, alp = alp2, track = track))
  } else {
    return(list(x_bar = x_mode, x_last = x2, x_sample = x_sample, xpar = xprior2, nab = n_ab2, alp = alp2))
  }
}

####============== functions used in fullSampler.R

###-------- sampleAlp
.sampleAlp <-function(x_bar, y_bar, z, alp0, n_ab) {
  ## sampling alpha using MH algorithm

  alp_new <- stats::rnorm(1, alp0, 0.005)
  while (alp_new < 0 | alp_new > 1) {
    alp_new <- stats::rnorm(1, alp0, 0.005)
  }

  a_new <- (x_bar*alp_new + y_bar*(1-alp_new))*n_ab
  b_new <- n_ab - a_new
  L_new <- sum(stats::dbeta(z, a_new+1, b_new+1, log = T), na.rm = T)

  a0 <- (x_bar*alp0 + y_bar*(1-alp0))*n_ab
  b0 <- n_ab - a0
  L0 <- sum(stats::dbeta(z, a0+1, b0+1, log = T), na.rm = T)

  acc <- exp(L_new - L0) #accept
  if (is.nan(acc)) {
    acc = 0.5
  }
  if (acc >= 1) {
    return(alp_new)
  } else {
    if (stats::runif(1) > acc) {
      return(alp0)
    } else {
      return(alp_new)
    }
  }
}

### ---------- sampleNab
.sampleNab <-function(x_bar, y_bar, z_all, alp_all, n_ab0) {
  # *_all means using the whole dataset
  # Assume uniform prior, log-normal proporsal
  # MH method
  if (is.null(n_ab0)) {
    return(n_ab0)
  } else {
    return(30)
  }
  return(30)
  n <- length(alp_all)
  m <- length(x_bar)

  propVar <- 0.01
  n_ab_new <- exp(stats::rnorm(1, log(n_ab0), propVar))

  L0 <- 0
  L_new <- 0
  a0 <- (x_bar %*% t(alp_all) + y_bar %*% t(1-alp_all))*n_ab0
  b0 <- n_ab0 - a0
  a_new <- (x_bar %*% t(alp_all) + y_bar %*% t(1-alp_all))*n_ab_new
  b_new <- n_ab_new - a_new

  L0 <- L0 + sum(stats::dbeta(z_all, a0+1, b0+1, log = T), na.rm = T)
  L_new <- L_new + sum(stats::dbeta(z_all, a_new+1, b_new+1, log = T), na.rm = T)

  acc <- exp(L_new - L0)
  if (is.nan(acc)) {
    acc = 0.5
  }

  if (acc >= 1) {
    return(n_ab_new)
  } else {
    if (stats::runif(1) > acc) {
      return(n_ab0)
    } else {
      return(n_ab_new)
    }
  }
}

### ---------sampleXbar
.sampleXbar <-function(y_bar, z, alp_all, n_ab, xprior = NULL, xbar0 = NULL) {
  # sample  xi (mode of beta value of each row in tumor)
  # using MH algorithm

  m <- length(y_bar)
  n <- length(alp_all)
  if (is.null(xprior)) {
    xprior <- c(1,1)
  }
  if (is.null(xbar0)) {
    return(stats::rbeta(m, xprior[1], xprior[2]))
  }

  x_bar <- stats::rnorm(m, xbar0, 0.005)
  temp <- which(x_bar < 0 | x_bar > 1)
  while(length(temp) > 0) {
    x_bar[temp] <- stats::rnorm(length(temp), xbar0[temp], 0.005)
    temp <- which(x_bar < 0 | x_bar > 1)
  }

  a0 <- (xbar0 %*% t(alp_all) + y_bar %*% t(1-alp_all))*n_ab
  b0 <- n_ab - a0

  a <- (x_bar %*% t(alp_all) + y_bar %*% t(1-alp_all))*n_ab
  b <- n_ab - a

  acc <- rowSums(stats::dbeta(z, a+1, b+1, log = T) - stats::dbeta(z, a0+1, b0+1, log = T), na.rm = T)
  acc <- exp(acc + stats::dbeta(x_bar, xprior[1], xprior[2], log = T) - stats::dbeta(xbar0, xprior[1], xprior[2], log = T))
  acc[is.nan(acc)] <- 0.5
  u <- stats::runif(m)
  x_bar[u > acc] <- xbar0[u > acc]

  return(x_bar)
}


### ---------- mv2ab
.mv2ab <-function(m1, v1) {
  ### matrix form
  a1 <- (m1*(1 - m1)/v1 - 1) * m1
  b1 <- (m1*(1 - m1)/v1 - 1) * (1 - m1)
  return(c(a1,b1))
}

###----------- sampleXprior  using mv2ab

.sampleXprior <-function (x_bar, xprior0 = NULL) {
  if (is.null(xprior0)) {
    m1 <- mean(x_bar, na.rm = T)
    v1 <- stats::var(as.vector(x_bar), na.rm = T)
    xprior <- .mv2ab(m1, v1)
    return(xprior)
  }

  # Assume xprior[1]/(xprior[1]+xprior[2]) is extreme
  xprior <- stats::rnorm(2, xprior0, 0.002)
  L_new <- sum(stats::dbeta(x_bar, xprior[1], xprior[2], log = T) +
                 stats::dbeta(xprior[1]/sum(xprior), 5, 5, log = T), na.rm = T)
  L0 <- sum(stats::dbeta(x_bar, xprior0[1], xprior0[2], log = T) +
              stats::dbeta(xprior0[1]/sum(xprior0), 5, 5, log = T), na.rm = T)
  acc <- exp(L_new - L0)

  if (is.nan(acc)) {
    acc = 0.5
  }

  if (acc >= 1) {
    return(xprior)
  } else {
    if (stats::runif(1) > acc) {
      return(xprior0)
    } else {
      return(xprior)
    }
  }
}

