require(rstan)
require(diffIRT)
require(sft)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#########################
###  Helper functions ###
#########################
getPr_ogival <- function(intensity, target, range, postSamps,
                    plotDensity=FALSE, ...) { 
    x <- with(postSamps, L * inv_logit(slope * (intensity - midpoint)))
    if(plotDensity==TRUE) { 
      plot(density(x), xlab=expression(alpha), ...)
      #lines(target + range*c(-1, 1),c(0,0), lwd=3, col='blue')
      abline(v=target, lwd=3, col='blue')
    }
    return( mean( x < target+range &  x > target-range))
}

######
getPr <- function(intensity, target, range, postSamps, polynomial_order, 
                    plotDensity=FALSE, ...) { 
    if (polynomial_order ==1) { 
      #x <- with(postSamps, intensity * alpha)
      x <- with(postSamps, intensity * alpha / mu)
    } else if (polynomial_order ==2) { 
      #x <- with(postSamps, intensity * alpha + intensity^2 * alpha2)
      x <- with(postSamps, (intensity * alpha + intensity^2 * alpha2) / mu)
    }
    if(plotDensity==TRUE) { 
      plot(density(x), xlab=expression(alpha), ...)
      lines(target + range*c(-1, 1),c(0,0), lwd=3, col='blue')
    }
    return( mean( x < target+range &  x > target-range))
}

######
get_log_pr <- function(intensity, target, range, postSamps, 
                        polynomial_order) { 
    if (polynomial_order ==1) { 
      #x <- with(postSamps, intensity * alpha)
      x <- with(postSamps, intensity * alpha / mu)
    } else if (polynomial_order ==2) { 
      #x <- with(postSamps, intensity * alpha + intensity^2 * alpha2)
      x <- with(postSamps, (intensity * alpha + intensity^2 * alpha2) / mu)
    }
    rval <- log(mean(x < target+range &  x > target-range))
    return( rval ) 
}

######
get_inv_log_pr <- function(intensity, target, range, postSamps, 
                            polynomial_order) {
   rval <- -get_log_pr(intensity, target, range, postSamps, 
                          polynomial_order)
   return(rval)
}


######
dlognormalrace <- function(x, m, psi, mu, sigmasq) { 
  sigma <- sqrt(sigmasq)
  g <- dlnorm(x-psi, mu[m], sigma[m], log=TRUE)
  #g <- dlnorm(x-psi, mu[m], sigma[m])
  G <- 0
  allchannels <- 1:length(mu)
  for ( i in allchannels[-m] ) {
     G <- G + plnorm(x-psi, mu[i], sigma[i], lower=FALSE, log=TRUE)
     #G <- G + (1-plnorm(x-psi, mu[i], sigma[i]))
  }
  rval <- exp(g + G)
  rval[x<psi] <- 0
  if(any(!is.finite(rval))){ 
    xerr <- x[!is.finite(rval)]
    print(c(xerr, m, psi, mu, sigmasq)) }
  return(rval)
}


######
plognormalrace <- function(x, m, psi, mu, sigmasq) { 
  px <- rep(NA, length(x))
  for(j in 1:length(x)) { 
     px[j] <- lnrm_adjusted_integral(x[j], m, psi, mu, sigmasq, x[j+1]-x[j])
  }
  return(px)
}


######
lnrm_adjusted_integral <- function(x, m, psi, mu, sigmasq, stepsize) { 
  tryCatch({
    f <- integrate(dlognormalrace, lower=0, upper=x, m=m, psi=psi,
                    mu=mu, sigmasq=sigmasq)$value
  }, error = function(e1) { 
    tryCatch({
      f <- integrate(dlognormalrace, lower=0, upper=x+stepsize/2, m=m,
                      psi=psi, mu=mu, sigmasq=sigmasqx)$value
    }, error = function(e2) { 
      if (dlognormalrace(x-stepsize, m, psi, mu, sigmasq) == 0) { 
        ff <- integrate(dlognormalrace, lower=0, upper=x+stepsize, m=2, 
                        psi=psi, mu=mu, sigmasq=sigmasq)$value
        ff <- ff/2
      } else { 
        ff <- NaN 
      }
    return(ff)
    })
  })
  return(f)
}
     

######
dfp_ddm <- function(N, drift.1, drift.2, a, ter, sdv, architecture, 
                    stopping.rule, pmix=.5) {
# Function to generate rt and accuracy from DDM in DFP

  if (architecture == "COA") { 
    channel12 <- simdiffT(N,a,drift.1+drift.2,sdv,ter)
    rt <- channel12$rt
    cr <- channel12$x
  } else { 
    channel1 <- simdiffT(N,a,drift.1,sdv,ter)
    channel2 <- simdiffT(N,a,drift.2,sdv,ter)
    if (architecture == "PAR") { 
      if (stopping.rule == "OR") { 
        rt <- pmin(channel1$rt, channel2$rt)
        cr <- channel2$x 
        cr[channel1$rt < channel2$rt] <- 
            channel1$x[channel1$rt < channel2$rt]
      } else if (stopping.rule == "AND") { 
        rt <- pmax(channel1$rt, channel2$rt)
        cr <- channel1$x & channel2$x 
      }
    } else if (architecture == "SER") { 
      if (stopping.rule == "OR") { 
        channel.samp <- runif(N) < pmix
        rt <- channel2$rt
        rt[channel.samp] <- channel1$rt[channel.samp]
        cr <- channel2$x 
        cr[channel.samp] <- channel1$x[channel.samp]
      } else if (stopping.rule == "AND") { 
        rt <- channel1$rt + channel2$rt
        cr <- channel1$x & channel2$x 
      }
    }
  }
  return(list(rt=rt, x=1*cr))
}

######
moc_ddm <- function(N, a, v, ter, sdv, intensity_levels) { 
# Function to generate method of constant stimuli data from DDM
  intensity <- c()
  correct <- c()
  rt <- c()
  for ( i in intensity_levels )  { 
    x <- simdiffT(N,a,i*v,sdv,ter)
    intensity <- c(intensity, rep(i, N))
    correct <- c(correct, x$x)
    rt <- c(rt, x$rt)
  }
  return(data.frame(intensity=intensity, rt=rt, correct=correct))
}

######
dataframe2stan <- function(dat) { 
# Reformat data for Stan
   standat <- with(dat, list(N=dim(dat)[1], intensity=intensity, 
                             correct=correct, minRT=min(rt), rt=rt) )
   return(standat)
}


#####
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) 1/(1+exp(-x))

find_salience_ogival <- function(dat, h_targ, l_targ, fitModel=NA) { 

  standatDiff <- dataframe2stan(dat)

  if (anyNA(fitModel)) { 
    fitModel <- stan(file="lnrm2a.stan", data=standatDiff, 
                     pars=c("slope", "midpoint", "mu", "varZ", "psi"),
                     open_progress=TRUE)
  }
  post.diff <- extract(fitModel, c("slope", "midpoint", "mu", "varZ","psi"))

  slope = post.diff$slope
  midpoint = post.diff$midpoint

  l_targ.dist = logit(l_targ / 10.) / slope + midpoint
  h_targ.dist = logit(h_targ / 10.) / slope + midpoint

  return(list(high=mean(h_targ.dist, na.rm=TRUE), 
               low=mean(l_targ.dist, na.rm=TRUE), fit=fitModel))
}


########
find_salience_polynomial <- function(dat, h_targ, l_targ, polynomial_order=2, 
                          fitModel=NA) { 

  standatDiff <- dataframe2stan(dat)

  if(polynomial_order==1) { 
    if (anyNA(fitModel)) { 
      fitModel <- stan(file="lnrm1.stan", data=standatDiff, 
                       pars=c("mu", "alpha", "varZ", "psi"))
    }
    post.diff <- extract(fitModel, c("mu", "alpha", "psi","varZ"))
  } else if(polynomial_order==2) { 
    if (anyNA(fitModel)) { 
      fitModel <- stan(file="lnrm2.stan", data=standatDiff, 
                       pars=c("mu", "alpha", "alpha2", "varZ", "psi"))
    }
    post.diff <- extract(fitModel, c("mu","alpha", "alpha2", "psi", "varZ"))
  }
  

  # If order=2
  # a2 i^2 +  a1 i - 1/2 h_targ
  # i^2 +  a1/a2 i - 1/(2a2)h_targ
  # -a1/a2 + sqrt((a1/a2)^2+2/a2 * h_targ) / 2

  if (post.diff$alpha2 <0) {
  l_targ.dist <- with(post.diff,  
        (-alpha/alpha2 - sqrt( (alpha/alpha2)^2 + 2 / alpha2 * l_targ)) / 2)
  h_targ.dist <- with(post.diff,  
        (-alpha/alpha2 - sqrt( (alpha/alpha2)^2 + 2 / alpha2 * h_targ)) / 2)
  }
  #l_targ.dist <- with(post.diff,  
  #      (-alpha + sqrt( alpha^2 + 4 * alpha2 * mu * l_targ)) / (2*alpha2) )
  #h_targ.dist <- with(post.diff,  
  #      (-alpha + sqrt( alpha^2 + 4 * alpha2 * mu * h_targ)) / (2*alpha2) )


  #nstarts <- 100
  #starts <- seq(0,5,length.out=nstarts)
  #
  ## Low Salience
  #currentmin <- Inf
  #x <- c()
  #x$value <- Inf
  #for (s in starts) { 
  #  try( 
  #    x <- optim(s, get_inv_log_pr, target=l_targ, range=.1, 
  #               postSamps=post.diff, polynomial_order=polynomial_order),
  #               silent=TRUE
  #  )
  #  if (x$value < currentmin) {
  #    opt.diff <- x
  #    currentmin <- x$value
  #  }
  #}
  #lowSalience <- opt.diff$par


  ## High Salience
  #x <- c()
  #x$value <- Inf
  #currentmin <- Inf
  #for (s in starts) { 
  #  try( 
  #    x <- optim(s, get_inv_log_pr, target=h_targ, range=.1, 
  #    postSamps=post.diff, polynomial_order=polynomial_order),
  #    silent=TRUE
  #  )
  #  if (x$value < currentmin) {
  #    opt.diff <- x
  #    currentmin <- x$value
  #  }
  #}
  #highSalience <- opt.diff$par
  

  return(list(high=mean(h_targ.dist, na.rm=TRUE), 
               low=mean(l_targ.dist, na.rm=TRUE), fit=fitModel))
}

