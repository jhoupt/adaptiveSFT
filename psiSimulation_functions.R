

pm.function <- function (x,a,b,d) .5 * d + (1-d) * pnorm(x,a,b) # What Psi uses to estimate response probabilities 

Est.Trial.Psi.Color <- function (nTrials) {
  
  prior <- NA
  
  
  sim.a <- 6
  sim.b <- 15
  sim.d <- .01
  
  #COLOR
  x.range <- c(-55,50); x.step <- 1
  a.range <- c(-25,45); a.step <- 1
  b.range <- c(1,50); b.step <- 1
  d <- .01
  
  
  
  #### Set up internal objects ####
  
  x <- seq(x.range[1],x.range[2],x.step)
  a <- seq(a.range[1],a.range[2],a.step)
  b <- seq(b.range[1],b.range[2],b.step)
  r <- c(0,1)
  
  
  # Initialize P(R|L,X) look up table
  pR.LX <- array(dim = c(length(r), length(a), length(b), length(x)), dimnames = list("response"=r,"alpha"=a,"beta"=b,"intensity"=x))
  for (i in 1:length(r)) {
    for (j in 1:length(a)) {
      for (k in 1:length(b)) {
        for (l in 1:length(x)) {
          pR.LX[i,j,k,l] <- (1-r[i]) + (2*r[i]-1) * pm.function(x[l],a[j],b[k],d)
        }
      }
    }
  }
  
  ret.lamb <- list("alpha"=numeric(0), "beta"=numeric(0))
  
  # Initialize P(Lambda)
  if (is.na(prior)) {
    pL <- array(data = 1/(length(a)*length(b)), dim = c(1,length(a),length(b),1), dimnames = list("response","alpha"=a,"beta"=b,"intensity"))
  } else {
    pL <- prior
  }
  
  # Initialize P(R|X)
  pR.X <- array(dim = c(length(r), 1, 1, length(x)), dimnames = list("response"=r, "alpha", "beta", "intensity"=x))
  for (i in 1:length(r)) {
    for (j in 1:length(x)) {
      pR.X[i,1,1,j] <- sum(pR.LX[i,,,j] * pL[1,,,1])
    }
  }
  
  # Initialize P(L|X,R)
  pL.XR <- array(dim = c(length(r), length(a), length(b), length(x)), dimnames = list("response"=r,"alpha"=a,"beta"=b,"intensity"=x))
  for (i in 1:length(r)) {
    for (j in 1:length(a)) {
      for (k in 1:length(b)) {
        for (l in 1:length(x)) {
          pL.XR[i,j,k,l] <- pL[1,j,k,1] * pR.LX[i,j,k,l] / pR.X[i,,,l]
        }
      }
    }
  }
  
  # Initialize H(X,R)
  entropy.XR <- array(dim = c(length(r), 1, 1, length(x)), dimnames = list("response"=r, "alpha", "beta", "intensity"=x))
  for (i in 1:length(r)) {
    for (j in 1:length(x)) {
      entropy.XR[i,1,1,j] <- -sum(pL.XR[i,,,j] * log10(pL.XR[i,,,j]))
    }
  }
  
  # Initialize E[H(X,R)]
  expected.entropy.X <- array(dim = c(1, 1, 1, length(x)), dimnames = list("response", "alpha", "beta", "intensity"=x))
  for (i in 1:length(x)) {
    expected.entropy.X[1,1,1,i] <- sum(entropy.XR[,1,1,i] * pR.X[,1,1,i])
  }
  
  # Calculate First Trial Intensity
  next.intensity.index <- which.min(expected.entropy.X[1,1,1,])
  next.intensity <- x[next.intensity.index]
  
  
  
  #estimated parameters from human data using psychometric function
  #thres50 = estimated intensity value to obtain 50% accuracy
  thres50 = 6
  
  # Ratcliff parameters
  threshold=1.45
  v= 1.6 
  ter=.1
  sdv=.25
  
  ### Simulate the trials
  for (trial in 1:nTrials) {
    
    scaled.intensity <- (next.intensity-thres50)/(x.range[2]-thres50)
    DDMresult <- simdiffT(1,threshold,scaled.intensity*v,sdv,ter)
    sim.response <- c(DDMresult$x)
    
    
    ### Update Internal Arrays
    
    # Update P(L)
    pL <- pL.XR[sim.response+1, , , next.intensity.index]
    dim(pL) <- c(1, dim(pL)[1], dim(pL)[2], 1)
    dimnames(pL) <- list("response","alpha"=a,"beta"=b,"intensity")
    
    # Update P(R|X)
    for (i in 1:length(r)) {
      for (j in 1:length(x)) {
        pR.X[i,1,1,j] <- sum(pR.LX[i,,,j] * pL[1,,,1])
      }
    }
    
    # Update P(L|X,R)
    for (i in 1:length(r)) {
      for (j in 1:length(a)) {
        for (k in 1:length(b)) {
          for (l in 1:length(x)) {
            pL.XR[i,j,k,l] <- pL[1,j,k,1] * pR.LX[i,j,k,l] / pR.X[i,,,l]
          }
        }
      }
    }
    
    # Update H(X,R)
    for (i in 1:length(r)) {
      for (j in 1:length(x)) {
        entropy.XR[i,1,1,j] <- -sum(pL.XR[i,,,j] * log10(pL.XR[i,,,j]))
      }
    }
    
    # Update E[H(X,R)]
    for (i in 1:length(x)) {
      expected.entropy.X[1,1,1,i] <- sum(entropy.XR[,1,1,i] * pR.X[,1,1,i])
    }
    
    
    ### Calculate Next Trial Intensity
    
    next.intensity.index <- which.min(expected.entropy.X[1,1,1,])
    next.intensity <- x[next.intensity.index]
    
    ### Estimate Lambda
    
    a.est <- 0; b.est <- 0
    for (i in 1:length(a)) {
      a.est <- a.est + sum(a[i] * pL[1,i,,1])
    }
    for (j in 1:length(b)) {
      b.est <- b.est + sum(b[j] * pL[1,,j,1])
    }
    
    ret.lamb$alpha <- append(ret.lamb$alpha, a.est)
    ret.lamb$beta <- append(ret.lamb$beta, b.est)
  }
  
  return(ret.lamb)  
  
}


inv.pm.function <- function (y,a,b,d) qnorm((y-.5*d)/(1-d), a, b) # Use this to calculate thresholds with estimated lambda (alpha, beta) parameters

psi_color_ddm <- function(result.color, nDFP) { 
  intensity <- c()
  rescaled.intensity <- c()
  correct <- c()
  rt <- c()
  sim.d = .01
  x.range <- c(-55,50); x.step <- 1
  
  
  #Find high (99% accuracy) and low (90% accuracy) intensity values
  highSalience.color <- inv.pm.function(.99, result.color$alpha[length(result.color$alpha)], result.color$beta[length(result.color$alpha)], sim.d)
  lowSalience.color <- inv.pm.function(.90, result.color$alpha[length(result.color$alpha)], result.color$beta[length(result.color$alpha)], sim.d)
  intensity <- c(highSalience.color, lowSalience.color)
  
  #rescale intensities for DDM
  thres50.color=6
  highSalience.color <- (highSalience.color-thres50.color)/(x.range[2]-thres50.color)
  lowSalience.color <- (lowSalience.color-thres50.color)/(x.range[2]-thres50.color)
  rescaled.intensity = c(highSalience.color, lowSalience.color)
  
  N=2000
  threshold=1.45
  v=1.6
  ter=.1
  sdv=.25
  
  #Simulate DDM data with high & low salience values
  x.low.color <- simdiffT(nDFP,threshold,lowSalience.color*v,sdv,ter)
  x.high.color <- simdiffT(nDFP,threshold,highSalience.color*v,sdv,ter)
  
  
  correct <- list(x.high.color$x, x.low.color$x)
  rt <- list(x.high.color$rt, x.low.color$rt)
  
  return(list(intensity=intensity, rescaled.intensity=rescaled.intensity, rt=rt, correct=correct))
}





Est.Trial.Psi.Orientation <- function (nTrials) {
  
  prior <- NA
  
  sim.a <- 63 
  sim.b <- 6 
  sim.d <- .01
  
  x.range <- c(45,90); x.step <- .5
  a.range <- c(45.5,75); a.step <- .5
  b.range <- c(1,10); b.step <- .5
  d <- .01
  
  
  #### Set up internal objects ####
  x <- seq(x.range[1],x.range[2],x.step)
  a <- seq(a.range[1],a.range[2],a.step)
  b <- seq(b.range[1],b.range[2],b.step)
  r <- c(0,1)
  
  # Initialize P(R|L,X) look up table
  pR.LX <- array(dim = c(length(r), length(a), length(b), length(x)), dimnames = list("response"=r,"alpha"=a,"beta"=b,"intensity"=x))
  for (i in 1:length(r)) {
    for (j in 1:length(a)) {
      for (k in 1:length(b)) {
        for (l in 1:length(x)) {
          pR.LX[i,j,k,l] <- (1-r[i]) + (2*r[i]-1) * pm.function(x[l],a[j],b[k],d)
        }
      }
    }
  }
  
  ret.lamb <- list("alpha"=numeric(0), "beta"=numeric(0))
  
  # Initialize P(Lambda)
  if (is.na(prior)) {
    pL <- array(data = 1/(length(a)*length(b)), dim = c(1,length(a),length(b),1), dimnames = list("response","alpha"=a,"beta"=b,"intensity"))
  } else {
    pL <- prior
  }
  
  # Initialize P(R|X)
  pR.X <- array(dim = c(length(r), 1, 1, length(x)), dimnames = list("response"=r, "alpha", "beta", "intensity"=x))
  for (i in 1:length(r)) {
    for (j in 1:length(x)) {
      pR.X[i,1,1,j] <- sum(pR.LX[i,,,j] * pL[1,,,1])
    }
  }
  
  # Initialize P(L|X,R)
  pL.XR <- array(dim = c(length(r), length(a), length(b), length(x)), dimnames = list("response"=r,"alpha"=a,"beta"=b,"intensity"=x))
  for (i in 1:length(r)) {
    for (j in 1:length(a)) {
      for (k in 1:length(b)) {
        for (l in 1:length(x)) {
          pL.XR[i,j,k,l] <- pL[1,j,k,1] * pR.LX[i,j,k,l] / pR.X[i,,,l]
        }
      }
    }
  }
  
  # Initialize H(X,R)
  entropy.XR <- array(dim = c(length(r), 1, 1, length(x)), dimnames = list("response"=r, "alpha", "beta", "intensity"=x))
  for (i in 1:length(r)) {
    for (j in 1:length(x)) {
      entropy.XR[i,1,1,j] <- -sum(pL.XR[i,,,j] * log10(pL.XR[i,,,j]))
    }
  }
  
  # Initialize E[H(X,R)]
  expected.entropy.X <- array(dim = c(1, 1, 1, length(x)), dimnames = list("response", "alpha", "beta", "intensity"=x))
  for (i in 1:length(x)) {
    expected.entropy.X[1,1,1,i] <- sum(entropy.XR[,1,1,i] * pR.X[,1,1,i])
  }
  
  # Calculate First Trial Intensity
  next.intensity.index <- which.min(expected.entropy.X[1,1,1,])
  next.intensity <- x[next.intensity.index]
  
  
  
  #estimated parameters from human data using psychometric function
  #thres50 = estimated intensity value to obtain 50% accuracy
  thres50 = 63
  
  # Ratcliff parameters
  threshold=1.45
  v= 1.6 
  ter=.1
  sdv=.25
  
  ### Simulate the trials
  for (trial in 1:nTrials) {
    
    
    scaled.intensity <- (next.intensity-thres50)/(x.range[2]-thres50)
    DDMresult <- simdiffT(1,threshold,scaled.intensity*v,sdv,ter)
    sim.response <- c(DDMresult$x)
    
    
    ### Update Internal Arrays
    
    # Update P(L)
    pL <- pL.XR[sim.response+1, , , next.intensity.index]
    dim(pL) <- c(1, dim(pL)[1], dim(pL)[2], 1)
    dimnames(pL) <- list("response","alpha"=a,"beta"=b,"intensity")
    
    # Update P(R|X)
    for (i in 1:length(r)) {
      for (j in 1:length(x)) {
        pR.X[i,1,1,j] <- sum(pR.LX[i,,,j] * pL[1,,,1])
      }
    }
    
    # Update P(L|X,R)
    for (i in 1:length(r)) {
      for (j in 1:length(a)) {
        for (k in 1:length(b)) {
          for (l in 1:length(x)) {
            pL.XR[i,j,k,l] <- pL[1,j,k,1] * pR.LX[i,j,k,l] / pR.X[i,,,l]
          }
        }
      }
    }
    
    # Update H(X,R)
    for (i in 1:length(r)) {
      for (j in 1:length(x)) {
        entropy.XR[i,1,1,j] <- -sum(pL.XR[i,,,j] * log10(pL.XR[i,,,j]))
      }
    }
    
    # Update E[H(X,R)]
    for (i in 1:length(x)) {
      expected.entropy.X[1,1,1,i] <- sum(entropy.XR[,1,1,i] * pR.X[,1,1,i])
    }
    
    
    ### Calculate Next Trial Intensity
    
    next.intensity.index <- which.min(expected.entropy.X[1,1,1,])
    next.intensity <- x[next.intensity.index]
    
    ### Estimate Lambda
    
    a.est <- 0; b.est <- 0
    for (i in 1:length(a)) {
      a.est <- a.est + sum(a[i] * pL[1,i,,1])
    }
    for (j in 1:length(b)) {
      b.est <- b.est + sum(b[j] * pL[1,,j,1])
    }
    
    ret.lamb$alpha <- append(ret.lamb$alpha, a.est)
    ret.lamb$beta <- append(ret.lamb$beta, b.est)
  }
  
  return(ret.lamb)  
  
}




psi_orientation_ddm <- function(result.orientation, nDFP) { 
  intensity <- c()
  rescaled.intensity <- c()
  correct <- c()
  rt <- c()
  sim.d = .01
  x.range <- c(45,90); x.step <- 1
  
  
  #Find high (99% accuracy) and low (90% accuracy) intensity values
  highSalience.orientation <- inv.pm.function(.99, result.orientation$alpha[length(result.orientation$alpha)], result.orientation$beta[length(result.orientation$alpha)], sim.d)
  lowSalience.orientation <- inv.pm.function(.90, result.orientation$alpha[length(result.orientation$alpha)], result.orientation$beta[length(result.orientation$alpha)], sim.d)
  intensity <- c(highSalience.orientation, lowSalience.orientation)
  
  #rescale intensities for DDM
  thres50.orientation=63
  highSalience.orientation <- (highSalience.orientation-thres50.orientation)/(x.range[2]-thres50.orientation)
  lowSalience.orientation <- (lowSalience.orientation-thres50.orientation)/(x.range[2]-thres50.orientation)
  rescaled.intensity = c(highSalience.orientation, lowSalience.orientation)
  
  N=2000
  threshold=1.45
  v=1.6
  ter=.1
  sdv=.25
  
  #Simulate DDM data with high & low salience values
  x.low.orientation <- simdiffT(nDFP,threshold,lowSalience.orientation*v,sdv,ter)
  x.high.orientation <- simdiffT(nDFP,threshold,highSalience.orientation*v,sdv,ter)
  
  
  correct <- list(x.high.orientation$x,x.low.orientation$x)
  rt <- list(x.high.orientation$rt, x.low.orientation$rt)
  
  return(list(intensity=intensity, rescaled.intensity=rescaled.intensity, rt=rt, correct=correct))
}



rmse_fun <- function (params) {
  sqrt(mean((pm.function(all.intensities, params[1], params[2], sim.d) - DDM.pCorrect)^2))
}






