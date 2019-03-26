source("adaptiveSFT_functions.R")
#########################################
### Generate data from Ratcliff Model ###
#########################################
PLOT <- FALSE
nDFP <- 1000

# Run parameter convergence tests?
runConvergence <- FALSE

# Run posterior predictive checks? 
posteriorPredictive <- FALSE
  # Fit separate LNRM models to each salience for comparison?
  fit.separate <- FALSE

# Simulate single clean DFP for each model? 
simulateSingleDFP <- FALSE

# Run full experiment
runFullExperiment <- FALSE


# Target DFP parameters in terms of LNRM
l_targ <- 1.3#.7
h_targ <- 8.0#1.4
L <- 10 # max separation 

# Times used for plotting
tvec <- seq(0, 5, length.out=1000)

# Ratcliff parameters
nPerLevel=100
nLevels=10
a=3
v=2
ter=.1
sdv=.20

#  MOC parameters
x.range <- c(45,90)
thres50 = 63

intensity_levels <- seq(x.range[1], x.range[2], length.out=nLevels) 
scaled.intensity <- (intensity_levels-thres50)/(x.range[2]-thres50)



# Simulate Data
datall <- moc_ddm(nPerLevel, a, v, ter, sdv, scaled.intensity)
psi <- rep(NA, nLevels)
psi.sem <- rep(NA, nLevels)
for (j in 1:nLevels) { 
  i <- scaled.intensity[j]
  psi[j] <- with(datall, mean(correct[intensity==i]))
  psi.sem[j] <- with(datall, sd(correct[intensity==i])) / sqrt(nPerLevel)
}

standatDiff <- dataframe2stan(datall)


# Fit LNRM
fitDiff <- stan(file="lnrm2a.stan", data=standatDiff, 
                     pars=c("mu", "midpoint", "slope", "varZ", "psi"))
post.diff <- extract(fitDiff, c("mu", "midpoint", "slope", "varZ", "psi"))

# Based on posterior distribution, determine the best intensity level
# to achieve target high and low salience values.
salience <- find_salience_ogival(datall, h_targ, l_targ, fitDiff)
highSalience <- salience$high
lowSalience <- salience$low

x.low <- simdiffT(nDFP,a,lowSalience*v,sdv,ter)
x.high <- simdiffT(nDFP,a,highSalience*v,sdv,ter)
standatHL <- list(N=2*nDFP,
                     intensity=rep(c(lowSalience,highSalience), each=nDFP),
                     correct=c(x.low$x, x.high$x),
                     minRT=min(x.low$rt, x.high$rt),
                     rt=c(x.low$rt, x.high$rt))




################################
####### Convergence Test #######
################################
Ns <- 1:300
if (runConvergence) { 
  nLevels <- 10
  intensity_levels <- seq(x.range[1], x.range[2], length.out=nLevels) 
  scaled.intensity <- (intensity_levels-thres50)/(x.range[2]-thres50)
  
  post.mean.midpoint <- matrix(NA, length(Ns))
  post.mean.slope <- matrix(NA, length(Ns))
  post.95.midpoint <- matrix(NA, length(Ns), 2)
  post.95.slope <- matrix(NA, length(Ns), 2)
  
  for (i in 1:length(Ns)) { 
    N <- Ns[i]
    datall <- moc_ddm(N, a, v, ter, sdv, scaled.intensity)
    standatDiff <- dataframe2stan(datall)
    fitDiff <- stan(fit=fitDiff, data=standatDiff, 
                    pars=c("mu", "midpoint", "slope", "varZ", "psi"))
    post.diff <- extract(fitDiff, c("midpoint", "slope"))
    post.mean.midpoint[i,] <- mean(post.diff$midpoint)
    post.95.midpoint[i,] <- quantile(post.diff$midpoint, c(0.05, 0.95))
    post.mean.slope[i,] <- mean(post.diff$slope)
    post.95.slope[i,] <- quantile(post.diff$slope, c(0.05, 0.95))
  }
  save(post.mean.midpoint, post.mean.slope, post.95.midpoint, post.95.slope, file="post95.Rdata") 
} else { 
  load("post95.Rdata")
}

setEPS()
postscript("LNRM_parameter-convergence.eps", width=7.4, height=3.2)
par(mfrow=c(1,2), mar=c(3.1, 3.6, 4.1, 1.1), oma=c(0,0,1,0))
matplot(1:300, post.95.midpoint, type='l', lty=2, col=1, ylab="",
        ylim=c(1, 1.6),
        #ylab="5/95% Posterior Quantiles", 
        main="Drift Difference Midpont", xlab="Trials per level")
mtext("5/95% Posterior Quantiles", 2, line=1.8)
lines(2:300, predict(loess(post.95.midpoint[2:300,1] ~ Ns[2:300], span=.5)), col=1)
lines(2:300, predict(loess(post.95.midpoint[2:300,2] ~ Ns[2:300], span=.5)), col=1)
#lines(1:300, post.mean.midpoint, col=1)

matplot(1:300, post.95.slope, type='l', lty=2, col=1, ylab="",
        #ylab="5/95% Posterior Quantiles", 
        main="Drift Difference Slope", xlab="Trials per level")
mtext("5/95% Posterior Quantiles", 2, line=1.8)
lines(2:300, predict(loess(post.95.slope[2:300,1] ~ Ns[2:300], span=.5)), col=1)
lines(2:300, predict(loess(post.95.slope[2:300,2] ~ Ns[2:300], span=.5)), col=1)
#lines(1:300, post.mean.slope, col=1)

title("LNRM Parameter Estimates", outer=TRUE, line=-1, cex=1.5)
dev.off()
################################
################################




############################
####### Fit Separate #######
############################
if (posteriorPredictive) { 
  ## Use this code to fit separate LNRM for each 
  ## intensity level.
  if (fit.separate) { 
    vmod <- vector("list", nLevels)
    drift_difference0 <- c()
    for ( ix in 1:nLevels)  { 
       i <- intensity_levels[ix]
       stanData0 <-  with(standatDiff, list(N=sum(intensity==i), 
                            correct=correct[intensity==i], 
                            minRT=min(rt[intensity==i]), 
                            rt=rt[intensity==i]))
       fitDiff0 <- stan(file="lnrm0.stan", data=stanData0)
       samps0 <- extract(fitDiff0, "mu", permute=TRUE)$mu
       vmod[[ix]] <- fitDiff0
       drift_difference0 <- c(drift_difference0, 
                              mean(samps0[,2] - samps0[,1]) / 2)
    }
    
    
    
    ## Plot individual LNRM fit versus LNRM-Rasch fit
    xx <- seq(0, max(intensity_levels)*1.1, length.out=1000)
    plot(intensity_levels, drift_difference0, xlab="Stimulus Intensity", 
          ylab="Difference between correct and incorrect mean in LNRM")
    if (polynomial_order==1) { 
       for (i in rsamp) { 
          lines(xx, xx * post.diff$alpha[i], col=grey(.8))
       }
    } else if (polynomial_order==2) { 
       for (i in rsamp) { 
          lines(xx, xx * post.diff$alpha[i] + xx^2 * 
                      post.diff$alpha2[i], col=grey(.8))
       }
    }
    legend("bottomright", c("Independent Fit", "Posterior LNRM-Rasch"), 
             pch=c(1, NA), lty=c(NA, 1))
    
  }
  
  
  ## Posterior predictive distribution
  if (PLOT) { 
  #dev.new()
  xx <- seq(0,5,length.out=100)
  png("posterior_predictive.png", 720*4, 360*4, res=72*4)
  par(mfrow=c(2,4))
  for (i in intensity_levels[1:8]) { 
     dat <- subset(datall, intensity==i)
     datecdf.cor <- ecdf(dat$rt[dat$correct==1])(tvec) * mean(dat$correct)
     if (sum(!dat$correct) > 0) { 
        datecdf.inc <- ecdf(dat$rt[!dat$correct])(tvec) * mean(!dat$correct)
     } else {
        datecdf.inc <- rep(0, length(tvec))
     }
     plot(tvec, datecdf.cor, type='l', main=paste("Intensity", round(i,2)),
              ylim=c(0,1), ylab="CDF", xlab="Time (s)")
     lines(tvec, datecdf.inc, col='red')
  
     sigmasqx <- c(1,1) * mean(post.diff$varZ)
  
     mux <- with(post.diff, t(
                 matrix(c(1,1), 2, 1) %*%mu[rsamp]  +  
                 matrix(c(-.5,.5)*L, 2, 1) %*%
                 inv_logit(slope[rsamp] * (i - midpoint[rsamp]))))
                 
     for (rn in 1:length(rsamp)) { 
        rs <- rsamp[rn]
        sigmasqx <- c(1,1) * post.diff$varZ[rs]
        px.cr <- plognormalrace(tvec, m=2, psi=post.diff$psi[rs], mu=mux[rn,],
                                  sigmasq=sigmasqx)
        px.in <- plognormalrace(tvec, m=1, psi=post.diff$psi[rs], mu=mux[rn,],
                                  sigmasq=sigmasqx)
  
        if(max(px.cr > 1)) { print(rs) }
        lines(tvec, px.cr, col=grey(rnorm(1, .8, .01)))
        lines(tvec, px.in, col=grey(rnorm(1, .8, .01)))
        #points(dat$rt, rnorm(length(dat$rt), 0, .0005), pch='.')
     }
     lines(tvec, datecdf.cor, col='green', lwd=3)
     lines(tvec, datecdf.inc, col='red', lwd=3)
  }
  dev.off()
  
  #dev.new()
  png("drift_difference.png", width=640, height=320)
  par(mfrow=c(1,2))
  getPr_ogival(lowSalience, target=l_targ, range=.01, postSamps=post.diff, 
           plotDensity=TRUE,
           ylab="Posterior Density", main="Low Salience")
  legend("topleft", c("Target Difference"), lwd=1, col='blue')
  getPr_ogival(highSalience, target=h_targ, range=.01, postSamps=post.diff, 
           plotDensity=TRUE,
           ylab="Posterior Density", main="High Salience")
  dev.off()
   
  
  mux.low <- with(post.diff, t(
                  matrix(c(1,1), 2, 1) %*%mu[rsamp]  +  
                  matrix(c(-.5,.5)*L, 2, 1) %*%
                  inv_logit(slope[rsamp] * 
                  (lowSalience - midpoint[rsamp]))))
  mux.high<- with(post.diff, t(
                  matrix(c(1,1), 2, 1) %*%mu[rsamp]  +  
                  matrix(c(-.5,.5)*L, 2, 1) %*%
                  inv_logit(slope[rsamp] * 
                  (highSalience - midpoint[rsamp]))))
  
  
  plot(tvec, ecdf(x.low$rt[x.low$x==1])(tvec) * mean(x.low$x==1), type='l', 
     ylim=c(0,1))
  for (rn in 1:length(rsamp)) { 
     rs <- rsamp[rn]
     sigmasqx <- c(1,1) * postOpt.diff$varZ[rs]
     px.cr <- plognormalrace(tvec, 2, postOpt.diff$psi[rs], mu=mux.low[rn,], 
                             sigmasq=sigmasqx)
     px.in <- plognormalrace(tvec, 1, postOpt.diff$psi[rs], mu=mux.low[rn,], 
                             sigmasq=sigmasqx)
     lines(tvec, px.cr, col=grey(rnorm(1, .8, .01)))
     lines(tvec, px.in, col=grey(rnorm(1, .8, .01)))
  
  }
  lines(tvec, ecdf(x.low$rt[x.low$x==1])(tvec) * mean(x.low$x==1), 
        col='green', lwd=3)
  if (sum(x.low$x==0)>0) { 
    lines(tvec, ecdf(x.low$rt[x.low$x==0])(tvec) * mean(x.low$x==0), 
          col='red', lwd=3)
  } else { 
    abline(h=0, col='red', lwd=3)
  }
  
  
  plot(tvec, ecdf(x.high$rt[x.high$x==1])(tvec) * mean(x.high$x==1), type='l',
     ylim=c(0,1))
  for (rn in 1:length(rsamp)) { 
     rs <- rsamp[rn]
     sigmasqx <- c(1,1) * postOpt.diff$varZ[rs]
     px.cr <- plognormalrace(tvec, 2, postOpt.diff$psi[rs], mu=mux.high[rn,], 
                             sigmasq=sigmasqx)
     px.in <- plognormalrace(tvec, 1, postOpt.diff$psi[rs], mu=mux.high[rn,], 
                             sigmasq=sigmasqx)
     lines(tvec, px.cr, col=grey(rnorm(1, .8, .01)))
     lines(tvec, px.in, col=grey(rnorm(1, .8, .01)))
  
  }
  lines(tvec, ecdf(x.high$rt[x.high$x==1])(tvec) * mean(x.high$x==1), 
        col='green', lwd=3)
  if (sum(x.high$x==0)>0) { 
    lines(tvec, ecdf(x.high$rt[x.high$x==0])(tvec) * mean(x.high$x==0), 
          col='red', lwd=3)
  } else { 
    abline(h=0, col='red', lwd=3)
  }
  
  }
}
############################
############################





####################
### Simulate DFP ###
####################
require(sft)
tvec.dense <- seq(0, 5, length.out=2000)

if (simulateSingleDFP) { 

# Coactive
HH <- dfp_ddm(250, highSalience*v, highSalience*v, a, ter, sdv, "COA", NA)
HL <- dfp_ddm(250, highSalience*v, lowSalience*v, a, ter, sdv, "COA", NA)
LH <- dfp_ddm(250, lowSalience*v, highSalience*v, a, ter, sdv, "COA", NA)
LL <- dfp_ddm(250, lowSalience*v, lowSalience*v, a, ter, sdv, "COA", NA)
p.or <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
            LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

#png("coactive.png", 640,320)
setEPS()
postscript("lnrm_coactive.eps", width=6.6, height=3.3)
par(mfrow=c(1,2), oma=c(0,0,.75,0), mar=c(3.6, 3.1, 3.1, 1.1))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="", ylab="")
mtext("Time (s)", side=1, line=1.9) 
mtext("S(t)", side=2, line=1.9)
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, p.or$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="", ylab="")
title("Coactive", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
mtext("Time (s)", side=1, line=1.9) 
mtext("SIC(t)", side=2, line=1.9)
dev.off()


# Parallel OR
HH <- dfp_ddm(250, highSalience*v, highSalience*v, a, ter, sdv, "PAR", "OR")
HL <- dfp_ddm(250, highSalience*v, lowSalience*v, a, ter, sdv, "PAR", "OR")
LH <- dfp_ddm(250, lowSalience*v, highSalience*v, a, ter, sdv, "PAR", "OR")
LL <- dfp_ddm(250, lowSalience*v, lowSalience*v, a, ter, sdv, "PAR", "OR")
p.or <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
            LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

#png("parallel_or.png", 640,320)
setEPS()
postscript("lnrm_parallel_or.eps", width=6.6, height=3.3)
par(mfrow=c(1,2), oma=c(0,0,.75,0), mar=c(3.6, 3.1, 3.1, 1.1))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="", ylab="")
mtext("Time (s)", side=1, line=1.9) 
mtext("S(t)", side=2, line=1.9)
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, p.or$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="", ylab="")
title("Parallel OR", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
mtext("Time (s)", side=1, line=1.9) 
mtext("SIC(t)", side=2, line=1.9)
dev.off()


# Parallel AND
HH <- dfp_ddm(250, highSalience*v, highSalience*v, a, ter, sdv, "PAR", "AND")
HL <- dfp_ddm(250, highSalience*v, lowSalience*v, a, ter, sdv, "PAR", "AND")
LH <- dfp_ddm(250, lowSalience*v, highSalience*v, a, ter, sdv, "PAR", "AND")
LL <- dfp_ddm(250, lowSalience*v, lowSalience*v, a, ter, sdv, "PAR", "AND")
p.and <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
             LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

#png("parallel_and.png", 640,320)
setEPS()
postscript("lnrm_parallel_and.eps", width=6.6, height=3.3)
par(mfrow=c(1,2), oma=c(0,0,.75,0), mar=c(3.6, 3.1, 3.1, 1.1))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="", ylab="")
mtext("Time (s)", side=1, line=1.9) 
mtext("S(t)", side=2, line=1.9)
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, p.and$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="", ylab="")
title("Parallel AND", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
mtext("Time (s)", side=1, line=1.9) 
mtext("SIC(t)", side=2, line=1.9)
dev.off()

# Serial OR
HH <- dfp_ddm(250, highSalience*v, highSalience*v, a, ter, sdv, "SER", "OR")
HL <- dfp_ddm(250, highSalience*v, lowSalience*v, a, ter, sdv, "SER", "OR")
LH <- dfp_ddm(250, lowSalience*v, highSalience*v, a, ter, sdv, "SER", "OR")
LL <- dfp_ddm(250, lowSalience*v, lowSalience*v, a, ter, sdv, "SER", "OR")
s.or <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
            LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

#png("serial_or.png", 640,320)
setEPS()
postscript("lnrm_serial_or.eps", width=6.6, height=3.3)
par(mfrow=c(1,2), oma=c(0,0,.75,0), mar=c(3.6, 3.1, 3.1, 1.1))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="", ylab="")
mtext("Time (s)", side=1, line=1.9) 
mtext("S(t)", side=2, line=1.9)
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, s.or$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="", ylab="")
title("Serial OR", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
mtext("Time (s)", side=1, line=1.9) 
mtext("SIC(t)", side=2, line=1.9)
dev.off()


# Serial AND
tvec.dense <- seq(0, 7, length.out=2000)
HH <- dfp_ddm(250, highSalience*v, highSalience*v, a, ter, sdv, "SER","AND")
HL <- dfp_ddm(250, highSalience*v, lowSalience*v, a, ter, sdv, "SER", "AND")
LH <- dfp_ddm(250, lowSalience*v, highSalience*v, a, ter, sdv, "SER", "AND")
LL <- dfp_ddm(250, lowSalience*v, lowSalience*v, a, ter, sdv, "SER", "AND")
s.and <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
             LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

#png("serial_and.png", 640,320)
postscript("lnrm_serial_and.eps", width=6.6, height=3.3)
par(mfrow=c(1,2), oma=c(0,0,.75,0), mar=c(3.6, 3.1, 3.1, 1.1))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="", ylab="")
mtext("Time (s)", side=1, line=1.9) 
mtext("S(t)", side=2, line=1.9)
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, s.and$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="", ylab="")
title("Serial AND", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
mtext("Time (s)", side=1, line=1.9) 
mtext("SIC(t)", side=2, line=1.9)
dev.off()
}

################################
### Simulate full experiment ###
################################
if (runFullExperiment) { 

  nLevels <- 10
  
  x.range <- c(45,90)
  thres50 = 63
  intensity_levels <- seq(x.range[1], x.range[2], length.out=nLevels) 
  orientation.intensity <- (intensity_levels-thres50)/(x.range[2]-thres50)
  
  x.range <- c(-55, 50)
  thres50 = 6
  intensity_levels <- seq(x.range[1], x.range[2], length.out=nLevels) 
  color.intensity <- (intensity_levels-thres50)/(x.range[2]-thres50)

  n.participants <- 10
  n.trials <- 100
  arch <- "PAR"
  srule <- "OR"
  sft.all <-c()
  sft.allx <-c()
  allfit.o <- vector("list", n.participants)
  allfit.c <- vector("list", n.participants)
  allpars <- array(NA, c(n.participants, 4))
  for (sn in 1:n.participants) { 
    a.p <- 0
    while(a.p <= 0) { 
      #a.p <- a
      a.p <- rnorm(1, mean=1.1*a, sd=a/6)
    }
    v.p <- 0
    while(v.p <= 0) { 
      v.p <- rnorm(1, mean=v/a*a.p, sd=v/6)
    }
    ter.p <- 0
    while(ter.p <= 0) { 
      ter.p <- rnorm(1, mean=ter, sd=ter/6)
    }
    sdv.p <- 0
    while(sdv.p <= 0) { 
      sdv.p <- sdv
      #sdv.p <- rnorm(1, mean=sdv, sd=sdv/10)
    }
    allpars[sn,] <- c(a.p, v.p, ter.p, sdv.p)
  
    #dat.p <- moc_ddm(N, a.p, v.p, ter.p, sdv.p, intensity_levels)
    dat.o <- moc_ddm(nPerLevel, a.p, v.p, ter.p, sdv.p, 
                     orientation.intensity)
    dat.c <- moc_ddm(nPerLevel, a.p, v.p, ter.p, sdv.p, color.intensity)
  
    allfit.o[[sn]] <- find_salience_ogival(dat.o, h_targ, l_targ)
    allfit.c[[sn]] <- find_salience_ogival(dat.c, h_targ, l_targ)
  
    high.o <- allfit.o[[sn]]$high
    low.o  <- allfit.o[[sn]]$low
  
    high.c <- allfit.c[[sn]]$high
    low.c  <- allfit.c[[sn]]$low
  
    HH <- dfp_ddm(n.trials,high.c*v.p, high.o*v.p, a.p, ter.p, sdv.p, arch, 
                  srule)
    HL <- dfp_ddm(n.trials,high.c*v.p,  low.o*v.p, a.p, ter.p, sdv.p, arch, 
                  srule)
    LH <- dfp_ddm(n.trials, low.c*v.p, high.o*v.p, a.p, ter.p, sdv.p, arch, 
                  srule)
    LL <- dfp_ddm(n.trials, low.c*v.p,  low.o*v.p, a.p, ter.p, sdv.p, arch, 
                  srule)
  
    sft.p <- data.frame(Subject=rep(sn, n.trials*4), 
                        Condition=rep(paste(arch,srule,sep="."), n.trials*4), 
                        Channel1=rep(c(2,1), each=2*n.trials),
                        Channel2=rep(rep(c(2,1), each=n.trials), 2),
                        Correct=c(HH$x, HL$x, LH$x, LL$x),
                        RT=c(HH$rt, HL$rt, LH$rt, LL$rt))
  
    sft.all <- rbind(sft.all, sft.p)
  
    HH.x <- dfp_ddm(n.trials, highSalience*v.p, highSalience*v.p, a.p, 
                    ter.p, sdv.p, arch, srule)
    HL.x <- dfp_ddm(n.trials, highSalience*v.p,  lowSalience*v.p, a.p,
                    ter.p, sdv.p, arch, srule)
    LH.x <- dfp_ddm(n.trials,  lowSalience*v.p, highSalience*v.p, a.p,
                    ter.p, sdv.p, arch, srule)
    LL.x <- dfp_ddm(n.trials,  lowSalience*v.p,  lowSalience*v.p, a.p,
                    ter.p, sdv.p, arch, srule)
  
    sft.x <- data.frame(Subject=rep(sn, n.trials*4), 
                        Condition=rep(paste(arch,srule,sep="."), n.trials*4), 
                        Channel1=rep(c(2,1), each=2*n.trials),
                        Channel2=rep(rep(c(2,1), each=n.trials), 2),
                        Correct=c(HH.x$x, HL$.xx, LH$.xx, LL.x$x),
                        RT=c(HH.x$rt, HL.x$rt, LH.x$rt, LL.x$rt))
  
    sft.allx <- rbind(sft.allx, sft.x)
  }
  colnames(allpars) <- c("a", "v", "Ter", "sdv")




}


printsft <- function() { 

  nLevels <- 10
  
  x.range <- c(45,90)
  thres50 = 63
  intensity_levels <- seq(x.range[1], x.range[2], length.out=nLevels) 
  orientation.intensity <- (intensity_levels-thres50)/(x.range[2]-thres50)
  
  x.range <- c(-55, 50)
  thres50 = 6
  intensity_levels <- seq(x.range[1], x.range[2], length.out=nLevels) 
  color.intensity <- (intensity_levels-thres50)/(x.range[2]-thres50)

  n.participants <- 20
  n.trials <- 100
  arch <- "PAR"
  srule <- "OR"
  sft.all <-c()
  sft.allx <-c()
  allfit.o <- vector("list", n.participants)
  allfit.c <- vector("list", n.participants)
  allpars <- array(NA, c(n.participants, 4))
  for (sn in 1:n.participants) { 
    a.p <- 0
    while(a.p <= 0) { 
      #a.p <- a
      a.p <- rnorm(1, mean=1.1*a, sd=a/8)
    }
    v.p <- 0
    while(v.p <= 0) { 
      v.p <- rnorm(1, mean=v/a*a.p, sd=v/8)
    }
    ter.p <- 0
    while(ter.p <= 0) { 
      ter.p <- rnorm(1, mean=ter, sd=ter/6)
    }
    sdv.p <- 0
    while(sdv.p <= 0) { 
      sdv.p <- sdv
      #sdv.p <- rnorm(1, mean=sdv, sd=sdv/10)
    }
    allpars[sn,] <- c(a.p, v.p, ter.p, sdv.p)
  
    #dat.p <- moc_ddm(N, a.p, v.p, ter.p, sdv.p, intensity_levels)
    dat.o <- moc_ddm(N, a.p, v.p, ter.p, sdv.p, orientation.intensity)
    dat.c <- moc_ddm(N, a.p, v.p, ter.p, sdv.p, color.intensity)
  
    allfit.o[[sn]] <- find_salience_ogival(dat.o, h_targ, l_targ)
    allfit.c[[sn]] <- find_salience_ogival(dat.c, h_targ, l_targ)
  
    high.o <- allfit.o[[sn]]$high
    low.o  <- allfit.o[[sn]]$low
  
    high.c <- allfit.c[[sn]]$high
    low.c  <- allfit.c[[sn]]$low
  
    HH <- dfp_ddm(n.trials,high.c*v.p, high.o*v.p, a.p, ter.p, sdv.p, arch, 
                  srule)
    HL <- dfp_ddm(n.trials,high.c*v.p,  low.o*v.p, a.p, ter.p, sdv.p, arch, 
                  srule)
    LH <- dfp_ddm(n.trials, low.c*v.p, high.o*v.p, a.p, ter.p, sdv.p, arch, 
                  srule)
    LL <- dfp_ddm(n.trials, low.c*v.p,  low.o*v.p, a.p, ter.p, sdv.p, arch, 
                  srule)
  
    sft.p <- data.frame(Subject=rep(sn, n.trials*4), 
                        Condition=rep(paste(arch,srule,sep="."), n.trials*4), 
                        Channel1=rep(c(2,1), each=2*n.trials),
                        Channel2=rep(rep(c(2,1), each=n.trials), 2),
                        Correct=c(HH$x, HL$x, LH$x, LL$x),
                        RT=c(HH$rt, HL$rt, LH$rt, LL$rt))
  
    sft.all <- rbind(sft.all, sft.p)
  
    HH.x <- dfp_ddm(n.trials, highSalience*v.p, highSalience*v.p, a.p, 
                    ter.p, sdv.p, arch, srule)
    HL.x <- dfp_ddm(n.trials, highSalience*v.p,  lowSalience*v.p, a.p,
                    ter.p, sdv.p, arch, srule)
    LH.x <- dfp_ddm(n.trials,  lowSalience*v.p, highSalience*v.p, a.p,
                    ter.p, sdv.p, arch, srule)
    LL.x <- dfp_ddm(n.trials,  lowSalience*v.p,  lowSalience*v.p, a.p,
                    ter.p, sdv.p, arch, srule)
  
    sft.x <- data.frame(Subject=rep(sn, n.trials*4), 
                        Condition=rep(paste(arch,srule,sep="."), n.trials*4), 
                        Channel1=rep(c(2,1), each=2*n.trials),
                        Channel2=rep(rep(c(2,1), each=n.trials), 2),
                        Correct=c(HH.x$x, HL$.xx, LH$.xx, LL.x$x),
                        RT=c(HH.x$rt, HL.x$rt, LH.x$rt, LL.x$rt))
  
    sft.allx <- rbind(sft.allx, sft.x)
  }
  colnames(allpars) <- c("a", "v", "Ter", "sdv")
}
