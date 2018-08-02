import numpy as np
import pystan
import os.path
import pickle
from scipy.stats import lognorm

#########################
###  Helper functions ###
#########################
def get_pr(intensity, target, range, posterior_samples, log=None) :
    x = ((posterior_samples['intensity'] * posterior_samples['alpha'] 
        + posterior_samples['intensity']^2 * posterior_samples['alpha2'])
        / posterior_samples['mu'])
    pr = np.mean(np.logical_and(x < target+range, x > target-range))

    if log is None or log is False : 
        return(pr)
    else :
        return(np.log(pr))


######
def lognormalrace_pdf(x, m, psi, mu, sigmasq, log=None) :
    sigma = np.sqrt(sigmasq)
    g = lognorm.logpdf(x-psi, loc=mu[m], s=sigma[m], scale=np.exp(mu[m]))
    G = 0
    for i in range(allchannels) :
        if i == m :
            continue
        else :
            G = G + lognorm.logsf(x-psi, loc=mu[i], s=sigma[i], 
                                  scale=np.exp(mu[i]))
    if log is None or log is False :
        rval = np.exp(g + G)
        rval[x<psi] = 0
    else : 
        rval = g + G
        rval[x < psi] = -np.inf
    return(rval)


#######
#plognormalrace <- function(x, m, psi, mu, sigmasq) { 
#  px <- rep(NA, length(x))
#  for(j in 1:length(x)) { 
#     px[j] <- lnrm_adjusted_integral(x[j], m, psi, mu, sigmasq, x[j+1]-x[j])
#  }
#  return(px)
#}


#######
#lnrm_adjusted_integral <- function(x, m, psi, mu, sigmasq, stepsize) { 
#  tryCatch({
#    f <- integrate(dlognormalrace, lower=0, upper=x, m=m, psi=psi,
#                    mu=mu, sigmasq=sigmasq)$value
#  }, error = function(e1) { 
#    tryCatch({
#      f <- integrate(dlognormalrace, lower=0, upper=x+stepsize/2, m=m,
#                      psi=psi, mu=mu, sigmasq=sigmasqx)$value
#    }, error = function(e2) { 
#      if (dlognormalrace(x-stepsize, m, psi, mu, sigmasq) == 0) { 
#        ff <- integrate(dlognormalrace, lower=0, upper=x+stepsize, m=2, 
#                        psi=psi, mu=mu, sigmasq=sigmasq)$value
#        ff <- ff/2
#      } else { 
#        ff <- NaN 
#      }
#    return(ff)
#    })
#  })
#  return(f)
#}
     

#######
#dfp_ddm <- function(N, drift.1, drift.2, a, ter, sdv, architecture, 
#                    stopping.rule, pmix=.5) {
## Function to generate rt and accuracy from DDM in DFP
#
#  if (architecture == "COA") { 
#    channel12 <- simdiffT(N,a,drift.1+drift.2,sdv,ter)
#    rt <- channel12$rt
#    cr <- channel12$x
#  } else { 
#    channel1 <- simdiffT(N,a,drift.1,sdv,ter)
#    channel2 <- simdiffT(N,a,drift.2,sdv,ter)
#    if (architecture == "PAR") { 
#      if (stopping.rule == "OR") { 
#        rt <- pmin(channel1$rt, channel2$rt)
#        cr <- channel2$x 
#        cr[channel1$rt < channel2$rt] <- 
#            channel1$x[channel1$rt < channel2$rt]
#      } else if (stopping.rule == "AND") { 
#        rt <- pmax(channel1$rt, channel2$rt)
#        cr <- channel1$x & channel2$x 
#      }
#    } else if (architecture == "SER") { 
#      if (stopping.rule == "OR") { 
#        channel.samp <- runif(N) < pmix
#        rt <- channel2$rt
#        rt[channel.samp] <- channel1$rt[channel.samp]
#        cr <- channel2$x 
#        cr[channel.samp] <- channel1$x[channel.samp]
#      } else if (stopping.rule == "AND") { 
#        rt <- channel1$rt + channel2$rt
#        cr <- channel1$x & channel2$x 
#      }
#    }
#  }
#  return(list(rt=rt, x=1*cr))
#}

######
#moc_ddm <- function(N, a, v, ter, sdv, intensity_levels) { 
## Function to generate method of constant stimuli data from DDM
#  intensity <- c()
#  correct <- c()
#  rt <- c()
#  for ( i in intensity_levels )  { 
#    x <- simdiffT(N,a,i*v,sdv,ter)
#    intensity <- c(intensity, rep(i, N))
#    correct <- c(correct, x$x)
#    rt <- c(rt, x$rt)
#  }
#  return(data.frame(intensity=intensity, rt=rt, correct=correct))
#}


######
#dataframe2stan <- function(dat) { 
## Reformat data for Stan
#   standat <- with(dat, list(N=dim(dat)[1], intensity=intensity, 
#                             correct=correct, minRT=min(rt), rt=rt) )
#   return(standat)
#}

#import pickle
#with open("temp_data.p", "rb") as f:
#  mydata = pickle.load(f)
#

#####
def find_salience(dat, h_targ, l_targ, fit_model = None):
# dat is a dictionary with:
#    'N':  total number of observations
#     intensity[]:  length N array-like containing stimulus intensity on 
#                   each trial
#     correct[]:  length N array-like containing indicator of correct
#                 on each trial
#     minRT:  smallest observed RT
#     rt[]:  length N array-like containing response time on each trial
#   
    from scipy.special import logit
    ML = False

    if fit_model is None :
        # Uncomment for quadratic
        #init_dict = {'alpha': -.1, 'alpha2': 0, 'mu': 1.5, 
        #             'psi': .1*dat['minRT'], 'varZ': 1}
        init_dict = {'slope': .1, 'midpoint': .5, 'mu': 1.5, 
                     'psi': .1*dat['minRT'], 'varZ': 1}
        if os.path.isfile("compiled_model.p"):
            with open("compiled_model.p", "rb") as f:
                sm = pickle.load(f)
        else: 
            sm = pystan.StanModel(file="lnrm2a.stan")
            #sm = pystan.StanModel(file="lnrm2.stan")
            with open("compiled_model.p", "wb") as f:
                pickle.dump(sm, f)
        if not ML : 
          fit_model = sm.sampling(data=dat, init=[init_dict, init_dict, 
                                                  init_dict, init_dict])

    if not ML: 
      post_diff = fit_model.extract(pars=["mu", "slope", "midpoint", "psi",
                                           "varZ"])
    else : 
      post_diff = sm.optimizing(data=dat, init=init_dict)

    slope = post_diff['slope']
    midpoint = post_diff['midpoint']


    l_targ_dist = logit(l_targ / 10.) / slope + midpoint
    h_targ_dist = logit(h_targ / 10.) / slope + midpoint
    
    rval = {}

    if not ML: 
      rval['high'] = np.nanmean(h_targ_dist)
      rval['high_var'] = np.var(h_targ_dist)
      rval['low'] = np.nanmean(l_targ_dist)
      rval['low_var'] = np.var(l_targ_dist)
      rval['fit'] = fit_model
    else : 
      rval['high'] = h_targ_dist
      rval['low'] = l_targ_dist
      rval['fit'] = post_diff

    return(rval)

