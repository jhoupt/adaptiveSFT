
setwd("C:/Users/w018elf/Google Drive/Publications/Adaptive SFT/Model")
source("adaptiveSFT_functions.R")
source("psiSimulation_functions.R")
require(diffIRT)

#############################PSI WITH RATCLIFF DDM SIMULATIONS FOR ORIENTATION#############################
########################################REPORT RESULTS IN MANUSCRIPT ######################################
###########################################################################################################

x.range <- c(45,90) #Range of orientations
axis.x <- seq(x.range[1], length.out=1000, x.range[2]) #Resolution of psychometric estimates

sim.d <- .01 #Lapse parameter
nsamps = 629 #want to run 629 simulations
trials = 300 #want to run 300 trials of psi for each simulated model

result1.a <- matrix(data=NA, nrow=nsamps, ncol=trials) #matrix with all alpha parameter estimates
result1.b <- matrix(data=NA, nrow=nsamps, ncol=trials) #matrix with all beta parameter estimates

#list of matrices with all psychometric functions 
pm.fun1 <- list() #Requires ~2.2 Gb of storage -- but can be regenerated anytime using only alpha & beta parameters

###########MODEL SIMULATIONS --> Runtime: ~10MINS FOR EACH SIMULATION (i.e., 10*nsamps)##################
for (sn in 1:nsamps){
  result = Est.Trial.Psi.Orientation(trials)
  result1.a[sn,] = result$alpha 
  result1.b[sn,] = result$beta
  pm.fun = matrix(data=NA, nrow=trials, ncol=length(axis.x))
  for (i in 1:trials){
    pm.fun[i,] = mapply(function(a,b) pm.function(axis.x, a, b, sim.d), a=result1.a[sn,i], b=result1.b[sn,i])
  }
  if(sn==1){pm.fun1 = list(pm.fun)}
  if(sn>1){pm.fun1 <- append(pm.fun1, list(pm.fun))}
}


###############SAVE ALPHA & BETA PARAMETERS TO CSV FILES###########################
#######CAN USE THESE TO READ-IN & REGENERATE PSYCHOMETRIC FUNCTIONS & PLOTS #######
write.csv(result1.a, file="orientation_alpha.csv", row.names=FALSE)
write.csv(result1.b, file="orientation_beta.csv", row.names=FALSE)


###USE COMMENTED CODE BELOW IF:
###ALPHA & BETA PARAMETERS ARE SAVED IN CSV & YOU WANT TO READ-IN & REGENERATE PSYCHOMETRIC FUNCTION DATA
###~2.2 GB OF STORAGE REQUIRED

# pm.fun1 <- list() 
# result1.a <- read.csv("orientation_alpha.csv")
# result1.b <- read.csv("orientation_beta.csv")
# for (sn in 1:nsamps){
#   pm.fun = matrix(data=NA, nrow=trials, ncol=length(axis.x))
#   for (i in 1:trials){
#     pm.fun[i,] = mapply(function(a,b) pm.function(axis.x, a, b, sim.d), a=result1.a[sn,i], b=result1.b[sn,i])
#   }
#   if(sn==1){pm.fun1 = list(pm.fun)}
#   if(sn>1){pm.fun1 <- append(pm.fun1, list(pm.fun))}
# }


#############CALCULATE THE MEAN PSYCHOMETRIC FUNCTION (length(x.axis)) FOR EACH TRIAL RUN (trials) ACROSS ALL MODEL SIMULATIONS (nsamps)##########
mean.pm.fun <- matrix(0, length(trials), length(x.axis))
for (ix in 1:length(nsamps)) {
  mean.pm.fun <- mean.pm.fun + pm.fun1[[ix]]
}
mean.pm.fun <- mean.pm.fun / length(nsamps)

#Plot estimated psychometric function across trials -- 
#Each is an average function across all model runs (nsamps)
#Red (Trial 1) --> Blue (Trial N)
windows()
plot(0,0,type='n', xlim=range(axis.x), ylim=c(0,1), xlab="Intensity", ylab="Resp Rate", main =paste("Psi: Psychometric Functions"))
for (i in 1:trials){
  lines(axis.x, mean.pm.fun[i,], col=rainbow(length(mean.pm.fun[,i]), start=0,end=4/6)[i])
}

#Estimate and plot the true DDM psychometric function that was fed into Psi
thres50 = 63 #50% threshold intensity value estimated from participant data
all.intensities <- seq(x.range[1], x.range[2], length.out=1000)
scaled.intensity <- (all.intensities-thres50)/(x.range[2]-thres50)



#DDM parameters: Fit to closely approximate parameter values of human psychometric data with the orientation task
threshold=1.45
v= 1.6 
ter=.1
sdv=.25

#Calculate P(Correct) for the DDM response on each psi trial
DDM.pCorrect = exp(threshold*scaled.intensity*v) / (exp(-1*threshold*scaled.intensity*v) + exp(threshold*scaled.intensity*v))

#Find the optimal alpha & beta parameters to fit the psychometric function of the DDM P(Correct) with each intensity value
ss = c(60, 6) #intensity search space for finding optimal parameters
fitted.params <- optim(ss, rmse_fun)

#Add line to plot TRUE psychometric function of the DDM model
lines(all.intensities, pm.function(all.intensities, fitted.params$par[1], fitted.params$par[2], sim.d), col='Black', lwd=3, lty=1)

####PLOT ALPHA CONVERGENCE ACROSS TRIALS
#Average (SOLID LINE), 5%-95% (DASHED LINES), and true alpha value (RED DASHED LINE) CI functions across all model runs (nsamps)
windows()
plot(apply(result1.a, 2, mean, na.rm=TRUE), lwd=2, type='l', main=paste("Psi: Alpha Convergence"), ylim = c(min(apply(result1.a, 2, quantile, probs=.05, na.rm=TRUE))-1, max(apply(result1.a, 2, quantile, probs=.95, na.rm=TRUE))+1))
lines(apply(result1.a, 2, quantile, probs=.05, na.rm=TRUE), lty=2)
lines(apply(result1.a, 2, quantile, probs=.95, na.rm=TRUE), lty=2)
abline(h=fitted.params$par[1], lty=2, col='red')

####PLOT BETA CONVERGENCE ACROSS TRIALS
#Average (SOLID LINE), 5%-95% (DASHED LINES), and true beta value (RED DASHED LINE) CI functions across all model runs (nsamps)
windows()
plot(apply(result1.b, 2, mean, na.rm=TRUE), lwd=2, type='l', main=paste("Psi: Beta Convergence"), ylim = c(min(apply(result1.b, 2, quantile, probs=.05, na.rm=TRUE))-1, max(apply(result1.b, 2, quantile, probs=.95, na.rm=TRUE))+1))
lines(apply(result1.b, 2, quantile, probs=.05, na.rm=TRUE), lty=2)
lines(apply(result1.b, 2, quantile, probs=.95, na.rm=TRUE), lty=2)
abline(h=fitted.params$par[2], lty=2, col='red')



##########################################PSI SIMULATIONS FOR COLOR########################################
#########################################NOT INCLUDED IN MANUSCRIPT########################################
###########################################################################################################
result2.a <- matrix(data=NA, nrow=nsamps, ncol=trials)
result2.b <- matrix(data=NA, nrow=nsamps, ncol=trials)
pm.fun2 <- list()

for (sn in 1:nsamps){
  result = Est.Trial.Psi.Color(trials)
  result2.a[sn,] = result$alpha
  result2.b[sn,] = result$beta
  pm.fun = matrix(data=NA, nrow=nsamps, ncol=length(x))
  for (i in 1:trials){
    pm.fun[i,] = mapply(function(a,b) pm.function(x, a, b, sim.d), a=result2.a[sn,i], b=result2.b[sn,i])
  }
  if(sn==1){pm.fun2 = list(pm.fun)}
  else if(sn>1){pm.fun2 <- append(pm.fun2, list(pm.fun))}
}

###########################################################################################################
###########################################################################################################
###########################################################################################################





#########################################
### Generate data from Ratcliff Model ###
#########################################
PLOT <- FALSE
nDFP <- 1000

# Fit separate LNRM models to each salience for comparison?
fit.separate <- FALSE

# Target DFP parameters in terms of LNRM
#l_targ <- 1.7#.7
#h_targ <- 3.5#1.4
#L <- 10 # max separation 

# Times used for plotting
tvec <- seq(0, 5, length.out=1000)
d <- .01
trials = 300

###############COLOR###############

#UNCOMMENT TO RESET HIGH/LOW VALUES FROM SIMULATION  
#NOTE: TAKES ~10 minutes
#result.color = Est.Trial.Psi.Color(trials) 
#highSalience.color <- inv.pm.function(.99, result.color$alpha[length(result.color$alpha)], result.color$beta[length(result.color$alpha)], sim.d)
#lowSalience.color <- inv.pm.function(.90, result.color$alpha[length(result.color$alpha)], result.color$beta[length(result.color$alpha)], sim.d)
x.range <- c(-55,50) #Color scale
thres50.color=6 #50% intensity estimate from human data

highSalience.color = 101.63970 #COMMENT THIS LINE IF YOU RERAN A SIMULATION ABOVE
highSalience.color <- (highSalience.color-thres50.color)/(x.range[2]-thres50.color)

lowSalience.color = 53.64586 #COMMENT THIS LINE IF YOU RERAN A SIMULATION ABOVE
lowSalience.color <- (lowSalience.color-thres50.color)/(x.range[2]-thres50.color)


x.low.color <- simdiffT(nDFP,threshold,lowSalience.color*v,sdv,ter)
x.high.color <- simdiffT(nDFP,threshold,highSalience.color*v,sdv,ter)


#########ORIENTATION###############

#UNCOMMENT TO RESET HIGH/LOW VALUES FROM SIMULATION  
#NOTE: TAKES ~10 minutes
#result.orientation = Est.Trial.Psi.Orientation(trials)
#highSalience.orientation <- inv.pm.function(.99, result.orientation$alpha[length(result.orientation$alpha)], result.orientation$beta[length(result.orientation$alpha)], sim.d)
#lowSalience.orientation <- inv.pm.function(.90, result.orientation$alpha[length(result.orientation$alpha)], result.orientation$beta[length(result.orientation$alpha)], sim.d)
x.range <- c(45,90) #orientation scale
thres50.orientation = 63 #50% intensity estimate from human data

highSalience.orientation = 87.11268 #COMMENT THIS LINE IF YOU RERAN A SIMULATION ABOVE
highSalience.orientation <- (highSalience.orientation-thres50.orientation)/(x.range[2]-thres50.orientation)

lowSalience.orientation = 75.06533 #COMMENT THIS LINE IF YOU RERAN A SIMULATION ABOVE
lowSalience.orientation <- (lowSalience.orientation-thres50.orientation)/(x.range[2]-thres50.orientation)


x.low.orientation<- simdiffT(nDFP,threshold,lowSalience.orientation*v,sdv,ter)
x.high.orientation<- simdiffT(nDFP,threshold,highSalience.orientation*v,sdv,ter)


###################################################################################################
############### SIMULATE A DDM MODEL WITH EACH ARCHITECTURE & STOPPING-RULE MODEL TYPE#############
##################PLOT THE SURVIVOR FUNCTIONS & SIC FOR EACH DFP MODEL PREDICTION##################
###################################################################################################
require(sft)
tvec.dense <- seq(0, 5, length.out=2000)

# Parallel OR
HH <- dfp_ddm(250, highSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "PAR", "OR")
HL <- dfp_ddm(250, highSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "PAR", "OR")
LH <- dfp_ddm(250, lowSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "PAR", "OR")
LL <- dfp_ddm(250, lowSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "PAR", "OR")
p.or <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
            LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

png("parallel_or.png", 640,320)
par(mfrow=c(1,2), oma=c(0,0,1,0))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="Time (s)", ylab="S(t)", xlim=c(0,max(LL$rt)))
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, p.or$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="Time (s)", ylab="SIC", xlim=c(0,max(LL$rt)))
title("Parallel OR", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
dev.off()


# Parallel AND
HH <- dfp_ddm(250, highSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "PAR", "AND")
HL <- dfp_ddm(250, highSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "PAR", "AND")
LH <- dfp_ddm(250, lowSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "PAR", "AND")
LL <- dfp_ddm(250, lowSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "PAR", "AND")
p.and <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
             LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

png("parallel_and.png", 640,320)
par(mfrow=c(1,2), oma=c(0,0,1,0))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="Time (s)", ylab="S(t)", xlim=c(0,max(LL$rt)))
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, p.and$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="Time (s)", ylab="SIC", xlim=c(0,max(LL$rt)))
title("Parallel AND", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
dev.off()

# Serial OR
HH <- dfp_ddm(250, highSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "SER", "OR")
HL <- dfp_ddm(250, highSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "SER", "OR")
LH <- dfp_ddm(250, lowSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "SER", "OR")
LL <- dfp_ddm(250, lowSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "SER", "OR")
s.or <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
            LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

png("serial_or.png", 640,320)
par(mfrow=c(1,2), oma=c(0,0,1,0))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="Time (s)", ylab="S(t)", xlim=c(0,max(LL$rt)))
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, s.or$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="Time (s)", ylab="SIC", xlim=c(0,max(LL$rt)))
title("Serial OR", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
dev.off()


# Serial AND
tvec.dense <- seq(0, 7, length.out=2000)
HH <- dfp_ddm(250, highSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "SER","AND")
HL <- dfp_ddm(250, highSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "SER", "AND")
LH <- dfp_ddm(250, lowSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "SER", "AND")
LL <- dfp_ddm(250, lowSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "SER", "AND")
s.and <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
             LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

png("serial_and.png", 640,320)
par(mfrow=c(1,2), oma=c(0,0,1,0))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="Time (s)", ylab="S(t)", xlim=c(0,max(LL$rt)))
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, s.and$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="Time (s)", ylab="SIC", xlim=c(0,max(LL$rt)))
title("Serial AND", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
dev.off()


###################################################################################################
############SIMULATE A FULL EXPERIMENT WITH 100 DDM SUBJECTS WITH PARALLEL-OR PROCESSES############
###################################################################################################
#FOR EACH SUBJECT: 
# 1) RUN DDM THROUGH PSI COLOR TASK
# 2) ESTIMATE COLOR INTENSITIES FOR HIGH & LOW SALIENCE 
# 3) RUN DDM THROUGH PSI ORIENTATION TASK 
# 4) ESTIMATE ORIENTATION INTENSITIES FOR HIGH & LOW SALIENCE 
# 5) RUN DDM THROUGH FULL DFP EXPERIMENT USING COLOR (CHANNEL 1) & ORIENTATION (CHANNEL 2) STIMULI
# 6) CALCULATE SIC  
# 7) SAVE PLOT OF SURVIVOR FUNCTIONS & SIC
###################################################################################################

nParticipants <- 100 #number of simulated participants
nTrials <- 200 #number of psi trials ---- NOTE: based on convergence of alpha/beta parameters in convergence study 


result.color <- list() #stores all psi color results
ddm.result.color <- list() #stores color ddm performance @ high/low values

result.orientation <- list() #stores all psi orientation results
ddm.result.orientation <- list() #stores orientation ddm performance @ high/low values

#H = high, L = low
#Channel 1 = Color , Channel 2 = Orientation
#List that stores each simulated DFP results with factorial design 
#HH = Channel 1: High salience & Channel 2: High salience
#LL = Channel 1: Low salience & Channel 2: Low salience
HH <- list()  
HL <- list()
LH <- list()
LL <- list()

#List that stores each simulated participants SIC results
p.or <- list()

for (sn in 1:nParticipants) { 
  nDFP = 1000 #number of DFP trials
  
 # result.color[[sn]] = Est.Trial.Psi.Color(nTrials) #ddm runs in psi color experiment
  ddm.result.color[[sn]] = psi_color_ddm(result.color[[sn]], nDFP) #ddm performance in 1000 trials with high & low salience color trials
   
  
 # result.orientation[[sn]] = Est.Trial.Psi.Orientation(nTrials) #ddm runs in psi orientation experiment
  ddm.result.orientation[[sn]] = psi_orientation_ddm(result.orientation[[sn]], nDFP)#ddm performance in 1000 trials with high & low salience orientation trials
  
  require(sft)
  tvec.dense <- seq(0, 5, length.out=2000)
  
  #Use same DDM paramaters from simulated convergence testing
  threshold=1.45
  v=1.6
  ter=.1
  sdv=.25
  
  #Simulate DDM performance for each trial type of DFP (assumes DDM has Parallel OR processes)
  HH[[sn]] <- dfp_ddm(250, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "PAR", "OR")
  HL[[sn]] <- dfp_ddm(250, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "PAR", "OR")
  LH[[sn]] <- dfp_ddm(250, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "PAR", "OR")
  LL[[sn]] <- dfp_ddm(250, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "PAR", "OR")
  
  #Estimates SIC from DDM DFP performance data
  p.or[[sn]] <- sic(HH=HH[[sn]]$rt[HH[[sn]]$x==1], HL=HL[[sn]]$rt[HL[[sn]]$x==1], 
                    LH=LH[[sn]]$rt[LH[[sn]]$x==1], LL=LL[[sn]]$rt[LL[[sn]]$x==1])
  
  #Plots & Saves the DDM survivor functions for each DFP trial type (HH, HL, LH, LL) & SIC Function
  png(paste("parallel_or Subject ", sn, ".png"), 640,320)
  par(mfrow=c(1,2), oma=c(0,0,1,0))
  plot(tvec.dense, 1-ecdf(HH[[sn]]$rt[HH[[sn]]$x==1])(tvec.dense), type='l', col='red',
       xlab="Time (s)", ylab="S(t)", xlim=c(0,max(LL[[sn]]$rt)))
  lines(tvec.dense, 1-ecdf(HL[[sn]]$rt[HL[[sn]]$x==1])(tvec.dense), col='orange')
  lines(tvec.dense, 1-ecdf(LH[[sn]]$rt[LH[[sn]]$x==1])(tvec.dense), col='purple')
  lines(tvec.dense, 1-ecdf(LL[[sn]]$rt[LL[[sn]]$x==1])(tvec.dense), col='blue')
  legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
         col=c("red", "orange", "purple", "blue"))
  
  plot(tvec.dense, p.or[[sn]]$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
       xlab="Time (s)", ylab="SIC", xlim=c(0,max(LL[[sn]]$rt)))
  title(paste("Parallel OR - Subject ", sn), outer=TRUE, line=-2, cex=1.5)
  abline(h=0)
  dev.off()
  
}

PsiParams = cbind("Subject" = rep(1, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[1]]$alpha, "Color_Beta" = result.color[[1]]$beta, "Orientation_Alpha" = result.orientation[[1]]$alpha, "Orientation_Beta" = result.orientation[[1]]$beta)
dominance = cbind("Subject"=rep(1, 8), p.or[[1]]$Dominance)
SFTresults = cbind("Subject"= 1, "H_Color Intensity" = ddm.result.color[[1]]$intensity[1],"L_Color Intensity" = ddm.result.color[[1]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[2],  "D+ Statistic" = p.or[[1]]$SICtest$positive$statistic, "D+ Pvalue" = p.or[[1]]$SICtest$positive$p.value,  "D- Statistic" = p.or[[1]]$SICtest$negative$statistic, "D- Pvalue" = p.or[[1]]$SICtest$negative$p.value, "MIC Statistic" = p.or[[1]]$MICtest$statistic, "MIC Pvalue" = p.or[[1]]$MICtest$p.value)
DDM.DFPperformance = cbind(cbind("Subject" = rep(1, 250), "Trial" = 1:250), "HH_RT" = HH[[1]]$rt[1:250], "HH_Correct" = HH[[1]]$x[1:250], "HL_RT" = HL[[1]]$rt[1:250], "HL_Correct" = HL[[1]]$x[1:250],  "LH_RT" = LH[[1]]$rt[1:250], "LH_Correct" = LH[[1]]$x[1:250],"LL_RT" = LL[[1]]$rt[1:250], "LL_Correct" = LL[[1]]$x[1:250])

for (sn in 2:nParticipants) { 
  PsiParams = rbind(PsiParams, cbind("Subject" = rep(sn, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[sn]]$alpha, "Color_Beta" = result.color[[sn]]$beta, "Orientation_Alpha" = result.orientation[[sn]]$alpha, "Orientation_Beta" = result.orientation[[sn]]$beta))
  SFTresults = rbind(SFTresults, cbind("Subject"= sn, "H_Color Intensity" = ddm.result.color[[sn]]$intensity[1],"L_Color Intensity" = ddm.result.color[[sn]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[2],  "D+ Statistic" = p.or[[sn]]$SICtest$positive$statistic, "D+ Pvalue" = p.or[[sn]]$SICtest$positive$p.value,  "D- Statistic" = p.or[[sn]]$SICtest$negative$statistic, "D- Pvalue" = p.or[[sn]]$SICtest$negative$p.value, "MIC Statistic" = p.or[[sn]]$MICtest$statistic, "MIC Pvalue" = p.or[[sn]]$MICtest$p.value))
  dominance = rbind(dominance, cbind("Subject"=rep(sn, 8), p.or[[sn]]$Dominance))
  DDM.DFPperformance = rbind(DDM.DFPperformance, cbind("Subject" = rep(sn, 250), "Trial" = 1:250, "HH_RT" = HH[[sn]]$rt[1:250], "HH_Correct" = HH[[1]]$x[1:250], "HL_RT" = HL[[1]]$rt[1:250], "HL_Correct" = HL[[sn]]$x[1:250],  "LH_RT" = LH[[sn]]$rt[1:250], "LH_Correct" = LH[[sn]]$x[1:250],"LL_RT" = LL[[sn]]$rt[1:250], "LL_Correct" = LL[[sn]]$x[1:250]))
  
  
}

write.csv(dominance, file="Psi_Simulation_KStests.csv", row.names=FALSE)
write.csv(SFTresults, file="Psi_Simulation_SFTresults.csv", row.names=FALSE)
write.csv(PsiParams, file="Psi_Simulation_ParameterEstimatesByTrial.csv", row.names=FALSE)
write.csv(DDM.DFPperformance, file="Psi_Simulation_DDM-DFP_RT&ACC.csv", row.names=FALSE)

