
setwd("C:/Users/w018elf/Google Drive/Projects/Adaptive SFT/Model")
#setwd("~/Desktop/Adaptive SFT/Model")

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

pm.fun1 <- list()
result1.a <- read.csv("orientation_alpha.csv")
result1.b <- read.csv("orientation_beta.csv")
for (sn in 1:nsamps){
  pm.fun = matrix(data=NA, nrow=trials, ncol=length(axis.x))
  for (i in 1:trials){
     pm.fun[i,] = mapply(function(a,b) pm.function(axis.x, a, b, sim.d), a=result1.a[sn,i], b=result1.b[sn,i])
  }
  if(sn==1){pm.fun1 = list(pm.fun)}
  if(sn>1){pm.fun1 <- append(pm.fun1, list(pm.fun))}
}


#############CALCULATE THE MEAN PSYCHOMETRIC FUNCTION (length(x.axis)) FOR EACH TRIAL RUN (trials) ACROSS ALL MODEL SIMULATIONS (nsamps)##########
mean.pm.fun <- matrix(0, length(trials), length(axis.x))
for (ix in 1:length(nsamps)) {
  mean.pm.fun <- mean.pm.fun + data.frame(pm.fun1[[ix]])
}
mean.pm.fun <- mean.pm.fun / length(nsamps)

#Plot estimated psychometric function across trials -- 
#Each is an average function across all model runs (nsamps)
#Red (Trial 1) --> Blue (Trial N)
#Plot estimated psychometric function across trials -- 
#Each is an average function across all model runs (nsamps)
#Red (Trial 1) --> Blue (Trial N)
postscript("Psi_psychometricConvergence.eps", horizontal=F, height = 6, width = 6)
#windows()
par(mar=c(5,5,5,0))
par(fig=c(0,8,0,10)/10)
plot(0,0,type='n', xlim=range(axis.x), ylim=c(0,1), xlab="Intensity", ylab="Detection Rate", main =paste("Psychometric Function Convergence"))
for (i in 1:trials){
  lines(axis.x, mean.pm.fun[i,], col=rainbow(length(mean.pm.fun[,i]), start=0,end=4/6)[i])
}

par(fig=c(8,10,0,10)/10)
par(mar=c(5,1,5,0))
par(new=T)

ys = seq(0,1,l=10)+.02
ys[10] = ys[10]-.02

#postscript("Psi_psychometricConvergence.eps", horizontal=F, height = 6, width = 6)

legend_image <- as.raster(matrix(rainbow(1000,start=0,end=4/6)), ncol=1)
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')#, main = 'Number of Trials')
rasterImage(legend_image, 0, 0, 1, 1)
text(x=1.5, y = ys, labels = floor(seq(300,0,l=10)))

dev.off()



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
dev.off()


####PLOT ALPHA CONVERGENCE ACROSS TRIALS
#Average (SOLID LINE), 5%-95% (DASHED LINES), and true alpha value (RED DASHED LINE) CI functions across all model runs (nsamps)
postscript("Psi_parsConvergence.eps", horizontal=F,width=7.4, height=3.2)
par(mfrow=c(1,2), oma=c(0,0,1,0))

plot(apply(result1.a, 2, mean, na.rm=TRUE), lwd=2, type='l', main=paste("Location"), xlab="Trial", ylab="Parameter Estimate", ylim = c(min(apply(result1.a, 2, quantile, probs=.05, na.rm=TRUE))-1, max(apply(result1.a, 2, quantile, probs=.95, na.rm=TRUE))+1))
lines(apply(result1.a, 2, quantile, probs=.05, na.rm=TRUE), lty=2)
lines(apply(result1.a, 2, quantile, probs=.95, na.rm=TRUE), lty=2)
abline(h=fitted.params$par[1], lty=2, col='red')

####PLOT BETA CONVERGENCE ACROSS TRIALS
#Average (SOLID LINE), 5%-95% (DASHED LINES), and true beta value (RED DASHED LINE) CI functions across all model runs (nsamps)

plot(apply(result1.b, 2, mean, na.rm=TRUE), lwd=2, type='l', main=paste("Slope"), xlab="Trial", ylab="Parameter Estimate", ylim = c(min(apply(result1.b, 2, quantile, probs=.05, na.rm=TRUE))-1, max(apply(result1.b, 2, quantile, probs=.95, na.rm=TRUE))+1))
lines(apply(result1.b, 2, quantile, probs=.05, na.rm=TRUE), lty=2)
lines(apply(result1.b, 2, quantile, probs=.95, na.rm=TRUE), lty=2)
abline(h=fitted.params$par[2], lty=2, col='red')
title("Psi Parameter Convergence", outer=TRUE, line=-1, cex=1.5)

dev.off()


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
#trials = 100
trials = 150 #PER REQUEST FOR REVISITION (21FEB2019)

###############COLOR###############

#UNCOMMENT TO RESET HIGH/LOW VALUES FROM SIMULATION  
#NOTE: TAKES ~10 minutes
#result.color = Est.Trial.Psi.Color(trials)
sim.d = .01
#highSalience.color <- inv.pm.function(.99, result.color$alpha[length(result.color$alpha)], result.color$beta[length(result.color$alpha)], sim.d)
#lowSalience.color <- inv.pm.function(.90, result.color$alpha[length(result.color$alpha)], result.color$beta[length(result.color$beta)], sim.d)
x.range <- c(-55,50) #Color scale
thres50.color=6 #50% intensity estimate from human data

highSalience.color = 50 #COMMENT THIS LINE IF YOU RERAN A SIMULATION ABOVE
highSalience.color <- (highSalience.color-thres50.color)/(x.range[2]-thres50.color)

lowSalience.color = 16.45 #DATA FROM SIMULATED PARTICIPANT 1 OF DFP #COMMENT THIS LINE IF YOU RERAN A SIMULATION ABOVE
lowSalience.color <- (lowSalience.color-thres50.color)/(x.range[2]-thres50.color)

#PARAMETERS FROM SIMULATED PARTICICPANT 1 FROM FULL DFP STUDY
threshold=3.8 
v= 2.05
ter=0.07
sdv=.2 #FIXED

x.low.color <- simdiffT(nDFP,threshold,lowSalience.color*v,sdv,ter)
x.high.color <- simdiffT(nDFP,threshold,highSalience.color*v,sdv,ter)


#########ORIENTATION###############

#UNCOMMENT TO RESET HIGH/LOW VALUES FROM SIMULATION  
#NOTE: TAKES ~10 minutes
#result.orientation = Est.Trial.Psi.Orientation(200)#trials)
#highSalience.orientation <- inv.pm.function(.99, result.orientation$alpha[length(result.orientation$alpha)], result.orientation$beta[length(result.orientation$alpha)], sim.d)
#lowSalience.orientation <- inv.pm.function(.90, result.orientation$alpha[length(result.orientation$alpha)], result.orientation$beta[length(result.orientation$alpha)], sim.d)
x.range <- c(45,90) #orientation scale
thres50.orientation = 63 #50% intensity estimate from human data

highSalience.orientation = 90 #75.50087 #COMMENT THIS LINE IF YOU RERAN A SIMULATION ABOVE
highSalience.orientation <- (highSalience.orientation-thres50.orientation)/(x.range[2]-thres50.orientation)

lowSalience.orientation = 70.81 #DATA FROM SIMULATED PARTICIPANT 1 OF DFP  #COMMENT THIS LINE IF YOU RERAN A SIMULATION ABOVE
lowSalience.orientation <- (lowSalience.orientation-thres50.orientation)/(x.range[2]-thres50.orientation)


x.low.orientation <- simdiffT(nDFP,threshold,lowSalience.orientation*v,sdv,ter)
x.high.orientation<- simdiffT(nDFP,threshold,highSalience.orientation*v,sdv,ter)


###################################################################################################
############### SIMULATE A DDM MODEL WITH EACH ARCHITECTURE & STOPPING-RULE MODEL TYPE#############
##################PLOT THE SURVIVOR FUNCTIONS & SIC FOR EACH DFP MODEL PREDICTION##################
###################################################################################################
require(sft)
tvec.dense <- seq(0, 20, length.out=2000)

dfptrials = 100 
##PER REQUEST FOR REVISITION (21FEB2019)

# Parallel OR
HH <- dfp_ddm(dfptrials, highSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "PAR", "OR")
HL <- dfp_ddm(dfptrials, highSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "PAR", "OR")
LH <- dfp_ddm(dfptrials, lowSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "PAR", "OR")
LL <- dfp_ddm(dfptrials, lowSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "PAR", "OR")
p.or <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
            LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

#HH = 0; HL = 0; LH = 0; LL = 0

p.or.results <- data.frame("HH.rt" = HH$rt, "HL.rt" = HL$rt, "LH.rt" = LH$rt, "LL.rt" = LL$rt, "HH.acc" = HH$x, "HL.acc" = HL$x, "LH.acc" = LH$x, "LL.acc" = LL$x)
write.csv(p.or.results, file="p.or.rts.accs.csv", row.names=FALSE)
p.or.results = read.csv("p.or.rts.accs.csv")
HH$rt = p.or.results$HH.rt
HL$rt = p.or.results$HL.rt
LH$rt = p.or.results$LH.rt
LL$rt = p.or.results$LL.rt
HH$x = p.or.results$HH.acc
HL$x = p.or.results$HL.acc
LH$x = p.or.results$LH.acc
LL$x = p.or.results$LL.acc

p.or.all.results = cbind("D+ Statistic" = p.or$SICtest$positive$statistic, "D+ Pvalue" = p.or$SICtest$positive$p.value,  "D- Statistic" = p.or$SICtest$negative$statistic, "D- Pvalue" = p.or$SICtest$negative$p.value, "MIC Statistic" = p.or$MICtest$statistic, "MIC Pvalue" = p.or$MICtest$p.value, "HH_RT" = mean(HH$rt[1:250]), "HL_RT" = mean(HL$rt[1:250]), "LH_RT" = mean(LH$rt[1:250]), "LL_RT" = mean(LL$rt[1:250]), "HH_Correct" = mean(HH$x[1:250]) , "HL_Correct" = mean(HL$x[1:250]), "LH_Correct" = mean(LH$x[1:250]), "LL_Correct" = mean(LL$x[1:250]))
write.csv(p.or.all.results, file="p.or.results.csv", row.names=FALSE)


setEPS()
postscript("Psi_parallel_or.eps", width=6.6, height=3.3)
par(mfrow=c(1,2), oma=c(0,0,.75,0), mar=c(3.6, 3.1, 3.1, 1.1))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="", ylab="", xlim=c(0,max(LL$rt[LL$x==1])))
mtext("Time (s)", side=1, line=1.9)
mtext("S(t)", side=2, line=1.9)
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, p.or$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="", ylab="", xlim=c(0,max(LL$rt[LL$x==1])))
title("Parallel OR", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
mtext("Time (s)", side=1, line=1.9)
mtext("SIC(t)", side=2, line=1.9)
dev.off()


# Parallel AND
HH <- dfp_ddm(dfptrials, highSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "PAR", "AND")
HL <- dfp_ddm(dfptrials, highSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "PAR", "AND")
LH <- dfp_ddm(dfptrials, lowSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "PAR", "AND")
LL <- dfp_ddm(dfptrials, lowSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "PAR", "AND")
p.and <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
             LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

p.and.results <- data.frame("HH.rt" = HH$rt, "HL.rt" = HL$rt, "LH.rt" = LH$rt, "LL.rt" = LL$rt, "HH.acc" = HH$x, "HL.acc" = HL$x, "LH.acc" = LH$x, "LL.acc" = LL$x)
write.csv(p.and.results, file="p.and.rts.accs.csv", row.names=FALSE)


HH = 0; HL = 0; LH = 0; LL = 0
p.and.results = read.csv("p.and.rts.accs.csv")


HH$rt = p.and.results$HH.rt
HL$rt = p.and.results$HL.rt
LH$rt = p.and.results$LH.rt
LL$rt = p.and.results$LL.rt
HH$x = p.and.results$HH.acc
HL$x = p.and.results$HL.acc
LH$x = p.and.results$LH.acc
LL$x = p.and.results$LL.acc

p.and.all.results = cbind("D+ Statistic" = p.and$SICtest$positive$statistic, "D+ Pvalue" = p.and$SICtest$positive$p.value,  "D- Statistic" = p.and$SICtest$negative$statistic, "D- Pvalue" = p.and$SICtest$negative$p.value, "MIC Statistic" = p.and$MICtest$statistic, "MIC Pvalue" = p.and$MICtest$p.value, "HH_RT" = mean(HH$rt[1:250]), "HL_RT" = mean(HL$rt[1:250]), "LH_RT" = mean(LH$rt[1:250]), "LL_RT" = mean(LL$rt[1:dfptrials]), "HH_Correct" = mean(HH$x[1:dfptrials]) , "HL_Correct" = mean(HL$x[1:dfptrials]), "LH_Correct" = mean(LH$x[1:dfptrials]), "LL_Correct" = mean(LL$x[1:dfptrials]))
write.csv(p.and.all.results, file="p.and.results.csv", row.names=FALSE)



#png("parallel_and.png", 640,320)
setEPS()
postscript("Psi_parallel_and.eps", width=6.6, height=3.3)
par(mfrow=c(1,2), oma=c(0,0,.75,0), mar=c(3.6, 3.1, 3.1, 1.1))

plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="", ylab="", xlim=c(0,max(LL$rt[LL$x==1])))
mtext("Time (s)", side=1, line=1.9)
mtext("S(t)", side=2, line=1.9)
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, p.and$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="", ylab="", xlim=c(0,max(LL$rt[LL$x==1])))
title("Parallel AND", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
mtext("Time (s)", side=1, line=1.9) 
mtext("SIC(t)", side=2, line=1.9)
dev.off()

# Serial OR
HH <- dfp_ddm(dfptrials, highSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "SER", "OR")
HL <- dfp_ddm(dfptrials, highSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "SER", "OR")
LH <- dfp_ddm(dfptrials, lowSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "SER", "OR")
LL <- dfp_ddm(dfptrials, lowSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "SER", "OR")
s.or <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
            LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

s.or.results <- data.frame("HH.rt" = HH$rt, "HL.rt" = HL$rt, "LH.rt" = LH$rt, "LL.rt" = LL$rt, "HH.acc" = HH$x, "HL.acc" = HL$x, "LH.acc" = LH$x, "LL.acc" = LL$x)
write.csv(s.or.results, file="s.or.rts.accs.csv", row.names=FALSE)


HH = 0; HL = 0; LH = 0; LL = 0
s.or.results = read.csv("s.or.rts.accs.csv")

HH$rt = s.or.results$HH.rt
HL$rt = s.or.results$HL.rt
LH$rt = s.or.results$LH.rt
LL$rt = s.or.results$LL.rt
HH$x = s.or.results$HH.acc
HL$x = s.or.results$HL.acc
LH$x = s.or.results$LH.acc
LL$x = s.or.results$LL.acc

s.or.all.results = cbind("D+ Statistic" = s.or$SICtest$positive$statistic, "D+ Pvalue" = s.or$SICtest$positive$p.value,  "D- Statistic" = s.or$SICtest$negative$statistic, "D- Pvalue" = s.or$SICtest$negative$p.value, "MIC Statistic" = s.or$MICtest$statistic, "MIC Pvalue" = s.or$MICtest$p.value, "HH_RT" = mean(HH$rt[1:250]), "HL_RT" = mean(HL$rt[1:250]), "LH_RT" = mean(LH$rt[1:250]), "LL_RT" = mean(LL$rt[1:250]), "HH_Correct" = mean(HH$x[1:250]) , "HL_Correct" = mean(HL$x[1:250]), "LH_Correct" = mean(LH$x[1:250]), "LL_Correct" = mean(LL$x[1:250]))
write.csv(s.or.all.results, file="s.or.results.csv", row.names=FALSE)



#png("serial_or.png", 640,320)
setEPS()
postscript("Psi_serial_or.eps", width=6.6, height=3.3)
par(mfrow=c(1,2), oma=c(0,0,.75,0), mar=c(3.6, 3.1, 3.1, 1.1))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="", ylab="", xlim=c(0,max(LL$rt[LL$x==1])))
mtext("Time (s)", side=1, line=1.9)
mtext("S(t)", side=2, line=1.9)
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, s.or$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="", ylab="", xlim=c(0,max(LL$rt[LL$x==1])))
title("Serial OR", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
mtext("Time (s)", side=1, line=1.9) 
mtext("SIC(t)", side=2, line=1.9)
dev.off()


# Serial AND
HH <- dfp_ddm(dfptrials, highSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "SER","AND")
HL <- dfp_ddm(dfptrials, highSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "SER", "AND")
LH <- dfp_ddm(dfptrials, lowSalience.color*v, highSalience.orientation*v, threshold, ter, sdv, "SER", "AND")
LL <- dfp_ddm(dfptrials, lowSalience.color*v, lowSalience.orientation*v, threshold, ter, sdv, "SER", "AND")
s.and <- sic(HH=HH$rt[HH$x==1], HL=HL$rt[HL$x==1], 
             LH=LH$rt[LH$x==1], LL=LL$rt[LL$x==1])

s.and.results <- data.frame("HH.rt" = HH$rt, "HL.rt" = HL$rt, "LH.rt" = LH$rt, "LL.rt" = LL$rt, "HH.acc" = HH$x, "HL.acc" = HL$x, "LH.acc" = LH$x, "LL.acc" = LL$x)
write.csv(s.and.results, file="s.and.rts.accs.csv", row.names=FALSE)

HH = 0; HL = 0; LH = 0; LL = 0
s.and.results = read.csv("s.and.rts.accs.csv")

HH$rt = s.and.results$HH.rt
HL$rt = s.and.results$HL.rt
LH$rt = s.and.results$LH.rt
LL$rt = s.and.results$LL.rt
HH$x = s.and.results$HH.acc
HL$x = s.and.results$HL.acc
LH$x = s.and.results$LH.acc
LL$x = s.and.results$LL.acc


s.and.all.results = cbind("D+ Statistic" = s.and$SICtest$positive$statistic, "D+ Pvalue" = s.and$SICtest$positive$p.value,  "D- Statistic" = s.and$SICtest$negative$statistic, "D- Pvalue" = s.and$SICtest$negative$p.value, "MIC Statistic" = s.and$MICtest$statistic, "MIC Pvalue" = s.and$MICtest$p.value, "HH_RT" = mean(HH$rt[1:250]), "HL_RT" = mean(HL$rt[1:250]), "LH_RT" = mean(LH$rt[1:250]), "LL_RT" = mean(LL$rt[1:250]), "HH_Correct" = mean(HH$x[1:250]) , "HL_Correct" = mean(HL$x[1:250]), "LH_Correct" = mean(LH$x[1:250]), "LL_Correct" = mean(LL$x[1:250]))
write.csv(s.and.all.results, file="s.and.results.csv", row.names=FALSE)



#png("serial_and.png", 640,320)
setEPS()
postscript("Psi_serial_and.eps", width=6.6, height=3.3)
par(mfrow=c(1,2), oma=c(0,0,.75,0), mar=c(3.6, 3.1, 3.1, 1.1))
plot(tvec.dense, 1-ecdf(HH$rt[HH$x==1])(tvec.dense), type='l', col='red',
     xlab="", ylab="", xlim=c(0,max(LL$rt[LL$x==1])))
mtext("Time (s)", side=1, line=1.9)
mtext("S(t)", side=2, line=1.9)
lines(tvec.dense, 1-ecdf(HL$rt[HL$x==1])(tvec.dense), col='orange')
lines(tvec.dense, 1-ecdf(LH$rt[LH$x==1])(tvec.dense), col='purple')
lines(tvec.dense, 1-ecdf(LL$rt[LL$x==1])(tvec.dense), col='blue')
legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
       col=c("red", "orange", "purple", "blue"))

plot(tvec.dense, s.and$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
     xlab="", ylab="", xlim=c(0,max(LL$rt[LL$x==1])))
title("Serial AND", outer=TRUE, line=-2, cex=1.5)
abline(h=0)
mtext("Time (s)", side=1, line=1.9) 
mtext("SIC(t)", side=2, line=1.9)
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

nParticipants <- 10 #number of simulated participants

sTrials = dfptrials #number of trials per salience combination
nDFP = sTrials*4 #number of DFP trials

result.color <- list()
result.color$alpha <- vector() #stores all psi color results
result.color$beta <- vector()
#result.color <- list(result.color, result.color)
ddm.result.color <- list() #stores color ddm performance @ high/low values

result.orientation <- list() #stores all psi orientation results
result.orientation$alpha <- vector() #stores all psi color results
result.orientation$beta <- vector()
#result.orientation <- list(result.orientation, result.orientation)
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

#Use same DDM paramaters from simulated convergence testing

allpars <- array(NA, c(nParticipants, 4))

sftresults = read.csv("Psi_Simulation_SFTresults.csv")  
allpars = cbind(sftresults$Threshold, sftresults$v, sftresults$ter, sftresults$sdv)

psiresults = read.csv("PsiDDM_Simulation_Pars.csv")  

for (sn in 1:nParticipants) { 
  #threshold=3#1.45
  #v=2#1.6
  #ter=.1
  #sdv=.2#.25
  #threshold.p <- 0
  #threshold.p <- rnorm(1, mean=1.1*threshold, sd=threshold/8)
  #v.p <- 0
  #v.p <- rnorm(1, mean=v/threshold*threshold.p, sd=v/8)
  #ter.p <- 0
  #ter.p <- rnorm(1, mean=ter, sd=ter/6)
  #sdv.p <- 0
  #sdv.p <- sdv
  
  #allpars[sn,] <- c(threshold.p, v.p, ter.p, sdv.p)
  
  result.color$alpha =  psiresults$Color_Alpha[psiresults$Subject==sn]
  result.color$beta = psiresults$Color_Beta[psiresults$Subject==sn]
  result.orientation$alpha = psiresults$Orientation_Alpha[psiresults$Subject==sn]
  result.orientation$beta = psiresults$Orientation_Beta[psiresults$Subject==sn]

  
  #result.color[[sn]] = Est.Trial.Psi.Color(nTrials, allpars[sn,]) #ddm runs in psi color experiment
  ddm.result.color[[sn]] = psi_color_ddm(result.color, nDFP, allpars[sn,]) #ddm performance in 1000 trials with high & low salience color trials
   
  
  #result.orientation[[sn]] = Est.Trial.Psi.Orientation(nTrials, allpars[sn,]) #ddm runs in psi orientation experiment
  ddm.result.orientation[[sn]] = psi_orientation_ddm(result.orientation, nDFP, allpars[sn,])#ddm performance in 1000 trials with high & low salience orientation trials
  
  require(sft)
  tvec.dense <- seq(0, 5, length.out=2000)
  

  threshold = allpars[sn,1]
  v = allpars[sn,2]
  ter = allpars[sn,3]
  sdv = allpars[sn,4]
  
  #Simulate DDM performance for each trial type of DFP (assumes DDM has Parallel OR processes)
  HH[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "PAR", "OR")
  HL[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "PAR", "OR")
  LH[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "PAR", "OR")
  LL[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "PAR", "OR")
  
  #Estimates SIC from DDM DFP performance data
  p.or[[sn]] <- sic(HH=HH[[sn]]$rt[HH[[sn]]$x==1], HL=HL[[sn]]$rt[HL[[sn]]$x==1], 
                    LH=LH[[sn]]$rt[LH[[sn]]$x==1], LL=LL[[sn]]$rt[LL[[sn]]$x==1])
  
  #Plots & Saves the DDM survivor functions for each DFP trial type (HH, HL, LH, LL) & SIC Function
  postscript("Psi_parallel_or.eps", horizontal=FALSE, width = 6, height = 4)
  
  #png(paste("parallel_or Subject ", sn, ".png"), 640,320)
  par(mfrow=c(1,2), oma=c(0,0,1,0))
  
  plot(tvec.dense, 1-ecdf(HH[[sn]]$rt[HH[[sn]]$x==1])(tvec.dense), type='l', col='red',
       xlab="Time (s)", ylab="S(t)", xlim=c(0,5))#max(LL[[sn]]$rt)))
  lines(tvec.dense, 1-ecdf(HL[[sn]]$rt[HL[[sn]]$x==1])(tvec.dense), col='orange')
  lines(tvec.dense, 1-ecdf(LH[[sn]]$rt[LH[[sn]]$x==1])(tvec.dense), col='purple')
  lines(tvec.dense, 1-ecdf(LL[[sn]]$rt[LL[[sn]]$x==1])(tvec.dense), col='blue')
  legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
         col=c("red", "orange", "purple", "blue"))
  
  plot(tvec.dense, p.or[[sn]]$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
       xlab="Time (s)", ylab="SIC", xlim=c(0,5))#max(LL[[sn]]$rt)))
  title("Parallel OR", outer=TRUE, line=-2, cex=1.5)
  abline(h=0)
  dev.off()
  
  print(sn)
}

#PsiParams = cbind("Subject" = rep(1, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[1]]$alpha, "Color_Beta" = result.color[[1]]$beta, "Orientation_Alpha" = result.orientation[[1]]$alpha, "Orientation_Beta" = result.orientation[[1]]$beta)
dominance = cbind("Subject"=rep(1, 8), p.or[[1]]$Dominance)
SFTresults = cbind("Subject"= 1,  "Threshold" = allpars[1,1], "v" = allpars[1,2], "ter" = allpars[1,3], "sdv" = allpars[1,4], "H_Color Intensity" = ddm.result.color[[1]]$intensity[1],"L_Color Intensity" = ddm.result.color[[1]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[2],  "D+ Statistic" = p.or[[1]]$SICtest$positive$statistic, "D+ Pvalue" = p.or[[1]]$SICtest$positive$p.value,  "D- Statistic" = p.or[[1]]$SICtest$negative$statistic, "D- Pvalue" = p.or[[1]]$SICtest$negative$p.value, "MIC Statistic" = p.or[[1]]$MICtest$statistic, "MIC Pvalue" = p.or[[1]]$MICtest$p.value, "HH_RT" = mean(HH[[1]]$rt[1:dfptrials]), "HL_RT" = mean(HL[[1]]$rt[1:dfptrials]), "LH_RT" = mean(LH[[1]]$rt[1:dfptrials]), "LL_RT" = mean(LL[[1]]$rt[1:dfptrials]), "HH_Correct" = mean(HH[[1]]$x[1:dfptrials]) , "HL_Correct" = mean(HL[[1]]$x[1:dfptrials]), "LH_Correct" = mean(LH[[1]]$x[1:dfptrials]), "LL_Correct" = mean(LL[[1]]$x[1:dfptrials]))
DDM.DFPperformance = cbind(cbind("Subject" = rep(1, dfptrials), "Trial" = 1:dfptrials), "HH_RT" = HH[[1]]$rt[1:dfptrials], "HH_Correct" = HH[[1]]$x[1:dfptrials], "HL_RT" = HL[[1]]$rt[1:dfptrials], "HL_Correct" = HL[[1]]$x[1:dfptrials],  "LH_RT" = LH[[1]]$rt[1:dfptrials], "LH_Correct" = LH[[1]]$x[1:dfptrials],"LL_RT" = LL[[1]]$rt[1:dfptrials], "LL_Correct" = LL[[1]]$x[1:dfptrials])

for (sn in 2:nParticipants) { 
  #PsiParams = rbind(PsiParams, cbind("Subject" = rep(sn, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[sn]]$alpha, "Color_Beta" = result.color[[sn]]$beta, "Orientation_Alpha" = result.orientation[[sn]]$alpha, "Orientation_Beta" = result.orientation[[sn]]$beta))
  SFTresults = rbind(SFTresults, cbind("Subject"= sn,  "Threshold" = allpars[sn,1], "v" = allpars[sn,2], "ter" = allpars[sn,3], "sdv" = allpars[sn,4], "H_Color Intensity" = ddm.result.color[[sn]]$intensity[1],"L_Color Intensity" = ddm.result.color[[sn]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[2],  "D+ Statistic" = p.or[[sn]]$SICtest$positive$statistic, "D+ Pvalue" = p.or[[sn]]$SICtest$positive$p.value,  "D- Statistic" = p.or[[sn]]$SICtest$negative$statistic, "D- Pvalue" = p.or[[sn]]$SICtest$negative$p.value, "MIC Statistic" = p.or[[sn]]$MICtest$statistic, "MIC Pvalue" = p.or[[sn]]$MICtest$p.value,  "HH_RT" = mean(HH[[sn]]$rt[1:dfptrials]), "HL_RT" = mean(HL[[sn]]$rt[1:dfptrials]), "LH_RT" = mean(LH[[sn]]$rt[1:dfptrials]), "LL_RT" = mean(LL[[sn]]$rt[1:dfptrials]), "HH_Correct" = mean(HH[[sn]]$x[1:dfptrials]) , "HL_Correct" = mean(HL[[sn]]$x[1:dfptrials]), "LH_Correct" = mean(LH[[sn]]$x[1:dfptrials]), "LL_Correct" = mean(LL[[sn]]$x[1:dfptrials])))
  dominance = rbind(dominance, cbind("Subject"=rep(sn, 8), p.or[[sn]]$Dominance))
  DDM.DFPperformance = rbind(DDM.DFPperformance, cbind("Subject" = rep(sn, dfptrials), "Trial" = 1:dfptrials, "HH_RT" = HH[[sn]]$rt[1:dfptrials], "HH_Correct" = HH[[1]]$x[1:dfptrials], "HL_RT" = HL[[1]]$rt[1:dfptrials], "HL_Correct" = HL[[sn]]$x[1:dfptrials],  "LH_RT" = LH[[sn]]$rt[1:dfptrials], "LH_Correct" = LH[[sn]]$x[1:dfptrials],"LL_RT" = LL[[sn]]$rt[1:dfptrials], "LL_Correct" = LL[[sn]]$x[1:dfptrials]))
}

write.csv(dominance, file="ParallelOR_Psi_Simulation_KStests.csv", row.names=FALSE)
write.csv(SFTresults, file="ParallelOR_Psi_Simulation_SFTresults.csv", row.names=FALSE)
#write.csv(PsiParams, file="ParallelOR_PsiDDM_Simulation_Pars.csv", row.names=FALSE)
write.csv(DDM.DFPperformance, file="ParallelOR_Psi_Simulation_DDM-DFP_RT&ACC.csv", row.names=FALSE)










###################################################################################################
############SIMULATE A FULL EXPERIMENT WITH 100 DDM SUBJECTS WITH PARALLEL-AND PROCESSES############
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

nParticipants <- 10 #number of simulated participants
nTrials <- dfptrials #100 #number of psi trials ---- NOTE: based on convergence of alpha/beta parameters in convergence study 

sTrials = dfptrials #number of trials per salience combination
nDFP = sTrials*4 #number of DFP trials

result.color <- list()
result.color$alpha <- vector() #stores all psi color results
result.color$beta <- vector()
#result.color <- list(result.color, result.color)
ddm.result.color <- list() #stores color ddm performance @ high/low values

result.orientation <- list() #stores all psi orientation results
result.orientation$alpha <- vector() #stores all psi color results
result.orientation$beta <- vector()
#result.orientation <- list(result.orientation, result.orientation)
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
p.and <- list()

#Use same DDM paramaters from simulated convergence testing

allpars <- array(NA, c(nParticipants, 4))

sftresults = read.csv("Psi_Simulation_SFTresults.csv")  
allpars = cbind(sftresults$Threshold, sftresults$v, sftresults$ter, sftresults$sdv)

psiresults = read.csv("PsiDDM_Simulation_Pars.csv")  

for (sn in 1:nParticipants) { 
  result.color$alpha =  psiresults$Color_Alpha[psiresults$Subject==sn]
  result.color$beta = psiresults$Color_Beta[psiresults$Subject==sn]
  result.orientation$alpha = psiresults$Orientation_Alpha[psiresults$Subject==sn]
  result.orientation$beta = psiresults$Orientation_Beta[psiresults$Subject==sn]
  
  
  #result.color[[sn]] = Est.Trial.Psi.Color(nTrials, allpars[sn,]) #ddm runs in psi color experiment
  ddm.result.color[[sn]] = psi_color_ddm(result.color, nDFP, allpars[sn,]) #ddm performance in 1000 trials with high & low salience color trials
  
  
  #result.orientation[[sn]] = Est.Trial.Psi.Orientation(nTrials, allpars[sn,]) #ddm runs in psi orientation experiment
  ddm.result.orientation[[sn]] = psi_orientation_ddm(result.orientation, nDFP, allpars[sn,])#ddm performance in 1000 trials with high & low salience orientation trials
  
  require(sft)
  tvec.dense <- seq(0, 5, length.out=2000)
  
  
  threshold = allpars[sn,1]
  v = allpars[sn,2]
  ter = allpars[sn,3]
  sdv = allpars[sn,4]
  
  #Simulate DDM performance for each trial type of DFP (assumes DDM has Parallel OR processes)
  HH[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "PAR", "AND")
  HL[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "PAR", "AND")
  LH[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "PAR", "AND")
  LL[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "PAR", "AND")
  
  
  #Estimates SIC from DDM DFP performance data
  p.and[[sn]] <- sic(HH=HH[[sn]]$rt[HH[[sn]]$x==1], HL=HL[[sn]]$rt[HL[[sn]]$x==1], 
                    LH=LH[[sn]]$rt[LH[[sn]]$x==1], LL=LL[[sn]]$rt[LL[[sn]]$x==1], mictest='art')
  
  #Plots & Saves the DDM survivor functions for each DFP trial type (HH, HL, LH, LL) & SIC Function
  #png(paste("Parallel_AND Subject ", sn, ".png"), 640,320)
  postscript("Psi_parallel_and.eps", horizontal=FALSE, width = 6, height = 4)
  par(mfrow=c(1,2), oma=c(0,0,1,0))
  plot(tvec.dense, 1-ecdf(HH[[sn]]$rt[HH[[sn]]$x==1])(tvec.dense), type='l', col='red',
       xlab="Time (s)", ylab="S(t)", xlim=c(0,5))#max(LL[[sn]]$rt)))
  lines(tvec.dense, 1-ecdf(HL[[sn]]$rt[HL[[sn]]$x==1])(tvec.dense), col='orange')
  lines(tvec.dense, 1-ecdf(LH[[sn]]$rt[LH[[sn]]$x==1])(tvec.dense), col='purple')
  lines(tvec.dense, 1-ecdf(LL[[sn]]$rt[LL[[sn]]$x==1])(tvec.dense), col='blue')
  legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
         col=c("red", "orange", "purple", "blue"))
  
  plot(tvec.dense, p.and[[sn]]$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
       xlab="Time (s)", ylab="SIC", xlim=c(0,5))#max(LL[[sn]]$rt)))
  title("Parallel AND", outer=TRUE, line=-2, cex=1.5)
  abline(h=0)
  dev.off()
  
  print(sn)
}

#PsiParams = cbind("Subject" = rep(1, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[1]]$alpha, "Color_Beta" = result.color[[1]]$beta, "Orientation_Alpha" = result.orientation[[1]]$alpha, "Orientation_Beta" = result.orientation[[1]]$beta)
dominance = cbind("Subject"=rep(1, 8), p.and[[1]]$Dominance)
SFTresults = cbind("Subject"= 1,  "Threshold" = allpars[1,1], "v" = allpars[1,2], "ter" = allpars[1,3], "sdv" = allpars[1,4], "H_Color Intensity" = ddm.result.color[[1]]$intensity[1],"L_Color Intensity" = ddm.result.color[[1]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[2],  "D+ Statistic" = p.and[[1]]$SICtest$positive$statistic, "D+ Pvalue" = p.and[[1]]$SICtest$positive$p.value,  "D- Statistic" = p.and[[1]]$SICtest$negative$statistic, "D- Pvalue" = p.and[[1]]$SICtest$negative$p.value, "MIC Statistic" = p.and[[1]]$MICtest$statistic, "MIC Pvalue" = p.and[[1]]$MICtest$p.value, "HH_RT" = mean(HH[[1]]$rt[1:dfptrials]), "HL_RT" = mean(HL[[1]]$rt[1:dfptrials]), "LH_RT" = mean(LH[[1]]$rt[1:dfptrials]), "LL_RT" = mean(LL[[1]]$rt[1:dfptrials]), "HH_Correct" = mean(HH[[1]]$x[1:dfptrials]) , "HL_Correct" = mean(HL[[1]]$x[1:dfptrials]), "LH_Correct" = mean(LH[[1]]$x[1:dfptrials]), "LL_Correct" = mean(LL[[1]]$x[1:dfptrials]))
DDM.DFPperformance = cbind(cbind("Subject" = rep(1, dfptrials), "Trial" = 1:dfptrials), "HH_RT" = HH[[1]]$rt[1:dfptrials], "HH_Correct" = HH[[1]]$x[1:dfptrials], "HL_RT" = HL[[1]]$rt[1:dfptrials], "HL_Correct" = HL[[1]]$x[1:dfptrials],  "LH_RT" = LH[[1]]$rt[1:dfptrials], "LH_Correct" = LH[[1]]$x[1:dfptrials],"LL_RT" = LL[[1]]$rt[1:dfptrials], "LL_Correct" = LL[[1]]$x[1:dfptrials])

for (sn in 2:nParticipants) { 
  # PsiParams = rbind(PsiParams, cbind("Subject" = rep(sn, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[sn]]$alpha, "Color_Beta" = result.color[[sn]]$beta, "Orientation_Alpha" = result.orientation[[sn]]$alpha, "Orientation_Beta" = result.orientation[[sn]]$beta))
  SFTresults = rbind(SFTresults, cbind("Subject"= sn,  "Threshold" = allpars[sn,1], "v" = allpars[sn,2], "ter" = allpars[sn,3], "sdv" = allpars[sn,4], "H_Color Intensity" = ddm.result.color[[sn]]$intensity[1],"L_Color Intensity" = ddm.result.color[[sn]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[2],  "D+ Statistic" = p.and[[sn]]$SICtest$positive$statistic, "D+ Pvalue" = p.and[[sn]]$SICtest$positive$p.value,  "D- Statistic" = p.and[[sn]]$SICtest$negative$statistic, "D- Pvalue" = p.and[[sn]]$SICtest$negative$p.value, "MIC Statistic" = p.and[[sn]]$MICtest$statistic, "MIC Pvalue" = p.and[[sn]]$MICtest$p.value,  "HH_RT" = mean(HH[[sn]]$rt[1:dfptrials]), "HL_RT" = mean(HL[[sn]]$rt[1:dfptrials]), "LH_RT" = mean(LH[[sn]]$rt[1:dfptrials]), "LL_RT" = mean(LL[[sn]]$rt[1:dfptrials]), "HH_Correct" = mean(HH[[sn]]$x[1:dfptrials]) , "HL_Correct" = mean(HL[[sn]]$x[1:dfptrials]), "LH_Correct" = mean(LH[[sn]]$x[1:dfptrials]), "LL_Correct" = mean(LL[[sn]]$x[1:dfptrials])))
  dominance = rbind(dominance, cbind("Subject"=rep(sn, 8), p.and[[sn]]$Dominance))
  DDM.DFPperformance = rbind(DDM.DFPperformance, cbind("Subject" = rep(sn, dfptrials), "Trial" = 1:dfptrials, "HH_RT" = HH[[sn]]$rt[1:dfptrials], "HH_Correct" = HH[[1]]$x[1:dfptrials], "HL_RT" = HL[[1]]$rt[1:dfptrials], "HL_Correct" = HL[[sn]]$x[1:dfptrials],  "LH_RT" = LH[[sn]]$rt[1:dfptrials], "LH_Correct" = LH[[sn]]$x[1:dfptrials],"LL_RT" = LL[[sn]]$rt[1:dfptrials], "LL_Correct" = LL[[sn]]$x[1:dfptrials]))
}

write.csv(dominance, file="ParallelAND_Psi_Simulation_KStests.csv", row.names=FALSE)
write.csv(SFTresults, file="ParallelAND_Psi_Simulation_SFTresults.csv", row.names=FALSE)
#write.csv(PsiParams, file="ParallelAND_PsiDDM_Simulation_Pars.csv", row.names=FALSE)
write.csv(DDM.DFPperformance, file="ParallelAND_Psi_Simulation_DDM-DFP_RT&ACC.csv", row.names=FALSE)










###################################################################################################
############SIMULATE A FULL EXPERIMENT WITH 100 DDM SUBJECTS WITH SERIAL-OR PROCESSES############
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

nParticipants <- 10 #number of simulated participants
nTrials <- dfptrials #100 #number of psi trials ---- NOTE: based on convergence of alpha/beta parameters in convergence study 

sTrials = dfptrials #number of trials per salience combination
nDFP = sTrials*4 #number of DFP trials

result.color <- list()
result.color$alpha <- vector() #stores all psi color results
result.color$beta <- vector()
#result.color <- list(result.color, result.color)
ddm.result.color <- list() #stores color ddm performance @ high/low values

result.orientation <- list() #stores all psi orientation results
result.orientation$alpha <- vector() #stores all psi color results
result.orientation$beta <- vector()
#result.orientation <- list(result.orientation, result.orientation)
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
s.or <- list()

#Use same DDM paramaters from simulated convergence testing

allpars <- array(NA, c(nParticipants, 4))

sftresults = read.csv("Psi_Simulation_SFTresults.csv")  
allpars = cbind(sftresults$Threshold, sftresults$v, sftresults$ter, sftresults$sdv)

psiresults = read.csv("PsiDDM_Simulation_Pars.csv")  

for (sn in 1:nParticipants) { 
  result.color$alpha =  psiresults$Color_Alpha[psiresults$Subject==sn]
  result.color$beta = psiresults$Color_Beta[psiresults$Subject==sn]
  result.orientation$alpha = psiresults$Orientation_Alpha[psiresults$Subject==sn]
  result.orientation$beta = psiresults$Orientation_Beta[psiresults$Subject==sn]
  
  
  #result.color[[sn]] = Est.Trial.Psi.Color(nTrials, allpars[sn,]) #ddm runs in psi color experiment
  ddm.result.color[[sn]] = psi_color_ddm(result.color, nDFP, allpars[sn,]) #ddm performance in 1000 trials with high & low salience color trials
  
  
  #result.orientation[[sn]] = Est.Trial.Psi.Orientation(nTrials, allpars[sn,]) #ddm runs in psi orientation experiment
  ddm.result.orientation[[sn]] = psi_orientation_ddm(result.orientation, nDFP, allpars[sn,])#ddm performance in 1000 trials with high & low salience orientation trials
  
  require(sft)
  tvec.dense <- seq(0, 5, length.out=2000)
  
  
  threshold = allpars[sn,1]
  v = allpars[sn,2]
  ter = allpars[sn,3]
  sdv = allpars[sn,4]
  
  #Simulate DDM performance for each trial type of DFP (assumes DDM has Parallel OR processes)
  HH[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "SER", "OR")
  HL[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "SER", "OR")
  LH[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "SER", "OR")
  LL[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "SER", "OR")
  
  
  #Estimates SIC from DDM DFP performance data
  s.or[[sn]] <- sic(HH=HH[[sn]]$rt[HH[[sn]]$x==1], HL=HL[[sn]]$rt[HL[[sn]]$x==1], 
                     LH=LH[[sn]]$rt[LH[[sn]]$x==1], LL=LL[[sn]]$rt[LL[[sn]]$x==1], mictest='art')
  
  #Plots & Saves the DDM survivor functions for each DFP trial type (HH, HL, LH, LL) & SIC Function
  #png(paste("Serial_OR Subject ", sn, ".png"), 640,320)
  postscript("Psi_serial_or.eps", horizontal=FALSE, width = 6, height = 4)
  
  par(mfrow=c(1,2), oma=c(0,0,1,0))
  plot(tvec.dense, 1-ecdf(HH[[sn]]$rt[HH[[sn]]$x==1])(tvec.dense), type='l', col='red',
       xlab="Time (s)", ylab="S(t)", xlim=c(0,5))#max(LL[[sn]]$rt)))
  lines(tvec.dense, 1-ecdf(HL[[sn]]$rt[HL[[sn]]$x==1])(tvec.dense), col='orange')
  lines(tvec.dense, 1-ecdf(LH[[sn]]$rt[LH[[sn]]$x==1])(tvec.dense), col='purple')
  lines(tvec.dense, 1-ecdf(LL[[sn]]$rt[LL[[sn]]$x==1])(tvec.dense), col='blue')
  legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
         col=c("red", "orange", "purple", "blue"))
  
  plot(tvec.dense, s.or[[sn]]$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
       xlab="Time (s)", ylab="SIC", xlim=c(0,5))#max(LL[[sn]]$rt)))
  title("Serial OR", outer=TRUE, line=-2, cex=1.5)
  abline(h=0)
  dev.off()
  
  print(sn)
}

#PsiParams = cbind("Subject" = rep(1, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[1]]$alpha, "Color_Beta" = result.color[[1]]$beta, "Orientation_Alpha" = result.orientation[[1]]$alpha, "Orientation_Beta" = result.orientation[[1]]$beta)
dominance = cbind("Subject"=rep(1, 8), s.or[[1]]$Dominance)
SFTresults = cbind("Subject"= 1,  "Threshold" = allpars[1,1], "v" = allpars[1,2], "ter" = allpars[1,3], "sdv" = allpars[1,4], "H_Color Intensity" = ddm.result.color[[1]]$intensity[1],"L_Color Intensity" = ddm.result.color[[1]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[2],  "D+ Statistic" = s.or[[1]]$SICtest$positive$statistic, "D+ Pvalue" = s.or[[1]]$SICtest$positive$p.value,  "D- Statistic" = s.or[[1]]$SICtest$negative$statistic, "D- Pvalue" = s.or[[1]]$SICtest$negative$p.value, "MIC Statistic" = s.or[[1]]$MICtest$statistic, "MIC Pvalue" = s.or[[1]]$MICtest$p.value, "HH_RT" = mean(HH[[1]]$rt[1:dfptrials]), "HL_RT" = mean(HL[[1]]$rt[1:dfptrials]), "LH_RT" = mean(LH[[1]]$rt[1:dfptrials]), "LL_RT" = mean(LL[[1]]$rt[1:dfptrials]), "HH_Correct" = mean(HH[[1]]$x[1:dfptrials]) , "HL_Correct" = mean(HL[[1]]$x[1:dfptrials]), "LH_Correct" = mean(LH[[1]]$x[1:dfptrials]), "LL_Correct" = mean(LL[[1]]$x[1:dfptrials]))
DDM.DFPperformance = cbind(cbind("Subject" = rep(1, dfptrials), "Trial" = 1:dfptrials), "HH_RT" = HH[[1]]$rt[1:dfptrials], "HH_Correct" = HH[[1]]$x[1:dfptrials], "HL_RT" = HL[[1]]$rt[1:dfptrials], "HL_Correct" = HL[[1]]$x[1:dfptrials],  "LH_RT" = LH[[1]]$rt[1:dfptrials], "LH_Correct" = LH[[1]]$x[1:dfptrials],"LL_RT" = LL[[1]]$rt[1:dfptrials], "LL_Correct" = LL[[1]]$x[1:dfptrials])

for (sn in 2:nParticipants) { 
  # PsiParams = rbind(PsiParams, cbind("Subject" = rep(sn, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[sn]]$alpha, "Color_Beta" = result.color[[sn]]$beta, "Orientation_Alpha" = result.orientation[[sn]]$alpha, "Orientation_Beta" = result.orientation[[sn]]$beta))
  SFTresults = rbind(SFTresults, cbind("Subject"= sn,  "Threshold" = allpars[sn,1], "v" = allpars[sn,2], "ter" = allpars[sn,3], "sdv" = allpars[sn,4], "H_Color Intensity" = ddm.result.color[[sn]]$intensity[1],"L_Color Intensity" = ddm.result.color[[sn]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[2],  "D+ Statistic" = s.or[[sn]]$SICtest$positive$statistic, "D+ Pvalue" = s.or[[sn]]$SICtest$positive$p.value,  "D- Statistic" = s.or[[sn]]$SICtest$negative$statistic, "D- Pvalue" = s.or[[sn]]$SICtest$negative$p.value, "MIC Statistic" = s.or[[sn]]$MICtest$statistic, "MIC Pvalue" = s.or[[sn]]$MICtest$p.value,  "HH_RT" = mean(HH[[sn]]$rt[1:dfptrials]), "HL_RT" = mean(HL[[sn]]$rt[1:dfptrials]), "LH_RT" = mean(LH[[sn]]$rt[1:dfptrials]), "LL_RT" = mean(LL[[sn]]$rt[1:dfptrials]), "HH_Correct" = mean(HH[[sn]]$x[1:dfptrials]) , "HL_Correct" = mean(HL[[sn]]$x[1:dfptrials]), "LH_Correct" = mean(LH[[sn]]$x[1:dfptrials]), "LL_Correct" = mean(LL[[sn]]$x[1:dfptrials])))
  dominance = rbind(dominance, cbind("Subject"=rep(sn, 8), s.or[[sn]]$Dominance))
  DDM.DFPperformance = rbind(DDM.DFPperformance, cbind("Subject" = rep(sn, dfptrials), "Trial" = 1:dfptrials, "HH_RT" = HH[[sn]]$rt[1:dfptrials], "HH_Correct" = HH[[1]]$x[1:dfptrials], "HL_RT" = HL[[1]]$rt[1:dfptrials], "HL_Correct" = HL[[sn]]$x[1:dfptrials],  "LH_RT" = LH[[sn]]$rt[1:dfptrials], "LH_Correct" = LH[[sn]]$x[1:dfptrials],"LL_RT" = LL[[sn]]$rt[1:dfptrials], "LL_Correct" = LL[[sn]]$x[1:dfptrials]))
}

write.csv(dominance, file="SerialOR_Psi_Simulation_KStests.csv", row.names=FALSE)
write.csv(SFTresults, file="SerialOR_Psi_Simulation_SFTresults.csv", row.names=FALSE)
#write.csv(PsiParams, file="SerialOR_PsiDDM_Simulation_Pars.csv", row.names=FALSE)
write.csv(DDM.DFPperformance, file="SerialOR_Psi_Simulation_DDM-DFP_RT&ACC.csv", row.names=FALSE)










###################################################################################################
############SIMULATE A FULL EXPERIMENT WITH 100 DDM SUBJECTS WITH SERIAL-AND PROCESSES############
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

nParticipants <- 10 #number of simulated participants
nTrials <- dfptrials #100 #number of psi trials ---- NOTE: based on convergence of alpha/beta parameters in convergence study 

sTrials = dfptrials #number of trials per salience combination
nDFP = sTrials*4 #number of DFP trials

result.color <- list()
result.color$alpha <- vector() #stores all psi color results
result.color$beta <- vector()
#result.color <- list(result.color, result.color)
ddm.result.color <- list() #stores color ddm performance @ high/low values

result.orientation <- list() #stores all psi orientation results
result.orientation$alpha <- vector() #stores all psi color results
result.orientation$beta <- vector()
#result.orientation <- list(result.orientation, result.orientation)
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
s.and <- list()

#Use same DDM paramaters from simulated convergence testing

allpars <- array(NA, c(nParticipants, 4))

sftresults = read.csv("Psi_Simulation_SFTresults.csv")  
allpars = cbind(sftresults$Threshold, sftresults$v, sftresults$ter, sftresults$sdv)

psiresults = read.csv("PsiDDM_Simulation_Pars.csv")  

for (sn in 1:nParticipants) { 
  result.color$alpha =  psiresults$Color_Alpha[psiresults$Subject==sn]
  result.color$beta = psiresults$Color_Beta[psiresults$Subject==sn]
  result.orientation$alpha = psiresults$Orientation_Alpha[psiresults$Subject==sn]
  result.orientation$beta = psiresults$Orientation_Beta[psiresults$Subject==sn]
  
  
  #result.color[[sn]] = Est.Trial.Psi.Color(nTrials, allpars[sn,]) #ddm runs in psi color experiment
  ddm.result.color[[sn]] = psi_color_ddm(result.color, nDFP, allpars[sn,]) #ddm performance in 1000 trials with high & low salience color trials
  
  
  #result.orientation[[sn]] = Est.Trial.Psi.Orientation(nTrials, allpars[sn,]) #ddm runs in psi orientation experiment
  ddm.result.orientation[[sn]] = psi_orientation_ddm(result.orientation, nDFP, allpars[sn,])#ddm performance in 1000 trials with high & low salience orientation trials
  
  require(sft)
  tvec.dense <- seq(0, 5, length.out=2000)
  
  
  threshold = allpars[sn,1]
  v = allpars[sn,2]
  ter = allpars[sn,3]
  sdv = allpars[sn,4]
  
  #Simulate DDM performance for each trial type of DFP (assumes DDM has Parallel OR processes)
  HH[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "SER", "AND")
  HL[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[1]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "SER", "AND")
  LH[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[1]]*v, threshold, ter, sdv, "SER", "AND")
  LL[[sn]] <- dfp_ddm(sTrials, ddm.result.color[[sn]]$rescaled.intensity[[2]]*v, ddm.result.orientation[[sn]]$rescaled.intensity[[2]]*v, threshold, ter, sdv, "SER", "AND")
  
  
  #Estimates SIC from DDM DFP performance data
  s.and[[sn]] <- sic(HH=HH[[sn]]$rt[HH[[sn]]$x==1], HL=HL[[sn]]$rt[HL[[sn]]$x==1], 
                    LH=LH[[sn]]$rt[LH[[sn]]$x==1], LL=LL[[sn]]$rt[LL[[sn]]$x==1])
  
  #Plots & Saves the DDM survivor functions for each DFP trial type (HH, HL, LH, LL) & SIC Function
  #png(paste("Serial_AND Subject ", sn, ".png"), 640,320)
  postscript("Psi_serial_and.eps", horizontal=FALSE, width = 6, height = 4)
  
  par(mfrow=c(1,2), oma=c(0,0,1,0))
  plot(tvec.dense, 1-ecdf(HH[[sn]]$rt[HH[[sn]]$x==1])(tvec.dense), type='l', col='red',
       xlab="Time (s)", ylab="S(t)", xlim=c(0,5))#max(LL[[sn]]$rt)))
  lines(tvec.dense, 1-ecdf(HL[[sn]]$rt[HL[[sn]]$x==1])(tvec.dense), col='orange')
  lines(tvec.dense, 1-ecdf(LH[[sn]]$rt[LH[[sn]]$x==1])(tvec.dense), col='purple')
  lines(tvec.dense, 1-ecdf(LL[[sn]]$rt[LL[[sn]]$x==1])(tvec.dense), col='blue')
  legend("topright", c("HH", "HL", "LH", "LL"), lty=1, 
         col=c("red", "orange", "purple", "blue"))
  
  plot(tvec.dense, s.and[[sn]]$SIC(tvec.dense), type='l', ylim=c(-.5, .5),
       xlab="Time (s)", ylab="SIC", xlim=c(0,5))#max(LL[[sn]]$rt)))
  title("Serial AND", outer=TRUE, line=-2, cex=1.5)
  abline(h=0)
  dev.off()
  
  print(sn)
}

#PsiParams = cbind("Subject" = rep(1, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[1]]$alpha, "Color_Beta" = result.color[[1]]$beta, "Orientation_Alpha" = result.orientation[[1]]$alpha, "Orientation_Beta" = result.orientation[[1]]$beta)
dominance = cbind("Subject"=rep(1, 8), s.and[[1]]$Dominance)
SFTresults = cbind("Subject"= 1,  "Threshold" = allpars[1,1], "v" = allpars[1,2], "ter" = allpars[1,3], "sdv" = allpars[1,4], "H_Color Intensity" = ddm.result.color[[1]]$intensity[1],"L_Color Intensity" = ddm.result.color[[1]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[1]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[1]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[1]]$rescaled.intensity[2],  "D+ Statistic" = s.and[[1]]$SICtest$positive$statistic, "D+ Pvalue" = s.and[[1]]$SICtest$positive$p.value,  "D- Statistic" = s.and[[1]]$SICtest$negative$statistic, "D- Pvalue" = s.and[[1]]$SICtest$negative$p.value, "MIC Statistic" = s.and[[1]]$MICtest$statistic, "MIC Pvalue" = s.and[[1]]$MICtest$p.value, "HH_RT" = mean(HH[[1]]$rt[1:dfptrials]), "HL_RT" = mean(HL[[1]]$rt[1:dfptrials]), "LH_RT" = mean(LH[[1]]$rt[1:dfptrials]), "LL_RT" = mean(LL[[1]]$rt[1:dfptrials]), "HH_Correct" = mean(HH[[1]]$x[1:dfptrials]) , "HL_Correct" = mean(HL[[1]]$x[1:dfptrials]), "LH_Correct" = mean(LH[[1]]$x[1:dfptrials]), "LL_Correct" = mean(LL[[1]]$x[1:dfptrials]))
DDM.DFPperformance = cbind(cbind("Subject" = rep(1, dfptrials), "Trial" = 1:dfptrials), "HH_RT" = HH[[1]]$rt[1:dfptrials], "HH_Correct" = HH[[1]]$x[1:dfptrials], "HL_RT" = HL[[1]]$rt[1:dfptrials], "HL_Correct" = HL[[1]]$x[1:dfptrials],  "LH_RT" = LH[[1]]$rt[1:dfptrials], "LH_Correct" = LH[[1]]$x[1:dfptrials],"LL_RT" = LL[[1]]$rt[1:dfptrials], "LL_Correct" = LL[[1]]$x[1:dfptrials])

for (sn in 2:nParticipants) { 
  # PsiParams = rbind(PsiParams, cbind("Subject" = rep(sn, nTrials), "Trial" = 1:nTrials, "Color_Alpha" = result.color[[sn]]$alpha, "Color_Beta" = result.color[[sn]]$beta, "Orientation_Alpha" = result.orientation[[sn]]$alpha, "Orientation_Beta" = result.orientation[[sn]]$beta))
  SFTresults = rbind(SFTresults, cbind("Subject"= sn,  "Threshold" = allpars[sn,1], "v" = allpars[sn,2], "ter" = allpars[sn,3], "sdv" = allpars[sn,4], "H_Color Intensity" = ddm.result.color[[sn]]$intensity[1],"L_Color Intensity" = ddm.result.color[[sn]]$intensity[2], "H_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[1],"L_Orientation Intensity" = ddm.result.orientation[[sn]]$intensity[2], "H_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[1],"L_Color Rescaled Intensity" = ddm.result.color[[sn]]$rescaled.intensity[2], "H_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[1],"L_Orientation Rescaled Intensity" = ddm.result.orientation[[sn]]$rescaled.intensity[2],  "D+ Statistic" = s.and[[sn]]$SICtest$positive$statistic, "D+ Pvalue" = s.and[[sn]]$SICtest$positive$p.value,  "D- Statistic" = s.and[[sn]]$SICtest$negative$statistic, "D- Pvalue" = s.and[[sn]]$SICtest$negative$p.value, "MIC Statistic" = s.and[[sn]]$MICtest$statistic, "MIC Pvalue" = s.and[[sn]]$MICtest$p.value,  "HH_RT" = mean(HH[[sn]]$rt[1:dfptrials]), "HL_RT" = mean(HL[[sn]]$rt[1:dfptrials]), "LH_RT" = mean(LH[[sn]]$rt[1:dfptrials]), "LL_RT" = mean(LL[[sn]]$rt[1:dfptrials]), "HH_Correct" = mean(HH[[sn]]$x[1:dfptrials]) , "HL_Correct" = mean(HL[[sn]]$x[1:dfptrials]), "LH_Correct" = mean(LH[[sn]]$x[1:dfptrials]), "LL_Correct" = mean(LL[[sn]]$x[1:dfptrials])))
  dominance = rbind(dominance, cbind("Subject"=rep(sn, 8), s.or[[sn]]$Dominance))
  DDM.DFPperformance = rbind(DDM.DFPperformance, cbind("Subject" = rep(sn, dfptrials), "Trial" = 1:dfptrials, "HH_RT" = HH[[sn]]$rt[1:dfptrials], "HH_Correct" = HH[[1]]$x[1:dfptrials], "HL_RT" = HL[[1]]$rt[1:dfptrials], "HL_Correct" = HL[[sn]]$x[1:dfptrials],  "LH_RT" = LH[[sn]]$rt[1:dfptrials], "LH_Correct" = LH[[sn]]$x[1:dfptrials],"LL_RT" = LL[[sn]]$rt[1:dfptrials], "LL_Correct" = LL[[sn]]$x[1:dfptrials]))
}

write.csv(dominance, file="SerialAND_Psi_Simulation_KStests.csv", row.names=FALSE)
write.csv(SFTresults, file="SerialAND_Psi_Simulation_SFTresults.csv", row.names=FALSE)
#write.csv(PsiParams, file="SerialAND_PsiDDM_Simulation_Pars.csv", row.names=FALSE)
write.csv(DDM.DFPperformance, file="SerialAND_Psi_Simulation_DDM-DFP_RT&ACC.csv", row.names=FALSE)

