
################################################################################
# This file contains the code for fitting the quantile curves reported in:     #
# 'Quantile regression for longitudinal within-race running data: the 2022 New #
# York City Marathon' - (submitted 2025) - D'Haen, Flórez, Molenberghs,        #
# Van Keilegom, Delecluse, Verhasselt.                                         #
#                                                                              #
# Code authors: Alvaro Flórez & Myrthe D'Haen                                  #
# Last revised on: 16/10/2025                                                  #
#                                                                              #
################################################################################



#### Part 0: initialisation (folders, functions, packages) #####################

# specify the correct working directory; the master directory should contain
# (i) NYClong.csv, (ii) NYCMarQu-main.R, (iii) NYCMarQu-functions.R
master.dir <- 'G:/My Drive/Marathon project'
setwd(master.dir)


MyCoreNum = 4 # specify your number of cores to be used for parallel computing


# loading all necessary functions and libraries
source(paste(master.dir, 'NYCMarQu-functions.R', sep = '/'))

library(quantreg)
library(rlist)
library(xtable)
library(ggplot2)
library(pbapply)
library(parallel)
library(maxLik)
library(gtools)
library(ks)
library(matrixStats)

library(patchwork)
library(RColorBrewer)
library(dplyr)


# folders for modelling output and plotting
output.dir <- paste(master.dir, 'Output files', sep = '/')
unstr.dir <- paste(output.dir, 'Unstructured models', sep = '/')
cubic.dir <- paste(output.dir, 'Cubic models', sep = '/')
uqr.dir <- paste(output.dir, 'UQR models', sep = '/')
plot.dir <- paste(output.dir, 'Model plots', sep = '/')

for (folder in c(output.dir, unstr.dir, cubic.dir, uqr.dir, plot.dir)){
  if (!dir.exists(folder)){
    dir.create(folder)
  }
}

################################################################################



##### Part 1: preparing the data ###############################################

### information on (scaled) distances, altitudes, slopes, etc.

Dist.vec <- c(5,10,15,20,21.0975,25,30,32.1869,35,40,42.195)

# extended with distance point 0
Dist.vec.ext <- c(0, Dist.vec)
alt.vec.ext <- c(40, 21, 15, 14, 5, 18, 44, 3, 11, 8, 24, 29)

# scaled distance vectors (d.scaled = D/42.195 in [0,1])
dist.scaled <- Dist.vec/42.195
dist.scaled.noD1 <- dist.scaled[-1] # omitting 5k point
dist.scaled.noDM <- dist.scaled[-11] # omitting MAR point

NYC.long = read.csv('NYClong.csv', header = T)

NYC.long$AgeGroup = as.factor(NYC.long$AgeGroup)
NYC.long$AgeScaled = (NYC.long$Age - 18)/10
NYC.long$Mtime = rep(NYC.long$time[NYC.long$distance=="MAR"],each=11)
NYC.long$dist = NYC.long$dis/42.195
NYC.long$revdist = 1-NYC.long$dist # for CCM model (splitVec)
NYC.long$Slope = diff(alt.vec.ext)/diff(Dist.vec.ext)
NYC.long$revSlope = NYC.long$Slope[11]-NYC.long$Slope


NYC.long$intSpeed = unlist(by(NYC.long,NYC.long$bib,IntervalSpeed))
NYC.long$cumSpeed = unlist(by(NYC.long,NYC.long$bib,CumulativeSpeed))
NYC.long$accel = unlist(by(NYC.long, NYC.long$bib, Acceleration))
NYC.long$speedDec = unlist(by(NYC.long,NYC.long$bib,SpeedDecline))
NYC.long$splitVec = unlist(by(NYC.long,NYC.long$bib,SplittingVector))



### Data with omission of some distance points ###
NYC.noD1 <- NYC.long[!(NYC.long$dis == 5),] # omit 5k point
NYC.noDM <- NYC.long[!(NYC.long$dis == 42.195),] # omit MAR point

################################################################################



##### Part 2: Copula-based PWE #################################################

# Note: the code for the PWE estimator was mostly written by Alvaro Flórez in
# the context of the paper
# Verhasselt, Flórez, Molenberghs, and Van Keilegom (2025) - Copula-based
# pairwise estimator for quantile regression with hierarchical missing data.
# Statistical Modelling, 25(2):129-149.


# For each variable, we fit the two PWE model options described in the paper:
#   option 1 = unstructured model (UM)
#   option 2 = cubic model (CM)
# Only for SplitVec, a third option is also considered:
#   option 3 = constrained cubic model (CCM) ("constrained" = 
#               "ensuring that splitting equals 1 at marathon distance")

# Each time, the quantile of level 0.1 is computed first. For the other quantile
# levels, the parameters of the previous quantile are used as initial values.

# We consider a version with less interaction terms (used for first submission
# on 28/1/2025), as well as a version with additional sex-age and age-distance
# interaction terms (used for submission on 17/10/2025, after the first review).

# (The value for epsilon can each time be modified.)

# All results are stored externally as .Rdata files, to enable running the code
# in chunks (since the procedure is highly time consuming).
# Only afterwards, all results are again put together and collected in tables
# and graphs (Part 4).



##### Part 2.1: fitting the PWE - unstructured model ###########################

### Choice of epsilon
eps.val <- 0.1

### Preparing the arguments for the pwe function
ID.full <- rep(1:(nrow(NYC.long)/11),each=11)
ID.noD1 <- rep(1:(nrow(NYC.noD1)/10),each=10)
ID.noDM <- rep(1:(nrow(NYC.noDM)/10),each=10)

# ## obtaining the model matrix X - first submission (28/01/2025)
# # (use arbitrary response variable for the lm)
# lm.full <- lm(intSpeed ~ Sex*as.factor(dis)+AgeScaled, data=NYC.long)
# lm.noD1 <- lm(intSpeed ~ Sex*as.factor(dis)+AgeScaled, data=NYC.noD1)
# lm.noDM <- lm(intSpeed ~ Sex*as.factor(dis)+AgeScaled, data=NYC.noDM)
# 
# X.full <- as.matrix(model.matrix(lm.full))
# X.noD1 <- as.matrix(model.matrix(lm.noD1))
# X.noDM <- as.matrix(model.matrix(lm.noDM))
# 
# colnames(X.full)[3:12] <- NYC.long$distance[2:11]
# colnames(X.full)[14:23] = paste('sexM.',NYC.long$distance[2:11],sep='')
# 
# colnames(X.noD1)[3:11] <- NYC.noD1$distance[2:10]
# colnames(X.noD1)[13:21] = paste('sexM.',NYC.noD1$distance[2:10],sep='')
# 
# colnames(X.noDM)[3:11] <- NYC.noDM$distance[2:10]
# colnames(X.noDM)[13:21] = paste('sexM.',NYC.noDM$distance[2:10],sep='')

## obtaining the model matrix X - submission after review (17/10/2025)
# (use arbitrary response variable for the lm. Intercept automatically included)
lm.full <- lm(intSpeed ~ Sex*AgeScaled + Sex*as.factor(dis) +
                AgeScaled*as.factor(dis), data = NYC.long)
lm.noD1 <- lm(intSpeed ~ Sex*AgeScaled + Sex*as.factor(dis)+
                AgeScaled*as.factor(dis), data = NYC.noD1)
lm.noDM <- lm(intSpeed ~ Sex*AgeScaled + Sex*as.factor(dis)+
                AgeScaled*as.factor(dis), data = NYC.noDM)

X.full <- as.matrix(model.matrix(lm.full))
X.noD1 <- as.matrix(model.matrix(lm.noD1))
X.noDM <- as.matrix(model.matrix(lm.noDM))

colnames(X.full)[3] <- 'AgeSc'
colnames(X.full)[4:13] <- NYC.long$distance[2:11]
colnames(X.full)[14] <- 'SexM.AgeSc'
colnames(X.full)[15:24] = paste('sexM.',NYC.long$distance[2:11],sep='')
colnames(X.full)[25:34] = paste('AgeSc.',NYC.long$distance[2:11],sep='')

colnames(X.noD1)[3] <- 'AgeSc'
colnames(X.noD1)[4:12] <- NYC.noD1$distance[2:10]
colnames(X.noD1)[13] <- 'SexM.AgeSc'
colnames(X.noD1)[14:22] = paste('sexM.',NYC.noD1$distance[2:10],sep='')
colnames(X.noD1)[23:31] = paste('AgeSc.',NYC.noD1$distance[2:10],sep='')


colnames(X.noDM)[3] <- 'AgeSc'
colnames(X.noDM)[4:12] <- NYC.noDM$distance[2:10]
colnames(X.noDM)[13] <- 'SexM.AgeSc'
colnames(X.noDM)[14:22] = paste('sexM.',NYC.noDM$distance[2:10],sep='')
colnames(X.noDM)[23:31] = paste('AgeSc.',NYC.noDM$distance[2:10],sep='')



### intSpeed (UM) ### ----------------------------------------------------------
est.10.UM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$intSpeed, tau=0.10,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = NULL,
                                         ncores = MyCoreNum)
list.save(est.10.UM.intSpeed,
          paste0(unstr.dir,paste0('/intSpeed.pwe.10.UM.eps.',eps.val,'.Rdata')))

# # Note: if you want to continue in another R session, use the following lines of
# # code at the appropriate place, with correct values for tau in the names:
# est.10.UM.intSpeed <- list.load(
#   paste0(unstr.dir, paste0('/intSpeed.pwe.10.UM.eps.',eps.val,'.Rdata')))


est.25.UM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$intSpeed, tau=0.25,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.10.UM.intSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.25.UM.intSpeed,
          paste0(unstr.dir,paste0('/intSpeed.pwe.25.UM.eps.',eps.val,'.Rdata')))

# est.25.UM.intSpeed <- list.load(
#   paste0(unstr.dir, paste0('/intSpeed.pwe.25.UM.eps.',eps.val,'.Rdata')))
# # etc., similar for variables and tau values below

est.50.UM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$intSpeed, tau=0.50,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.25.UM.intSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.50.UM.intSpeed,
          paste0(unstr.dir,paste0('/intSpeed.pwe.50.UM.eps.',eps.val,'.Rdata')))

est.75.UM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$intSpeed, tau=0.75,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.50.UM.intSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.75.UM.intSpeed,
          paste0(unstr.dir,paste0('/intSpeed.pwe.75.UM.eps.',eps.val,'.Rdata')))

est.90.UM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$intSpeed, tau=0.90,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.75.UM.intSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.90.UM.intSpeed,
          paste0(unstr.dir,paste0('/intSpeed.pwe.90.UM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------

### cumSpeed (UM) ### ----------------------------------------------------------
est.10.UM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$cumSpeed, tau=0.10,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = NULL,
                                         ncores = MyCoreNum)
list.save(est.10.UM.cumSpeed,
          paste0(unstr.dir,paste0('/cumSpeed.pwe.10.UM.eps.',eps.val,'.Rdata')))

est.25.UM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$cumSpeed, tau=0.25,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.10.UM.cumSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.25.UM.cumSpeed,
          paste0(unstr.dir,paste0('/cumSpeed.pwe.25.UM.eps.',eps.val,'.Rdata')))

est.50.UM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$cumSpeed, tau=0.50,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.25.UM.cumSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.50.UM.cumSpeed,
          paste0(unstr.dir,paste0('/cumSpeed.pwe.50.UM.eps.',eps.val,'.Rdata')))

est.75.UM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$cumSpeed, tau=0.75,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.50.UM.cumSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.75.UM.cumSpeed,
          paste0(unstr.dir,paste0('/cumSpeed.pwe.75.UM.eps.',eps.val,'.Rdata')))

est.90.UM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = X.full,
                                         y = NYC.long$cumSpeed, tau=0.90,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.75.UM.cumSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.90.UM.cumSpeed,
          paste0(unstr.dir,paste0('/cumSpeed.pwe.90.UM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------

### accel (UM) ### ----------------------------------------------------------
est.10.UM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                         y = NYC.noD1$accel, tau=0.10,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = NULL,
                                         ncores = MyCoreNum)
list.save(est.10.UM.accel,
          paste0(unstr.dir,paste0('/accel.pwe.10.UM.eps.',eps.val,'.Rdata')))

est.25.UM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                         y = NYC.noD1$accel, tau=0.25,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.10.UM.accel$parm,
                                         ncores = MyCoreNum)
list.save(est.25.UM.accel,
          paste0(unstr.dir,paste0('/accel.pwe.25.UM.eps.',eps.val,'.Rdata')))

est.50.UM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                         y = NYC.noD1$accel, tau=0.50,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.25.UM.accel$parm,
                                         ncores = MyCoreNum)
list.save(est.50.UM.accel,
          paste0(unstr.dir,paste0('/accel.pwe.50.UM.eps.',eps.val,'.Rdata')))

est.75.UM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                         y = NYC.noD1$accel, tau=0.75,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.50.UM.accel$parm,
                                         ncores = MyCoreNum)
list.save(est.75.UM.accel,
          paste0(unstr.dir,paste0('/accel.pwe.75.UM.eps.',eps.val,'.Rdata')))

est.90.UM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                         y = NYC.noD1$accel, tau=0.90,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.75.UM.accel$parm,
                                         ncores = MyCoreNum)
list.save(est.90.UM.accel,
          paste0(unstr.dir,paste0('/accel.pwe.90.UM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------

### speedDec (UM) ### ----------------------------------------------------------
est.10.UM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                      y = NYC.noD1$speedDec, tau=0.10,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = NULL,
                                      ncores = MyCoreNum)
list.save(est.10.UM.speedDec,
          paste0(unstr.dir,paste0('/speedDec.pwe.10.UM.eps.',eps.val,'.Rdata')))

est.25.UM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                      y = NYC.noD1$speedDec, tau=0.25,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = est.10.UM.speedDec$parm,
                                      ncores = MyCoreNum)
list.save(est.25.UM.speedDec,
          paste0(unstr.dir,paste0('/speedDec.pwe.25.UM.eps.',eps.val,'.Rdata')))

est.50.UM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                      y = NYC.noD1$speedDec, tau=0.50,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = est.25.UM.speedDec$parm,
                                      ncores = MyCoreNum)
list.save(est.50.UM.speedDec,
          paste0(unstr.dir,paste0('/speedDec.pwe.50.UM.eps.',eps.val,'.Rdata')))

est.75.UM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                      y = NYC.noD1$speedDec, tau=0.75,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = est.50.UM.speedDec$parm,
                                      ncores = MyCoreNum)
list.save(est.75.UM.speedDec,
          paste0(unstr.dir,paste0('/speedDec.pwe.75.UM.eps.',eps.val,'.Rdata')))

est.90.UM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = X.noD1,
                                      y = NYC.noD1$speedDec, tau=0.90,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = est.75.UM.speedDec$parm,
                                      ncores = MyCoreNum)
list.save(est.90.UM.speedDec,
          paste0(unstr.dir,paste0('/speedDec.pwe.90.UM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------

### splitVec (UM) ### ----------------------------------------------------------
est.10.UM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = X.noDM,
                                         y = NYC.noDM$splitVec, tau=0.10,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = NULL,
                                         ncores = MyCoreNum)
list.save(est.10.UM.splitVec,
          paste0(unstr.dir,paste0('/splitVec.pwe.10.UM.eps.',eps.val,'.Rdata')))

est.25.UM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = X.noDM,
                                         y = NYC.noDM$splitVec, tau=0.25,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.10.UM.splitVec$parm,
                                         ncores = MyCoreNum)
list.save(est.25.UM.splitVec,
          paste0(unstr.dir,paste0('/splitVec.pwe.25.UM.eps.',eps.val,'.Rdata')))

est.50.UM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = X.noDM,
                                         y = NYC.noDM$splitVec, tau=0.50,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.25.UM.splitVec$parm,
                                         ncores = MyCoreNum)
list.save(est.50.UM.splitVec,
          paste0(unstr.dir,paste0('/splitVec.pwe.50.UM.eps.',eps.val,'.Rdata')))

est.75.UM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = X.noDM,
                                         y = NYC.noDM$splitVec, tau=0.75,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.50.UM.splitVec$parm,
                                         ncores = MyCoreNum)
list.save(est.75.UM.splitVec,
          paste0(unstr.dir,paste0('/splitVec.pwe.75.UM.eps.',eps.val,'.Rdata')))

# est.75.UM.splitVec <- list.load(
#   paste0(unstr.dir, paste0('/splitVec.pwe.75.UM.eps.',eps.val,'.Rdata')))



est.90.UM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = X.noDM,
                                         y = NYC.noDM$splitVec, tau=0.90,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.75.UM.splitVec$parm,
                                         ncores = MyCoreNum)
list.save(est.90.UM.splitVec,
          paste0(unstr.dir,paste0('/splitVec.pwe.90.UM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------

################################################################################


##### Part 2.2: fitting the PWE - cubic model ##################################

### Choice of epsilon
eps.val <- 0.1

### Preparing the arguments for the pwe function
ID.full <- rep(1:(nrow(NYC.long)/11),each=11)
ID.noD1 <- rep(1:(nrow(NYC.noD1)/10),each=10)
ID.noDM <- rep(1:(nrow(NYC.noDM)/10),each=10)

## obtaining the model matrix X
# full data
AgeSc.full <- NYC.long$AgeScaled
Sex.full = as.double(NYC.long$Sex == 'M')
dist.full = NYC.long$dist
Slope.full = NYC.long$Slope

# omitting 5k point
AgeSc.noD1 <- NYC.noD1$AgeScaled
Sex.noD1 = as.double(NYC.noD1$Sex == 'M')
dist.noD1 = NYC.noD1$dist
Slope.noD1 = NYC.noD1$Slope

# omitting MAR point - note: different variables, since CCM vs. CM!
AgeSc.noDM <- NYC.noDM$AgeScaled
Sex.noDM = as.double(NYC.noDM$Sex == 'M')
revdist = NYC.noDM$revdist
revSlope = NYC.noDM$revSlope


# ## model matrix X - first submission (28/01/2025)
# Xcub.full = cbind(1,Sex.full,dist.full,dist.full^2,dist.full^3,
#                   dist.full*Sex.full,dist.full^2*Sex.full,dist.full^3*Sex.full,
#                   AgeSc.full,Slope.full)
# colnames(Xcub.full) = c('int','sex','d','d2','d3','ds','d2s','d3s','age','slope')
# 
# Xcub.noD1 = cbind(1,Sex.noD1,dist.noD1,dist.noD1^2,dist.noD1^3,
#                   dist.noD1*Sex.noD1,dist.noD1^2*Sex.noD1,dist.noD1^3*Sex.noD1,
#                   AgeSc.noD1,Slope.noD1)
# colnames(Xcub.noD1) = c('int','sex','d','d2','d3','ds','d2s','d3s','age','slope')
# 
# Xccub.noDM = cbind(revSlope, revdist, revdist^2, revdist^3,
#                    revdist*Sex.noDM, revdist^2*Sex.noDM, revdist^3*Sex.noDM,
#                    revdist*AgeSc.noDM)
# colnames(Xccub.noDM) = c('rSlope','rd','rd2','rd3','rds','rd2s','rd3s','rdAge')


## model matrix X - submission after review (17/10/2025)
Xcub.full = cbind(1,Sex.full,dist.full,dist.full^2,dist.full^3,
                  dist.full*Sex.full,dist.full^2*Sex.full,dist.full^3*Sex.full,
                  AgeSc.full,Slope.full,Sex.full*AgeSc.full,
                  dist.full*AgeSc.full,dist.full^2*AgeSc.full,dist.full^3*AgeSc.full)
colnames(Xcub.full) = c('int','sex','d','d2','d3','ds','d2s','d3s','age','slope',
                        'sAge', 'dAge', 'd2Age', 'd3Age')

Xcub.noD1 = cbind(1,Sex.noD1,dist.noD1,dist.noD1^2,dist.noD1^3,
                  dist.noD1*Sex.noD1,dist.noD1^2*Sex.noD1,dist.noD1^3*Sex.noD1,
                  AgeSc.noD1,Slope.noD1,Sex.noD1*AgeSc.noD1,
                  dist.noD1*AgeSc.noD1,dist.noD1^2*AgeSc.noD1,dist.noD1^3*AgeSc.noD1)
colnames(Xcub.noD1) = c('int','sex','d','d2','d3','ds','d2s','d3s','age','slope',
                        'sAge', 'dAge', 'd2Age', 'd3Age')

Xccub.noDM = cbind(revSlope, revdist, revdist^2, revdist^3,
                   revdist*Sex.noDM, revdist^2*Sex.noDM, revdist^3*Sex.noDM,
                   revdist*AgeSc.noDM, revdist^2*AgeSc.noDM, revdist^3*AgeSc.noDM,
                   revdist*Sex.noDM*AgeSc.noDM)
colnames(Xccub.noDM) = c('rSlope','rd','rd2','rd3','rds','rd2s','rd3s','rdAge',
                         'rd2Age', 'rd3Age', 'rdsAge')



### intSpeed (CM) ### ----------------------------------------------------------
est.10.CM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$intSpeed, tau=0.10,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = NULL,
                                         ncores = MyCoreNum)
list.save(est.10.CM.intSpeed,
          paste0(cubic.dir,paste0('/intSpeed.pwe.10.CM.eps.',eps.val,'.Rdata')))

# # Note: if you want to continue in another R session, use the following lines of
# # code at the appropriate place, with correct values for tau in the names:
# est.10.CM.intSpeed <- list.load(
#   paste0(cubic.dir, paste0('/intSpeed.pwe.10.CM.eps.',eps.val,'.Rdata')))


est.25.CM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$intSpeed, tau=0.25,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.10.CM.intSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.25.CM.intSpeed,
          paste0(cubic.dir,paste0('/intSpeed.pwe.25.CM.eps.',eps.val,'.Rdata')))

# est.25.CM.intSpeed <- list.load(
#   paste0(cubic.dir, paste0('/intSpeed.pwe.25.CM.eps.',eps.val,'.Rdata')))
# # etc., similar for variables and tau values below

est.50.CM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$intSpeed, tau=0.50,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.25.CM.intSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.50.CM.intSpeed,
          paste0(cubic.dir,paste0('/intSpeed.pwe.50.CM.eps.',eps.val,'.Rdata')))

est.75.CM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$intSpeed, tau=0.75,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.50.CM.intSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.75.CM.intSpeed,
          paste0(cubic.dir,paste0('/intSpeed.pwe.75.CM.eps.',eps.val,'.Rdata')))

est.90.CM.intSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$intSpeed, tau=0.90,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.75.CM.intSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.90.CM.intSpeed,
          paste0(cubic.dir,paste0('/intSpeed.pwe.90.CM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------

### cumSpeed (CM) ### ----------------------------------------------------------
est.10.CM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$cumSpeed, tau=0.10,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = NULL,
                                         ncores = MyCoreNum)
list.save(est.10.CM.cumSpeed,
          paste0(cubic.dir,paste0('/cumSpeed.pwe.10.CM.eps.',eps.val,'.Rdata')))

est.25.CM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$cumSpeed, tau=0.25,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.10.CM.cumSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.25.CM.cumSpeed,
          paste0(cubic.dir,paste0('/cumSpeed.pwe.25.CM.eps.',eps.val,'.Rdata')))

est.50.CM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$cumSpeed, tau=0.50,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.25.CM.cumSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.50.CM.cumSpeed,
          paste0(cubic.dir,paste0('/cumSpeed.pwe.50.CM.eps.',eps.val,'.Rdata')))

est.75.CM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$cumSpeed, tau=0.75,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.50.CM.cumSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.75.CM.cumSpeed,
          paste0(cubic.dir,paste0('/cumSpeed.pwe.75.CM.eps.',eps.val,'.Rdata')))

est.90.CM.cumSpeed = pwe.ALDcop.parallel(ID = ID.full, X = Xcub.full,
                                         y = NYC.long$cumSpeed, tau=0.90,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.75.CM.cumSpeed$parm,
                                         ncores = MyCoreNum)
list.save(est.90.CM.cumSpeed,
          paste0(cubic.dir,paste0('/cumSpeed.pwe.90.CM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------

### accel (CM) ### ----------------------------------------------------------
est.10.CM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                      y = NYC.noD1$accel, tau=0.10,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = NULL,
                                      ncores = MyCoreNum)
list.save(est.10.CM.accel,
          paste0(cubic.dir,paste0('/accel.pwe.10.CM.eps.',eps.val,'.Rdata')))

est.25.CM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                      y = NYC.noD1$accel, tau=0.25,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = est.10.CM.accel$parm,
                                      ncores = MyCoreNum)
list.save(est.25.CM.accel,
          paste0(cubic.dir,paste0('/accel.pwe.25.CM.eps.',eps.val,'.Rdata')))

est.50.CM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                      y = NYC.noD1$accel, tau=0.50,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = est.25.CM.accel$parm,
                                      ncores = MyCoreNum)
list.save(est.50.CM.accel,
          paste0(cubic.dir,paste0('/accel.pwe.50.CM.eps.',eps.val,'.Rdata')))

est.75.CM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                      y = NYC.noD1$accel, tau=0.75,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = est.50.CM.accel$parm,
                                      ncores = MyCoreNum)
list.save(est.75.CM.accel,
          paste0(cubic.dir,paste0('/accel.pwe.75.CM.eps.',eps.val,'.Rdata')))

est.90.CM.accel = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                      y = NYC.noD1$accel, tau=0.90,
                                      EPS=eps.val, Rstr='AR',
                                      mtimes=dist.scaled.noD1,
                                      miss.data.method= 'AC', var.estimate=T,
                                      parm.ini = est.75.CM.accel$parm,
                                      ncores = MyCoreNum)
list.save(est.90.CM.accel,
          paste0(cubic.dir,paste0('/accel.pwe.90.CM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------

### speedDec (CM) ### ----------------------------------------------------------
est.10.CM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                         y = NYC.noD1$speedDec, tau=0.10,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = NULL,
                                         ncores = MyCoreNum)
list.save(est.10.CM.speedDec,
          paste0(cubic.dir,paste0('/speedDec.pwe.10.CM.eps.',eps.val,'.Rdata')))

est.25.CM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                         y = NYC.noD1$speedDec, tau=0.25,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.10.CM.speedDec$parm,
                                         ncores = MyCoreNum)
list.save(est.25.CM.speedDec,
          paste0(cubic.dir,paste0('/speedDec.pwe.25.CM.eps.',eps.val,'.Rdata')))

est.50.CM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                         y = NYC.noD1$speedDec, tau=0.50,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.25.CM.speedDec$parm,
                                         ncores = MyCoreNum)
list.save(est.50.CM.speedDec,
          paste0(cubic.dir,paste0('/speedDec.pwe.50.CM.eps.',eps.val,'.Rdata')))

est.75.CM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                         y = NYC.noD1$speedDec, tau=0.75,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.50.CM.speedDec$parm,
                                         ncores = MyCoreNum)
list.save(est.75.CM.speedDec,
          paste0(cubic.dir,paste0('/speedDec.pwe.75.CM.eps.',eps.val,'.Rdata')))

est.90.CM.speedDec = pwe.ALDcop.parallel(ID = ID.noD1, X = Xcub.noD1,
                                         y = NYC.noD1$speedDec, tau=0.90,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noD1,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.75.CM.speedDec$parm,
                                         ncores = MyCoreNum)
list.save(est.90.CM.speedDec,
          paste0(cubic.dir,paste0('/speedDec.pwe.90.CM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------

### splitVec (CCM) ### ---------------------------------------------------------
## Note: y = NYC.noDM$splitVec - 1 is modelled here, including shift -1!
est.10.CCM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = Xccub.noDM,
                                         y = NYC.noDM$splitVec - 1, tau=0.10,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = NULL,
                                         ncores = MyCoreNum)
list.save(est.10.CCM.splitVec,
          paste0(cubic.dir,paste0('/splitVec.pwe.10.CCM.eps.',eps.val,'.Rdata')))

est.25.CCM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = Xccub.noDM,
                                         y = NYC.noDM$splitVec - 1, tau=0.25,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.10.CCM.splitVec$parm,
                                         ncores = MyCoreNum)
list.save(est.25.CCM.splitVec,
          paste0(cubic.dir,paste0('/splitVec.pwe.25.CCM.eps.',eps.val,'.Rdata')))

est.50.CCM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = Xccub.noDM,
                                         y = NYC.noDM$splitVec - 1, tau=0.50,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.25.CCM.splitVec$parm,
                                         ncores = MyCoreNum)
list.save(est.50.CCM.splitVec,
          paste0(cubic.dir,paste0('/splitVec.pwe.50.CCM.eps.',eps.val,'.Rdata')))

est.75.CCM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = Xccub.noDM,
                                         y = NYC.noDM$splitVec - 1, tau=0.75,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.50.CCM.splitVec$parm,
                                         ncores = MyCoreNum)
list.save(est.75.CCM.splitVec,
          paste0(cubic.dir,paste0('/splitVec.pwe.75.CCM.eps.',eps.val,'.Rdata')))

est.90.CCM.splitVec = pwe.ALDcop.parallel(ID = ID.noDM, X = Xccub.noDM,
                                         y = NYC.noDM$splitVec - 1, tau=0.90,
                                         EPS=eps.val, Rstr='AR',
                                         mtimes=dist.scaled.noDM,
                                         miss.data.method= 'AC', var.estimate=T,
                                         parm.ini = est.75.CCM.splitVec$parm,
                                         ncores = MyCoreNum)
list.save(est.90.CCM.splitVec,
          paste0(cubic.dir,paste0('/splitVec.pwe.90.CCM.eps.',eps.val,'.Rdata')))
# ------------------------------------------------------------------------------


##### Part 2.3: example of t-copula instead of Gaussian; all code above (parts
# 2.1 and 2.2) can be modified accordingly.
est.10.UM.intSpeed.tCop =
  pwe.ALDcop.t.parallel(ID = ID.full,X = X.full , y = NYC.long$intSpeed,
                        tau=0.5,v=4,EPS=eps.val,parm.ini=NULL,
                        weighted=F,weights=NULL,Rstr='AR',
                        mtimes=dist.scaled, miss.data.method='AC',var.estimate=T)

################################################################################



##### Part 3: fitting the UQR for comparison ###################################

## Before review: without Sex * AgeSc interaction term -------------------------
uqr.SexAge.intSpeed <- lapply(1:11, function(index){
  summary(rq(intSpeed ~ Sex + AgeScaled, data = NYC.long,
             subset = (dist == dist.scaled[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

uqr.SexAge.cumSpeed <- lapply(1:11, function(index){
  summary(rq(cumSpeed ~ Sex + AgeScaled, data = NYC.long,
             subset = (dist == dist.scaled[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

uqr.SexAge.accel <- lapply(1:10, function(index){
  summary(rq(accel ~ Sex + AgeScaled, data = NYC.noD1,
             subset = (dist == dist.scaled.noD1[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

uqr.SexAge.speedDec <- lapply(1:10, function(index){
  summary(rq(speedDec ~ Sex + AgeScaled, data = NYC.noD1,
             subset = (dist == dist.scaled.noD1[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

uqr.SexAge.splitVec <- lapply(1:10, function(index){
  summary(rq(splitVec ~ Sex + AgeScaled, data = NYC.noDM,
             subset = (dist == dist.scaled.noDM[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

# Storing summaries externally in .Rdata file
list.save(uqr.SexAge.intSpeed, paste0(uqr.dir, '/UQR.SexAge.intSpeed.Rdata'))
list.save(uqr.SexAge.cumSpeed, paste0(uqr.dir, '/UQR.SexAge.cumSpeed.Rdata'))
list.save(uqr.SexAge.accel, paste0(uqr.dir, '/UQR.SexAge.accel.Rdata'))
list.save(uqr.SexAge.speedDec, paste0(uqr.dir, '/UQR.SexAge.speedDec.Rdata'))
list.save(uqr.SexAge.splitVec, paste0(uqr.dir, '/UQR.SexAge.splitVec.Rdata'))
# ------------------------------------------------------------------------------

## After review: with Sex * AgeSc interaction term -----------------------------
uqr.SexAgeInter.intSpeed <- lapply(1:11, function(index){
  summary(rq(intSpeed ~ Sex*AgeScaled, data = NYC.long,
             subset = (dist == dist.scaled[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

uqr.SexAgeInter.cumSpeed <- lapply(1:11, function(index){
  summary(rq(cumSpeed ~ Sex*AgeScaled, data = NYC.long,
             subset = (dist == dist.scaled[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

uqr.SexAgeInter.accel <- lapply(1:10, function(index){
  summary(rq(accel ~ Sex*AgeScaled, data = NYC.noD1,
             subset = (dist == dist.scaled.noD1[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

uqr.SexAgeInter.speedDec <- lapply(1:10, function(index){
  summary(rq(speedDec ~ Sex*AgeScaled, data = NYC.noD1,
             subset = (dist == dist.scaled.noD1[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

uqr.SexAgeInter.splitVec <- lapply(1:10, function(index){
  summary(rq(splitVec ~ Sex*AgeScaled, data = NYC.noDM,
             subset = (dist == dist.scaled.noDM[index]), tau = c(0.1,0.25,0.5,0.75,0.9)))
})

# Storing summaries externally in .Rdata file
list.save(uqr.SexAgeInter.intSpeed, paste0(uqr.dir, '/UQR.SexAgeInter.intSpeed.Rdata'))
list.save(uqr.SexAgeInter.cumSpeed, paste0(uqr.dir, '/UQR.SexAgeInter.cumSpeed.Rdata'))
list.save(uqr.SexAgeInter.accel, paste0(uqr.dir, '/UQR.SexAgeInter.accel.Rdata'))
list.save(uqr.SexAgeInter.speedDec, paste0(uqr.dir, '/UQR.SexAgeInter.speedDec.Rdata'))
list.save(uqr.SexAgeInter.splitVec, paste0(uqr.dir, '/UQR.SexAgeInter.splitVec.Rdata'))
# ------------------------------------------------------------------------------

################################################################################



#### Part 4: obtaining coefficient estimates ###################################

# verify that the appropriate values for epsilon and corresponding folder paths
# are specified

# Note: function summary_LQC assumes a model specification belonging to one of the
# options above ('SexAge' before review, 'SexAgeInter' after first review).
# In case another one is to be considered, update summary_LQC by adding another
# model formula option.

# ! Note: due to specificities of the summary.rq function (package quantreg),
# the UQR method returns estimates with their confidence intervals rather than
# standard errors in case the sample size (per measurement time) does not exceed
# 1001. (https://cran.r-project.org/web/packages/quantreg/quantreg.pdf).
# The code below then produces lower bd rather than SE and needs to be modified.

### SexAge models (before review: without interaction) ### ---------------------

# Note: {UM,CM}_infix and suffix are just for compatibility with results from
# older versions of the code, where files were stored under different names.

# intSpeed # -------------------------------------------------------------------
res.intSpeed <- summary_LQC(varName = 'intSpeed',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.01, CM_eps = 0.05,
                            UM_path = paste0(master.dir, "/Old UM models/"),
                            CM_path = paste0(master.dir, "/Old CM models/"),
                            UM_formula = "SexAge", CM_formula = "SexAge",
                            UM_infix = "unstr", CM_infix = "cubic",
                            UM_suffix = "NR.Rdata", CM_suffix = "NR.Rdata")

intSpeed.UM <- res.intSpeed$outUns
intSpeed.CM <- res.intSpeed$outCub

# Unstructured model results
print(xtable(cbind(intSpeed.UM[[1]][,c(2,3,4,7)],NA,intSpeed.UM[[2]][,c(3,4,7)],NA,
                   intSpeed.UM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(intSpeed.UM[[4]][,c(2,3,4,7)],NA,intSpeed.UM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# Cubic model results
print(xtable(cbind(intSpeed.CM[[1]][,c(2,3,4,7)],NA,intSpeed.CM[[2]][,c(3,4,7)],NA,
                   intSpeed.CM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(intSpeed.CM[[4]][,c(2,3,4,7)],NA,intSpeed.CM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# UQR results
uqr.SexAge.intSpeed <- list.load(paste0(uqr.dir,'/UQR.SexAge.intSpeed.Rdata'))
uqr.SexAge.intSpeed.df <- get_all_est_SE(summary_list = uqr.SexAge.intSpeed,
                                  taus = c(0.1,0.25,0.5,0.75,0.9))
xtable(uqr.SexAge.intSpeed.df, digits = 4)
# ------------------------------------------------------------------------------

# cumSpeed # -------------------------------------------------------------------
res.cumSpeed <- summary_LQC(varName = 'cumSpeed',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.01, CM_eps = 0.01,
                            UM_path = paste0(master.dir, "/Old UM models/"),
                            CM_path = paste0(master.dir, "/Old CM models/"),
                            UM_formula = "SexAge", CM_formula = "SexAge",
                            UM_infix = "unstr", CM_infix = "cubic",
                            UM_suffix = "NR.Rdata", CM_suffix = "NR.Rdata")

cumSpeed.UM <- res.cumSpeed$outUns
cumSpeed.CM <- res.cumSpeed$outCub

# Unstructured model results
print(xtable(cbind(cumSpeed.UM[[1]][,c(2,3,4,7)],NA,cumSpeed.UM[[2]][,c(3,4,7)],NA,
                   cumSpeed.UM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(cumSpeed.UM[[4]][,c(2,3,4,7)],NA,cumSpeed.UM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# Cubic model results
print(xtable(cbind(cumSpeed.CM[[1]][,c(2,3,4,7)],NA,cumSpeed.CM[[2]][,c(3,4,7)],NA,
                   cumSpeed.CM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(cumSpeed.CM[[4]][,c(2,3,4,7)],NA,cumSpeed.CM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# UQR results
uqr.SexAge.cumSpeed <- list.load(paste0(uqr.dir,'/UQR.SexAge.cumSpeed.Rdata'))
uqr.SexAge.cumSpeed.df <- get_all_est_SE(summary_list = uqr.SexAge.cumSpeed,
                                         taus = c(0.1,0.25,0.5,0.75,0.9))
xtable(uqr.SexAge.cumSpeed.df, digits = 4)
# ------------------------------------------------------------------------------

# accel # ----------------------------------------------------------------------
res.accel <- summary_LQC(varName = 'accel',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.1, CM_eps = 0.1,
                            UM_path = paste0(master.dir, "/Old UM models/"),
                            CM_path = paste0(master.dir, "/Old CM models/"),
                            UM_formula = "SexAge", CM_formula = "SexAge", 
                            UM_infix = "unstr", CM_infix = "cubic",
                            UM_suffix = "NR.Rdata", CM_suffix = "NR.Rdata")

accel.UM <- res.accel$outUns
accel.CM <- res.accel$outCub

# Unstructured model results
print(xtable(cbind(accel.UM[[1]][,c(2,3,4,7)],NA,accel.UM[[2]][,c(3,4,7)],NA,
                   accel.UM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(accel.UM[[4]][,c(2,3,4,7)],NA,accel.UM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# Cubic model results
print(xtable(cbind(accel.CM[[1]][,c(2,3,4,7)],NA,accel.CM[[2]][,c(3,4,7)],NA,
                   accel.CM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(accel.CM[[4]][,c(2,3,4,7)],NA,accel.CM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# UQR results
uqr.SexAge.accel <- list.load(paste0(uqr.dir,'/UQR.SexAge.accel.Rdata'))
uqr.SexAge.accel.df <- get_all_est_SE(summary_list = uqr.SexAge.accel,
                                         taus = c(0.1,0.25,0.5,0.75,0.9))
xtable(uqr.SexAge.accel.df, digits = 4)
# ------------------------------------------------------------------------------

# speedDec # -------------------------------------------------------------------
res.speedDec <- summary_LQC(varName = 'speedDec',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.01, CM_eps = 0.1,
                            UM_path = paste0(master.dir, "/Old UM models/"),
                            CM_path = paste0(master.dir, "/Old CM models/"),
                            UM_formula = "SexAge", CM_formula = "SexAge",
                            UM_infix = "unstr", CM_infix = "cubic",
                            UM_suffix = "NR.Rdata", CM_suffix = "NR.Rdata")

speedDec.UM <- res.speedDec$outUns
speedDec.CM <- res.speedDec$outCub

# Unstructured model results
print(xtable(cbind(speedDec.UM[[1]][,c(2,3,4,7)],NA,speedDec.UM[[2]][,c(3,4,7)],NA,
                   speedDec.UM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(speedDec.UM[[4]][,c(2,3,4,7)],NA,speedDec.UM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# Cubic model results
print(xtable(cbind(speedDec.CM[[1]][,c(2,3,4,7)],NA,speedDec.CM[[2]][,c(3,4,7)],NA,
                   speedDec.CM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(speedDec.CM[[4]][,c(2,3,4,7)],NA,speedDec.CM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# UQR results
uqr.SexAge.speedDec <- list.load(paste0(uqr.dir,'/UQR.SexAge.speedDec.Rdata'))
uqr.SexAge.speedDec.df <- get_all_est_SE(summary_list = uqr.SexAge.speedDec,
                                         taus = c(0.1,0.25,0.5,0.75,0.9))
xtable(uqr.SexAge.speedDec.df, digits = 4)
# ------------------------------------------------------------------------------

# splitVec # -------------------------------------------------------------------
res.splitVec <- summary_LQC(varName = 'splitVec',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.05, CM_eps = 0.01,
                            UM_path = paste0(master.dir, "/Old UM models/"),
                            CM_path = paste0(master.dir, "/Old CM models/"),
                            UM_formula = "SexAge", CM_formula = "SexAge", 
                            UM_infix = "unstr", CM_infix = "CC",
                            UM_suffix = "NR.Rdata", CM_suffix = "NR.Rdata")

splitVec.UM <- res.splitVec$outUns
splitVec.CM <- res.splitVec$outCub

# Unstructured model results
print(xtable(cbind(splitVec.UM[[1]][,c(2,3,4,7)],NA,splitVec.UM[[2]][,c(3,4,7)],NA,
                   splitVec.UM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(splitVec.UM[[4]][,c(2,3,4,7)],NA,splitVec.UM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# Constrained cubic model results
print(xtable(cbind(splitVec.CM[[1]][,c(2,3,4,7)],NA,splitVec.CM[[2]][,c(3,4,7)],NA,
                   splitVec.CM[[3]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

print(xtable(cbind(splitVec.CM[[4]][,c(2,3,4,7)],NA,splitVec.CM[[5]][,c(3,4,7)]),
             digits = c(1,1,3,5,5,1,3,5,5)),
      include.rownames = FALSE)

# UQR results
uqr.SexAge.splitVec <- list.load(paste0(uqr.dir,'/UQR.SexAge.splitVec.Rdata'))
uqr.SexAge.splitVec.df <- get_all_est_SE(summary_list = uqr.SexAge.splitVec,
                                         taus = c(0.1,0.25,0.5,0.75,0.9))
xtable(uqr.SexAge.splitVec.df, digits = 4)
# ------------------------------------------------------------------------------


### SexAgeInter models (after review: with interaction terms age-sex and
# age-distance) ### -----

# intSpeed # -------------------------------------------------------------------
res.intSpeed <- summary_LQC(varName = 'intSpeed',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.01, CM_eps = 0.1,
                            UM_path = paste0(unstr.dir, "/"),
                            CM_path = paste0(cubic.dir, "/"),
                            UM_formula = "SexAgeInter",
                            CM_formula = "SexAgeInter",
                            UM_infix = "UM", CM_infix = "CM",
                            UM_suffix = "Rdata", CM_suffix = "Rdata")

intSpeed.UM <- res.intSpeed$outUns
intSpeed.CM <- res.intSpeed$outCub

## Unstructured model results

# with est, SE, gradient...
print(xtable(cbind(intSpeed.UM[[1]][,c(2,3,4,7)],NA,intSpeed.UM[[2]][,c(3,4,7)],NA,
                   intSpeed.UM[[3]][,c(3,4,7)],NA,
                   intSpeed.UM[[4]][,c(3,4,7)],NA,intSpeed.UM[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)

# # ... or with est, SE, p-values
# print(xtable(cbind(intSpeed.UM[[1]][,c(2,3,4,6)],NA,intSpeed.UM[[2]][,c(3,4,6)],NA,
#                    intSpeed.UM[[3]][,c(3,4,6)],NA,
#                    intSpeed.UM[[4]][,c(3,4,6)],NA,intSpeed.UM[[5]][,c(3,4,6)]),
#              digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
#       include.rownames = FALSE)


## Cubic model results

# with est, SE, gradient...
print(xtable(cbind(intSpeed.CM[[1]][,c(2,3,4,7)],NA,intSpeed.CM[[2]][,c(3,4,7)],NA,
                   intSpeed.CM[[3]][,c(3,4,7)],NA,
                   intSpeed.CM[[4]][,c(3,4,7)],NA,intSpeed.CM[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)

## UQR model results
uqr.SexAgeInter.intSpeed <- list.load(paste0(uqr.dir,'/UQR.SexAgeInter.intSpeed.Rdata'))
uqr.SexAgeInter.intSpeed.df <- get_all_est_SE(summary_list = uqr.SexAgeInter.intSpeed,
                                              taus = c(0.1,0.25,0.5,0.75,0.9))
xtable(uqr.SexAgeInter.intSpeed.df, digits = 4)




## Optimal epsilon models with est and SE only
print(xtable(cbind(intSpeed.UM[[1]][,c(2,3,4)],NA,intSpeed.UM[[2]][,c(3,4)],NA,
                   intSpeed.UM[[3]][,c(3,4)],NA,
                   intSpeed.UM[[4]][,c(3,4)],NA,intSpeed.UM[[5]][,c(3,4)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)

print(xtable(cbind(intSpeed.CM[[1]][,c(2,3,4)],NA,intSpeed.CM[[2]][,c(3,4)],NA,
                   intSpeed.CM[[3]][,c(3,4)],NA,
                   intSpeed.CM[[4]][,c(3,4)],NA,intSpeed.CM[[5]][,c(3,4)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)

## Optimal epsilon models with est and p only
print(xtable(cbind(intSpeed.UM[[1]][,c(2,3,6)],NA,intSpeed.UM[[2]][,c(3,6)],NA,
                   intSpeed.UM[[3]][,c(3,6)],NA,
                   intSpeed.UM[[4]][,c(3,6)],NA,intSpeed.UM[[5]][,c(3,6)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)

print(xtable(cbind(intSpeed.CM[[1]][,c(2,3,6)],NA,intSpeed.CM[[2]][,c(3,6)],NA,
                   intSpeed.CM[[3]][,c(3,6)],NA,
                   intSpeed.CM[[4]][,c(3,6)],NA,intSpeed.CM[[5]][,c(3,6)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)


# ## To also get p-values and confidence interval bounds ##
# intSpeed.UM.tau10 <- intSpeed.UM[[1]]
# intSpeed.UM.tau10$lbound <- intSpeed.UM.tau10$est - 1.96*intSpeed.UM.tau10$SE
# intSpeed.UM.tau10$rbound <- intSpeed.UM.tau10$est + 1.96*intSpeed.UM.tau10$SE
# 
# intSpeed.UM.tau25 <- intSpeed.UM[[2]]
# intSpeed.UM.tau25$lbound <- intSpeed.UM.tau25$est - 1.96*intSpeed.UM.tau25$SE
# intSpeed.UM.tau25$rbound <- intSpeed.UM.tau25$est + 1.96*intSpeed.UM.tau25$SE
# 
# intSpeed.UM.tau50 <- intSpeed.UM[[3]]
# intSpeed.UM.tau50$lbound <- intSpeed.UM.tau50$est - 1.96*intSpeed.UM.tau50$SE
# intSpeed.UM.tau50$rbound <- intSpeed.UM.tau50$est + 1.96*intSpeed.UM.tau50$SE
# 
# intSpeed.UM.tau75 <- intSpeed.UM[[4]]
# intSpeed.UM.tau75$lbound <- intSpeed.UM.tau75$est - 1.96*intSpeed.UM.tau75$SE
# intSpeed.UM.tau75$rbound <- intSpeed.UM.tau75$est + 1.96*intSpeed.UM.tau75$SE
# 
# intSpeed.UM.tau90 <- intSpeed.UM[[5]]
# intSpeed.UM.tau90$lbound <- intSpeed.UM.tau90$est - 1.96*intSpeed.UM.tau90$SE
# intSpeed.UM.tau90$rbound <- intSpeed.UM.tau90$est + 1.96*intSpeed.UM.tau90$SE

## Unstructured model with Student t-copula
res.intSpeed <- summary_LQC(varName = 'intSpeed',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.01, CM_eps = 0.1,
                            UM_path = paste0(unstr.dir, "/"),
                            CM_path = paste0(cubic.dir, "/"),
                            UM_formula = "SexAgeInter",
                            CM_formula = "SexAgeInter",
                            UM_infix = "UM.t", CM_infix = "CM",
                            UM_suffix = "Rdata", CM_suffix = "Rdata", noCM = TRUE)

intSpeed.UM.t <- res.intSpeed$outUns

print(xtable(cbind(intSpeed.UM.t[[1]][,c(2,3,4,7)],NA,intSpeed.UM.t[[2]][,c(3,4,7)],NA,
                   intSpeed.UM.t[[3]][,c(3,4,7)],NA,
                   intSpeed.UM.t[[4]][,c(3,4,7)],NA,intSpeed.UM.t[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)

# ------------------------------------------------------------------------------

# cumSpeed # -------------------------------------------------------------------
res.cumSpeed <- summary_LQC(varName = 'cumSpeed',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.01, CM_eps = 0.01,
                            UM_path = paste0(unstr.dir, "/"),
                            CM_path = paste0(cubic.dir, "/"),
                            UM_formula = "SexAgeInter",
                            CM_formula = "SexAgeInter",
                            UM_infix = "UM", CM_infix = "CM",
                            UM_suffix = "Rdata", CM_suffix = "Rdata")

cumSpeed.UM <- res.cumSpeed$outUns
cumSpeed.CM <- res.cumSpeed$outCub


## Unstructured model results

# with est, SE, gradient...
print(xtable(cbind(cumSpeed.UM[[1]][,c(2,3,4,7)],NA,cumSpeed.UM[[2]][,c(3,4,7)],NA,
                   cumSpeed.UM[[3]][,c(3,4,7)],NA,
                   cumSpeed.UM[[4]][,c(3,4,7)],NA,cumSpeed.UM[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)

# # ... or with est, SE, p-values
# print(xtable(cbind(cumSpeed.UM[[1]][,c(2,3,4,6)],NA,cumSpeed.UM[[2]][,c(3,4,6)],NA,
#                    cumSpeed.UM[[3]][,c(3,4,6)],NA,
#                    cumSpeed.UM[[4]][,c(3,4,6)],NA,cumSpeed.UM[[5]][,c(3,4,6)]),
#              digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
#       include.rownames = FALSE)


## Cubic model results

# with est, SE, gradient...
print(xtable(cbind(cumSpeed.CM[[1]][,c(2,3,4,7)],NA,cumSpeed.CM[[2]][,c(3,4,7)],NA,
                   cumSpeed.CM[[3]][,c(3,4,7)],NA,
                   cumSpeed.CM[[4]][,c(3,4,7)],NA,cumSpeed.CM[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)


## UQR model results
uqr.SexAgeInter.cumSpeed <- list.load(paste0(uqr.dir,'/UQR.SexAgeInter.cumSpeed.Rdata'))
uqr.SexAgeInter.cumSpeed.df <- get_all_est_SE(summary_list = uqr.SexAgeInter.cumSpeed,
                                              taus = c(0.1,0.25,0.5,0.75,0.9))
xtable(uqr.SexAgeInter.cumSpeed.df, digits = 4)



## Optimal epsilon models with est and SE only
print(xtable(cbind(cumSpeed.UM[[1]][,c(2,3,4)],NA,cumSpeed.UM[[2]][,c(3,4)],NA,
                   cumSpeed.UM[[3]][,c(3,4)],NA,
                   cumSpeed.UM[[4]][,c(3,4)],NA,cumSpeed.UM[[5]][,c(3,4)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)

print(xtable(cbind(cumSpeed.CM[[1]][,c(2,3,4)],NA,cumSpeed.CM[[2]][,c(3,4)],NA,
                   cumSpeed.CM[[3]][,c(3,4)],NA,
                   cumSpeed.CM[[4]][,c(3,4)],NA,cumSpeed.CM[[5]][,c(3,4)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)



# ------------------------------------------------------------------------------

# speedDec # -------------------------------------------------------------------
res.speedDec <- summary_LQC(varName = 'speedDec',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.1, CM_eps = 0.1,
                            UM_path = paste0(unstr.dir, "/"),
                            CM_path = paste0(cubic.dir, "/"),
                            UM_formula = "SexAgeInter",
                            CM_formula = "SexAgeInter",
                            UM_infix = "UM", CM_infix = "CM",
                            UM_suffix = "Rdata", CM_suffix = "Rdata")

speedDec.UM <- res.speedDec$outUns
speedDec.CM <- res.speedDec$outCub

## Unstructured model results

# with est, SE, gradient...
print(xtable(cbind(speedDec.UM[[1]][,c(2,3,4,7)],NA,speedDec.UM[[2]][,c(3,4,7)],NA,
                   speedDec.UM[[3]][,c(3,4,7)],NA,
                   speedDec.UM[[4]][,c(3,4,7)],NA,speedDec.UM[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)

# # ... or with est, SE, p-values
# print(xtable(cbind(speedDec.UM[[1]][,c(2,3,4,6)],NA,speedDec.UM[[2]][,c(3,4,6)],NA,
#                    speedDec.UM[[3]][,c(3,4,6)],NA,
#                    speedDec.UM[[4]][,c(3,4,6)],NA,speedDec.UM[[5]][,c(3,4,6)]),
#              digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
#       include.rownames = FALSE)


## Cubic model results

# with est, SE, gradient...
print(xtable(cbind(speedDec.CM[[1]][,c(2,3,4,7)],NA,speedDec.CM[[2]][,c(3,4,7)],NA,
                   speedDec.CM[[3]][,c(3,4,7)],NA,
                   speedDec.CM[[4]][,c(3,4,7)],NA,speedDec.CM[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)


## UQR model results
uqr.SexAgeInter.speedDec <- list.load(paste0(uqr.dir,'/UQR.SexAgeInter.speedDec.Rdata'))
uqr.SexAgeInter.speedDec.df <- get_all_est_SE(summary_list = uqr.SexAgeInter.speedDec,
                                              taus = c(0.1,0.25,0.5,0.75,0.9))
xtable(uqr.SexAgeInter.speedDec.df, digits = 4)


## Optimal epsilon models with est and SE only
print(xtable(cbind(speedDec.UM[[1]][,c(2,3,4)],NA,speedDec.UM[[2]][,c(3,4)],NA,
                   speedDec.UM[[3]][,c(3,4)],NA,
                   speedDec.UM[[4]][,c(3,4)],NA,speedDec.UM[[5]][,c(3,4)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)

print(xtable(cbind(speedDec.CM[[1]][,c(2,3,4)],NA,speedDec.CM[[2]][,c(3,4)],NA,
                   speedDec.CM[[3]][,c(3,4)],NA,
                   speedDec.CM[[4]][,c(3,4)],NA,speedDec.CM[[5]][,c(3,4)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)

# ------------------------------------------------------------------------------

# splitVec # -------------------------------------------------------------------
res.splitVec <- summary_LQC(varName = 'splitVec',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.01, CM_eps = 0.05,
                            UM_path = paste0(unstr.dir, "/"),
                            CM_path = paste0(cubic.dir, "/"),
                            UM_formula = "SexAgeInter",
                            CM_formula = "SexAgeInter", 
                            UM_infix = "UM", CM_infix = "CCM",
                            UM_suffix = "Rdata", CM_suffix = "Rdata")

splitVec.UM <- res.splitVec$outUns
splitVec.CM <- res.splitVec$outCub

## Unstructured model results

# with est, SE, gradient...
print(xtable(cbind(splitVec.UM[[1]][,c(2,3,4,7)],NA,splitVec.UM[[2]][,c(3,4,7)],NA,
                   splitVec.UM[[3]][,c(3,4,7)],NA,
                   splitVec.UM[[4]][,c(3,4,7)],NA,splitVec.UM[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)

# # ... or with est, SE, p-values
# print(xtable(cbind(splitVec.UM[[1]][,c(2,3,4,6)],NA,splitVec.UM[[2]][,c(3,4,6)],NA,
#                    splitVec.UM[[3]][,c(3,4,6)],NA,
#                    splitVec.UM[[4]][,c(3,4,6)],NA,splitVec.UM[[5]][,c(3,4,6)]),
#              digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
#       include.rownames = FALSE)


## Cubic model results

# with est, SE, gradient...
print(xtable(cbind(splitVec.CM[[1]][,c(2,3,4,7)],NA,splitVec.CM[[2]][,c(3,4,7)],NA,
                   splitVec.CM[[3]][,c(3,4,7)],NA,
                   splitVec.CM[[4]][,c(3,4,7)],NA,splitVec.CM[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)


# UQR results
uqr.SexAgeInter.splitVec <- list.load(paste0(uqr.dir,'/UQR.SexAgeInter.splitVec.Rdata'))
uqr.SexAgeInter.splitVec.df <- get_all_est_SE(summary_list = uqr.SexAgeInter.splitVec,
                                              taus = c(0.1,0.25,0.5,0.75,0.9))
xtable(uqr.SexAgeInter.splitVec.df, digits = 4)


## Optimal epsilon models with est and SE only
print(xtable(cbind(splitVec.UM[[1]][,c(2,3,4)],NA,splitVec.UM[[2]][,c(3,4)],NA,
                   splitVec.UM[[3]][,c(3,4)],NA,
                   splitVec.UM[[4]][,c(3,4)],NA,splitVec.UM[[5]][,c(3,4)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)

print(xtable(cbind(splitVec.CM[[1]][,c(2,3,4)],NA,splitVec.CM[[2]][,c(3,4)],NA,
                   splitVec.CM[[3]][,c(3,4)],NA,
                   splitVec.CM[[4]][,c(3,4)],NA,splitVec.CM[[5]][,c(3,4)]),
             digits = c(1,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4)),
      include.rownames = FALSE)



## Cubic model with Student t-copula
res.splitVec <- summary_LQC(varName = 'splitVec',
                            qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                            UM_eps = 0.1, CM_eps = 0.05,
                            UM_path = paste0(unstr.dir, "/"),
                            CM_path = paste0(cubic.dir, "/"),
                            UM_formula = "SexAgeInter",
                            CM_formula = "SexAgeInter", 
                            UM_infix = "UM", CM_infix = "CCM.t",
                            UM_suffix = "Rdata", CM_suffix = "Rdata", noUM = TRUE)

splitVec.CM.t <- res.splitVec$outCub

print(xtable(cbind(splitVec.CM.t[[1]][,c(2,3,4,7)],NA,splitVec.CM.t[[2]][,c(3,4,7)],NA,
                   splitVec.CM.t[[3]][,c(3,4,7)],NA,
                   splitVec.CM.t[[4]][,c(3,4,7)],NA,splitVec.CM.t[[5]][,c(3,4,7)]),
             digits = c(1,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4,1,4,4,4)),
      include.rownames = FALSE)
# ------------------------------------------------------------------------------

################################################################################



##### Part 5: obtaining plots ##################################################

# verify that the appropriate values for epsilon and corresponding folder paths
# are specified

# Note: function plot_LQC assumes a model specification belonging to one of the
# options above ('SexAge' before review, 'SexAgeInter' after first review).
# In case another one is to be considered, update plot_LQC by adding another
# model formula option; idem for plot_LQC_SplitVec.

# in this part, always use full data NYC.long (not NYC.noD1 or NYC.noDM; index
# selection is already taken care of in the functions)


### SexAge models (before review: without interaction) ### ---------------------
plot_LQC(Data = NYC.long, y_var_name = 'intSpeed',
         y_label = 'Interval speed (km/h)',
         y_range = c(3,17.5),
         y_ticks = c(5,7.5,10,12.5,15,17.5),
         female_box_location = c(39.5,17.2),
         male_box_location = c(40.5,17.2),
         qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
         plot_title = 'LQC interval speed',
         storage_path = plot.dir,
         storage_suffix = "SexAge",
         UM_path = paste0(master.dir, "/Old UM models/"),
         CM_path = paste0(master.dir, "/Old CM models/"),
         UQR_path = paste0(master.dir, "/Old UQR models/"),
         UM_formula = "SexAge", CM_formula = "SexAge", UQR_formula = "SexAge",
         UM_infix = "unstr", CM_infix = "cubic",
         UM_suffix = "NR.Rdata", CM_suffix = "NR.Rdata",
         UM_eps = 0.01, CM_eps = 0.05, Age = 35, SampleSize = 200,
         omit_indices = c(),
         save_separate_panels = FALSE)


plot_LQC(Data = NYC.long, y_var_name = 'cumSpeed',
         y_label = 'Cumulative speed (km/h)',
         y_range = c(3,17.5),
         y_ticks = c(5,7.5,10,12.5,15,17.5),
         female_box_location = c(39.5,17.2),
         male_box_location = c(40.5,17.2),
         qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
         plot_title = 'LQC cumulative speed',
         storage_path = plot.dir,
         storage_suffix = "SexAge",
         UM_path = paste0(master.dir, "/Old UM models/"),
         CM_path = paste0(master.dir, "/Old CM models/"),
         UQR_path = paste0(master.dir, "/Old UQR models/"),
         UM_formula = "SexAge", CM_formula = "SexAge", UQR_formula = "SexAge",
         UM_infix = "unstr", CM_infix = "cubic",
         UM_suffix = "NR.Rdata", CM_suffix = "NR.Rdata",
         UM_eps = 0.01, CM_eps = 0.01, Age = 35, SampleSize = 200,
         omit_indices = c(),
         save_separate_panels = FALSE)

plot_LQC(Data = NYC.long, y_var_name = 'speedDec',
         y_label = 'Speed variation',
         y_range = c(0.5,1.25),
         y_ticks = c(0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25),
         female_box_location = c(39.5,1.23),
         male_box_location = c(40.5,1.23),
         qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
         plot_title = 'LQC speed variation',
         storage_path = plot.dir,
         storage_suffix = "SexAge",
         UM_path = paste0(master.dir, "/Old UM models/"),
         CM_path = paste0(master.dir, "/Old CM models/"),
         UQR_path = paste0(master.dir, "/Old UQR models/"),
         UM_formula = "SexAge", CM_formula = "SexAge", UQR_formula = "SexAge",
         UM_infix = "unstr", CM_infix = "cubic",
         UM_suffix = "NR.Rdata", CM_suffix = "NR.Rdata",
         UM_eps = 0.01, CM_eps = 0.1, Age = 35, SampleSize = 200,
         omit_indices = c(1),
         save_separate_panels = FALSE)


plot_LQC_splitVec(Data = NYC.long, y_var_name = 'splitVec',
                  y_label = 'Splitting indicator',
                  y_range = c(0.9,1.6),
                  y_ticks = c(0.9,1,1.1,1.2,1.3,1.4,1.5,1.6),
                  female_box_location = c(37.5,1.59),
                  male_box_location = c(38.5,1.59),
                  qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                  plot_title = 'LQC splitting indicator',
                  storage_path = plot.dir,
                  storage_suffix = "SexAge",
                  UM_path = paste0(master.dir, "/Old UM models/"),
                  CM_path = paste0(master.dir, "/Old CM models/"),
                  UQR_path = paste0(master.dir, "/Old UQR models/"),
                  UM_formula = "SexAge", CCM_formula = "SexAge", UQR_formula = "SexAge",
                  UM_infix = "unstr", CCM_infix = "CC",
                  UM_suffix = "NR.Rdata", CCM_suffix = "NR.Rdata",
                  UM_eps = 0.05, CCM_eps = 0.01, Age = 35, SampleSize = 200,
                  omit_indices = c(11),
                  save_separate_panels = FALSE)
# ------------------------------------------------------------------------------


### SexAgeInter models (after review: with interaction for UM and UQR) ### -----
plot_LQC(Data = NYC.long, y_var_name = 'intSpeed',
         y_label = 'Interval speed (km/h)',
         y_range = c(3,17.5),
         y_ticks = c(5,7.5,10,12.5,15,17.5),
         female_box_location = c(39.5,17.2),
         male_box_location = c(40.5,17.2),
         qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
         plot_title = 'LQC interval speed',
         storage_path = plot.dir, storage_suffix = "SexAgeInter",
         UM_path = paste0(unstr.dir,"/"),
         CM_path = paste0(cubic.dir,"/"),
         UQR_path = paste0(uqr.dir,"/"),
         UM_formula = "SexAgeInter", CM_formula = "SexAgeInter",
         UQR_formula = "SexAgeInter",
         UM_infix = "UM", CM_infix = "CM",
         UM_suffix = "Rdata", CM_suffix = "Rdata",
         UM_eps = 0.01, CM_eps = 0.1, Age = 35, SampleSize = 200,
         omit_indices = c(),
         save_separate_panels = FALSE)


plot_LQC(Data = NYC.long, y_var_name = 'cumSpeed',
         y_label = 'Cumulative speed (km/h)',
         y_range = c(3,17.5),
         y_ticks = c(5,7.5,10,12.5,15,17.5),
         female_box_location = c(39.5,17.2),
         male_box_location = c(40.5,17.2),
         qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
         plot_title = 'LQC cumulative speed',
         storage_path = plot.dir, storage_suffix = "SexAgeInter",
         UM_path = paste0(unstr.dir,"/"),
         CM_path = paste0(cubic.dir,"/"),
         UQR_path = paste0(uqr.dir,"/"),
         UM_formula = "SexAgeInter", CM_formula = "SexAgeInter",
         UQR_formula = "SexAgeInter",
         UM_infix = "UM", CM_infix = "CM",
         UM_suffix = "Rdata", CM_suffix = "Rdata",
         UM_eps = 0.01, CM_eps = 0.01, Age = 35, SampleSize = 200,
         omit_indices = c(),
         save_separate_panels = FALSE)

plot_LQC(Data = NYC.long, y_var_name = 'speedDec',
         y_label = 'Speed variation',
         y_range = c(0.5,1.25),
         y_ticks = c(0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25),
         female_box_location = c(39.5,1.23),
         male_box_location = c(40.5,1.23),
         qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
         plot_title = 'LQC speed variation',
         storage_path = plot.dir, storage_suffix = "SexAgeInter",
         UM_path = paste0(unstr.dir,"/"),
         CM_path = paste0(cubic.dir,"/"),
         UQR_path = paste0(uqr.dir,"/"),
         UM_formula = "SexAgeInter", CM_formula = "SexAgeInter",
         UQR_formula = "SexAgeInter",
         UM_infix = "UM", CM_infix = "CM",
         UM_suffix = "Rdata", CM_suffix = "Rdata",
         UM_eps = 0.05, CM_eps = 0.1, Age = 35, SampleSize = 200,
         omit_indices = c(1),
         save_separate_panels = FALSE)


plot_LQC_splitVec(Data = NYC.long, y_var_name = 'splitVec',
                  y_label = 'Splitting indicator',
                  y_range = c(0.9,1.6),
                  y_ticks = c(0.9,1,1.1,1.2,1.3,1.4,1.5,1.6),
                  female_box_location = c(37.5,1.59),
                  male_box_location = c(38.5,1.59),
                  qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                  plot_title = 'LQC splitting indicator',
                  storage_path = plot.dir, storage_suffix = "SexAgeInter",
                  UM_path = paste0(unstr.dir,"/"),
                  CCM_path = paste0(cubic.dir,"/"),
                  UQR_path = paste0(uqr.dir,"/"),
                  UM_formula = "SexAgeInter", CCM_formula = "SexAgeInter",
                  UQR_formula = "SexAgeInter",
                  UM_infix = "UM", CCM_infix = "CCM",
                  UM_suffix = "Rdata", CCM_suffix = "Rdata",
                  UM_eps = 0.1, CCM_eps = 0.05, Age = 35, SampleSize = 200,
                  omit_indices = c(11),
                  save_separate_panels = FALSE)
# ------------------------------------------------------------------------------

################################################################################



##### Part 6: comparison to classical methods ##################################
library(lme4)
library(tidyr)
library(emmeans)
library(afex)
library(performance)
# library(see)
library(qqplotr)
library(lavaan) # Structural equation modelling: https://lavaan.ugent.be/tutorial/sem.html

##### Part 6.1: Comparison to Nikolaidis & Knechtle (2018),
# partly ignoring longitudinality ##############################################

NYC.short <- NYC.long[,c("bib", "Sex", "AgeGroup", "AgeScaled", "time",
                         "distance", "dis", "Slope", "intSpeed", "cumSpeed",
                         "accel", "speedDec", "splitVec")]
NYC.short$Sex <- as.factor(NYC.short$Sex)
NYC.short$bib <- as.factor(NYC.short$bib)

# extra variables for linear regression...
# NYC.short$distIndex <- seq_along(unique(NYC.short$dis)) # this makes sense only
# when the measurements are equidistant, possibly the initial or final segment
# aside, as they can be taken care of in the initSegm and finalSegm variables
NYC.short$distIndex <- NYC.short$dis/5
NYC.short$initSegm <- c(1, rep(0, length(Dist.vec)-1))
NYC.short$finalSegm <- c(rep(0, length(Dist.vec)-1),1)

# ... and for ANOVA: categorical for distance
NYC.short$distCat <- as.factor(c("D01", "D02", "D03", "D04", "D05", "D06",
                                 "D07", "D08","D09", "D10","D11"))

#### (1) Repeated measures ANOVA #### ------------------------------------------

### split (within-subjects) and sex (between-subjects) on speed:
aov.sex <- aov_ez(id = "bib", dv = "intSpeed", data = NYC.short,
                  between = c("Sex"), within = c("distCat"),
                  anova_table = list(es = "pes"))

# post-hoc pairwise comparisons (Bonferroni corrected)
NK.sex.posthoc <- emmeans(aov.sex, specs = list(pairwise ~ distCat),
                          adjust = "bonferroni")
pw.diff.sex <- as.data.frame(NK.sex.posthoc$`pairwise differences of distCat`)
matrix.signif.sex <- get.signif.mat(pw.diff.sex)
pw.consec.sex <- pw.diff.sex[c(1,11,20,28,35,41,46,50,53,55),] # only consecutive race segments


# verification of assumptions (none satisfied)
check_homogeneity(aov.sex)
check_sphericity(aov.sex)
is_norm.sex <- check_normality(aov.sex)
plot(is_norm.sex)



### per sex: split (within-subjects) and AgeGroup (between-subjects) on speed:
NYC.short.M <- NYC.short[NYC.short$Sex == "M", ]
NYC.short.F <- NYC.short[NYC.short$Sex == "F", ]

aov.age.M <- aov_ez(id = "bib", dv = "intSpeed", data = NYC.short.M,
                    between = c("AgeGroup"), within = c("distCat"),
                    anova_table = list(es = "pes"))
aov.age.F <- aov_ez(id = "bib", dv = "intSpeed", data = NYC.short.F,
                    between = c("AgeGroup"), within = c("distCat"),
                    anova_table = list(es = "pes"))

# post-hoc pairwise comparisons (Bonferroni corrected)
NK.age.M.posthoc <- emmeans(aov.age.M, specs = list(pairwise ~ distCat),
                            adjust = "bonferroni")
NK.age.F.posthoc <- emmeans(aov.age.F, specs = list(pairwise ~ distCat),
                            adjust = "bonferroni")
pw.diff.age.M <- as.data.frame(NK.age.M.posthoc$`pairwise differences of distCat`)
pw.diff.age.F <- as.data.frame(NK.age.F.posthoc$`pairwise differences of distCat`)
matrix.signif.age.M <- get.signif.mat(pw.diff.age.M)
matrix.signif.age.F <- get.signif.mat(pw.diff.age.F)
pw.consec.age.M <- pw.diff.age.M[c(1,11,20,28,35,41,46,50,53,55),] # only consecutive race segments
pw.consec.age.F <- pw.diff.age.F[c(1,11,20,28,35,41,46,50,53,55),] # only consecutive race segments


# verification of assumptions (none satisfied)
check_homogeneity(aov.age.M)
check_sphericity(aov.age.M)
is_norm.age.M <- check_normality(aov.age.M)
plot(is_norm.age.M)

check_homogeneity(aov.age.F)
check_sphericity(aov.age.F)
is_norm.age.F <- check_normality(aov.age.F)
plot(is_norm.age.F)



# Combine output for all post-hoc tests (only consecutive race segments)
pw.consec.all <- cbind(pw.consec.sex[,c("estimate", "SE", "p.value")], NA,
                       pw.consec.age.M[,c("estimate", "SE", "p.value")], NA,
                       pw.consec.age.F[,c("estimate", "SE", "p.value")])

xtable(pw.consec.all, digits = 4)
# ------------------------------------------------------------------------------


#### (2) Linear regression #### ------------------------------------------------
# NOTE: what NK (2018) describe as "multivariate" comes down to several multiple
# regression models, i.e., one response variable with different groups underlying the
# analysis, with multiple explanatory variables each time

ageGroups <- levels(NYC.short$AgeGroup)

# omit ageGroup 1 and 8 due to small sample sizes

modelsM <- lapply(2:7, function(k){
  lm(intSpeed ~ distIndex + as.factor(initSegm) +
       as.factor(finalSegm) + Slope,
     data = NYC.short[NYC.short$Sex=="M" &
                        (NYC.short$AgeGroup == ageGroups[k]),])
})
modelsF <- lapply(2:7, function(k){
  lm(intSpeed ~ distIndex + as.factor(initSegm) +
       as.factor(finalSegm) + Slope,
     data = NYC.short[NYC.short$Sex=="F" &
                        (NYC.short$AgeGroup == ageGroups[k]),])
})


## get output summary; each time row of estimates, row of SE and row of p-values
# male output
df.M <- do.call('rbind', lapply(modelsM, get.summ))
colnames(df.M) <- c("Intercept", "distIndex", "initSegm", "finSegm", "Slope")
rownames(df.M) <- do.call("c", lapply(
  list("M 20-29", "M 30-39", "M 40-49", "M 50-59", "M 60-69", "M 70-79"),
  FUN = function(strInp){
    paste(strInp, c("est", "SE", "p"))
  }))
print(xtable(df.M, digits = c(0,4,4,4,4,4)), include.rownames = TRUE)

# female output
df.F <- do.call('rbind', lapply(modelsF, get.summ))
colnames(df.F) <- c("Intercept", "distIndex", "initSegm", "finSegm", "Slope")
rownames(df.F) <- do.call("c", lapply(
  list("F 20-29", "F 30-39", "F 40-49", "F 50-59", "F 60-69", "F 70-79"),
  FUN = function(strInp){
    paste(strInp, c("est", "SE", "p"))
  }))
print(xtable(df.F, digits = c(0,4,4,4,4,4)), include.rownames = TRUE)


## verifying linear regression assumptions
check_model(modelsM[[1]]) # just one example
# looks generally okay, patterns just seem to be a bit weird due to the 'continuous'
# distIndex variable taking on only 11 different values.



### Reducing to non-longitudinal data first, to avoid wrong precision estimates
# (should still be unbiased)

## Create new data frame, containing only one (random) observation per subject
NYC.red <- NYC.short
NYC.red$row <- 1:nrow(NYC.short)

set.seed(123)
red.indices <- sapply(unique(NYC.red$bib),
                      FUN = function(k){
                        sample(NYC.red$row[NYC.red$bib == k], size = 1)
                      })
NYC.red <- NYC.red[NYC.red$row %in% red.indices, ]


modelsM.red <- lapply(2:7, function(k){
  lm(intSpeed ~ distIndex + as.factor(initSegm) +
       as.factor(finalSegm) + Slope,
     data = NYC.red[NYC.red$Sex=="M" &
                      (NYC.red$AgeGroup == ageGroups[k]),])
})
modelsF.red <- lapply(2:7, function(k){
  lm(intSpeed ~ distIndex + as.factor(initSegm) +
       as.factor(finalSegm) + Slope,
     data = NYC.red[NYC.red$Sex=="F" &
                      (NYC.red$AgeGroup == ageGroups[k]),])
})


## get output summary; each time row of estimates, row of SE and row of p-values
# male output
df.M.red <- do.call('rbind', lapply(modelsM.red, get.summ))
colnames(df.M.red) <- c("Intercept", "distIndex", "initSegm", "finSegm", "Slope")
rownames(df.M.red) <- do.call("c", lapply(
  list("M 20-29", "M 30-39", "M 40-49", "M 50-59", "M 60-69", "M 70-79"),
  FUN = function(strInp){
    paste(strInp, c("est", "SE", "p"))
  }))
print(xtable(df.M.red, digits = c(0,4,4,4,4,4)), include.rownames = TRUE)

# female output
df.F.red <- do.call('rbind', lapply(modelsF.red, get.summ))
colnames(df.F.red) <- c("Intercept", "distIndex", "initSegm", "finSegm", "Slope")
rownames(df.F.red) <- do.call("c", lapply(
  list("F 20-29", "F 30-39", "F 40-49", "F 50-59", "F 60-69", "F 70-79"),
  FUN = function(strInp){
    paste(strInp, c("est", "SE", "p"))
  }))
print(xtable(df.F.red, digits = c(0,4,4,4,4,4)), include.rownames = TRUE)


# Comparison of estimates (ratios should be close to 1)
(df.M.red/df.M)[c(1,4,7,10,13,16), ]
(df.F.red/df.F)[c(1,4,7,10,13,16), ]

# or differences should be close to 0
(df.M.red-df.M)[c(1,4,7,10,13,16), ]
(df.F.red-df.F)[c(1,4,7,10,13,16), ]

print(xtable((df.M.red-df.M)[c(1,4,7,10,13,16), ],
             digits = c(0,4,4,4,4,4)), include.rownames = TRUE)

print(xtable((df.F.red-df.F)[c(1,4,7,10,13,16), ],
             digits = c(0,4,4,4,4,4)), include.rownames = TRUE)


# reduction of sample size by 1/11 --> expected increase of SE: by a factor
# between 2 sqrt(2) and 4, roughly between (3,3.5), as halving sample size comes
# with factor sqrt(2) SE inflation
# more precisely, factor (sqrt(2))^log2(11) = 3.316625 is expected

(df.M.red/df.M)[c(2,5,8,11,14,17), ]
(df.F.red/df.F)[c(2,5,8,11,14,17), ]

print(xtable((df.M.red/df.M)[c(2,5,8,11,14,17), ],
             digits = c(0,4,4,4,4,4)), include.rownames = TRUE)

print(xtable((df.F.red/df.F)[c(2,5,8,11,14,17), ],
             digits = c(0,4,4,4,4,4)), include.rownames = TRUE)
# ------------------------------------------------------------------------------

################################################################################


##### Part 6.2: With longitudinality: multivariate linear regression ###########

# Note: AgeScaled now rather than AgeGroup (factorial) in ANOVA

# get data in wide format
wide.NYC.intSpeed <- pivot_wider(
  data = NYC.short[, c("bib", "Sex", "AgeScaled", "distCat", "intSpeed")],
  id_cols = c("bib", "Sex", "AgeScaled"), names_from = distCat,
  values_from = c(intSpeed))

## Models with unstructured mean and unstructured covariance matrix
# (unstructured mean in the sense that both intercept and slopes per distance
# can vary in any way, not necessarily parallel / same intercepts etc.)

multiv.lm <- lm(formula =
                  cbind(D01,D02,D03,D04,D05,D06,D07,D08,D09,D10,D11) ~ Sex + AgeScaled,
                data = wide.NYC.intSpeed)
# Slope not included, as it is identical for each subject, hence automatically
# accounted for in the coefficients per distance.
multiv.lm
# seems that AgeScaled coefficient is very similar in all cases, maybe we can
# use a common slope there?


# check whether the impact of age is different per sex
multiv.inter.lm <- lm(formula =
                        cbind(D01,D02,D03,D04,D05,D06,D07,D08,D09,D10,D11) ~ Sex*AgeScaled,
                      data = wide.NYC.intSpeed)
anova(multiv.lm, multiv.inter.lm) # highly significant, i.e., there is indeed a
# sex*Age interaction for the joint multivariate model


## Testing for equality of the coefficients: apparently not so easy in R, run lm
# in a different way (not with lm() function; using lavaan package + constraints)

wide.NYC.extended <- wide.NYC.intSpeed
# different encoding (to make sense of interactions): 0 for females, 1 for males
wide.NYC.extended$Sex <- as.numeric(1*(wide.NYC.extended$Sex == "M"))




## Basis model
model.unconstr <- '
# regressions
D01 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D02 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D03 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D04 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D05 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D06 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D07 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D08 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D09 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D10 ~ 1 + Sex + AgeScaled + Sex:AgeScaled
D11 ~ 1 + Sex + AgeScaled + Sex:AgeScaled

# allow residual correlations (multivariate modelling) between all outcomes
# note: unnecessary if no constraints are to be imposed; assumed anyway
# could be used to compare to the model where e.g. D10 ~~ 0*D11, to see whether
# there is a correlation. Or use D01 ~~ a* D02, D10 ~~ b*D11 and a == b to test
# whether performance is significantly different for the model with different
# covariances, to test whether these covariances are significantly different, etc.
D01 ~~ D02 + D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D02 ~~ D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D03 ~~ D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D04 ~~ D05 + D06 + D07 + D08 + D09 + D10 + D11
D05 ~~ D06 + D07 + D08 + D09 + D10 + D11
D06 ~~ D07 + D08 + D09 + D10 + D11
D07 ~~ D08 + D09 + D10 + D11
D08 ~~ D09 + D10 + D11
D09 ~~ D10 + D11
D10 ~~ D11
'

# model fitting and output
fit.unconstr <- sem(model.unconstr, data = wide.NYC.extended)
summary(fit.unconstr)
# more output can be obtained through  summary(fit.unconstr)$pe

coef.unconstr <- matrix(coef(fit.unconstr)[1:44], nrow = 4, ncol = 11, byrow = FALSE)
se.unconstr <- matrix(summary(fit.unconstr)$pe$se[1:44], nrow = 4, ncol = 11,
                      byrow = FALSE)

df.coef.unconstr <- as.data.frame(coef.unconstr)
rownames(df.coef.unconstr) <- c("Intercept", "SexM", "AgeSc", "SexM:AgeSc")
colnames(df.coef.unconstr) <- c("D01", "D02", "D03", "D04", "D05", "D06",
                                "D07", "D08","D09", "D10","D11")

df.all.unconstr <- data.frame(Intercept = coef.unconstr[1,],
                              Intercept.SE = se.unconstr[1,],
                              SexM = coef.unconstr[2,],
                              SexM.SE = se.unconstr[2,],
                              AgeSc = coef.unconstr[3,],
                              AgeSc.SE = se.unconstr[3,],
                              SexM.AgeSc = coef.unconstr[4,],
                              SexM.AgeSc.SE = se.unconstr[4,])
rownames(df.all.unconstr) <- c("K5", "K10", "K15", "K20", "Half", "K25",
                               "K30", "M20","K35", "K40","MAR")

xtable(cbind(df.all.unconstr[,c(1,2)], NA,
             df.all.unconstr[,c(3,4)], NA,
             df.all.unconstr[,c(5,6)], NA,
             df.all.unconstr[,c(7,8)]), digits = 4)
# note: exactly the same output as multiv.inter.lm above

# to also get the covariance matrix
mySumm <- summary(fit.unconstr)$pe

covMat.unconstr <- matrix(NA,11,11)
covMat.unconstr[1,2:11] <- mySumm$est[45:54]
covMat.unconstr[2,3:11] <- mySumm$est[55:63]
covMat.unconstr[3,4:11] <- mySumm$est[64:71]
covMat.unconstr[4,5:11] <- mySumm$est[72:78]
covMat.unconstr[5,6:11] <- mySumm$est[79:84]
covMat.unconstr[6,7:11] <- mySumm$est[85:89]
covMat.unconstr[7,8:11] <- mySumm$est[90:93]
covMat.unconstr[8,9:11] <- mySumm$est[94:96]
covMat.unconstr[9,10:11] <- mySumm$est[97:98]
covMat.unconstr[10,11] <- mySumm$est[99]
diag(covMat.unconstr) <- mySumm$est[100:110]

xtable(covMat.unconstr, digits = 4)


covMat.SE.unconstr <- matrix(NA,11,11)
covMat.SE.unconstr[1,2:11] <- mySumm$se[45:54]
covMat.SE.unconstr[2,3:11] <- mySumm$se[55:63]
covMat.SE.unconstr[3,4:11] <- mySumm$se[64:71]
covMat.SE.unconstr[4,5:11] <- mySumm$se[72:78]
covMat.SE.unconstr[5,6:11] <- mySumm$se[79:84]
covMat.SE.unconstr[6,7:11] <- mySumm$se[85:89]
covMat.SE.unconstr[7,8:11] <- mySumm$se[90:93]
covMat.SE.unconstr[8,9:11] <- mySumm$se[94:96]
covMat.SE.unconstr[9,10:11] <- mySumm$se[97:98]
covMat.SE.unconstr[10,11] <- mySumm$se[99]
diag(covMat.SE.unconstr) <- mySumm$se[100:110]

# to avoid lengthy printing in paper, note that all standard errors are extremely
# similar; reporting one range should suffice
range(mySumm$se[45:110])



# assumption verification (one example dependent variable)
check_model(lm(D01 ~ Sex*AgeScaled, data = wide.NYC.intSpeed))



## Basis model without interaction terms
model.unconstr.nointer <- '
# regressions
D01 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D02 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D03 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D04 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D05 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D06 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D07 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D08 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D09 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D10 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled
D11 ~ 1 + Sex + AgeScaled + 0*Sex:AgeScaled

# allow residual correlations (multivariate modelling) between all outcomes
D01 ~~ D02 + D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D02 ~~ D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D03 ~~ D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D04 ~~ D05 + D06 + D07 + D08 + D09 + D10 + D11
D05 ~~ D06 + D07 + D08 + D09 + D10 + D11
D06 ~~ D07 + D08 + D09 + D10 + D11
D07 ~~ D08 + D09 + D10 + D11
D08 ~~ D09 + D10 + D11
D09 ~~ D10 + D11
D10 ~~ D11
'

# model fitting and output
fit.unconstr.nointer <- sem(model.unconstr.nointer, data = wide.NYC.extended)
summary(fit.unconstr.nointer)


# verification: including interaction terms indeed improves the model
anova(fit.unconstr, fit.unconstr.nointer)




## Basis model without interaction terms (implemented differently)
model.unconstr.nointer2 <- '
# regressions
D01 ~ 1 + Sex + AgeScaled
D02 ~ 1 + Sex + AgeScaled
D03 ~ 1 + Sex + AgeScaled
D04 ~ 1 + Sex + AgeScaled
D05 ~ 1 + Sex + AgeScaled
D06 ~ 1 + Sex + AgeScaled
D07 ~ 1 + Sex + AgeScaled
D08 ~ 1 + Sex + AgeScaled
D09 ~ 1 + Sex + AgeScaled
D10 ~ 1 + Sex + AgeScaled
D11 ~ 1 + Sex + AgeScaled

# allow residual correlations (multivariate modelling) between all outcomes
D01 ~~ D02 + D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D02 ~~ D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D03 ~~ D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D04 ~~ D05 + D06 + D07 + D08 + D09 + D10 + D11
D05 ~~ D06 + D07 + D08 + D09 + D10 + D11
D06 ~~ D07 + D08 + D09 + D10 + D11
D07 ~~ D08 + D09 + D10 + D11
D08 ~~ D09 + D10 + D11
D09 ~~ D10 + D11
D10 ~~ D11
'

# model fitting and output
fit.unconstr.nointer2 <- sem(model.unconstr.nointer2, data = wide.NYC.extended)
summary(fit.unconstr.nointer2)


# verification: including interaction terms indeed improves the model
anova(fit.unconstr, fit.unconstr.nointer2)
# note: slightly different testing output from the nointer implementation above,
# but the same conclusion + exactly the same estimates, SE etc. in summary().







## Basis model + imposed equality of Sex effect

model.constr.Sex <- '
# regressions
D01 ~ 1 + a*Sex + AgeScaled + l*Sex:AgeScaled
D02 ~ 1 + b*Sex + AgeScaled + m*Sex:AgeScaled
D03 ~ 1 + c*Sex + AgeScaled + n*Sex:AgeScaled
D04 ~ 1 + d*Sex + AgeScaled + o*Sex:AgeScaled
D05 ~ 1 + e*Sex + AgeScaled + p*Sex:AgeScaled
D06 ~ 1 + f*Sex + AgeScaled + q*Sex:AgeScaled
D07 ~ 1 + g*Sex + AgeScaled + r*Sex:AgeScaled
D08 ~ 1 + h*Sex + AgeScaled + s*Sex:AgeScaled
D09 ~ 1 + i*Sex + AgeScaled + t*Sex:AgeScaled
D10 ~ 1 + j*Sex + AgeScaled + u*Sex:AgeScaled
D11 ~ 1 + k*Sex + AgeScaled + v*Sex:AgeScaled

# allow residual correlations (multivariate modelling) between all outcomes
D01 ~~ D02 + D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D02 ~~ D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D03 ~~ D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D04 ~~ D05 + D06 + D07 + D08 + D09 + D10 + D11
D05 ~~ D06 + D07 + D08 + D09 + D10 + D11
D06 ~~ D07 + D08 + D09 + D10 + D11
D07 ~~ D08 + D09 + D10 + D11
D08 ~~ D09 + D10 + D11
D09 ~~ D10 + D11
D10 ~~ D11

# constrain paths
a == b
b == c
c == d
d == e
e == f
f == g
g == h
h == i
i == j
j == k

# also constrain interaction terms, since otherwise we still can get different
# contributions per sex through the Sex:AgeScaled terms
l == m
m == n
n == o
o == p
p == q
q == r
r == s
s == t
t == u
u == v
'


fit.constr.Sex <- sem(model.constr.Sex, data = wide.NYC.extended)
# summary(fit.constr.Sex)

# check whether significant difference between both models
anova(fit.unconstr, fit.constr.Sex)
# minimal AIC preferred => with non-common Sex (+ interactions) significantly better







## Basis model + imposed equality of AgeScaled effect

model.constr.Age <- '
# regressions
D01 ~ 1 + Sex + a*AgeScaled + l*Sex:AgeScaled
D02 ~ 1 + Sex + b*AgeScaled + m*Sex:AgeScaled
D03 ~ 1 + Sex + c*AgeScaled + n*Sex:AgeScaled
D04 ~ 1 + Sex + d*AgeScaled + o*Sex:AgeScaled
D05 ~ 1 + Sex + e*AgeScaled + p*Sex:AgeScaled
D06 ~ 1 + Sex + f*AgeScaled + q*Sex:AgeScaled
D07 ~ 1 + Sex + g*AgeScaled + r*Sex:AgeScaled
D08 ~ 1 + Sex + h*AgeScaled + s*Sex:AgeScaled
D09 ~ 1 + Sex + i*AgeScaled + t*Sex:AgeScaled
D10 ~ 1 + Sex + j*AgeScaled + u*Sex:AgeScaled
D11 ~ 1 + Sex + k*AgeScaled + v*Sex:AgeScaled

# allow residual correlations (multivariate modelling) between all outcomes
D01 ~~ D02 + D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D02 ~~ D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D03 ~~ D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D04 ~~ D05 + D06 + D07 + D08 + D09 + D10 + D11
D05 ~~ D06 + D07 + D08 + D09 + D10 + D11
D06 ~~ D07 + D08 + D09 + D10 + D11
D07 ~~ D08 + D09 + D10 + D11
D08 ~~ D09 + D10 + D11
D09 ~~ D10 + D11
D10 ~~ D11

# constrain paths
a == b
b == c
c == d
d == e
e == f
f == g
g == h
h == i
i == j
j == k

# also constrain interaction terms, since otherwise we still can get different
# contributions per age through the Sex:AgeScaled terms
l == m
m == n
n == o
o == p
p == q
q == r
r == s
s == t
t == u
u == v
'


fit.constr.Age <- sem(model.constr.Age, data = wide.NYC.extended)
# summary(fit.constr.Age)


# check whether significant difference between both models
anova(fit.unconstr, fit.constr.Age)
# minimal AIC preferred => with variable Age is significantly better






## Basis model + imposed equality of interaction effect Sex*AgeScaled

model.constr.AgeSex <- '
# regressions
D01 ~ 1 + Sex + AgeScaled + l*Sex:AgeScaled
D02 ~ 1 + Sex + AgeScaled + m*Sex:AgeScaled
D03 ~ 1 + Sex + AgeScaled + n*Sex:AgeScaled
D04 ~ 1 + Sex + AgeScaled + o*Sex:AgeScaled
D05 ~ 1 + Sex + AgeScaled + p*Sex:AgeScaled
D06 ~ 1 + Sex + AgeScaled + q*Sex:AgeScaled
D07 ~ 1 + Sex + AgeScaled + r*Sex:AgeScaled
D08 ~ 1 + Sex + AgeScaled + s*Sex:AgeScaled
D09 ~ 1 + Sex + AgeScaled + t*Sex:AgeScaled
D10 ~ 1 + Sex + AgeScaled + u*Sex:AgeScaled
D11 ~ 1 + Sex + AgeScaled + v*Sex:AgeScaled

# allow residual correlations (multivariate modelling) between all outcomes
D01 ~~ D02 + D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D02 ~~ D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D03 ~~ D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11
D04 ~~ D05 + D06 + D07 + D08 + D09 + D10 + D11
D05 ~~ D06 + D07 + D08 + D09 + D10 + D11
D06 ~~ D07 + D08 + D09 + D10 + D11
D07 ~~ D08 + D09 + D10 + D11
D08 ~~ D09 + D10 + D11
D09 ~~ D10 + D11
D10 ~~ D11

# constrain paths
l == m
m == n
n == o
o == p
p == q
q == r
r == s
s == t
t == u
u == v
'


fit.constr.AgeSex <- sem(model.constr.AgeSex, data = wide.NYC.extended)
mySumm <- summary(fit.constr.AgeSex)



# check whether significant difference between both models
anova(fit.unconstr, fit.constr.AgeSex)
# minimal AIC preferred => with variable Age is significantly better
################################################################################



##### Part 7: Comparison with vs. without interaction term #####################
### PWE - UM models

# intSpeed ---------------------------------------------------------------------
comp.UM.intSpeed <- compare_LQC(Data = NYC.long,
                             y_var_name = "intSpeed",
                             y_label = 'Interval speed (km/h)',
                             y_range = c(3,17.5),
                             y_ticks = c(5,7.5,10,12.5,15,17.5),
                             female_box_location = c(39.5,17.2),
                             male_box_location = c(40.5,17.2),
                             qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                             plot_title = 'LQC interval speed UM (comp.)',
                             storage_path = plot.dir,
                             modelType = "UM",
                             res_formula_vec = c("SexAge", "SexAgeInter"),
                             res_path_vec = c(paste0(master.dir, "/Old UM models/"),
                                              paste0(unstr.dir,"/")),
                             res_eps_vec = c(0.01,0.05),
                             res_infix_vec = c("unstr", "UM"),
                             res_suffix_vec = c("NR.Rdata","Rdata"),
                             Age = 35, SampleSize = 200,
                             omit_indices = c())

comp.UM.intSpeed$AIC
comp.UM.intSpeed$BIC
# ------------------------------------------------------------------------------


# cumSpeed ---------------------------------------------------------------------
comp.UM.cumSpeed <- compare_LQC(Data = NYC.long,
                             y_var_name = "cumSpeed",
                             y_label = 'Cumulative speed (km/h)',
                             y_range = c(3,17.5),
                             y_ticks = c(5,7.5,10,12.5,15,17.5),
                             female_box_location = c(39.5,17.2),
                             male_box_location = c(40.5,17.2),
                             qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                             plot_title = 'LQC cumulative speed UM (comp.)',
                             storage_path = plot.dir,
                             modelType = "UM",
                             res_formula_vec = c("SexAge", "SexAgeInter"),
                             res_path_vec = c(paste0(master.dir, "/Old UM models/"),
                                              paste0(unstr.dir,"/")),
                             res_eps_vec = c(0.01,0.01),
                             res_infix_vec = c("unstr", "UM"),
                             res_suffix_vec = c("NR.Rdata","Rdata"),
                             Age = 35, SampleSize = 200,
                             omit_indices = c())

comp.UM.cumSpeed$AIC
comp.UM.cumSpeed$BIC
# ------------------------------------------------------------------------------


# SplitVec ---------------------------------------------------------------------
comp.UM.splitVec <- compare_LQC(Data = NYC.long,
                             y_var_name = "splitVec",
                             y_label = 'Splitting indicator',
                             y_range = c(0.9,1.6),
                             y_ticks = c(0.9,1,1.1,1.2,1.3,1.4,1.5,1.6),
                             female_box_location = c(37.5,1.59),
                             male_box_location = c(38.5,1.59),
                             qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                             plot_title = 'LQC splitting indicator UM (comp.)',
                             storage_path = plot.dir,
                             modelType = "UM",
                             res_formula_vec = c("SexAge", "SexAgeInter"),
                             res_path_vec = c(paste0(master.dir, "/Old UM models/"),
                                              paste0(unstr.dir,"/")),
                             res_eps_vec = c(0.05,0.1),
                             res_infix_vec = c("unstr", "UM"),
                             res_suffix_vec = c("NR.Rdata","Rdata"),
                             Age = 35, SampleSize = 200,
                             omit_indices = c(11))

comp.UM.splitVec$AIC
comp.UM.splitVec$BIC
# ------------------------------------------------------------------------------


### PWE - CM models

# intSpeed ---------------------------------------------------------------------
comp.CM.intSpeed <- compare_LQC(Data = NYC.long,
                                y_var_name = "intSpeed",
                                y_label = 'Interval speed (km/h)',
                                y_range = c(3,17.5),
                                y_ticks = c(5,7.5,10,12.5,15,17.5),
                                female_box_location = c(39.5,17.2),
                                male_box_location = c(40.5,17.2),
                                qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                plot_title = 'LQC interval speed CM (comp.)',
                                storage_path = plot.dir,
                                modelType = "CM",
                                res_formula_vec = c("SexAge", "SexAgeInter"),
                                res_path_vec = c(paste0(master.dir, "/Old CM models/"),
                                                 paste0(cubic.dir,"/")),
                                res_eps_vec = c(0.05,0.1),
                                res_infix_vec = c("cubic", "CM"),
                                res_suffix_vec = c("NR.Rdata","Rdata"),
                                Age = 35, SampleSize = 200,
                                omit_indices = c())

comp.CM.intSpeed$AIC
comp.CM.intSpeed$BIC
# ------------------------------------------------------------------------------


# SplitVec ---------------------------------------------------------------------
comp.CCM.splitVec <- compare_LQC(Data = NYC.long,
                             y_var_name = "splitVec",
                             y_label = 'Splitting indicator',
                             y_range = c(0.9,1.6),
                             y_ticks = c(0.9,1,1.1,1.2,1.3,1.4,1.5,1.6),
                             female_box_location = c(37.5,1.59),
                             male_box_location = c(38.5,1.59),
                             qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                             plot_title = 'LQC splitting indicator CCM (comp.)',
                             storage_path = plot.dir,
                             modelType = "CCM",
                             res_formula_vec = c("SexAge", "SexAgeInter"),
                             res_path_vec = c(paste0(master.dir, "/Old CM models/"),
                                              paste0(cubic.dir,"/")),
                             res_eps_vec = c(0.01,0.01),
                             res_infix_vec = c("CC", "CCM"),
                             res_suffix_vec = c("NR.Rdata","Rdata"),
                             Age = 35, SampleSize = 200,
                             omit_indices = c(11))

comp.CCM.splitVec$AIC
comp.CCM.splitVec$BIC
# ------------------------------------------------------------------------------

################################################################################



##### Part 8: Comparison to UQR: construction of CI ############################
library(data.table)


## SexAge (before review)

UM_formula <- "SexAge"

# intSpeed ---------------------------------------------------------------------
# for UM
tab.intSpeed <- get.UM.quartiles.CI.df(varName = 'intSpeed',
                                       UM_eps = 0.01,
                                       UM_path = paste0(master.dir,
                                                        "/Old UM models/"),
                                       UM_infix = 'unstr',
                                       UM_suffix = 'NR.Rdata',
                                       UM_formula = 'SexAge')
tab.intSpeed <- as.data.frame(tab.intSpeed)
colnames(tab.intSpeed) <- c("interc25", "sex25","age25", "-",
                            "interc50", "sex50","age50", "--",
                            "interc75", "sex75","age75")
tab.intSpeed$distCol <- 1:11
  
# for UQR
tab.intSpeed.uqr <- get.UQR.quartiles.CI.df(varName = 'intSpeed',
                                            UM_formula = 'SexAge')
tab.intSpeed.uqr <- as.data.frame(tab.intSpeed.uqr)
colnames(tab.intSpeed.uqr) <- c("interc25", "sex25","age25", "-",
                                "interc50", "sex50","age50", "--",
                                "interc75", "sex75","age75")
tab.intSpeed.uqr$distCol <- 1:11


# combining them in one frame
tab.intSpeed.combi <- rbindlist(list(tab.intSpeed.uqr, tab.intSpeed))[order(distCol)]

# printing
xtable(tab.intSpeed.combi[,-c('distCol')])
# ------------------------------------------------------------------------------

# splitVec ---------------------------------------------------------------------

# turn off scientific notation to get correct number of digits instead of e-02 etc
options(scipen=999)

# for UM
tab.splitVec <- get.UM.quartiles.CI.df(varName = 'splitVec',
                                       UM_eps = 0.05,
                                       UM_path = paste0(master.dir,
                                                        "/Old UM models/"),
                                       UM_infix = 'unstr',
                                       UM_suffix = 'NR.Rdata',
                                       UM_formula = 'SexAge')
tab.splitVec <- as.data.frame(tab.splitVec)
colnames(tab.splitVec) <- c("interc25", "sex25","age25", "-",
                            "interc50", "sex50","age50", "--",
                            "interc75", "sex75","age75")
tab.splitVec$distCol <- 1:10

# for UQR
tab.splitVec.uqr <- get.UQR.quartiles.CI.df(varName = 'splitVec',
                                            UM_formula = 'SexAge')
tab.splitVec.uqr <- as.data.frame(tab.splitVec.uqr)
colnames(tab.splitVec.uqr) <- c("interc25", "sex25","age25", "-",
                                "interc50", "sex50","age50", "--",
                                "interc75", "sex75","age75")
tab.splitVec.uqr$distCol <- 1:10


# combining them in one frame
tab.splitVec.combi <- rbindlist(list(tab.splitVec.uqr, tab.splitVec))[order(distCol)]

# printing
xtable(tab.splitVec.combi[,-c('distCol')])

# back to default
options(scipen=0)
# ------------------------------------------------------------------------------



## SexAgeInter (after review)
UM_formula <- "SexAgeInter"

# intSpeed ---------------------------------------------------------------------
# for UM
tab.intSpeed <- get.UM.quartiles.CI.df(varName = 'intSpeed',
                                       UM_eps = 0.01,
                                       UM_path = paste0(unstr.dir, '/'),
                                       UM_infix = 'UM',
                                       UM_suffix = 'Rdata',
                                       UM_formula = 'SexAgeInter')
tab.intSpeed <- as.data.frame(tab.intSpeed)
colnames(tab.intSpeed) <- c("interc25", "sex25","age25", "sexage25","-",
                            "interc50", "sex50","age50", "sexage50", "--",
                            "interc75", "sex75","age75", "sexage75")
tab.intSpeed$distCol <- 1:11

# for UQR
tab.intSpeed.uqr <- get.UQR.quartiles.CI.df(varName = 'intSpeed',
                                            UM_formula = 'SexAgeInter')
tab.intSpeed.uqr <- as.data.frame(tab.intSpeed.uqr)
colnames(tab.intSpeed.uqr) <- c("interc25", "sex25","age25", "sexage25", "-",
                                "interc50", "sex50","age50", "sexage50", "--",
                                "interc75", "sex75","age75", "sexage75")
tab.intSpeed.uqr$distCol <- 1:11


# combining them in one frame
tab.intSpeed.combi <- rbindlist(list(tab.intSpeed.uqr, tab.intSpeed))[order(distCol)]

# printing
xtable(tab.intSpeed.combi[,-c('distCol')])
# ------------------------------------------------------------------------------


# splitVec ---------------------------------------------------------------------

# turn off scientific notation to get correct number of digits instead of e-02 etc
options(scipen=999)

# for UM
tab.splitVec <- get.UM.quartiles.CI.df(varName = 'splitVec',
                                       UM_eps = 0.1,
                                       UM_path = paste0(unstr.dir, '/'),
                                       UM_infix = 'UM',
                                       UM_suffix = 'Rdata',
                                       UM_formula = 'SexAgeInter')
tab.splitVec <- as.data.frame(tab.splitVec)
colnames(tab.splitVec) <- c("interc25", "sex25","age25", "sexage25","-",
                            "interc50", "sex50","age50", "sexage50", "--",
                            "interc75", "sex75","age75", "sexage75")
tab.splitVec$distCol <- 1:10

# for UQR
tab.splitVec.uqr <- get.UQR.quartiles.CI.df(varName = 'splitVec',
                                            UM_formula = 'SexAgeInter')
tab.splitVec.uqr <- as.data.frame(tab.splitVec.uqr)
colnames(tab.splitVec.uqr) <- c("interc25", "sex25","age25", "sexage25", "-",
                                "interc50", "sex50","age50", "sexage50", "--",
                                "interc75", "sex75","age75", "sexage75")
tab.splitVec.uqr$distCol <- 1:10


# combining them in one frame
tab.splitVec.combi <- rbindlist(list(tab.splitVec.uqr, tab.splitVec))[order(distCol)]

# printing
xtable(tab.splitVec.combi[,-c('distCol')])

# back to default
options(scipen=0)
# ------------------------------------------------------------------------------

################################################################################



##### Part 9: Quantile exploration (some sample code only) #####################


explor.dir <- paste(master.dir, 'Exploratory plots', sep = '/')
class.dir <- paste(explor.dir, 'Class curves', sep = '/')
summary.dir <- paste(explor.dir, 'Summary curves', sep = '/')



for (folder in c(explor.dir, class.dir, summary.dir)){
  if (!dir.exists(folder)){
    dir.create(folder)
  }
}


NYC.longM = NYC.long[NYC.long$Sex=='M',]
NYC.longF = NYC.long[NYC.long$Sex=='F',]



##########################################
### QUANTILE-CONDITIONED CURVES (EQCC) ###
##########################################


### EQCC Summary curves --------------------------------------------------------

# Quantile determining positions: d2 = 10k, d5 = MAR/2, d11 = MAR

# Variable: intSpeed
plot_EQCC(data_female = NYC.longF, data_male = NYC.longM,
          y_var_name = 'intSpeed', y_label = 'Interval speed (km/h)',
          qu_var_name = "time",
          qu_lvls = c(0.01, 0.05, 0.1, 0.25,0.5,0.75, 0.9, 0.95, 0.99),
          qu_determ_pos = c(2,5,11),
          plot_title = 'EQCC interval speed',
          file_path = summary.dir,
          female_box_location = c(36.5, 16.9),
          male_box_location = c(37.5,16.9),
          omit_indices = c(),
          save_separate_panels = FALSE)

# Variable: cumSpeed
plot_EQCC(data_female = NYC.longF, data_male = NYC.longM,
          y_var_name = 'cumSpeed', y_label = 'Cumulative speed (km/h)',
          qu_var_name = "time",
          qu_lvls = c(0.01, 0.05, 0.1, 0.25,0.5,0.75, 0.9, 0.95, 0.99),
          qu_determ_pos = c(2,5,11),
          plot_title = 'EQCC cumulative speed',
          file_path = summary.dir,
          female_box_location = c(36.5, 16.9),
          male_box_location = c(37.5,16.9),
          omit_indices = c(),
          save_separate_panels = FALSE)

# variable: SpeedDec
plot_EQCC(data_female = NYC.longF, data_male = NYC.longM,
          y_var_name = 'speedDec', y_label = 'Speed variation',
          qu_var_name = "time",
          qu_lvls = c(0.01, 0.05, 0.1, 0.25,0.5,0.75, 0.9, 0.95, 0.99),
          qu_determ_pos = c(2,5,11),
          plot_title = 'EQCC speed variation',
          file_path = summary.dir,
          female_box_location = c(14.6, 0.81),
          male_box_location = c(14.6,0.81),
          omit_indices = c(1),
          save_separate_panels = FALSE)

# variable: SplitVec
plot_EQCC(data_female = NYC.longF, data_male = NYC.longM,
          y_var_name = 'splitVec', y_label = 'Splitting indicator',
          qu_var_name = "time",
          qu_lvls = c(0.01, 0.05, 0.1, 0.25,0.5,0.75, 0.9, 0.95, 0.99),
          qu_determ_pos = c(2,5,11),
          plot_title = 'EQCC splitting indicator',
          file_path = summary.dir,
          female_box_location = c(34.4, 1.255),
          male_box_location = c(35.4,1.255),
          omit_indices = c(11),
          save_separate_panels = FALSE)
# ------------------------------------------------------------------------------


### EQCC Class curves ----------------------------------------------------------

### intSpeed ###
## men ## ----
panels_M_d11 <- get_panels_EQCC_class(data_sub = NYC.longM, y_var_name = 'intSpeed',
                                      y_label = 'Interval speed (km/h)',
                                      qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                      qu_determ_pos = 11,
                                      box_location = c(34, 17),
                                      omit_indices = c(),
                                      overrule_yrange = c(2.5,17.5))
panels_M_d5 <- get_panels_EQCC_class(data_sub = NYC.longM, y_var_name = 'intSpeed',
                                     y_label = 'Interval speed (km/h)',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 5,
                                     box_location = c(34, 17),
                                     omit_indices = c(),
                                     overrule_yrange = c(2.5,17.5))
panels_M_d2 <- get_panels_EQCC_class(data_sub = NYC.longM, y_var_name = 'intSpeed',
                                     y_label = 'Interval speed (km/h)',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 2,
                                     box_location = c(34, 17),
                                     omit_indices = c(),
                                     overrule_yrange = c(2.5,17.5))

all_panels <- c(panels_M_d2, panels_M_d5, panels_M_d11)
pCombi <- wrap_plots(all_panels, ncol = 3, byrow = FALSE) +
  plot_layout(axes = 'collect_x', axis_titles = 'collect') +
  plot_annotation(title = 'EQCC (class) interval speed: men',
                  theme = theme(plot.title = element_text(hjust = 0.55),
                                text = element_text(size=15)))

ggsave(path = class.dir, filename = paste0("EQCC_class_intSpeed_M_all.png"),
       width = 10, height = 15, device='png')
# ----

## women ## ----
panels_F_d11 <- get_panels_EQCC_class(data_sub = NYC.longF, y_var_name = 'intSpeed',
                                      y_label = 'Interval speed (km/h)',
                                      qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                      qu_determ_pos = 11,
                                      box_location = c(34, 17),
                                      omit_indices = c(),
                                      overrule_yrange = c(2.5,17.5))
panels_F_d5 <- get_panels_EQCC_class(data_sub = NYC.longF, y_var_name = 'intSpeed',
                                     y_label = 'Interval speed (km/h)',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 5,
                                     box_location = c(34, 17),
                                     omit_indices = c(),
                                     overrule_yrange = c(2.5,17.5))
panels_F_d2 <- get_panels_EQCC_class(data_sub = NYC.longF, y_var_name = 'intSpeed',
                                     y_label = 'Interval speed (km/h)',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 2,
                                     box_location = c(34, 17),
                                     omit_indices = c(),
                                     overrule_yrange = c(2.5,17.5))

all_panels <- c(panels_F_d2, panels_F_d5, panels_F_d11)
pCombi <- wrap_plots(all_panels, ncol = 3, byrow = FALSE) +
  plot_layout(axes = 'collect_x', axis_titles = 'collect') +
  plot_annotation(title = 'EQCC (class) interval speed: women',
                  theme = theme(plot.title = element_text(hjust = 0.55),
                                text = element_text(size=15)))

ggsave(path = class.dir, filename = paste0("EQCC_class_intSpeed_F_all.png"),
       width = 10, height = 15, device='png')
# ----


### SplitVec ###

# here, we need to overrule the y-range, otherwise too wide

## men ## ----
panels_M_d11 <- get_panels_EQCC_class(data_sub = NYC.longM, y_var_name = 'SplitVec',
                                      y_label = 'Splitting indicator',
                                      qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                      qu_determ_pos = 11,
                                      box_location = c(34, 1.57),
                                      omit_indices = c(),
                                      overrule_yrange = c(0.8,1.6))
panels_M_d5 <- get_panels_EQCC_class(data_sub = NYC.longM, y_var_name = 'SplitVec',
                                     y_label = 'Splitting indicator',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 5,
                                     box_location = c(34, 1.57),
                                     omit_indices = c(),
                                     overrule_yrange = c(0.8,1.6))
panels_M_d2 <- get_panels_EQCC_class(data_sub = NYC.longM, y_var_name = 'SplitVec',
                                     y_label = 'Splitting indicator',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 2,
                                     box_location = c(34, 1.57),
                                     omit_indices = c(),
                                     overrule_yrange = c(0.8,1.6))

all_panels <- c(panels_M_d2, panels_M_d5, panels_M_d11)
pCombi <- wrap_plots(all_panels, ncol = 3, byrow = FALSE) +
  plot_layout(axes = 'collect_x', axis_titles = 'collect') +
  plot_annotation(title = 'EQCC (class) splitting indicator: men',
                  theme = theme(plot.title = element_text(hjust = 0.55),
                                text = element_text(size=15)))

ggsave(path = class.dir, filename = paste0("EQCC_class_SplitVec_M_all.png"),
       width = 10, height = 15, device='png')
# ----

## women ## ----
panels_F_d11 <- get_panels_EQCC_class(data_sub = NYC.longF, y_var_name = 'SplitVec',
                                      y_label = 'Splitting indicator',
                                      qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                      qu_determ_pos = 11,
                                      box_location = c(34, 1.57),
                                      omit_indices = c(),
                                      overrule_yrange = c(0.8,1.6))
panels_F_d5 <- get_panels_EQCC_class(data_sub = NYC.longF, y_var_name = 'SplitVec',
                                     y_label = 'Splitting indicator',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 5,
                                     box_location = c(34, 1.57),
                                     omit_indices = c(),
                                     overrule_yrange = c(0.8,1.6))
panels_F_d2 <- get_panels_EQCC_class(data_sub = NYC.longF, y_var_name = 'SplitVec',
                                     y_label = 'Splitting indicator',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 2,
                                     box_location = c(34, 1.57),
                                     omit_indices = c(),
                                     overrule_yrange = c(0.8,1.6))

all_panels <- c(panels_F_d2, panels_F_d5, panels_F_d11)
pCombi <- wrap_plots(all_panels, ncol = 3, byrow = FALSE) +
  plot_layout(axes = 'collect_x', axis_titles = 'collect') +
  plot_annotation(title = 'EQCC (class) splitting indicator: women',
                  theme = theme(plot.title = element_text(hjust = 0.55),
                                text = element_text(size=15)))

ggsave(path = class.dir, filename = paste0("EQCC_class_SplitVec_F_all.png"),
       width = 10, height = 15, device='png')
# ----

# ------------------------------------------------------------------------------




###################################
### META-QUANTILE CURVES (EMQC) ###
###################################

### EMQC Summary curves --------------------------------------------------------

# Quantile determining positions: d11 = MAR only

plot_EMQC(data_female = NYC.longF, data_male = NYC.longM,
          y_var_name = 'qu.short', y_label = 'Quantile level',
          qu_var_name = "time",
          qu_lvls = c(0.01, 0.05, 0.1, 0.25,0.5,0.75, 0.9, 0.95, 0.99),
          qu_determ_pos = c(11),
          plot_title = 'Meta-quantile curves',
          file_path = summary.dir,
          female_box_location = c(36.5, 1.08),
          male_box_location = c(37.5,1.08),
          omit_indices = c(),
          save_separate_panels = FALSE)
# ------------------------------------------------------------------------------


### EMQC Class curves ----------------------------------------------------------

## men ## ----
panels_M_d11 <- get_panels_EQCC_class(data_sub = NYC.longM, y_var_name = 'qu.short',
                                      y_label = 'Quantile level',
                                      qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                      qu_determ_pos = 11,
                                      box_location = c(34, 1.08),
                                      omit_indices = c(),
                                      overrule_yrange = c(0,1.1))
panels_M_d5 <- get_panels_EQCC_class(data_sub = NYC.longM, y_var_name = 'qu.short',
                                     y_label = 'Quantile level',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 5,
                                     box_location = c(34, 1.08),
                                     omit_indices = c(),
                                     overrule_yrange = c(0,1.1))
panels_M_d2 <- get_panels_EQCC_class(data_sub = NYC.longM, y_var_name = 'qu.short',
                                     y_label = 'Quantile level',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 2,
                                     box_location = c(34, 1.08),
                                     omit_indices = c(),
                                     overrule_yrange = c(0,1.1))

all_panels <- c(panels_M_d2, panels_M_d5, panels_M_d11)
pCombi <- wrap_plots(all_panels, ncol = 3, byrow = FALSE) +
  plot_layout(axes = 'collect_x', axis_titles = 'collect') +
  plot_annotation(title = 'Meta-quantile (class) curves: men',
                  theme = theme(plot.title = element_text(hjust = 0.55),
                                text = element_text(size=15)))

ggsave(path = 'Thesis', filename = paste0("EMQC_class_M_all.png"),
       width = 10, height = 15, device='png')
# ----


## women ## ----
panels_F_d11 <- get_panels_EQCC_class(data_sub = NYC.longF, y_var_name = 'qu.short',
                                      y_label = 'Quantile level',
                                      qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                      qu_determ_pos = 11,
                                      box_location = c(34, 1.08),
                                      omit_indices = c(),
                                      overrule_yrange = c(0,1.1))
panels_F_d5 <- get_panels_EQCC_class(data_sub = NYC.longF, y_var_name = 'qu.short',
                                     y_label = 'Quantile level',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 5,
                                     box_location = c(34, 1.08),
                                     omit_indices = c(),
                                     overrule_yrange = c(0,1.1))
panels_F_d2 <- get_panels_EQCC_class(data_sub = NYC.longF, y_var_name = 'qu.short',
                                     y_label = 'Quantile level',
                                     qu_var_name = 'time', qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                                     qu_determ_pos = 2,
                                     box_location = c(34, 1.08),
                                     omit_indices = c(),
                                     overrule_yrange = c(0,1.1))

all_panels <- c(panels_F_d2, panels_F_d5, panels_F_d11)
pCombi <- wrap_plots(all_panels, ncol = 3, byrow = FALSE) +
  plot_layout(axes = 'collect_x', axis_titles = 'collect') +
  plot_annotation(title = 'Meta-quantile (class) curves: women',
                  theme = theme(plot.title = element_text(hjust = 0.55),
                                text = element_text(size=15)))

ggsave(path = class.dir, filename = paste0("EMQC_class_F_all.png"),
       width = 10, height = 15, device='png')
# ----

# ------------------------------------------------------------------------------


################################################################################




