##########################################################################################################################################
##########################################################################################################################################
##### Code for Parasite species co-occurrence patterns on North American red squirrels (Tamiasciurus hudsonicus), Veitch et al. 2024 #####
##########################################################################################################################################
##########################################################################################################################################

###Boral analysis
library(nlme)
library(tidyverse)
library(boral)
library(corrplot)
library(ggboral)

RS  <- read.table("RS-subsetted data.csv", header=T, na.strings=c("","NA"), sep = ",", fill=TRUE)
RS <- RS %>%
  mutate(Year = as.factor(Year), Sex = as.factor(Sex), ID = as.factor(ID), ReproCond3 = as.factor(ReproCond3),
         TagLeft = as.factor(TagLeft), Date = as.numeric(Date), Weight = as.numeric(Weight))

RS2 <- RS[!is.na(RS$Last.Capture),]

m <- lme(Orchopeas.caedensP ~ Last.Capture + Date + Year, random = ~1|ID, data = RS2)
summary(m)

m <- lme(Ceratophyllus.visonP ~ Last.Capture + Date + Year, random = ~1|ID, data = RS2)
summary(m)

m <- lme(Orange.Mite ~ Last.Capture + Date + Year, random = ~1|ID, data = RS2)
summary(m)

###Using Warten et al. "so many variables.." starter code to fit multivariate GLMM
RSprev <- RS %>% select(Orange.Mite, Orchopeas.caedensP, Ceratophyllus.visonP)
RScovX <- cbind(scale(RS[,1]), RS[,3])
colnames(RScovX) <- c("2Date", "1Sex")
#Creating matrix for ID to use as a random effect
random <- RS %>% select(ID)
random <- data.matrix(random, rownames.force = FALSE)
#RScovX2 <- cbind(RS[,14])
#colnames(RScovX2) <- "ID"
#Date, sex, weight, PD, repro, Year, ID
#RStraits <- as.matrix(RS[,c(2,14)])
#RS_which_traits <- vector("list",ncol(RScovX))
#for(i in 1:length(RS_which_traits)) 
#  RS_which_traits[[i]] <- 1:3
##Unconstrained model
fit.lvm2 <- boral(y = RSprev, family = "binomial", num.lv = NULL,
                  lv.control = list(num.lv = 2), save.model = TRUE, row.eff = "random",
                  ranef.ids = random,
                  mcmc.control = list(n.burnin = 200000, n.iteration = 300000, n.thin = 100,
                                      seed = 123))
plot(fit.lvm2) ##Checking diagnostics
summary(fit.lvm2) # To look at estimated parameter values
fit.lvm2$hpdintervals # 95% credible intervals for model parameters.
##Constrained model
fit.lvm <- boral(y = RSprev, X = RScovX, family = "binomial", num.lv = NULL,
                 lv.control = list(num.lv = 2), save.model = TRUE, row.eff = "random",
                 ranef.ids = random,
                 mcmc.control = list(n.burnin = 200000, n.iteration = 300000, n.thin = 100,
                                     seed = 123))
plot(fit.lvm) ##Checking diagnostics
summary(fit.lvm) # To look at estimated parameter values
fit.lvm$hpdintervals # 95% credible intervals for model parameters.

###SSVS###
priors <- list(ssvs.index = c(0,0))
fit.lvm.test <- boral(y = RSprev, X = RScovX, family = "binomial", num.lv = NULL,
                      lv.control = list(num.lv = 2), row.eff = "random", ranef.ids = random,
                      prior.control = priors,
                      mcmc.control = list(n.burnin = 200000, n.iteration = 300000, n.thin = 100,
                                          seed = 123))
summary(fit.lvm.test$ssvs.indcoefs.mean)
#Based on SSVS, I need to remove:
#(1) body mass
#(3) year
#(4) repro
#This is because their mean is <0.5

##Correlations between species due to environmental responses
envcors <- get.enviro.cor(fit.lvm)
rescors <- get.residual.cor(fit.lvm)
#estimates & CI
rescors[["cor"]]
rescors[["cor.lower"]]
rescors[["cor.upper"]]
envcors[["cor"]]
envcors[["cor.lower"]]
envcors[["cor.upper"]]
#plots of correlations
corrplot(envcors$sig.cor, type = "lower", diag = FALSE, mar = c(3,0.5,2,1), tl.srt = 45)
corrplot(rescors$sig.cor, type = "lower", diag = FALSE, mar = c(3,0.5,2,1), tl.srt = 45)
#lvsplot(fit.lvm2)
#lvsplot(fit.lvm)
summary(rescors$trace)
rescors2 <- get.residual.cor(fit.lvm2)
summary(rescors2$trace)
##Env covariates accounted for ~14% of covariation between species
gg_coefsplot(fit.lvm, palette = NULL, single.colour = "Black") + theme_test(base_size=15)
gg_varpart(fit.lvm) + theme_test()

##Figures
pdf("Fig 1.pdf", width = 8, height = 6)
dev.off()

#Diagnostics (geweke, Dunn-Smyth residual plots)
fit.lvm$geweke.diag
par(mfrow=c(2,2))
plot.boral(fit.lvm)