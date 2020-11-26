
# ----------------------------------------------------------------------------#
#    routine to modeling fish and benthos SR and FD  relative to environment
#                          using GAM
# ----------------------------------------------------------------------------#

## chamar os pacotes
source("R/packages.R")
source("R/functions.R")

# ------------------------------------------ #
# Load rarefied data of fishes and benthos
# sites x transect / video x species

load (here ("output","random_composition_bentos.RData"))
load (here ("output","random_composition_fish.RData"))

# ------------------------------------------ #
#
#     Load results of the rarefaction

load (here ("output","FD_fish_ASSi.RData"))
load (here ("output","FD_fish_MSS.RData"))
load (here ("output","FD_benthos_ASSi.RData"))
load (here ("output","FD_benthos_MSS.RData"))

# ------------------------------------------ #
#
#     Load environment data

load (here ("output","env_data.RData"))

# minimum sample size
### fish 
FRic_fish_MSS <- do.call (cbind ,sapply (
   sapply (FD_fish_MSS, "[", "Fdindexes"),"[", "FRic")
   )
# get the average
av_FRic_fish_MSS <- apply (FRic_fish_MSS,1,mean,na.rm=T)

### benthos
FRic_benthos_MSS <- do.call (cbind, sapply (
   sapply (FD_benthos_MSS, "[", "Fdindexes"),"[", "FRic")
   )
# get the average
av_FRic_benthos_MSS <- apply (FRic_benthos_MSS,1,mean,na.rm=T)

# asyntotic sample size
### fish
FRic_fish_ASSi <- do.call (cbind, sapply (
   sapply (FD_fish_ASSi, "[", "Fdindexes"),"[", "FRic")
   )
# get the average
av_FRic_fish_ASSi <- apply (FRic_fish_ASSi,1,mean,na.rm=T)

### benthos
FRic_benthos_ASSi <- do.call (cbind, sapply (
   sapply (FD_benthos_ASSi, "[", "Fdindexes"),"[", "FRic")
   )
# get the average
av_FRic_benthos_ASSi <- apply (FRic_benthos_ASSi,1,mean,na.rm=T)

# ----------------------------------------- #
#            connect with covariates
# ------------------------------------------#
#---------------
# fish under MMS
#---------------
site_fish <- unlist(lapply (strsplit (rownames(covariates_site$sea_data[[1]]), "\\."), function (i)
   i[2]))

cov_fish_MMS <- covariates_site$sea_data[[1]][match (res_table_samples$Site,site_fish),]
# bind richness estimate
cov_fish_MMS <- data.frame(cov_fish_MMS, EST.rich = res_table_samples$EST.Rich)
# bind FD estimate
cov_fish_MMS <- cbind(cov_fish_MMS, FD = av_FRic_fish_MSS)

#---------------
# fish under ASSi
#---------------
cov_fish_ASSi <- covariates_site$sea_data[[1]][match (res_sp_accum_fish_asymptote$Site[which(is.na(res_sp_accum_fish_asymptote$EST.Rich)!= T)],
                                                      site_fish),]
# bind richness estimate
cov_fish_ASSi <- data.frame(cov_fish_ASSi, 
                       EST.rich = res_sp_accum_fish_asymptote$EST.Rich[which(is.na(res_sp_accum_fish_asymptote$EST.Rich)!= T)])
# bind FD estimate
cov_fish_ASSi <- cbind(cov_fish_ASSi, FD = av_FRic_fish_ASSi)

#-------------------
# benthos under MMS
#-------------------
site_benthos <- unlist(lapply (strsplit (rownames(covariates_site$sea_data[[2]]), "\\."), function (i)
   i[2]))

cov_benthos_MMS <- covariates_site$sea_data[[2]][match (res_table_samples_bentos$Site,site_benthos),]
# bind richness estimate
cov_benthos_MMS <- data.frame(cov_benthos_MMS, EST.rich = res_table_samples_bentos$EST.Rich)
# bind FD estimate
cov_benthos_MMS <- cbind(cov_benthos_MMS, FD = av_FRic_benthos_MSS)

#-------------------
# benthos under ASSi
#-------------------
cov_benthos_ASSi <- covariates_site$sea_data[[2]][match (res_sp_accum_bentos_asymptote$Site[which(is.na(res_sp_accum_bentos_asymptote$EST.Rich)!= T)],
                                                      site_benthos),]
# bind richness estimate
cov_benthos_ASSi <- data.frame(cov_benthos_ASSi, 
                            EST.rich = res_sp_accum_bentos_asymptote$EST.Rich[which(is.na(res_sp_accum_bentos_asymptote$EST.Rich)!= T)])
# bind FD estimate
cov_benthos_ASSi <- cbind(cov_benthos_ASSi, FD = av_FRic_benthos_ASSi)

# ------------------------
# inverse of temperature, standardize covarites, insert region and locality
# ------------------------

boltzmann_factor <- 8.62e-5

# FISH
## MMS
temp_kelvin<- cov_fish_MMS$BO2_tempmean_ss + 273.15
inv_temp <- 1 / boltzmann_factor * (1 / mean(temp_kelvin) - 1 / temp_kelvin)
cov_fish_MMS <- cbind(cov_fish_MMS,
                      inv_temp=inv_temp)

## standardize other variables
cov_fish_MMS$BO2_tempmean_ss <- (cov_fish_MMS$BO2_tempmean_ss-mean(cov_fish_MMS$BO2_tempmean_ss))/sd(cov_fish_MMS$BO2_tempmean_ss)
cov_fish_MMS$BO2_ppmean_ss <- (cov_fish_MMS$BO2_ppmean_ss-mean(cov_fish_MMS$BO2_ppmean_ss))/sd(cov_fish_MMS$BO2_ppmean_ss)
cov_fish_MMS$BO2_salinitymean_ss <- (cov_fish_MMS$BO2_salinitymean_ss-mean(cov_fish_MMS$BO2_salinitymean_ss))/sd(cov_fish_MMS$BO2_salinitymean_ss)
cov_fish_MMS$BO_damean <- (cov_fish_MMS$BO_damean-mean(cov_fish_MMS$BO_damean))/sd(cov_fish_MMS$BO_damean)
cov_fish_MMS$inv_temp <- (cov_fish_MMS$inv_temp-mean(cov_fish_MMS$inv_temp))/sd(cov_fish_MMS$inv_temp)

## region
cov_fish_MMS <- cbind(cov_fish_MMS,
                      Region=covariates_site$region [match(rownames(cov_fish_MMS),rownames(covariates_site$region))])
cov_fish_MMS$Region <- as.factor (cov_fish_MMS$Region)

## locality
locality <- unlist(lapply (strsplit (covariates_site$site_names,"\\."), function (i) i[1]))
cov_fish_MMS <- cbind(cov_fish_MMS,
                      Locality=locality [match(rownames(cov_fish_MMS),rownames(covariates_site$region))])
cov_fish_MMS$Locality <- as.factor (cov_fish_MMS$Locality)


## ASSi
temp_kelvin_ASSi<- cov_fish_ASSi$BO2_tempmean_ss + 273.15
inv_temp_ASSi <- 1 / boltzmann_factor * (1 / mean(temp_kelvin_ASSi) - 1 / temp_kelvin_ASSi)
cov_fish_ASSi <- cbind(cov_fish_ASSi,
                      inv_temp=inv_temp_ASSi)

## standardize other variables
cov_fish_ASSi$BO2_tempmean_ss <- (cov_fish_ASSi$BO2_tempmean_ss-mean(cov_fish_ASSi$BO2_tempmean_ss))/sd(cov_fish_ASSi$BO2_tempmean_ss)
cov_fish_ASSi$BO2_ppmean_ss <- (cov_fish_ASSi$BO2_ppmean_ss-mean(cov_fish_ASSi$BO2_ppmean_ss))/sd(cov_fish_ASSi$BO2_ppmean_ss)
cov_fish_ASSi$BO2_salinitymean_ss <- (cov_fish_ASSi$BO2_salinitymean_ss-mean(cov_fish_ASSi$BO2_salinitymean_ss))/sd(cov_fish_ASSi$BO2_salinitymean_ss)
cov_fish_ASSi$BO_damean <- (cov_fish_ASSi$BO_damean-mean(cov_fish_ASSi$BO_damean))/sd(cov_fish_ASSi$BO_damean)
cov_fish_ASSi$inv_temp <- (cov_fish_ASSi$inv_temp-mean(cov_fish_ASSi$inv_temp))/sd(cov_fish_ASSi$inv_temp)

## region
cov_fish_ASSi <- cbind(cov_fish_ASSi,
                      Region=covariates_site$region [match(rownames(cov_fish_ASSi),rownames(covariates_site$region))])
cov_fish_ASSi$Region <- as.factor (cov_fish_ASSi$Region)

## locality
locality <- unlist(lapply (strsplit (covariates_site$site_names,"\\."), function (i) i[1]))
cov_fish_ASSi <- cbind(cov_fish_ASSi,
                      Locality=locality [match(rownames(cov_fish_ASSi),rownames(covariates_site$region))])
cov_fish_ASSi$Locality <- as.factor (cov_fish_ASSi$Locality)

# benthos
## MMS
temp_kelvinB <- cov_benthos_MMS$BO2_tempmean_ss + 273.15
inv_tempB <- 1 / boltzmann_factor * (1 / mean(temp_kelvinB) - 1 / temp_kelvinB)
cov_benthos_MMS <- cbind(cov_benthos_MMS,
                      inv_temp=inv_tempB)

## standardize other variables
cov_benthos_MMS$BO2_tempmean_ss <- (cov_benthos_MMS$BO2_tempmean_ss-mean(cov_benthos_MMS$BO2_tempmean_ss))/sd(cov_benthos_MMS$BO2_tempmean_ss)
cov_benthos_MMS$BO2_ppmean_ss <- (cov_benthos_MMS$BO2_ppmean_ss-mean(cov_benthos_MMS$BO2_ppmean_ss))/sd(cov_benthos_MMS$BO2_ppmean_ss)
cov_benthos_MMS$BO2_salinitymean_ss <- (cov_benthos_MMS$BO2_salinitymean_ss-mean(cov_benthos_MMS$BO2_salinitymean_ss))/sd(cov_benthos_MMS$BO2_salinitymean_ss)
cov_benthos_MMS$BO_damean <- (cov_benthos_MMS$BO_damean-mean(cov_benthos_MMS$BO_damean))/sd(cov_benthos_MMS$BO_damean)
cov_benthos_MMS$inv_temp <- (cov_benthos_MMS$inv_temp-mean(cov_benthos_MMS$inv_temp))/sd(cov_benthos_MMS$inv_temp)

## region
cov_benthos_MMS <- cbind(cov_benthos_MMS,
                       Region=covariates_site$region [match(rownames(cov_benthos_MMS),rownames(covariates_site$region))])
cov_benthos_MMS$Region <- as.factor (cov_benthos_MMS$Region)

## locality
locality <- unlist(lapply (strsplit (covariates_site$site_names,"\\."), function (i) i[1]))
cov_benthos_MMS <- cbind(cov_benthos_MMS,
                       Locality=locality [match(rownames(cov_benthos_MMS),rownames(covariates_site$region))])
cov_benthos_MMS$Locality <- as.factor (cov_benthos_MMS$Locality)

## ASSi
temp_kelvinB_ASSi <- cov_benthos_ASSi$BO2_tempmean_ss + 273.15
inv_tempB_ASSi <- 1 / boltzmann_factor * (1 / mean(temp_kelvinB_ASSi) - 1 / temp_kelvinB_ASSi)
cov_benthos_ASSi <- cbind(cov_benthos_ASSi,
                         inv_temp=inv_tempB_ASSi)

## standardize other variables
cov_benthos_ASSi$BO2_tempmean_ss <- (cov_benthos_ASSi$BO2_tempmean_ss-mean(cov_benthos_ASSi$BO2_tempmean_ss))/sd(cov_benthos_ASSi$BO2_tempmean_ss)
cov_benthos_ASSi$BO2_ppmean_ss <- (cov_benthos_ASSi$BO2_ppmean_ss-mean(cov_benthos_ASSi$BO2_ppmean_ss))/sd(cov_benthos_ASSi$BO2_ppmean_ss)
cov_benthos_ASSi$BO2_salinitymean_ss <- (cov_benthos_ASSi$BO2_salinitymean_ss-mean(cov_benthos_ASSi$BO2_salinitymean_ss))/sd(cov_benthos_ASSi$BO2_salinitymean_ss)
cov_benthos_ASSi$BO_damean <- (cov_benthos_ASSi$BO_damean-mean(cov_benthos_ASSi$BO_damean))/sd(cov_benthos_ASSi$BO_damean)
cov_benthos_ASSi$inv_temp <- (cov_benthos_ASSi$inv_temp-mean(cov_benthos_ASSi$inv_temp))/sd(cov_benthos_ASSi$inv_temp)

## region
cov_benthos_ASSi <- cbind(cov_benthos_ASSi,
                      Region=covariates_site$region [match(rownames(cov_benthos_ASSi),rownames(covariates_site$region))])
cov_benthos_ASSi$Region <- as.factor (cov_benthos_ASSi$Region)
                      
## locality
locality <- unlist(lapply (strsplit (covariates_site$site_names,"\\."), function (i) i[1]))
cov_benthos_ASSi <- cbind(cov_benthos_ASSi,
                      Locality=locality [match(rownames(cov_benthos_ASSi),rownames(covariates_site$region))])
cov_benthos_ASSi$Locality <- as.factor (cov_benthos_ASSi$Locality)

############################################################################
# -------------------------------------------------------------------------
#          Non-linear models
#
#           RICHNESS
#
# -------------------------------------------------------------------------
############################################################################

# ----------------
#      FISH
#-----------------

# MMS

## Richness
model1 <- glm (EST.rich ~ 
        poly(BO2_ppmean_ss, 2)+
        poly(BO2_salinitymean_ss, 2)+
        inv_temp +
        poly(BO_damean, 2),
     family= Gamma (link='log'),
     data = cov_fish_MMS,
     na.action = na.exclude)

# 
summary (model1)

# simplify
model2 <- glm (EST.rich ~ 
                  poly(BO2_salinitymean_ss, 2)+
                  inv_temp +
                  poly(BO_damean, 2),
               family= Gamma (link='log'),
               data = cov_fish_MMS,
               na.action = na.exclude)

summary(model2)

# compare them
anova (model1,model2,test = "Chisq")


# simplify again
model3 <- glm (EST.rich ~ 
                  poly(BO2_salinitymean_ss, 2)+
                  inv_temp +
                  BO_damean,
               family= Gamma (link='log'),
               data = cov_fish_MMS,
               na.action = na.exclude)

summary(model3)

#
anova (model2,model3, test = "Chisq")

# simplify again
model4 <- glm (EST.rich ~ 
                  poly(BO2_salinitymean_ss, 1)+
                  inv_temp+
                  BO_damean,
               family= Gamma (link='log'),
               data = cov_fish_MMS,
               na.action = na.exclude)

summary(model4)

#
anova (model3,model4, test = "Chisq")

## we can't simplify anymore

## diagnose model fit
plot(model3) 

## check coeffs and other statistics
summary (model3)

## unadjusted R2
# help here: https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r
round (with(summary(model3), 1 - deviance/null.deviance),2)

## newd
newd <- data.frame (BO2_salinitymean_ss= 0,
                    inv_temp=0,
                    BO_damean=seq (range(cov_fish_MMS$BO_damean)[1],
                                   range(cov_fish_MMS$BO_damean)[2],
                                   0.05))

## data to predict (based on the model)
pred.vals <- predict (model3,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$BO_damean, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,37),
     ylab = "Estimated fish richness",
     xlab = "Turbidity")
lines (exp(upr) ~  newd$BO_damean,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$BO_damean,lwd=2,col="gray50")

#

## newd
newd <- data.frame (BO2_salinitymean_ss= 0,
                    inv_temp=seq (range(cov_fish_MMS$inv_temp)[1],
                                  range(cov_fish_MMS$inv_temp)[2],
                                  0.05),
                    BO_damean=0)

## data to predict (based on the model)
pred.vals <- predict (model3,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$inv_temp, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,37),
     ylab = "Estimated fish richness",
     xlab = "Temperature")
lines (exp(upr) ~  newd$inv_temp,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$inv_temp,lwd=2,col="gray50")


## newd
newd <- data.frame (BO2_salinitymean_ss= seq (range(cov_fish_MMS$BO2_salinitymean_ss)[1],
                                              range(cov_fish_MMS$BO2_salinitymean_ss)[2],
                                              0.05),
                    inv_temp=0,
                    BO_damean=0)

## data to predict (based on the model)
pred.vals <- predict (model3,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$BO2_salinitymean_ss, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,37),
     ylab = "Estimated fish richness",
     xlab = "Salinity")
lines (exp(upr) ~  newd$BO2_salinitymean_ss,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$BO2_salinitymean_ss,lwd=2,col="gray50")


# -------------------------
#           ASSi

## Richness
model1 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss, 2)+
                  poly(BO2_salinitymean_ss, 2)+
                  inv_temp +
                  poly(BO_damean, 2),
               family= Gamma (link='log'),
               data = cov_fish_ASSi,
               na.action = na.exclude)

# 
summary (model1)

# simplify
model2 <- glm (EST.rich ~ poly(BO2_ppmean_ss, 2)+
                  poly(BO2_salinitymean_ss, 2)+
                  inv_temp,
               family= Gamma (link='log'),
               data = cov_fish_ASSi,
               na.action = na.exclude)

summary(model2)

# compare them
anova (model1,model2,test = "Chisq")

# simplify again
model3 <- glm (EST.rich ~ BO2_ppmean_ss +
                   poly(BO2_salinitymean_ss, 2)+
                   inv_temp,
                family= Gamma (link='log'),
                data = cov_fish_ASSi,
                na.action = na.exclude)

summary(model3)

#
anova (model2,model3, test = "Chisq")

# simplify again
model4 <- glm (EST.rich ~ BO2_ppmean_ss +
                  poly(BO2_salinitymean_ss, 1)+
                  inv_temp,
               family= Gamma (link='log'),
               data = cov_fish_ASSi,
               na.action = na.exclude)

summary(model4)

#
anova (model3,model4, test = "Chisq")

## we can't simplify anymore

## diagnose model fit
plot(model3) 

## check coeffs and other statistics
summary (model3)

## unadjusted R2
# help here: https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r
round (with(summary(model3), 1 - deviance/null.deviance),2)

## newd
newd <- data.frame (BO2_salinitymean_ss= 0,
                    inv_temp=0,
                    BO2_ppmean_ss=seq (range(cov_fish_MMS$BO_damean)[1],
                                   range(cov_fish_MMS$BO_damean)[2],
                                   0.05))

## data to predict (based on the model)
pred.vals <- predict (model3,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$BO2_ppmean_ss, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,90),
     ylab = "Estimated fish richness",
     xlab = "Productivity")
lines (exp(upr) ~  newd$BO2_ppmean_ss,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$BO2_ppmean_ss,lwd=2,col="gray50")

#

## newd
newd <- data.frame (BO2_salinitymean_ss= 0,
                    inv_temp=seq (range(cov_fish_MMS$inv_temp)[1],
                                  range(cov_fish_MMS$inv_temp)[2],
                                  0.05),
                    BO2_ppmean_ss=0)

## data to predict (based on the model)
pred.vals <- predict (model3,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$inv_temp, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,150),
     ylab = "Estimated fish richness",
     xlab = "Temperature")
lines (exp(upr) ~  newd$inv_temp,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$inv_temp,lwd=2,col="gray50")


## newd
newd <- data.frame (BO2_salinitymean_ss= seq (range(cov_fish_MMS$BO2_salinitymean_ss)[1],
                                              range(cov_fish_MMS$BO2_salinitymean_ss)[2],
                                              0.05),
                    inv_temp=0,
                    BO2_ppmean_ss=0)

## data to predict (based on the model)
pred.vals <- predict (model3,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$BO2_salinitymean_ss, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,70),
     ylab = "Estimated fish richness",
     xlab = "Salinity")
lines (exp(upr) ~  newd$BO2_salinitymean_ss,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$BO2_salinitymean_ss,lwd=2,col="gray50")

# ----------------
#      BENTHOS
#-----------------

# MMS

## Richness
model1 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss, 2)+
                  poly(BO2_salinitymean_ss, 2)+
                  inv_temp +
                  poly(BO_damean, 2),
               family= Gamma (link='log'),
               data = cov_benthos_MMS,
               na.action = na.exclude)

# 
summary (model1)

# simplify
model2 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss, 2)+
                  poly(BO2_salinitymean_ss, 2)+
                  inv_temp,
               family= Gamma (link='log'),
               data = cov_benthos_MMS,
               na.action = na.exclude)

summary(model2)

# compare them
anova (model1,model2,test = "Chisq")


# simplify again
model3 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss, 2)+
                  poly(BO2_salinitymean_ss, 2),
               family= Gamma (link='log'),
               data = cov_benthos_MMS,
               na.action = na.exclude)

summary(model3)

#
anova (model2,model3, test = "Chisq")

# simplify again
model4 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss, 2)+
                  I(BO2_salinitymean_ss^2),
               family= Gamma (link='log'),
               data = cov_benthos_MMS,
               na.action = na.exclude)

summary(model4)

#
anova (model3,model4, test = "Chisq")

## we can't simplify anymore

## diagnose model fit
plot(model4) 

## check coeffs and other statistics
summary (model4)

## unadjusted R2
# help here: https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r
round (with(summary(model4), 1 - deviance/null.deviance),2)

## newd
newd <- data.frame (BO2_salinitymean_ss= 0,
                    BO2_ppmean_ss=seq (range(cov_benthos_MMS$BO2_ppmean_ss)[1],
                                   range(cov_benthos_MMS$BO2_ppmean_ss)[2],
                                   0.05))

## data to predict (based on the model)
pred.vals <- predict (model4,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$BO2_ppmean_ss, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(5,25),
     ylab = "Estimated benthos richness",
     xlab = "Productivity")
lines (exp(upr) ~  newd$BO2_ppmean_ss,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$BO2_ppmean_ss,lwd=2,col="gray50")

#

## newd
newd <- data.frame (BO2_salinitymean_ss= seq (range(cov_benthos_MMS$BO2_salinitymean_ss)[1],
                                              range(cov_benthos_MMS$BO2_salinitymean_ss)[2],
                                              0.05),
                    BO2_ppmean_ss=0)

## data to predict (based on the model)
pred.vals <- predict (model4,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$BO2_salinitymean_ss, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(5,25),
     ylab = "Estimated fish richness",
     xlab = "Salinity")
lines (exp(upr) ~  newd$BO2_salinitymean_ss,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$BO2_salinitymean_ss,lwd=2,col="gray50")

# -------------------------
#      ASSi

## Richness
model1 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss, 2)+
                  poly(BO2_salinitymean_ss, 2)+
                  inv_temp +
                  poly(BO_damean, 2),
               family= Gamma (link='log'),
               data = cov_benthos_ASSi,
               na.action = na.exclude)

# 
summary (model1)

# simplify
model2 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss, 2)+
                 inv_temp +
                  poly(BO_damean, 2),
               family= Gamma (link='log'),
               data = cov_benthos_ASSi,
               na.action = na.exclude)

summary(model2)

# compare them
anova (model1,model2,test = "Chisq")


# simplify again
model3 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss, 2)+
                  inv_temp ,
               family= Gamma (link='log'),
               data = cov_benthos_ASSi,
               na.action = na.exclude)

summary(model3)

#
anova (model2,model3, test = "Chisq")

# simplify again
model4 <- glm (EST.rich ~ 
                  I(BO2_ppmean_ss^2)+
                  inv_temp ,
               family= Gamma (link='log'),
               data = cov_benthos_ASSi,
               na.action = na.exclude)

summary(model4)

#
anova (model3,model4, test = "Chisq")

# simplify again
model5 <- glm (EST.rich ~ 
                  inv_temp ,
               family= Gamma (link='log'),
               data = cov_benthos_ASSi,
               na.action = na.exclude)

summary(model5)

#
anova (model4,model5, test = "Chisq")

## we can't simplify anymore

## diagnose model fit
plot(model5) 

## check coeffs and other statistics
summary (model5)

## unadjusted R2
# help here: https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r
round (with(summary(model5), 1 - deviance/null.deviance),2)

## newd
newd <- data.frame (inv_temp= seq (range(cov_benthos_ASSi$inv_temp)[1],
                                       range(cov_benthos_ASSi$inv_temp)[2],
                                       0.05))

## data to predict (based on the model)
pred.vals <- predict (model5,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$inv_temp, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(5,35),
     ylab = "Estimated benthos richness",
     xlab = "Temperature")
lines (exp(upr) ~  newd$inv_temp,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$inv_temp,lwd=2,col="gray50")

############################################################################
# -------------------------------------------------------------------------
#          Non-linear models
#
#           FUNCTIONAL DIVERSITY
#
# -------------------------------------------------------------------------
############################################################################

# ----------------
#      FISH
#-----------------

# MMS

## FD
model1 <- glm (FD ~ poly(BO2_ppmean_ss, 2)+
                  poly(BO2_salinitymean_ss,2)+
                  inv_temp +
                  poly(BO_damean, 2),
               family = gaussian(link="identity"),
               data = cov_fish_MMS,
               na.action = na.exclude)
#
summary (model1)

# simplify
model2 <- glm (FD ~ poly(BO2_salinitymean_ss,2)+
                  inv_temp +
                  poly(BO_damean, 2),
               family = gaussian(link="identity"),
               data = cov_fish_MMS,
               na.action = na.exclude)

summary(model2)
#
anova (model1,model2, test = "Chisq")

# simplify
model3 <- glm (FD ~ inv_temp +
                  poly(BO_damean, 2),
               family = gaussian(link="identity"),
               data = cov_fish_MMS,
               na.action = na.exclude)

summary(model3)

# compare them
anova(model2,
      model3,test="Chisq")

# simplify again
model4 <- glm (FD ~ inv_temp +
                  BO_damean,
               family = gaussian(link="identity"),
               data = cov_fish_MMS,
               na.action = na.exclude)

summary(model4)

#
anova(model3, model4,test="Chisq")

## we can't simplify anymore

## diagnose model fit
plot(model3,pages=1) 

## check coeffs and other statistics
## pseudo R2
summary (model3)

## newd
newd <- data.frame (inv_temp=0,
                    BO_damean=seq (range(cov_fish_MMS$BO_damean)[1],
                                   range(cov_fish_MMS$BO_damean)[2],
                                   0.05))

## data to predict (based on the model)
## help here : https://stats.stackexchange.com/questions/230501/variance-vs-standard-deviation-in-beta-regression/230681#230681
pred.resp <- predict (model3,
                      newdata = newd,
                      type= "link",
                      se.fit=T)

upr <- pred.resp$fit+(2*pred.resp$se.fit)
lwr <- pred.resp$fit-(2*pred.resp$se.fit)

# plotting 

plot(pred.resp$fit ~ newd$BO_damean, ## plogis is because the inverse of log-link in its plogisonential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(0,0.3),
     ylab = "Estimated functional diversity",
     xlab = "Turbidity")
lines (upr ~  newd$BO_damean,lwd=2,col="gray50")
lines (lwr ~  newd$BO_damean,lwd=2,col="gray50")
points (cov_fish_MMS$FD ~ cov_fish_MMS$BO_damean,
        col=rgb(0,0,0.2,alpha=0.1),pch=19)

#

## newd
newd <- data.frame (BO_damean=0,
                    inv_temp=seq (range(cov_fish_MMS$inv_temp)[1],
                                   range(cov_fish_MMS$inv_temp)[2],
                                   0.05))

## data to predict (based on the model)
## help here : https://stats.stackexchange.com/questions/230501/variance-vs-standard-deviation-in-beta-regression/230681#230681
pred.resp <- predict (model3,
                      newdata = newd,
                      type= "link",
                      se.fit=T)

upr <- pred.resp$fit+(2*pred.resp$se.fit)
lwr <- pred.resp$fit-(2*pred.resp$se.fit)

# plotting 

plot(pred.resp$fit ~ newd$inv_temp, ## plogis is because the inverse of log-link in its plogisonential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(0,0.3),
     ylab = "Estimated functional diversity",
     xlab = "Temperature")
lines (upr ~  newd$inv_temp,lwd=2,col="gray50")
lines (lwr ~  newd$inv_temp,lwd=2,col="gray50")

points (cov_fish_MMS$FD ~ cov_fish_MMS$inv_temp,
        col=rgb(0,0,0.2,alpha=0.2),pch=19)

## newd
newd <- data.frame (BO2_ppmean_ss=seq (range(cov_fish_MMS$BO2_ppmean_ss)[1],
                                       range(cov_fish_MMS$BO2_ppmean_ss)[2],
                                       0.05),
                    BO2_salinitymean_ss=0,
                    BO_damean=0)

## data to predict (based on the model)
## help here : https://stats.stackexchange.com/questions/230501/variance-vs-standard-deviation-in-beta-regression/230681#230681
pred.resp <- predict (model2,
                      newdata = newd,
                      type= "link",
                      se.fit=T)

upr <- pred.resp$fit+(2*pred.resp$se.fit)
lwr <- pred.resp$fit-(2*pred.resp$se.fit)

# plotting 

plot(plogis(pred.resp$fit) ~ newd$BO2_ppmean_ss, ## plogis is because the inverse of log-link in its plogisonential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(0,0.3),
     ylab = "Estimated fish functional diversity",
     xlab = "Productivity")
lines (plogis(upr) ~  newd$BO2_ppmean_ss,lwd=2,col="gray50")
lines (plogis(lwr) ~  newd$BO2_ppmean_ss,lwd=2,col="gray50")
