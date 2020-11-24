
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
#          Generalized Additive Models
# -------------------------------------------------------------------------
############################################################################

# ----------------
#      FISH
#-----------------

### MMS

## Richness

# alternatives for the number of knots
kts <- seq (3,10)

# richness
GAM_fish_rich_MMS <- lapply (kts, function (kts)
   
   gam (EST.rich ~ 
           s(BO2_ppmean_ss, k=kts,bs="cr")+
           s(BO2_salinitymean_ss, k=kts,bs="cr")+
           inv_temp +
           s(BO_damean, k=kts,bs="cr"),
        family= Gamma (link='log'),
        data = cov_fish_MMS,
        na.action = na.exclude,
        drop.unused.levels=TRUE,
        method="REML"
   )
)

# check 
par(mfrow=c(2,2))
lapply (GAM_fish_rich_MMS, function (i)
   gam.check(i, pch=19,cex=.5))

# refit
kts.fit <- 5

GAM_fish_rich_MMS <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                           s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                           inv_temp +
                           s(BO_damean, k=kts.fit,bs="cr"),
                        family= Gamma (link='log'),
                        data = cov_fish_MMS,
                        na.action = na.exclude,
                        drop.unused.levels=TRUE,
                        method="REML"
)

summary(GAM_fish_rich_MMS)

# simplify
s1_GAM_fish_rich_MMS <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                              s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                                inv_temp,
                           family= Gamma (link='log'),
                           data = cov_fish_MMS,
                           na.action = na.exclude,
                           drop.unused.levels=TRUE,
                           method="REML")
summary(s1_GAM_fish_rich_MMS)

# compare
anova(GAM_fish_rich_MMS, s1_GAM_fish_rich_MMS,test = "Chisq")

# keep on simplifying
s2_GAM_fish_rich_MMS <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                              s(BO2_salinitymean_ss, k=kts.fit,bs="cr"),
                           family= Gamma (link='log'),
                           data = cov_fish_MMS,
                           na.action = na.exclude,
                           drop.unused.levels=TRUE,
                           method="REML")

# compare
anova(s1_GAM_fish_rich_MMS, s2_GAM_fish_rich_MMS,test = "Chisq")

# we can't simplify more
plot(s1_GAM_fish_rich_MMS,pages=1)
gam.vcomp(s1_GAM_fish_rich_MMS)


### MMS

## Functional diversity

GAM_fish_FD_MMS <- lapply (kts, function (kts)
   
   gam (FD ~ 
             s(BO2_ppmean_ss, k=kts,bs="cr")+
             s(BO2_salinitymean_ss, k=kts,bs="cr")+
             inv_temp +
             s(BO_damean, k=kts,bs="cr"),
     family= betar (link='logit'),
     data = cov_fish_MMS,
     na.action = na.exclude,
     drop.unused.levels=TRUE,
     method="REML"
     )
)

# check 
par(mfrow=c(2,2))
lapply (GAM_fish_FD_MMS, function (i)
   gam.check(i, pch=19,cex=.5))

# refit
kts.fit <- 9

GAM_fish_FD_MMS <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
           s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
           inv_temp +
           s(BO_damean, k=kts.fit,bs="cr"),
           family=betar (link='logit'),
           data = cov_fish_MMS,
           na.action = na.exclude,
           drop.unused.levels=TRUE,
           method="REML"
   )

summary(GAM_fish_FD_MMS)

# simplify
s1_GAM_fish_FD_MMS <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                                s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                                s(BO_damean, k=kts.fit,bs="cr"),
                             family= betar (link='logit'),
                             data = cov_fish_MMS,
                             na.action = na.exclude,
                             drop.unused.levels=TRUE,
                             method="REML")
summary(s1_GAM_fish_FD_MMS)

# compare
anova(GAM_fish_FD_MMS, s1_GAM_fish_FD_MMS,test = "Chisq")

# keep on simplifying
s2_GAM_fish_FD_MMS <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                                s(BO2_salinitymean_ss, k=kts.fit,bs="cr"),
                             family= betar (link='logit'),
                             data = cov_fish_MMS,
                             na.action = na.exclude,
                             drop.unused.levels=TRUE,
                             method="REML")

# compare
anova(s1_GAM_fish_FD_MMS, s2_GAM_fish_FD_MMS,test = "Chisq")

# we can't simplify more
plot(s1_GAM_fish_FD_MMS,pages=1)
gam.vcomp(s1_GAM_fish_FD_MMS)

# -----------------------------
###         ASSi
# -----------------------------

## Richness

# alternatives for the number of knots
kts <- seq (3,10)

# richness
GAM_fish_rich_ASSi <- lapply (kts, function (kts)
   
   gam (EST.rich ~ 
           s(BO2_ppmean_ss, k=kts,bs="cr")+
           s(BO2_salinitymean_ss, k=kts,bs="cr")+
           inv_temp +
           s(BO_damean, k=kts,bs="cr"),
        family= Gamma (link='log'),
        data = cov_fish_ASSi,
        na.action = na.exclude,
        drop.unused.levels=TRUE,
        method="REML"
   )
)

# check 
par(mfrow=c(2,2))
lapply (GAM_fish_rich_ASSi, function (i)
   gam.check(i, pch=19,cex=.5))

# refit
kts.fit <- 5

GAM_fish_rich_ASSi <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                             s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                             inv_temp +
                             s(BO_damean, k=kts.fit,bs="cr"),
                          family= Gamma (link='log'),
                          data = cov_fish_ASSi,
                          na.action = na.exclude,
                          drop.unused.levels=TRUE,
                          method="REML"
)

summary(GAM_fish_rich_ASSi)

# simplify
s1_GAM_fish_rich_ASSi <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                                s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                                inv_temp,
                             family= Gamma (link='log'),
                             data = cov_fish_ASSi,
                             na.action = na.exclude,
                             drop.unused.levels=TRUE,
                             method="REML")
summary(s1_GAM_fish_rich_ASSi)

# compare
anova(GAM_fish_rich_ASSi, s1_GAM_fish_rich_ASSi,test = "Chisq")

# keep on simplifying
s2_GAM_fish_rich_ASSi <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                                s(BO2_salinitymean_ss, k=kts.fit,bs="cr"),
                             family= Gamma (link='log'),
                             data = cov_fish_ASSi,
                             na.action = na.exclude,
                             drop.unused.levels=TRUE,
                             method="REML")

# compare
anova(s1_GAM_fish_rich_ASSi, s2_GAM_fish_rich_ASSi,test = "Chisq")

# we can't simplify more
plot(GAM_fish_rich_ASSi,pages=1,shade.col="gray90",shade=T)
gam.vcomp(s1_GAM_fish_rich_ASSi)


### ASSi

## Functional diversity

GAM_fish_FD_ASSi <- lapply (kts, function (kts)
   
   gam (FD ~ 
           s(BO2_ppmean_ss, k=kts,bs="cr")+
           s(BO2_salinitymean_ss, k=kts,bs="cr")+
           inv_temp +
           s(BO_damean, k=kts,bs="cr"),
        family= betar (link='logit'),
        data = cov_fish_ASSi,
        na.action = na.exclude,
        drop.unused.levels=TRUE,
        method="REML"
   )
)

# check 
par(mfrow=c(2,2))
lapply (GAM_fish_FD_ASSi, function (i)
   gam.check(i, pch=19,cex=.5))

# refit
kts.fit <- 6

GAM_fish_FD_ASSi <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                           s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                           inv_temp +
                           s(BO_damean, k=kts.fit,bs="cr"),
                        family=betar (link='logit'),
                        data = cov_fish_ASSi,
                        na.action = na.exclude,
                        drop.unused.levels=TRUE,
                        method="REML"
)

summary(GAM_fish_FD_ASSi)

# simplify
s1_GAM_fish_FD_ASSi <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                              s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                               inv_temp,
                           family= betar (link='logit'),
                           data = cov_fish_ASSi,
                           na.action = na.exclude,
                           drop.unused.levels=TRUE,
                           method="REML")
summary(s1_GAM_fish_FD_ASSi)

# compare
anova(GAM_fish_FD_ASSi, s1_GAM_fish_FD_ASSi,test = "Chisq")

# keep on simplifying
s2_GAM_fish_FD_ASSi <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                              s(BO2_salinitymean_ss, k=kts.fit,bs="cr"),
                           family= betar (link='logit'),
                           data = cov_fish_ASSi,
                           na.action = na.exclude,
                           drop.unused.levels=TRUE,
                           method="REML")

summary(s2_GAM_fish_FD_ASSi)

# compare
anova(s1_GAM_fish_FD_ASSi, s2_GAM_fish_FD_ASSi,test = "Chisq")

# we can't simplify more
plot(s1_GAM_fish_FD_ASSi,pages=1,shade.col="gray90",shade=T)
gam.vcomp(s1_GAM_fish_FD_ASSi)

# ----------------
#      BENTHOS
#-----------------

### MMS

## Richness

# alternatives for the number of knots
kts <- seq (3,10)

# richness
GAM_benthos_rich_MMS <- lapply (kts, function (kts)
   
   gam (EST.rich ~ 
           s(BO2_ppmean_ss, k=kts,bs="cr")+
           s(BO2_salinitymean_ss, k=kts,bs="cr")+
           inv_temp +
           s(BO_damean, k=kts,bs="cr"),
        family= Gamma (link='log'),
        data = cov_benthos_MMS,
        na.action = na.exclude,
        drop.unused.levels=TRUE,
        method="REML"
   )
)

# check 
par(mfrow=c(2,2))
lapply (GAM_benthos_rich_MMS, function (i)
   gam.check(i, pch=19,cex=.5))

# refit
kts.fit <- 8

GAM_benthos_rich_MMS <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                             s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                             inv_temp +
                             s(BO_damean, k=kts.fit,bs="cr"),
                          family= Gamma (link='log'),
                          data = cov_benthos_MMS,
                          na.action = na.exclude,
                          drop.unused.levels=TRUE,
                          method="REML"
)

summary(GAM_benthos_rich_MMS)

# simplify
s1_GAM_benthos_rich_MMS <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                                s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                                inv_temp,
                             family= Gamma (link='log'),
                             data = cov_benthos_MMS,
                             na.action = na.exclude,
                             drop.unused.levels=TRUE,
                             method="REML")
summary(s1_GAM_benthos_rich_MMS)

# compare
anova(GAM_benthos_rich_MMS, s1_GAM_benthos_rich_MMS,test = "Chisq")

# keep on simplifying
s2_GAM_benthos_rich_MMS <- gam (EST.rich ~ inv_temp+
                                s(BO2_salinitymean_ss, k=kts.fit,bs="cr"),
                             family= Gamma (link='log'),
                             data = cov_benthos_MMS,
                             na.action = na.exclude,
                             drop.unused.levels=TRUE,
                             method="REML")

# compare
anova(s1_GAM_benthos_rich_MMS, s2_GAM_benthos_rich_MMS,test = "Chisq")

# we can't simplify more
plot(s1_GAM_benthos_rich_MMS,pages=1)
gam.vcomp(s1_GAM_benthos_rich_MMS)


### MMS

## Functional diversity

GAM_benthos_FD_MMS <- lapply (kts, function (kts)
   
   gam (FD ~ 
           s(BO2_ppmean_ss, k=kts,bs="cr")+
           s(BO2_salinitymean_ss, k=kts,bs="cr")+
           inv_temp +
           s(BO_damean, k=kts,bs="cr"),
        family= betar (link='logit'),
        data = cov_benthos_MMS,
        na.action = na.exclude,
        drop.unused.levels=TRUE,
        method="REML"
   )
)

# check 
par(mfrow=c(2,2))
lapply (GAM_benthos_FD_MMS, function (i)
   gam.check(i, pch=19,cex=.5))

# refit
kts.fit <- 9

GAM_benthos_FD_MMS <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                           s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                           inv_temp +
                           s(BO_damean, k=kts.fit,bs="cr"),
                        family=betar (link='logit'),
                        data = cov_benthos_MMS,
                        na.action = na.exclude,
                        drop.unused.levels=TRUE,
                        method="REML"
)

summary(GAM_benthos_FD_MMS)

# simplify
s1_GAM_benthos_FD_MMS <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                              s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                              s(BO_damean, k=kts.fit,bs="cr"),
                           family= betar (link='logit'),
                           data = cov_benthos_MMS,
                           na.action = na.exclude,
                           drop.unused.levels=TRUE,
                           method="REML")
summary(s1_GAM_benthos_FD_MMS)

# compare
anova(GAM_benthos_FD_MMS, s1_GAM_benthos_FD_MMS,test = "Chisq")

# keep on simplifying
s2_GAM_benthos_FD_MMS <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                              s(BO2_salinitymean_ss, k=kts.fit,bs="cr"),
                           family= betar (link='logit'),
                           data = cov_benthos_MMS,
                           na.action = na.exclude,
                           drop.unused.levels=TRUE,
                           method="REML")

# compare
anova(s1_GAM_benthos_FD_MMS, s2_GAM_benthos_FD_MMS,test = "Chisq")

# we can't simplify more
plot(s1_GAM_benthos_FD_MMS,pages=1)
gam.vcomp(s1_GAM_benthos_FD_MMS)

# -----------------------------
###         ASSi
# -----------------------------

## Richness

# alternatives for the number of knots
kts <- seq (3,10)

# richness
GAM_benthos_rich_ASSi <- lapply (kts, function (kts)
   
   gam (EST.rich ~ 
           s(BO2_ppmean_ss, k=kts,bs="cr")+
           s(BO2_salinitymean_ss, k=kts,bs="cr")+
           inv_temp +
           s(BO_damean, k=kts,bs="cr"),
        family= Gamma (link='log'),
        data = cov_benthos_ASSi,
        na.action = na.exclude,
        drop.unused.levels=TRUE,
        method="REML"
   )
)

# check 
par(mfrow=c(2,2))
lapply (GAM_benthos_rich_ASSi, function (i)
   gam.check(i, pch=19,cex=.5))

# refit
kts.fit <- 5

GAM_benthos_rich_ASSi <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                              s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                              inv_temp +
                              s(BO_damean, k=kts.fit,bs="cr"),
                           family= Gamma (link='log'),
                           data = cov_benthos_ASSi,
                           na.action = na.exclude,
                           drop.unused.levels=TRUE,
                           method="REML"
)

summary(GAM_benthos_rich_ASSi)

# simplify
s1_GAM_benthos_rich_ASSi <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                                 s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                                 inv_temp,
                              family= Gamma (link='log'),
                              data = cov_benthos_ASSi,
                              na.action = na.exclude,
                              drop.unused.levels=TRUE,
                              method="REML")
summary(s1_GAM_benthos_rich_ASSi)

# compare
anova(GAM_benthos_rich_ASSi, s1_GAM_benthos_rich_ASSi,test = "Chisq")

# keep on simplifying
s2_GAM_benthos_rich_ASSi <- gam (EST.rich ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                                 s(BO2_salinitymean_ss, k=kts.fit,bs="cr"),
                              family= Gamma (link='log'),
                              data = cov_benthos_ASSi,
                              na.action = na.exclude,
                              drop.unused.levels=TRUE,
                              method="REML")

# compare
anova(s1_GAM_benthos_rich_ASSi, s2_GAM_benthos_rich_ASSi,test = "Chisq")

# we can't simplify more
plot(GAM_benthos_rich_ASSi,pages=1,shade.col="gray90",shade=T)
gam.vcomp(s1_GAM_benthos_rich_ASSi)


### ASSi

## Functional diversity

GAM_benthos_FD_ASSi <- lapply (kts, function (kts)
   
   gam (FD ~ 
           s(BO2_ppmean_ss, k=kts,bs="cr")+
           s(BO2_salinitymean_ss, k=kts,bs="cr")+
           inv_temp +
           s(BO_damean, k=kts,bs="cr"),
        family= betar (link='logit'),
        data = cov_benthos_ASSi,
        na.action = na.exclude,
        drop.unused.levels=TRUE,
        method="REML"
   )
)

# check 
par(mfrow=c(2,2))
lapply (GAM_benthos_FD_ASSi, function (i)
   gam.check(i, pch=19,cex=.5))

# refit
kts.fit <- 6

GAM_benthos_FD_ASSi <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                            s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                            inv_temp +
                            s(BO_damean, k=kts.fit,bs="cr"),
                         family=betar (link='logit'),
                         data = cov_benthos_ASSi,
                         na.action = na.exclude,
                         drop.unused.levels=TRUE,
                         method="REML"
)

summary(GAM_benthos_FD_ASSi)

# simplify
s1_GAM_benthos_FD_ASSi <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                               s(BO2_salinitymean_ss, k=kts.fit,bs="cr")+
                               inv_temp,
                            family= betar (link='logit'),
                            data = cov_benthos_ASSi,
                            na.action = na.exclude,
                            drop.unused.levels=TRUE,
                            method="REML")
summary(s1_GAM_benthos_FD_ASSi)

# compare
anova(GAM_benthos_FD_ASSi, s1_GAM_benthos_FD_ASSi,test = "Chisq")

# keep on simplifying
s2_GAM_benthos_FD_ASSi <- gam (FD ~ s(BO2_ppmean_ss, k=kts.fit,bs="cr")+
                               s(BO2_salinitymean_ss, k=kts.fit,bs="cr"),
                            family= betar (link='logit'),
                            data = cov_benthos_ASSi,
                            na.action = na.exclude,
                            drop.unused.levels=TRUE,
                            method="REML")

summary(s2_GAM_benthos_FD_ASSi)

# compare
anova(s1_GAM_benthos_FD_ASSi, s2_GAM_benthos_FD_ASSi,test = "Chisq")

# we can't simplify more
plot(s1_GAM_benthos_FD_ASSi,pages=1,shade.col="gray90",shade=T)
gam.vcomp(s1_GAM_benthos_FD_ASSi)


#















## newdata to predict
newd <- data.frame (BO2_ppmean_ss = seq (-2,2,0.1),
                    BO2_salinitymean_ss=0)
# predict
p <- predict.gam(s2_GAM_fish_FD_ASSi,
                 newdata=newd,
                 type="link",
                 se.fit=T)
upr <- p$fit + (1.96 * p$fit)
lwr <- p$fit - (1.96 * p$fit)

plot(newd$BO2_ppmean_ss, 
     p$fit,
     lwd=2,type="l",
     ylim = c(-10,4))
lines (newd$BO2_ppmean_ss,
       upr, lwd=1,col="gray")
lines (newd$BO2_ppmean_ss,
       lwr, lwd=1,col="gray")















## newdata to predict
newd <- data.frame (BO2_ppmean_ss = seq (-2,2,0.1),
                    BO2_salinitymean_ss=0)
# predict
p <- predict.gam(s2_GAM_fish_FD_ASSi,
                 newdata=newd,
                 type="link",
                 se.fit=T)
upr <- p$fit + (1.96 * p$fit)
lwr <- p$fit - (1.96 * p$fit)

plot(newd$BO2_ppmean_ss, 
     p$fit,
     lwd=2,type="l",
     ylim = c(-10,4))
lines (newd$BO2_ppmean_ss,
       upr, lwd=1,col="gray")
lines (newd$BO2_ppmean_ss,
       lwr, lwd=1,col="gray")
