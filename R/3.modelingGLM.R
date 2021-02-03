
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

cov_fish_MMS <- covariates_site$sea_data[[1]][match (res_table_samples$Site,sites_fish_complete),]
# bind richness estimate
cov_fish_MMS <- data.frame(cov_fish_MMS, EST.rich = res_table_samples$EST.Rich)
# bind FD estimate
cov_fish_MMS <- cbind(cov_fish_MMS, FD = av_FRic_fish_MSS)
# bind coordinates
cov_fish_MMS <- cbind(cov_fish_MMS,covariates_site$coord$coord_peixes [match(rownames(cov_fish_MMS),covariates_site$coord$coord_peixes$Group.1),2:3])

#---------------
# fish under ASSi
#---------------
cov_fish_ASSi <- covariates_site$sea_data[[1]][match (res_sp_accum_fish_asymptote$Site[which(is.na(res_sp_accum_fish_asymptote$EST.Rich)!= T)],
                                                      sites_fish_complete),]
# bind richness estimate
cov_fish_ASSi <- data.frame(cov_fish_ASSi, 
                       EST.rich = res_sp_accum_fish_asymptote$EST.Rich[which(is.na(res_sp_accum_fish_asymptote$EST.Rich)!= T)])
# bind FD estimate
cov_fish_ASSi <- cbind(cov_fish_ASSi, FD = av_FRic_fish_ASSi)
# bind coordinates
cov_fish_ASSi <- cbind(cov_fish_ASSi,covariates_site$coord$coord_peixes [match(rownames(cov_fish_ASSi),covariates_site$coord$coord_peixes$Group.1),2:3])

#-------------------
# benthos under MMS
#-------------------

cov_benthos_MMS <- covariates_site$sea_data[[2]][match (res_table_samples_bentos$Site,sites_bentos_complete),]
# bind richness estimate
cov_benthos_MMS <- data.frame(cov_benthos_MMS, EST.rich = res_table_samples_bentos$EST.Rich)
# bind FD estimate
cov_benthos_MMS <- cbind(cov_benthos_MMS, FD = av_FRic_benthos_MSS)
# bind coordinates
cov_benthos_MMS <- cbind(cov_benthos_MMS,
                         covariates_site$coord$coord_bentos [match(rownames(cov_benthos_MMS),covariates_site$coord$coord_bentos$Group.1),2:3])

#-------------------
# benthos under ASSi
#-------------------
cov_benthos_ASSi <- covariates_site$sea_data[[2]][match (res_sp_accum_bentos_asymptote$Site[which(is.na(res_sp_accum_bentos_asymptote$EST.Rich)!= T)],
                                                         sites_bentos_complete),]
# bind richness estimate
cov_benthos_ASSi <- data.frame(cov_benthos_ASSi, 
                            EST.rich = res_sp_accum_bentos_asymptote$EST.Rich[which(is.na(res_sp_accum_bentos_asymptote$EST.Rich)!= T)])
# bind FD estimate
cov_benthos_ASSi <- cbind(cov_benthos_ASSi, FD = av_FRic_benthos_ASSi)
# bind coordienates
cov_benthos_ASSi <- cbind(cov_benthos_ASSi,
                         covariates_site$coord$coord_bentos [match(rownames(cov_benthos_ASSi),covariates_site$coord$coord_bentos$Group.1),2:3])

# ------------------------
# inverse of temperature, standardize covarites, insert region and locality
# ------------------------

boltzmann_factor <- 8.62e-5
a<-20:25 + 273.15
plot(a, (1/boltzmann_factor * (1/mean (a) - 1/(a))))

# FISH
## MMS
temp_kelvin<- cov_fish_MMS$BO2_tempmean_ss + 273.15
inv_temp <- 1 / boltzmann_factor * (1 / mean(temp_kelvin) - 1 / temp_kelvin)
cov_fish_MMS <- cbind(cov_fish_MMS,
                      inv_temp=inv_temp)

## standardize other variables
cov_fish_MMS$BO2_tempmean_ss_std <- (cov_fish_MMS$BO2_tempmean_ss-mean(cov_fish_MMS$BO2_tempmean_ss))/sd(cov_fish_MMS$BO2_tempmean_ss)
cov_fish_MMS$BO2_ppmean_ss_std <- (cov_fish_MMS$BO2_ppmean_ss-mean(cov_fish_MMS$BO2_ppmean_ss))/sd(cov_fish_MMS$BO2_ppmean_ss)
cov_fish_MMS$BO2_salinitymean_ss_std <- (cov_fish_MMS$BO2_salinitymean_ss-mean(cov_fish_MMS$BO2_salinitymean_ss))/sd(cov_fish_MMS$BO2_salinitymean_ss)
cov_fish_MMS$BO_damean_std <- (cov_fish_MMS$BO_damean-mean(cov_fish_MMS$BO_damean))/sd(cov_fish_MMS$BO_damean)
cov_fish_MMS$inv_temp_std <- (cov_fish_MMS$inv_temp-mean(cov_fish_MMS$inv_temp))/sd(cov_fish_MMS$inv_temp)

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
cov_fish_ASSi$BO2_tempmean_ss_std <- (cov_fish_ASSi$BO2_tempmean_ss-mean(cov_fish_ASSi$BO2_tempmean_ss))/sd(cov_fish_ASSi$BO2_tempmean_ss)
cov_fish_ASSi$BO2_ppmean_ss_std <- (cov_fish_ASSi$BO2_ppmean_ss-mean(cov_fish_ASSi$BO2_ppmean_ss))/sd(cov_fish_ASSi$BO2_ppmean_ss)
cov_fish_ASSi$BO2_salinitymean_ss_std <- (cov_fish_ASSi$BO2_salinitymean_ss-mean(cov_fish_ASSi$BO2_salinitymean_ss))/sd(cov_fish_ASSi$BO2_salinitymean_ss)
cov_fish_ASSi$BO_damean_std <- (cov_fish_ASSi$BO_damean-mean(cov_fish_ASSi$BO_damean))/sd(cov_fish_ASSi$BO_damean)
cov_fish_ASSi$inv_temp_std <- (cov_fish_ASSi$inv_temp-mean(cov_fish_ASSi$inv_temp))/sd(cov_fish_ASSi$inv_temp)

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
cov_benthos_MMS$BO2_tempmean_ss_std <- (cov_benthos_MMS$BO2_tempmean_ss-mean(cov_benthos_MMS$BO2_tempmean_ss))/sd(cov_benthos_MMS$BO2_tempmean_ss)
cov_benthos_MMS$BO2_ppmean_ss_std <- (cov_benthos_MMS$BO2_ppmean_ss-mean(cov_benthos_MMS$BO2_ppmean_ss))/sd(cov_benthos_MMS$BO2_ppmean_ss)
cov_benthos_MMS$BO2_salinitymean_ss_std <- (cov_benthos_MMS$BO2_salinitymean_ss-mean(cov_benthos_MMS$BO2_salinitymean_ss))/sd(cov_benthos_MMS$BO2_salinitymean_ss)
cov_benthos_MMS$BO_damean_std <- (cov_benthos_MMS$BO_damean-mean(cov_benthos_MMS$BO_damean))/sd(cov_benthos_MMS$BO_damean)
cov_benthos_MMS$inv_temp_std <- (cov_benthos_MMS$inv_temp-mean(cov_benthos_MMS$inv_temp))/sd(cov_benthos_MMS$inv_temp)

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
cov_benthos_ASSi$BO2_tempmean_ss_std <- (cov_benthos_ASSi$BO2_tempmean_ss-mean(cov_benthos_ASSi$BO2_tempmean_ss))/sd(cov_benthos_ASSi$BO2_tempmean_ss)
cov_benthos_ASSi$BO2_ppmean_ss_std <- (cov_benthos_ASSi$BO2_ppmean_ss-mean(cov_benthos_ASSi$BO2_ppmean_ss))/sd(cov_benthos_ASSi$BO2_ppmean_ss)
cov_benthos_ASSi$BO2_salinitymean_ss_std <- (cov_benthos_ASSi$BO2_salinitymean_ss-mean(cov_benthos_ASSi$BO2_salinitymean_ss))/sd(cov_benthos_ASSi$BO2_salinitymean_ss)
cov_benthos_ASSi$BO_damean_std <- (cov_benthos_ASSi$BO_damean-mean(cov_benthos_ASSi$BO_damean))/sd(cov_benthos_ASSi$BO_damean)
cov_benthos_ASSi$inv_temp_std <- (cov_benthos_ASSi$inv_temp-mean(cov_benthos_ASSi$inv_temp))/sd(cov_benthos_ASSi$inv_temp)

## region
cov_benthos_ASSi <- cbind(cov_benthos_ASSi,
                      Region=covariates_site$region [match(rownames(cov_benthos_ASSi),rownames(covariates_site$region))])
cov_benthos_ASSi$Region <- as.factor (cov_benthos_ASSi$Region)
                      
## locality
locality <- unlist(lapply (strsplit (covariates_site$site_names,"\\."), function (i) i[1]))
cov_benthos_ASSi <- cbind(cov_benthos_ASSi,
                      Locality=locality [match(rownames(cov_benthos_ASSi),rownames(covariates_site$region))])
cov_benthos_ASSi$Locality <- as.factor (cov_benthos_ASSi$Locality)

## correlation between variables
cor(cov_fish_MMS[,16:19])
cor(cov_fish_ASSi[,16:19])
cor(cov_benthos_MMS[,16:19])
cor(cov_benthos_ASSi[,16:19])

## aveerage of covariates
apply(cov_fish_MMS [,c(1,3,5,8,14)],2,mean)
apply(cov_fish_MMS [,c(1,3,5,8,14)],2,sd)

#
apply(cov_fish_ASSi [,c(1,3,5,8,14)],2,mean)
apply(cov_fish_ASSi [,c(1,3,5,8,14)],2,sd)
#
apply(cov_benthos_MMS [,c(1,3,5,8,14)],2,mean)
apply(cov_benthos_MMS [,c(1,3,5,8,14)],2,sd)
#
apply(cov_benthos_ASSi [,c(1,3,5,8,14)],2,mean)
apply(cov_benthos_ASSi [,c(1,3,5,8,14)],2,sd)

#############################################################################
## maps of estimated values

agg_fish_MMS <- aggregate (cov_fish_MMS, by = list(cov_fish_MMS$Locality), FUN="mean")

## maps
# mapa mundi
world <- ne_countries(scale = "medium", returnclass = "sf")

# cortar o mapa para ver a america do Sul e parte da central
wm <- ggplot() + 
   geom_sf (data=world, size = 0.1, 
            fill= "gray90",colour="gray90") +
   coord_sf (xlim = c(-50, -25),  ylim = c(-30, 4), expand = FALSE) +
   theme_bw() + #xlab ("Longitude")  + ylab ("Latitude") +
   theme(panel.border = element_blank(), 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "aliceblue",#darkslategray1
                                         colour = "aliceblue"),
         axis.text.x = element_text(size=6),
         axis.ticks.x=element_line(size=1),
         axis.text.y = element_text(size=6),
         axis.ticks.y=element_line(size=1),
         axis.title.x = element_text(size=10),
         axis.title.y = element_text(size=10),
         title = element_blank(),
         plot.margin = unit(c(0,-0.8,0,0.3), "cm")) +
   xlab("Longitude") + ylab("Latitude")
                            

fish_MMS_map <- wm + geom_point(data = cov_fish_MMS, aes (x=Lon,y=Lat),
                size=3, col = "red",alpha=0.2)
                #size=cov_fish_MMS$EST.rich*0.1)

fish_MMS_ASSi_map <- fish_MMS_map + 
   geom_point(data = cov_fish_ASSi, aes (x=Lon+1,y=Lat+0.1),
              size=3, col = "darkred",alpha=0.5)
                          #size=cov_fish_ASSi$EST.rich*0.05)

benthos_MMS_map<- fish_MMS_ASSi_map + 
   geom_point(data = cov_benthos_MMS, aes (x=Lon,y=Lat+1),
              size=3, col = "green",alpha=0.5)

benthos_MMS_ASSi_map<- benthos_MMS_map + 
   geom_point(data = cov_benthos_ASSi, aes (x=Lon+1,y=Lat+1),
              size=3, col = "darkgreen",alpha=0.5)

#

ggsave (here ("output","vectorized","benthos_MMS_ASSi_map.pdf"))

############################################################################
# -------------------------------------------------------------------------
#          Non-linear models
#
#        Environmental drivers
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
               poly(BO2_ppmean_ss_std,2)+
               poly(BO2_salinitymean_ss_std,2)+
               inv_temp_std +
               poly(BO_damean_std,2),
     family= Gamma (link='log'),
     data = cov_fish_MMS,
     na.action = na.exclude)

# 
summary (model1)

# simplify
model2 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss_std,2)+
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std ,
               family= Gamma (link='log'),
               data = cov_fish_MMS,
               na.action = na.exclude)

summary(model2)

# compare them
anova (model1,model2,test = "Chisq")
AIC(model1,
    model2)

# simplify again
model3 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss_std,2)+
                  inv_temp_std ,
               family= Gamma (link='log'),
               data = cov_fish_MMS,
               na.action = na.exclude)

summary(model3)

#
anova (model2,model3, test = "Chisq")
AIC(model2,model3)

# simplify again
model4 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss_std, 2),
               family= Gamma (link='log'),
               data = cov_fish_MMS,
               na.action = na.exclude)

summary(model4)

#
anova (model3,model4, test = "Chisq")
AIC(model3,model4)
## we can't simplify anymore

## diagnose model fitpar(mfrow=c(2,2))
par(mfrow=c(2,2))
plot(model3) 

## check coeffs and other statistics
summary (model3)

## unadjusted R2
# help here: https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r
round (with(summary(model3), 1 - deviance/null.deviance),2)

## newd
pdf (here("output","vectorized", "fish_MMS.pdf"),height=4,width=7)

plot(allEffects(model3))

dev.off()

## newd
newd <- data.frame (BO2_ppmean_ss_std= 0,
                    inv_temp_std= seq (range(cov_fish_MMS$inv_temp_std)[1],
                                      range(cov_fish_MMS$inv_temp_std)[2],
                                      0.01))

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
     ylim=c(10,40),
     ylab = "",
     xlab = "Inverse of temperature",
     xaxt='n')
axis(1,
     at=cov_fish_MMS$inv_temp_std,
     round(cov_fish_MMS$inv_temp,2))

lines (exp(upr) ~  newd$inv_temp_std,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$inv_temp_std,lwd=2,col="gray50")
points (cov_fish_MMS$EST.rich ~ cov_fish_MMS$inv_temp_std,
        col=rgb(0,0,0.2,alpha=0.1),pch=19)


## newd
newd <- data.frame (BO2_ppmean_ss_std= seq (range(cov_fish_MMS$BO2_ppmean_ss_std)[1],
                                                  range(cov_fish_MMS$BO2_ppmean_ss_std)[2],
                                                  0.0001),
                    inv_temp_std=0)

## data to predict (based on the model)
pred.vals <- predict (model3,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$BO2_ppmean_ss_std, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,70),
     ylab = "",
     xlab = "Salinity",
     xaxt='n')
axis(1,
     at=cov_fish_MMS$BO2_ppmean_ss_std,
     round(cov_fish_MMS$BO2_ppmean_ss,2))

lines (exp(upr) ~  newd$BO2_ppmean_ss_std,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$BO2_ppmean_ss_std,lwd=2,col="gray50")
points (cov_fish_MMS$EST.rich ~ cov_fish_MMS$BO2_ppmean_ss_std,
        col=rgb(0,0,0.2,alpha=0.1),pch=19)


# -------------------------
#           ASSi

## Richness
model1 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss_std,2)+
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std +
                  poly(BO_damean_std,2),
               family= Gamma (link='log'),
               data = cov_fish_ASSi,
               na.action = na.exclude)

# 
summary (model1)

# simplify
model2 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss_std,2)+
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std ,
               family= Gamma (link='log'),
               data = cov_fish_ASSi,
               na.action = na.exclude)

summary(model2)

# compare them
anova (model1,model2,test = "Chisq")
AIC(model1,model2)

# simplify again
model3 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss_std, 2)+
                  inv_temp_std,
               family= Gamma (link='log'),
               data = cov_fish_ASSi,
               na.action = na.exclude)

summary(model3)

#
anova (model2,model3, test = "Chisq")
AIC(model2,model3)

# simplify again
model4 <- glm (EST.rich ~ poly(BO2_ppmean_ss_std, 2),
               family= Gamma (link='log'),
               data = cov_fish_ASSi,
               na.action = na.exclude)
#
summary(model4)

#
anova (model3,model4, test = "Chisq")
AIC(model3,model4)
## we can't simplify anymore

## diagnose model fit
par(mfrow=c(2,2))
plot(model4) 

## check coeffs and other statistics
summary (model4)

## unadjusted R2
# help here: https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r
round (with(summary(model4), 1 - deviance/null.deviance),2)

## newd
pdf (here("output","vectorized", "fish_ASSi.pdf"),height=4,width=4)

plot(allEffects(model4))

dev.off()


# ----------------
#      BENTHOS
#-----------------

# MMS

## Richness
model1 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss_std,2)+
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std +
                  poly(BO_damean_std,2),
               family= Gamma (link='log'),
               data = cov_benthos_MMS,
               na.action = na.exclude)

# 
summary (model1)

# simplify
model2 <-  glm (EST.rich ~ 
                   poly(BO2_ppmean_ss_std,2)+
                   poly(BO2_salinitymean_ss_std,2)+
                   poly(BO_damean_std,2),
                family= Gamma (link='log'),
                data = cov_benthos_MMS,
                na.action = na.exclude)

summary(model2)

# compare them
anova (model1,model2,test = "Chisq")
AIC(model1,model2)

# simplify again
model3 <-glm (EST.rich ~ 
                 poly(BO2_ppmean_ss_std,2)+
                 poly(BO2_salinitymean_ss_std,2),
              family= Gamma (link='log'),
              data = cov_benthos_MMS,
              na.action = na.exclude)

summary(model3)

#
anova (model2,model3, test = "Chisq")
AIC(model2,model3)

# simplify again
model4 <- glm (EST.rich ~ 
                  poly(BO2_salinitymean_ss_std,2),
               family= Gamma (link='log'),
               data = cov_benthos_MMS,
               na.action = na.exclude)

summary(model4)

#
anova (model3,model4, test = "Chisq")
AIC(model3,model4)

## we can't simplify anymore

## diagnose model fit
par(mfrow=c(2,2))
plot(model4) 

## check coeffs and other statistics
summary (model4)

## unadjusted R2
# help here: https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r
round (with(summary(model4), 1 - deviance/null.deviance),2)

## newd
pdf (here ("output","vectorized","benthos_MMS.pdf"),height=4,width=4)

plot(allEffects(model4))

dev.off()

# -------------------------
#      ASSi

## Richness
model1 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss_std,2)+
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std +
                  poly(BO_damean_std,2),
               family= Gamma (link='log'),
               data = cov_benthos_ASSi,
               na.action = na.exclude)

# 
summary (model1)


# simplify
model2 <- glm (EST.rich ~ 
                  poly(BO2_ppmean_ss_std,2)+
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std ,
               family= Gamma (link='log'),
               data = cov_benthos_ASSi,
               na.action = na.exclude)

summary(model2)

# compare them
anova (model1,model2,test = "Chisq")
AIC(model1,model2)

# simplify again
model3 <- glm (EST.rich ~ 
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std ,
               family= Gamma (link='log'),
               data = cov_benthos_ASSi,
               na.action = na.exclude)

summary(model3)

#
anova (model2,model3, test = "Chisq")
AIC(model2,model3)

# simplify again
model4 <- glm (EST.rich ~ 
                  inv_temp_std,
               family= Gamma (link='log'),
               data = cov_benthos_ASSi,
               na.action = na.exclude)

summary(model4)

#
anova (model3,model4, test = "Chisq")
AIC(model3,model4)

## we can't simplify anymore

## diagnose model fit
plot(model4) 

## check coeffs and other statistics
summary (model4)

## unadjusted R2
# help here: https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r
round (with(summary(model4), 1 - deviance/null.deviance),2)

## newd
pdf (here("output","vectorized","benthos_ASSi.pdf"),height=4,width=4)

plot(allEffects(model4))

dev.off()

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
model1 <- betareg ((FD) ~ poly(BO2_ppmean_ss_std, 2)+
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std +
                  poly(BO_damean_std, 2),
               link = 'logit',
               data = cov_fish_MMS,
               na.action = na.exclude)
#
summary (model1)

# simplify
model2 <- betareg ((FD) ~ poly(BO2_ppmean_ss_std, 2)+
                      inv_temp_std +
                      poly(BO2_salinitymean_ss_std, 2),
                   link = 'logit',
                   data = cov_fish_MMS,
                   na.action = na.exclude)

summary(model2)
#
AIC (model1,model2)

# simplify
model3 <- betareg ((FD) ~ inv_temp_std +
                      poly(BO2_salinitymean_ss_std, 2),
                   link = 'logit',
                   data = cov_fish_MMS,
                   na.action = na.exclude)


summary(model3)

# compare them
AIC(model2,
      model3)

# simplify again
model4 <-  betareg ((FD) ~ poly(BO2_salinitymean_ss_std, 2),
                    link = 'logit',
                    data = cov_fish_MMS,
                    na.action = na.exclude)
summary(model4)

#
AIC(model3, model4)

## we can't simplify anymore

## diagnose model fit
plot(model3) 

## check coeffs and other statistics
## pseudo R2
summary(model3)

## newd
pdf(here ("output","vectorized", "fish_MMS_FD.pdf"),heigh=4,width=7)

plot(allEffects(model3))

dev.off()

########################
# ASSi

## FD
model1 <- betareg ((FD) ~ poly(BO2_ppmean_ss_std, 2)+
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std +
                  poly(BO_damean_std, 2),
               link = "logit",
               data = cov_fish_ASSi,
               na.action = na.exclude)
#
summary (model1)

# simplify
model2 <- betareg ((FD) ~  poly(BO2_salinitymean_ss_std,2)+
                      inv_temp_std +
                      poly(BO_damean_std, 2),
                   link = "logit",
                   data = cov_fish_ASSi,
                   na.action = na.exclude)

summary(model2)

#
anova (model1,model2, test = "Chisq")
AIC(model1,model2)

# simplify
model3 <- betareg ((FD) ~  poly(BO_damean_std, 2)+
                      inv_temp_std,
                   link = "logit",
                   data = cov_fish_ASSi,
                   na.action = na.exclude)

summary(model3)

# compare them
AIC(model2,
      model3)

# simplify again
model4 <- betareg ((FD) ~  inv_temp_std,
                   link = "logit",
                   data = cov_fish_ASSi,
                   na.action = na.exclude)

summary(model4)

#
AIC(model3, model4)


## we can't simplify anymore

## diagnose model fit
plot(model4,pages=1) 

## check coeffs and other statistics
## pseudo R2
round (with(summary(model4), 1 - deviance/null.deviance),2)
summary(model4)

## newd
pdf(here ("output","vectorized","fish_ASSi_FD.pdf"),heigh=4,width=4)

plot(allEffects(model4))

dev.off()

########################
# benthos

# MMS

model1 <- betareg ((FD) ~ poly(BO2_ppmean_ss_std, 2)+
                      poly(BO2_salinitymean_ss_std,2)+
                      inv_temp_std +
                      poly(BO_damean_std, 2),
                   link = "logit",
                   data = cov_benthos_MMS,
                   na.action = na.exclude)
summary (model1)

# simplify

model2 <- betareg ((FD) ~ poly(BO2_ppmean_ss_std, 2)+
                      poly(BO2_salinitymean_ss_std,2)+
                      inv_temp_std ,
                   link = "logit",
                   data = cov_benthos_MMS,
                   na.action = na.exclude)
summary(model2)

#
anova (model1,model2, test = "Chisq")
AIC(model1,model2)


# simplify

model3 <- betareg ((FD) ~ poly(BO2_salinitymean_ss_std,2)+
                      inv_temp_std ,
                   link = "logit",
                   data = cov_benthos_MMS,
                   na.action = na.exclude)
summary(model3)

# compared
AIC(model2,model3)

plot(model2)

pdf(here ("output","vectorized","benthos_MMS_FD_betareg.pdf"))

plot(allEffects(model2))

dev.off()

###############
## ASSi

model1 <- betareg ((FD) ~ poly(BO2_ppmean_ss_std, 2)+
                  poly(BO2_salinitymean_ss_std,2)+
                  inv_temp_std +
                  poly(BO_damean_std, 2),
               link = "logit",
               data = cov_benthos_ASSi,
               na.action = na.exclude)

summary(model1)

# simplify
model2 <- betareg ((FD) ~ poly(BO_damean_std, 2)+
                      poly(BO2_salinitymean_ss_std,2)+
                      inv_temp_std,
                   link = "logit",
                   data = cov_benthos_ASSi,
                   na.action = na.exclude)
summary(model2)
#
anova (model1,model2, test = "Chisq")
AIC(model1,model2)

# simplify
model3 <- betareg ((FD) ~ poly(BO2_salinitymean_ss_std,2)+
                      poly(BO_damean_std, 2),
                   link = "logit",
                   data = cov_benthos_ASSi,
                   na.action = na.exclude)
summary(model3)

AIC(model2,model3)

# simplify

model4 <- betareg ((FD) ~ poly(BO_damean_std, 2),
                   link = "logit",
                   data = cov_benthos_ASSi,
                   na.action = na.exclude)
summary(model4)

AIC(model3,model4)

## we can't simplify anymore

## diagnose model fit
plot(model4) 

## check coeffs and other statistics
## pseudo R2
summary(model4)

pdf(here ("output","vectorized","benthos_ASSi_FD_betareg.pdf"),heigh=4,width=4)

plot(allEffects(model4))

dev.off()

##------------------------------------------------- #
# 
# Relationship between fish fd and bethos fd
#
# --------------------------------------------------#

# help here
# https://stats.stackexchange.com/questions/33013/what-test-can-i-use-to-compare-slopes-from-two-or-more-regression-models
# MSS
data.ancova.MSS <- rbind(data.frame (EST.rich=cov_benthos_MMS$EST.rich,
            FD=cov_benthos_MMS$FD,
            Organism="Benthos"),
data.frame (EST.rich=cov_fish_MMS$EST.rich,
            FD=cov_fish_MMS$FD,
            Organism="Fishes"))

model.ancova <-glm (FD ~ EST.rich*Organism,
                    data.ancova.MSS,
     family = gaussian (link="identity"))
summary(model.ancova)
anova(model.ancova)

m.lst <- emtrends(model.ancova, "Organism", var="EST.rich")
m.lst
# ASSi
data.ancova.ASSi <- rbind(data.frame (EST.rich=cov_benthos_ASSi$EST.rich,
                                     FD=cov_benthos_ASSi$FD,
                                     Organism="Benthos"),
                         data.frame (EST.rich=cov_fish_ASSi$EST.rich,
                                     FD=cov_fish_ASSi$FD,
                                     Organism="Fishes"))

model.ancova.ASSi <-glm (FD ~ EST.rich*Organism,
                    data.ancova.ASSi,
                    family = gaussian (link="identity"))
summary(model.ancova.ASSi)
anova(model.ancova.ASSi)

m.lst.ASSi <- emtrends(model.ancova.ASSi, "Organism", var="EST.rich")
m.lst.ASSi

## plottign results
pdf(here("output","vectorized","rel_fish_benthos.pdf"),height = 6,width=8)
par (mfrow = c(1,2),mar = c(7,4,7,1))
## MSS 
m1 <- glm (FD ~ EST.rich,
                data = cov_benthos_MMS,
          family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (EST.rich = seq(5.5,25,0.5))
pred.vals <- predict (m1,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 
plot(NA,
     ylim=c(-0.03,0.8),
     xlim=c(0,50),
     ylab = "Functional diversity",
     xlab = "Species richness",
     main = "Minimum Sample Size (MSS)")
lines((pred.vals$fit) ~ newd$EST.rich, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "green",
     lwd =3)

#
lines (upr ~  newd$EST.rich,lwd=2,col="gray50")
lines (lwr ~  newd$EST.rich,lwd=2,col="gray50")
points (cov_benthos_MMS$FD ~ cov_benthos_MMS$EST.rich,
        col=rgb(0,0.2,0,alpha=0.4),pch=19)

m2 <- glm (FD ~ poly(EST.rich,2),
           data = cov_fish_MMS,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (EST.rich = seq(12,50,0.5))
pred.vals <- predict (m2,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

lines((pred.vals$fit) ~ newd$EST.rich, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "red",
     type="l",
     lwd =3,
     ylab = "Functional diversity",
     xlab = "Species richness")

#
lines (upr ~  newd$EST.rich,lwd=2,col="gray50")
lines (lwr ~  newd$EST.rich,lwd=2,col="gray50")
points (cov_fish_MMS$FD ~ cov_fish_MMS$EST.rich,
        col=rgb(0.2,0,0,alpha=0.4),pch=19)

# ASSi

m1 <- glm (FD ~ EST.rich,
           data = cov_benthos_ASSi,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (EST.rich = seq(7,26,0.5))
pred.vals <- predict (m1,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 
plot(NA,
     ylim=c(-0.03,0.8),
     xlim=c(0,90),
     ylab = "",
     xlab = "Species richness",
     main = "Asymptotic Sample Size (ASSi)")
lines((pred.vals$fit) ~ newd$EST.rich, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "darkgreen",
      lwd =3)

#
lines (upr ~  newd$EST.rich,lwd=2,col="gray50")
lines (lwr ~  newd$EST.rich,lwd=2,col="gray50")
points (cov_benthos_ASSi$FD ~ cov_benthos_ASSi$EST.rich,
        col=rgb(0,0.2,0,alpha=0.4),pch=19)

m2 <- glm (FD ~ poly(EST.rich,2),
           data = cov_fish_ASSi,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (EST.rich = seq(14,90,0.5))
pred.vals <- predict (m2,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

lines((pred.vals$fit) ~ newd$EST.rich, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "darkred",
      type="l",
      lwd =3)

#
lines (upr ~  newd$EST.rich,lwd=2,col="gray50")
lines (lwr ~  newd$EST.rich,lwd=2,col="gray50")
points (cov_fish_ASSi$FD ~ cov_fish_ASSi$EST.rich,
        col=rgb(0.2,0,0,alpha=0.4),pch=19)


legend ("topleft", legend = c("Benthos", "Fishes"),
        lwd = 3, col =c("darkgreen","darkred"),
        bty="n")

dev.off()

# ---------------------------------------------------
# latitudinal variation

# bind richness data from MSS and ASSi
# fish
bind.rich.data <- rbind (data.frame (Lat= cov_fish_MMS$Lat,
                                     EST.rich = cov_fish_MMS$EST.rich,
                                     Organism="Fishes"),
                         data.frame (Lat= cov_benthos_MMS$Lat,
                                     EST.rich = cov_benthos_MMS$EST.rich,
                                     Organism="Benthos"))

## test of slopes

model.SR.MMS <-glm (EST.rich ~ poly(Lat,2)*Organism,
                       bind.rich.data,
                         family = gaussian (link="identity"))
summary(model.SR.MMS)
anova(model.SR.MMS)

m.lst.SR.MMS <- emtrends(model.SR.MMS, "Organism", var="Lat")
m.lst.SR.MMS

## plottign results
pdf(here("output","vectorized","rel_lat_fish_benthos.pdf"),height = 6,width=8)
par (mfrow = c(1,2),mar = c(7,4,7,1))
## MSS 
m1 <- glm (EST.rich ~ poly(Lat,2),
           data = cov_benthos_MMS,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (Lat = seq(-27.85,1,0.5))
pred.vals <- predict (m1,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 
plot(NA,
     ylim=c(0,50),
     xlim=c(-28,1),
     ylab = "Species richness",
     xlab = "Latitude",
     main = "Minimum Sample Size (MSS)")
lines((pred.vals$fit) ~ newd$Lat, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "green",
      lwd =3)

#
lines (upr ~  newd$Lat,lwd=2,col="gray50")
lines (lwr ~  newd$Lat,lwd=2,col="gray50")
points (cov_benthos_MMS$EST.rich ~ cov_benthos_MMS$Lat,
        col=rgb(0,0.2,0,alpha=0.4),pch=19)

m2 <- glm (EST.rich ~ poly(Lat,2),
           data = cov_fish_MMS,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (Lat = seq(-27.85,1,0.5))
pred.vals <- predict (m2,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

lines((pred.vals$fit) ~ newd$Lat, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "red",
      type="l",
      lwd =3)

#
lines (upr ~  newd$Lat,lwd=2,col="gray50")
lines (lwr ~  newd$Lat,lwd=2,col="gray50")
points (cov_fish_MMS$EST.rich ~ cov_fish_MMS$Lat,
        col=rgb(0.2,0,0,alpha=0.4),pch=19)

# ASSi
bind.rich.data.ASSi <- rbind (data.frame (Lat= cov_fish_ASSi$Lat,
                                     EST.rich = cov_fish_ASSi$EST.rich,
                                     Organism="Fishes"),
                         data.frame (Lat= cov_benthos_ASSi$Lat,
                                     EST.rich = cov_benthos_ASSi$EST.rich,
                                     Organism="Benthos"))

model.SR.ASSi <-glm (EST.rich ~ Lat*Organism,
                       bind.rich.data.ASSi,
                      family = gaussian (link="identity"))
summary(model.SR.ASSi)
anova(model.SR.ASSi)

m.lst.SR.ASSi <- emtrends(model.SR.ASSi, "Organism", var="Lat")
m.lst.SR.ASSi

## model
m1 <- glm (EST.rich ~ poly(Lat,2),
           data = cov_benthos_ASSi,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (Lat = seq(-27.85,1,0.5))
pred.vals <- predict (m1,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 
plot(NA,
     ylim=c(0,90),
     xlim=c(-28,1),
     ylab = "",
     xlab = "Latitude",
     main = "Asymptotic Sample Size (ASSi)")
lines((pred.vals$fit) ~ newd$Lat, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "darkgreen",
      lwd =3)

#
lines (upr ~  newd$Lat,lwd=2,col="gray50")
lines (lwr ~  newd$Lat,lwd=2,col="gray50")
points (cov_benthos_ASSi$EST.rich ~ cov_benthos_ASSi$Lat,
        col=rgb(0,0.2,0,alpha=0.4),pch=19)

m2 <- glm (EST.rich ~ poly(Lat,2),
           data = cov_fish_ASSi,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (Lat = seq(-27.85,1,0.5))
pred.vals <- predict (m2,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

lines((pred.vals$fit) ~ newd$Lat, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "darkred",
      type="l",
      lwd =3)

#
lines (upr ~  newd$Lat,lwd=2,col="gray50")
lines (lwr ~  newd$Lat,lwd=2,col="gray50")
points (cov_fish_ASSi$EST.rich ~ cov_fish_ASSi$Lat,
        col=rgb(0.2,0,0,alpha=0.4),pch=19)


legend ("topright", legend = c("Benthos", "Fishes"),
        lwd = 3, col =c("darkgreen","darkred"),
        bty="n")

dev.off()

############################################
### FD

# bind FD data from MSS and ASSi
# fish
bind.FD.data <- rbind (data.frame (Lat= cov_fish_MMS$Lat,
                                     FD = cov_fish_MMS$FD,
                                     Organism="Fishes"),
                         data.frame (Lat= cov_benthos_MMS$Lat,
                                     FD = cov_benthos_MMS$FD,
                                     Organism="Benthos"))

## test of slopes

model.FD.MMS <-glm (FD ~ poly(Lat,2)*Organism,
                    bind.FD.data,
                      family = gaussian (link="identity"))
summary(model.FD.MMS)
anova(model.FD.MMS)

m.lst.FD.MMS <- emtrends(model.FD.MMS, "Organism", var="Lat")
m.lst.FD.MMS

## plottign results
pdf(here("output","vectorized","rel_lat_FD_fish_benthos.pdf"),height = 6,width=8)
par (mfrow = c(1,2),mar = c(7,4,7,1))
## MSS 
m1 <- glm (FD ~ poly(Lat,2),
           data = cov_benthos_MMS,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (Lat = seq(-27.85,1,0.5))
pred.vals <- predict (m1,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 
plot(NA,
     ylim=c(0,0.7),
     xlim=c(-28,1),
     ylab = "Functional diversity",
     xlab = "Latitude",
     main = "Minimum Sample Size (MSS)")
lines((pred.vals$fit) ~ newd$Lat, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "green",
      lwd =3)

#
lines (upr ~  newd$Lat,lwd=2,col="gray50")
lines (lwr ~  newd$Lat,lwd=2,col="gray50")
points (cov_benthos_MMS$FD ~ cov_benthos_MMS$Lat,
        col=rgb(0,0.2,0,alpha=0.4),pch=19)

m2 <- glm (FD ~ poly(Lat,2),
           data = cov_fish_MMS,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (Lat = seq(-27.85,1,0.5))
pred.vals <- predict (m2,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

lines((pred.vals$fit) ~ newd$Lat, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "red",
      type="l",
      lwd =3)

#
lines (upr ~  newd$Lat,lwd=2,col="gray50")
lines (lwr ~  newd$Lat,lwd=2,col="gray50")
points (cov_fish_MMS$FD ~ cov_fish_MMS$Lat,
        col=rgb(0.2,0,0,alpha=0.4),pch=19)

# ASSi
bind.FD.data.ASSi <- rbind (data.frame (Lat= cov_fish_ASSi$Lat,
                                          FD = cov_fish_ASSi$FD,
                                          Organism="Fishes"),
                              data.frame (Lat= cov_benthos_ASSi$Lat,
                                          FD = cov_benthos_ASSi$FD,
                                          Organism="Benthos"))

model.FD.ASSi <-glm (FD ~ poly(Lat,2)*Organism,
                       bind.FD.data.ASSi,
                       family = gaussian (link="identity"))
summary(model.FD.ASSi)
anova(model.FD.ASSi)

m.lst.FD.ASSi <- emtrends(model.FD.ASSi, "Organism", var="Lat")
m.lst.FD.ASSi

## model
m1 <- glm (FD ~ poly(Lat,2),
           data = cov_benthos_ASSi,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (Lat = seq(-27.85,1,0.5))
pred.vals <- predict (m1,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 
plot(NA,
     ylim=c(0,0.7),
     xlim=c(-28,1),
     ylab = "",
     xlab = "Latitude",
     main = "Asymptotic Sample Size (ASSi)")
lines((pred.vals$fit) ~ newd$Lat, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "darkgreen",
      lwd =3)

#
lines (upr ~  newd$Lat,lwd=2,col="gray50")
lines (lwr ~  newd$Lat,lwd=2,col="gray50")
points (cov_benthos_ASSi$FD ~ cov_benthos_ASSi$Lat,
        col=rgb(0,0.2,0,alpha=0.4),pch=19)

m2 <- glm (FD ~ poly(Lat,2),
           data = cov_fish_ASSi,
           family = gaussian (link="identity"))

# data to predict (based on the model)
newd <- data.frame (Lat = seq(-27.85,1,0.5))
pred.vals <- predict (m2,
                      newdata = newd,
                      type="response", # predictions in the link-function scale
                      se.fit = T)

# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

lines((pred.vals$fit) ~ newd$Lat, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
      col = "darkred",
      type="l",
      lwd =3)

#
lines (upr ~  newd$Lat,lwd=2,col="gray50")
lines (lwr ~  newd$Lat,lwd=2,col="gray50")
points (cov_fish_ASSi$FD ~ cov_fish_ASSi$Lat,
        col=rgb(0.2,0,0,alpha=0.4),pch=19)


legend ("topright", legend = c("Benthos", "Fishes"),
        lwd = 3, col =c("darkgreen","darkred"),
        bty="n")

dev.off()







#############################

par(mfrow=c(2,2))

newd <- data.frame (BO2_salinitymean_ss_std= 0,
                    inv_temp_std=0,
                    BO2_ppmean_ss_std=seq (range(cov_fish_MMS$BO_damean_std)[1],
                                           range(cov_fish_MMS$BO_damean_std)[2],
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

plot(exp(pred.vals$fit) ~ newd$BO2_ppmean_ss_std, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,90),
     ylab = "Estimated fish richness",
     xlab = "Productivity",
     xaxt='n')
axis(1,
     at=cov_fish_ASSi$BO2_ppmean_ss_std,
     round(cov_fish_ASSi$BO2_ppmean_ss,3))

lines (exp(upr) ~  newd$BO2_ppmean_ss_std,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$BO2_ppmean_ss_std,lwd=2,col="gray50")
points (cov_fish_ASSi$EST.rich ~ cov_fish_ASSi$BO2_ppmean_ss_std,
        col=rgb(0,0,0.2,alpha=0.1),pch=19)

#

## newd
newd <- data.frame (BO2_salinitymean_ss_std= 0,
                    inv_temp_std=seq (range(cov_fish_MMS$inv_temp_std)[1],
                                      range(cov_fish_MMS$inv_temp_std)[2],
                                      0.05),
                    BO2_ppmean_ss_std=0)

## data to predict (based on the model)
pred.vals <- predict (model3,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$inv_temp_std, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,150),
     ylab = "",
     xlab = "Inverse of temperature",
     xaxt='n')
axis(1,
     at=cov_fish_ASSi$inv_temp_std,
     round(cov_fish_ASSi$inv_temp,2))

lines (exp(upr) ~  newd$inv_temp_std,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$inv_temp_std,lwd=2,col="gray50")
points (cov_fish_ASSi$EST.rich ~ cov_fish_ASSi$inv_temp_std,
        col=rgb(0,0,0.2,alpha=0.1),pch=19)


## newd
newd <- data.frame (BO2_salinitymean_ss_std= seq (range(cov_fish_MMS$BO2_salinitymean_ss_std)[1],
                                                  range(cov_fish_MMS$BO2_salinitymean_ss_std)[2],
                                                  0.05),
                    inv_temp_std=0,
                    BO2_ppmean_ss_std=0)

## data to predict (based on the model)
pred.vals <- predict (model3,
                      newdata = newd,
                      type="link", # predictions in the link-function scale
                      se.fit = T)
# get confidence interval (in the link-function scale)
upr <- pred.vals$fit + (1.96 * pred.vals$se.fit)
lwr <- pred.vals$fit - (1.96 * pred.vals$se.fit)

# plotting 

plot(exp(pred.vals$fit) ~ newd$BO2_salinitymean_ss_std, ## exp is because the inverse of log-link in its exponential - see basics of GLM Poisson 
     col = "darkred",
     type="l",
     lwd =3,
     ylim=c(10,70),
     ylab = "",
     xlab = "Salinity",
     xaxt='n')
axis(1,
     at=cov_fish_ASSi$BO2_salinitymean_ss_std,
     round(cov_fish_ASSi$BO2_salinitymean_ss,2))

lines (exp(upr) ~  newd$BO2_salinitymean_ss_std,lwd=2,col="gray50")
lines (exp(lwr) ~  newd$BO2_salinitymean_ss_std,lwd=2,col="gray50")
points (cov_fish_ASSi$EST.rich ~ cov_fish_ASSi$BO2_salinitymean_ss_std,
        col=rgb(0,0,0.2,alpha=0.1),pch=19)

