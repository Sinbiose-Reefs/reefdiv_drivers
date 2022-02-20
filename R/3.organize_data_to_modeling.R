
# ----------------------------------------------------------------------------#
#    routine to modeling fish and benthos SR and FD  relative to environment
#                          using GLM
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

load (here ("output","FD_fish_MSS.RData"))
load (here ("output","FD_benthos_MSS.RData"))

# ------------------------------------------ #
#
#     Load environment data

load (here ("output","env_data.RData"))

############# -
# checking quality of Functional space

# fishes

naxes_fishes <- data.frame (meanAxesMSS = mean(unlist(sapply (FD_fish_MSS, "[","chosenAxes")),na.rm=T),
                            # sd
                            sdAxesMSS = sd(unlist(sapply (FD_fish_MSS, "[","chosenAxes")),na.rm=T),
                            # quality
                            qualityMSS= mean (unlist(lapply (sapply (sapply (FD_fish_MSS, "[","Fdindexes"), "[", "qual.FRic"),unique)),na.rm=T),
                            # sd quality
                            sdqualityMSS= sd (unlist(lapply (sapply (sapply (FD_fish_MSS, "[","Fdindexes"), "[", "qual.FRic"),unique)),na.rm=T)
                            )                            


# benthos

naxes_benthos <- data.frame (meanAxesMSS = mean(unlist(sapply (FD_benthos_MSS, "[","chosenAxes")),na.rm=T),
                             # sd
                             sdAxesMSS = sd(unlist(sapply (FD_benthos_MSS, "[","chosenAxes")),na.rm=T),
                             # quality
                             qualityMSS= mean (unlist(lapply (sapply (sapply (FD_benthos_MSS, "[","Fdindexes"), "[", "qual.FRic"),unique)),na.rm=T),
                             # sd quality
                             sdqualityMSS= sd (unlist(lapply (sapply (sapply (FD_benthos_MSS, "[","Fdindexes"), "[", "qual.FRic"),unique)),na.rm=T)
                             )                            


save (naxes_fishes,
      naxes_benthos, file = here("output","FE_quality.RData") )


#################################################################
## PROCESSING DATA FOR MODELING

# ------------------------------- #
# results to organize (average, connect with covariates)

## list of indexes
FDindexes <- c("FRic", "FEve", "FDiv")

# apply org
FRic_fish <- lapply (FDindexes, function (index) {
  
  pre_res <- sapply (
                sapply (FD_fish_MSS, "[", "Fdindexes"),"[", index)
  
  # rm estimates that did not work
  pre_res <- do.call(cbind,
                     pre_res [(which(unlist(lapply (lapply (pre_res,dim),is.null)) != TRUE))])
  # get the average
  av_FRic <- apply (pre_res,1,mean,na.rm=T)
  ; # return
  av_FRic
  
  }
)

## bind estimated richness and  site names 
# bind
FRic_fish<-  data.frame (do.call(cbind, FRic_fish), 
              res_table_samples[,"EST.Rich"],
              res_table_samples[,"Site"])

# adjust colnames
colnames(FRic_fish) <- c("FRic","FEve","FDiv","EstRich","Sites")

# -------------------------------------------------------------------------- #
# Benthos

# apply org
FRic_benthos <- lapply (FDindexes, function (index) {
  
  pre_res <- sapply (
                  sapply (FD_benthos_MSS, "[", "Fdindexes"),"[", index)
  
  # rm estimates that did not work
  pre_res <- do.call(cbind,
                     pre_res [(which(unlist(lapply (lapply (pre_res,dim),is.null)) != TRUE))])
  # get the average
  av_FRic <- apply (pre_res,1,mean,na.rm=T)
  ; # return
  av_FRic
  
}
)


## bind estimated richness and  site names 

FRic_benthos <-  data.frame (do.call(cbind,FRic_benthos), 
                            res_table_samples_bentos[,"EST.Rich"],
                            res_table_samples_bentos[,"Site"])

# adjust colnames
colnames(FRic_benthos) <- c("FRic","FEve","FDiv","EstRich","Sites")


# ----------------------------------------------- #
# organize data to modeling

# organizing data
# fish
cov_fish <- cbind(FRic_fish,
                    Organism="Fishes")
  
# environmental data
cov_fish <- cbind(cov_fish,
                       covariates_site$sea_data[[1]])

# spatial data
cov_fish <- cbind(cov_fish,
                       covariates_site$coord$coord_peixes[,c("Group.1","Lon","Lat","distance")])

# ------------------------
# organize and standardize covariates
# inverse of temperature, standardize covariates, insert region and locality
# ------------------------

## standardize other variables
cov_fish$BO2_tempmean_ss_std <- (cov_fish$BO2_tempmean_ss-mean(cov_fish$BO2_tempmean_ss))/sd(cov_fish$BO2_tempmean_ss)
cov_fish$BO2_ppmean_ss_std <- (cov_fish$BO2_ppmean_ss-mean(cov_fish$BO2_ppmean_ss))/sd(cov_fish$BO2_ppmean_ss)
cov_fish$BO2_salinitymean_ss_std <- (cov_fish$BO2_salinitymean_ss-mean(cov_fish$BO2_salinitymean_ss))/sd(cov_fish$BO2_salinitymean_ss)
cov_fish$BO_damean_std <- (cov_fish$BO_damean-mean(cov_fish$BO_damean))/sd(cov_fish$BO_damean)
cov_fish$distanceLog <- log(cov_fish$distance)
cov_fish$distance_std <- (cov_fish$distanceLog-mean(cov_fish$distanceLog))/sd(cov_fish$distanceLog)
cov_fish$area <- (covariates_site$reef_area-mean(covariates_site$reef_area))/sd(covariates_site$reef_area)

# depth
cov_fish <- cbind(cov_fish,
                  Depth=covariates_site$prof)
cov_fish$Depth <- as.factor (cov_fish$Depth)
  
# ----------------------------- #
#            Benthos
# ----------------------------- #

# organizing data
# benthos
cov_benthos <- cbind(FRic_benthos,
                  Organism="Benthos")

# environmental data
cov_benthos <- cbind(cov_benthos,
                  covariates_site$sea_data[[1]])

# spatial data
cov_benthos <- cbind(cov_benthos,
                  covariates_site$coord$coord_peixes[,c("Group.1","Lon","Lat","distance")])

# ------------------------
# organize and standardize covariates
# inverse of temperature, standardize covariates, insert region and locality
# ------------------------

## standardize other variables
cov_benthos$BO2_tempmean_ss_std <- (cov_benthos$BO2_tempmean_ss-mean(cov_benthos$BO2_tempmean_ss))/sd(cov_benthos$BO2_tempmean_ss)
cov_benthos$BO2_ppmean_ss_std <- (cov_benthos$BO2_ppmean_ss-mean(cov_benthos$BO2_ppmean_ss))/sd(cov_benthos$BO2_ppmean_ss)
cov_benthos$BO2_salinitymean_ss_std <- (cov_benthos$BO2_salinitymean_ss-mean(cov_benthos$BO2_salinitymean_ss))/sd(cov_benthos$BO2_salinitymean_ss)
cov_benthos$BO_damean_std <- (cov_benthos$BO_damean-mean(cov_benthos$BO_damean))/sd(cov_benthos$BO_damean)
cov_benthos$distanceLog <- log(cov_benthos$distance)
cov_benthos$distance_std <- (cov_benthos$distanceLog-mean(cov_benthos$distanceLog))/sd(cov_benthos$distanceLog)
cov_benthos$area <- (covariates_site$reef_area-mean(covariates_site$reef_area))/sd(covariates_site$reef_area)

# depth
cov_benthos <- cbind(cov_benthos,
                  Depth=covariates_site$prof)
cov_benthos$Depth <- as.factor (cov_benthos$Depth)

## aggregate list of results across organisms by rarefaction method

data_to_modeling_GLM <- rbind(cov_fish, 
                              cov_benthos)
  
# save
save (data_to_modeling_GLM, 
      cov_fish,
      cov_benthos,
      file=here("output", "data_to_modeling_GLM.RData"))

rm(list=ls())
