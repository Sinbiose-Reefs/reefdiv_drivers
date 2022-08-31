
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

load (here ("output_no_resampling","FD_obs_data_fish.RData"))
load (here ("output_no_resampling","FD_obs_data_benthos.RData"))

# ------------------------------------------ #
#
#     Load environment data

load (here ("output","env_data.RData"))

############# -
# checking quality of Functional space
FD_obs_data_fish$chosenAxes # fish
FD_obs_data_benthos$chosenAxes # benthos

#################################################################
## PROCESSING DATA FOR MODELING

# ------------------------------- #
# results to organize (average, connect with covariates)

## bind FD, estimated richness and  site names 
# bind
FRic_fish<-  data.frame ("FRic"= FD_obs_data_fish$Fdindexes$FRic,
                         "FEve"= FD_obs_data_fish$Fdindexes$FEve,
                         "FDiv"= FD_obs_data_fish$Fdindexes$FDiv,
                         "SR" = rowSums(ifelse (comp_peixes>0,1,0)),
                         "Ntrans" = res_table_samples$n_total,
                         "SamplingArea" = res_table_samples$n_total * 40) # 40m2 per transect

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
cov_fish$SamplingArea_std <- log(cov_fish$SamplingArea)
cov_fish$SamplingArea_std <- (cov_fish$SamplingArea_std-mean(cov_fish$SamplingArea_std))/sd(cov_fish$SamplingArea_std)

# depth
cov_fish <- cbind(cov_fish,
                  Depth=covariates_site$prof)
cov_fish$Depth <- as.factor (cov_fish$Depth)
  
# ----------------------------- #
#            Benthos
# ----------------------------- #

# -------------------------------------------------------------------------- #
# Benthos

# bind
FRic_benthos<-  data.frame ("FRic"= FD_obs_data_benthos$Fdindexes$FRic,
                            "FEve"= FD_obs_data_benthos$Fdindexes$FEve,
                            "FDiv"= FD_obs_data_benthos$Fdindexes$FDiv,
                            "SR" = rowSums(ifelse (comp_bentos>0,1,0)),
                            "Ntrans" = res_table_samples_bentos$n_total,
                            "SamplingArea" = res_table_samples_bentos$n_total * 2) # 40m2 per transect
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
cov_benthos$SamplingArea_std <- log(cov_benthos$SamplingArea)
cov_benthos$SamplingArea_std <- (cov_benthos$SamplingArea_std-mean(cov_benthos$SamplingArea_std))/sd(cov_benthos$SamplingArea_std)

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
      file=here("output_no_resampling", "data_to_modeling_GLM.RData"))

rm(list=ls())
