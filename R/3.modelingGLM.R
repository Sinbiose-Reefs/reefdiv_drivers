
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

load (here ("output","FD_fish_ASSi.RData"))
load (here ("output","FD_fish_MSS.RData"))
load (here ("output","FD_fish_lomolino.RData"))

load (here ("output","FD_benthos_ASSi.RData"))
load (here ("output","FD_benthos_MSS.RData"))
load (here ("output","FD_benthos_lomolino.RData"))

# ------------------------------------------ #
#
#     Load environment data

load (here ("output","env_data.RData"))


#################################################################
## PROCESSING DATA FOR MODELING

# ------------------------------- #
# minimum sample size

## FRIC

### fish 
FRic_fish_MSS <- sapply (
   sapply (FD_fish_MSS, "[", "Fdindexes"),"[", "FRic")
# rm estimates that did not work
FRic_fish_MSS <- do.call(cbind,
                         FRic_fish_MSS [(which(unlist(lapply (lapply (FRic_fish_MSS,dim),is.null)) != TRUE))])

# get the average
av_FRic_fish_MSS <- apply (FRic_fish_MSS,1,mean,na.rm=T)

### benthos
FRic_benthos_MSS <- sapply (
   sapply (FD_benthos_MSS, "[", "Fdindexes"),"[", "FRic")
# rm estimates that did not work
FRic_benthos_MSS <- do.call(cbind,
                            FRic_benthos_MSS [(which(unlist(lapply (lapply (FRic_benthos_MSS,length),is.null)) != TRUE))])

# get the average
av_FRic_benthos_MSS <- apply (FRic_benthos_MSS,1,mean,na.rm=T)

# Feve

### fish 
FEve_fish_MSS <- sapply (
   sapply (FD_fish_MSS, "[", "Fdindexes"),"[", "FEve")
# rm estimates that did not work
FEve_fish_MSS <- do.call(cbind,
                         FEve_fish_MSS [(which(unlist(lapply (lapply (FEve_fish_MSS,dim),is.null)) != TRUE))])

# get the average
av_FEve_fish_MSS <- apply (FEve_fish_MSS,1,mean,na.rm=T)

### benthos
FEve_benthos_MSS <- sapply (
   sapply (FD_benthos_MSS, "[", "Fdindexes"),"[", "FEve")
# rm estimates that did not work
FEve_benthos_MSS <- do.call(cbind,
                            FEve_benthos_MSS [(which(unlist(lapply (lapply (FEve_benthos_MSS,length),is.null)) != TRUE))])

# get the average
av_FEve_benthos_MSS <- apply (FEve_benthos_MSS,1,mean,na.rm=T)

## FDIV

### fish 
FDiv_fish_MSS <- sapply (
   sapply (FD_fish_MSS, "[", "Fdindexes"),"[", "FDiv")
# rm estimates that did not work
FDiv_fish_MSS <- do.call(cbind,
                         FDiv_fish_MSS [(which(unlist(lapply (lapply (FDiv_fish_MSS,dim),is.null)) != TRUE))])

# get the average
av_FDiv_fish_MSS <- apply (FDiv_fish_MSS,1,mean,na.rm=T)

### benthos
FDiv_benthos_MSS <- sapply (
   sapply (FD_benthos_MSS, "[", "Fdindexes"),"[", "FDiv")
# rm estimates that did not work
FDiv_benthos_MSS <- do.call(cbind,
                            FDiv_benthos_MSS [(which(unlist(lapply (lapply (FDiv_benthos_MSS,length),is.null)) != TRUE))])

# get the average
av_FDiv_benthos_MSS <- apply (FDiv_benthos_MSS,1,mean,na.rm=T)

# ------------------------------- #
# Precision Based sample size

# FRic
### fish 
FRic_fish_PSS <- sapply (
   sapply (FD_fish_ASSi, "[", "Fdindexes"),"[", "FRic")
# rm estimates that did not work
FRic_fish_PSS <- do.call(cbind,
                         FRic_fish_PSS [(which(unlist(lapply (lapply (FRic_fish_PSS,dim),is.null)) != TRUE))])

# get the average
av_FRic_fish_PSS <- apply (FRic_fish_PSS,1,mean,na.rm=T)

### benthos
FRic_benthos_PSS <- sapply (
   sapply (FD_benthos_ASSi, "[", "Fdindexes"),"[", "FRic")
# rm estimates that did not work
FRic_benthos_PSS <- do.call(cbind,
                            FRic_benthos_PSS [(which(unlist(lapply (lapply (FRic_benthos_PSS,length),is.null)) != TRUE))])

# get the average
av_FRic_benthos_PSS <- apply (FRic_benthos_PSS,1,mean,na.rm=T)

# FEve

### fish 
FEve_fish_PSS <- sapply (
   sapply (FD_fish_ASSi, "[", "Fdindexes"),"[", "FEve")
# rm estimates that did not work
FEve_fish_PSS <- do.call(cbind,
                         FEve_fish_PSS [(which(unlist(lapply (lapply (FEve_fish_PSS,dim),is.null)) != TRUE))])

# get the average
av_FEve_fish_PSS <- apply (FEve_fish_PSS,1,mean,na.rm=T)

### benthos
FEve_benthos_PSS <- sapply (
   sapply (FD_benthos_ASSi, "[", "Fdindexes"),"[", "FEve")
# rm estimates that did not work
FEve_benthos_PSS <- do.call(cbind,
                            FEve_benthos_PSS [(which(unlist(lapply (lapply (FEve_benthos_PSS,length),is.null)) != TRUE))])

# get the average
av_FEve_benthos_PSS <- apply (FEve_benthos_PSS,1,mean,na.rm=T)

# FDiv

### fish 
FDiv_fish_PSS <- sapply (
   sapply (FD_fish_ASSi, "[", "Fdindexes"),"[", "FDiv")
# rm estimates that did not work
FDiv_fish_PSS <- do.call(cbind,
                         FDiv_fish_PSS [(which(unlist(lapply (lapply (FDiv_fish_PSS,dim),is.null)) != TRUE))])

# get the average
av_FDiv_fish_PSS <- apply (FDiv_fish_PSS,1,mean,na.rm=T)

### benthos
FDiv_benthos_PSS <- sapply (
   sapply (FD_benthos_ASSi, "[", "Fdindexes"),"[", "FDiv")
# rm estimates that did not work
FDiv_benthos_PSS <- do.call(cbind,
                            FDiv_benthos_PSS [(which(unlist(lapply (lapply (FDiv_benthos_PSS,length),is.null)) != TRUE))])

# get the average
av_FDiv_benthos_PSS <- apply (FDiv_benthos_PSS,1,mean,na.rm=T)

# ------------------------------- #
# LOMOLINO'S FUNCTION

# FRic
### fish 
FRic_fish_lomolino <- sapply (
   sapply (FD_fish_lomolino, "[", "Fdindexes"),"[", "FRic")
# rm estimates that did not work
FRic_fish_lomolino <- do.call(cbind,
                              FRic_fish_lomolino [(which(unlist(lapply (lapply (FRic_fish_lomolino,length),is.null)) != TRUE))])

# get the average
av_FRic_fish_lomolino <- apply (FRic_fish_lomolino,1,mean,na.rm=T)

### benthos
FRic_benthos_lomolino <- sapply (
   sapply (FD_benthos_lomolino, "[", "Fdindexes"),"[", "FRic")
# rm estimates that did not work
FRic_benthos_lomolino <- do.call(cbind,
                                 FRic_benthos_lomolino [(which(unlist(lapply (lapply (FRic_benthos_lomolino,dim),is.null)) != TRUE))])

# get the average
av_FRic_benthos_lomolino <- apply (FRic_benthos_lomolino,1,mean,na.rm=T)

## FEve

### fish 
FEve_fish_lomolino <- sapply (
   sapply (FD_fish_lomolino, "[", "Fdindexes"),"[", "FEve")
# rm estimates that did not work
FEve_fish_lomolino <- do.call(cbind,
                              FEve_fish_lomolino [(which(unlist(lapply (lapply (FEve_fish_lomolino,length),is.null)) != TRUE))])

# get the average
av_FEve_fish_lomolino <- apply (FEve_fish_lomolino,1,mean,na.rm=T)

### benthos
FEve_benthos_lomolino <- sapply (
   sapply (FD_benthos_lomolino, "[", "Fdindexes"),"[", "FEve")
# rm estimates that did not work
FEve_benthos_lomolino <- do.call(cbind,
                                 FEve_benthos_lomolino [(which(unlist(lapply (lapply (FEve_benthos_lomolino,dim),is.null)) != TRUE))])

# get the average
av_FEve_benthos_lomolino <- apply (FEve_benthos_lomolino,1,mean,na.rm=T)

## FDIV

### fish 
FDiv_fish_lomolino <- sapply (
   sapply (FD_fish_lomolino, "[", "Fdindexes"),"[", "FDiv")
# rm estimates that did not work
FDiv_fish_lomolino <- do.call(cbind,
                              FDiv_fish_lomolino [(which(unlist(lapply (lapply (FDiv_fish_lomolino,length),is.null)) != TRUE))])

# get the average
av_FDiv_fish_lomolino <- apply (FDiv_fish_lomolino,1,mean,na.rm=T)

### benthos
FDiv_benthos_lomolino <- sapply (
   sapply (FD_benthos_lomolino, "[", "Fdindexes"),"[", "FDiv")
# rm estimates that did not work
FDiv_benthos_lomolino <- do.call(cbind,
                                 FDiv_benthos_lomolino [(which(unlist(lapply (lapply (FDiv_benthos_lomolino,dim),is.null)) != TRUE))])

# get the average
av_FDiv_benthos_lomolino <- apply (FDiv_benthos_lomolino,1,mean,na.rm=T)


# ------------------------------------------------------- #
## bind FISH data to test correlation

# MSS

res_table_samples <- cbind (res_table_samples,
       FRic=av_FRic_fish_MSS,
       FEve=av_FEve_fish_MSS,
       FDiv=av_FDiv_fish_MSS)
       
# PSSi

res_sp_accum_fish_asymptote <- cbind (res_sp_accum_fish_asymptote,
       FRic=ifelse (is.na (res_sp_accum_fish_asymptote$EST.Rich) == F,
                  av_FRic_fish_PSS, NA),
       FEve=ifelse (is.na (res_sp_accum_fish_asymptote$EST.Rich) == F,
                    av_FEve_fish_PSS, NA),
       FDiv=ifelse (is.na (res_sp_accum_fish_asymptote$EST.Rich) == F,
                    av_FDiv_fish_PSS, NA))

# lomolino

res_lomolino <- cbind (res_lomolino, 
         FRic = ifelse ( res_lomolino[,1] < res_lomolino[,2], 
                  av_FRic_fish_lomolino, NA),
       
       FEve = ifelse ( res_lomolino[,1] < res_lomolino[,2], 
                       av_FEve_fish_lomolino, NA),
       FDiv = ifelse ( res_lomolino[,1] < res_lomolino[,2], 
                       av_FDiv_fish_lomolino, NA))
# ------------------------------ #       
# correlation between indexes - FISH

## MSS 

cormat_MSS_fish <-  (cbind (Richness=res_table_samples$EST.Rich,
                            FRic=res_table_samples$FRic,
                            FEve=res_table_samples$FEve,
                            FDiv=res_table_samples$FDiv)
)

res1 <- cor.mtest(cormat_MSS_fish, conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_MSS_fish.pdf"),width=6,heigh=5)
corrplot(cor(cormat_MSS_fish), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

## PSSi
cormat_PSS_fish <-   cbind (Richness= res_sp_accum_fish_asymptote$EST.Rich,
                            FRic = res_sp_accum_fish_asymptote$FRic,
                            FEve = res_sp_accum_fish_asymptote$FEve,
                            FDiv = res_sp_accum_fish_asymptote$FDiv)[which(is.na(res_sp_accum_fish_asymptote$EST.Rich)!= T),]

res1 <- cor.mtest(cormat_PSS_fish, conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_PSSi_fish.pdf"),width=6,heigh=5)
corrplot(cor(cormat_PSS_fish), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()


## lomolino's S50
cormat_S50_fish <- cbind (Richness = res_lomolino[,"asymp"],
                          FRic = res_lomolino [,"FRic"],
                          FEve = res_lomolino [,"FEve"],
                          FDiv = res_lomolino [,"FDiv"])[which(is.na(res_lomolino[,"FRic"])!= T),]

res1 <- cor.mtest(cormat_S50_fish, conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_S50_fish.pdf"),width=6,heigh=5)
corrplot(cor(cormat_S50_fish), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

# ------------------------------ #       
# correlation between algorithms - FISH

## richness 
estimates_alg <- cbind (MSS =res_table_samples$EST.Rich,
                        PSSi=res_sp_accum_fish_asymptote$EST.Rich,
                        S50=res_lomolino[,"asymp"])

res1 <- cor.mtest(estimates_alg[which(is.na(estimates_alg[,2] )!= T),], conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_rich_algorithm_fish.pdf"),width=6,heigh=5)
corrplot(cor(estimates_alg[which(is.na(estimates_alg[,2] )!= T),]), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

# FRic

estimates_alg_FRic <- cbind (MSS=res_table_samples$FRic,
                             PSSi=res_sp_accum_fish_asymptote$FRic,
                             S50=res_lomolino[,"FRic"])

res1 <- cor.mtest(estimates_alg_FRic[which(is.na(estimates_alg_FRic[,2]) !=T & 
                                              (is.na(estimates_alg_FRic[,3]) !=T)),], 
                  conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_FRic_algorithm_fish.pdf"),width=6,heigh=5)
corrplot(cor(estimates_alg_FRic[which(is.na(estimates_alg_FRic[,2]) !=T & 
                                         (is.na(estimates_alg_FRic[,3]) !=T)),]), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

# FEve

estimates_alg_FEve <- cbind (MSS=res_table_samples$FEve,
                             PSSi=res_sp_accum_fish_asymptote$FEve,
                             S50=res_lomolino[,"FEve"])

res1 <- cor.mtest(estimates_alg_FEve[which(is.na(estimates_alg_FEve[,2]) !=T & 
                                              (is.na(estimates_alg_FEve[,3]) !=T)),], 
                  conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_FEve_algorithm_fish.pdf"),width=6,heigh=5)
corrplot(cor(estimates_alg_FEve[which(is.na(estimates_alg_FEve[,2]) !=T & 
                                         (is.na(estimates_alg_FEve[,3]) !=T)),]), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()


# FDiv

estimates_alg_FDiv <- cbind (MSS=res_table_samples$FDiv,
                             PSSi=res_sp_accum_fish_asymptote$FDiv,
                             S50=res_lomolino[,"FDiv"])

res1 <- cor.mtest(estimates_alg_FDiv[which(is.na(estimates_alg_FDiv[,2]) !=T & 
                                              (is.na(estimates_alg_FDiv[,3]) !=T)),], 
                  conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_FDiv_algorithm_fish.pdf"),width=6,heigh=5)
corrplot(cor(estimates_alg_FDiv[which(is.na(estimates_alg_FDiv[,2]) !=T & 
                                         (is.na(estimates_alg_FDiv[,3]) !=T)),]), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()



## average values of indexes

rbind (apply(estimates_alg,2,mean,na.rm=T),
       apply(estimates_alg_FRic,2,mean,na.rm=T),
       apply(estimates_alg_FEve,2,mean,na.rm=T),
       apply(estimates_alg_FDiv,2,mean,na.rm=T)
)

## sd
rbind (apply(estimates_alg,2,sd,na.rm=T),
       apply(estimates_alg_FRic,2,sd,na.rm=T),
       apply(estimates_alg_FEve,2,sd,na.rm=T),
       apply(estimates_alg_FDiv,2,sd,na.rm=T)
)

# ---------------- #
# Benthos
# bind benthos data to test correlation

# MSS

res_table_samples_bentos <- cbind (res_table_samples_bentos,
                            FRic=av_FRic_benthos_MSS,
                            FEve=av_FEve_benthos_MSS,
                            FDiv=av_FDiv_benthos_MSS)

# PSSi

res_sp_accum_bentos_asymptote <- cbind (res_sp_accum_bentos_asymptote,
                                      FRic=ifelse (is.na (res_sp_accum_bentos_asymptote$EST.Rich) == F,
                                                   av_FRic_benthos_PSS, NA),
                                      FEve=ifelse (is.na (res_sp_accum_bentos_asymptote$EST.Rich) == F,
                                                   av_FEve_benthos_PSS, NA),
                                      FDiv=ifelse (is.na (res_sp_accum_bentos_asymptote$EST.Rich) == F,
                                                   av_FDiv_benthos_PSS, NA))

# lomolino

res_lomolino_bentos <- cbind (res_lomolino_bentos, 
                       FRic = ifelse ( res_lomolino_bentos[,1] < res_lomolino_bentos[,2], 
                                       av_FRic_benthos_lomolino, NA),
                       FEve = ifelse ( res_lomolino_bentos[,1] < res_lomolino_bentos[,2], 
                                       av_FEve_benthos_lomolino, NA),
                       FDiv = ifelse ( res_lomolino_bentos[,1] < res_lomolino_bentos[,2], 
                                       av_FDiv_benthos_lomolino, NA))

# ------------------------------ #       
# correlation between indexes - benthos

## MSS 

cormat_MSS_benthos <-  (cbind (Richness=res_table_samples_bentos$EST.Rich,
                            FRic=res_table_samples_bentos$FRic,
                            FEve=res_table_samples_bentos$FEve,
                            FDiv=res_table_samples_bentos$FDiv)
)

res1 <- cor.mtest(cormat_MSS_benthos, conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_MSS_benthos.pdf"),width=6,heigh=5)
corrplot(cor(cormat_MSS_benthos), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

## PSSi
cormat_PSS_benthos <-   cbind (Richness= res_sp_accum_bentos_asymptote$EST.Rich,
                            FRic = res_sp_accum_bentos_asymptote$FRic,
                            FEve = res_sp_accum_bentos_asymptote$FEve,
                            FDiv = res_sp_accum_bentos_asymptote$FDiv)[which(is.na(res_sp_accum_bentos_asymptote$EST.Rich)!= T),]

res1 <- cor.mtest(cormat_PSS_benthos, conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_PSSi_benthos.pdf"),width=6,heigh=5)
corrplot(cor(cormat_PSS_benthos), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()


## lomolino's S50
cormat_S50_benthos <- cbind (Richness = res_lomolino_bentos[,"asymp"],
                          FRic = res_lomolino_bentos [,"FRic"],
                          FEve = res_lomolino_bentos [,"FEve"],
                          FDiv = res_lomolino_bentos [,"FDiv"])[which(is.na(res_lomolino_bentos[,"FRic"])!= T),]

res1 <- cor.mtest(cormat_S50_benthos, conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_S50_benthos.pdf"),width=6,heigh=5)
corrplot(cor(cormat_S50_benthos), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

# ------------------------------ #       
# correlation between algorithms - benthos

## richness 
estimates_alg_benthos <- cbind (MSS =res_table_samples_bentos$EST.Rich,
                        PSSi=res_sp_accum_bentos_asymptote$EST.Rich,
                        S50=res_lomolino_bentos[,"asymp"])

res1 <- cor.mtest(estimates_alg_benthos[which(is.na(estimates_alg_benthos[,2] )!= T),], conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_rich_algorithm_benthos.pdf"),width=6,heigh=5)
corrplot(cor(estimates_alg_benthos[which(is.na(estimates_alg_benthos[,2] )!= T),]), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

# FRic

estimates_alg_FRic_benthos <- cbind (MSS=res_table_samples_bentos$FRic,
                             PSSi=res_sp_accum_bentos_asymptote$FRic,
                             S50=res_lomolino_bentos[,"FRic"])

res1 <- cor.mtest(estimates_alg_FRic_benthos[which(is.na(estimates_alg_FRic_benthos[,2]) !=T & 
                                              (is.na(estimates_alg_FRic_benthos[,3]) !=T)),], 
                  conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_FRic_algorithm_benthos.pdf"),width=6,heigh=5)
corrplot(cor(estimates_alg_FRic_benthos[which(is.na(estimates_alg_FRic_benthos[,2]) !=T & 
                                         (is.na(estimates_alg_FRic_benthos[,3]) !=T)),]), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

# FEve

estimates_alg_FEve_benthos <- cbind (MSS=res_table_samples_bentos$FEve,
                             PSSi=res_sp_accum_bentos_asymptote$FEve,
                             S50=res_lomolino_bentos[,"FEve"])

res1 <- cor.mtest(estimates_alg_FEve_benthos[which(is.na(estimates_alg_FEve_benthos[,2]) !=T & 
                                              (is.na(estimates_alg_FEve_benthos[,3]) !=T)),], 
                  conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_FEve_algorithm_benthos.pdf"),width=6,heigh=5)
corrplot(cor(estimates_alg_FEve_benthos[which(is.na(estimates_alg_FEve_benthos[,2]) !=T & 
                                         (is.na(estimates_alg_FEve_benthos[,3]) !=T)),]), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

# FDiv

estimates_alg_FDiv_benthos <- cbind (MSS=res_table_samples_bentos$FDiv,
                             PSSi=res_sp_accum_bentos_asymptote$FDiv,
                             S50=res_lomolino_bentos[,"FDiv"])

res1 <- cor.mtest(estimates_alg_FDiv_benthos[which(is.na(estimates_alg_FDiv_benthos[,2]) !=T & 
                                              (is.na(estimates_alg_FDiv_benthos[,3]) !=T)),], 
                  conf.level = .99)

## specialized the insignificant value according to the significant level

pdf (here ("output", "vectorized", "corr_FDiv_algorithm_benthos.pdf"),width=6,heigh=5)
corrplot(cor(estimates_alg_FDiv_benthos[which(is.na(estimates_alg_FDiv_benthos[,2]) !=T & 
                                         (is.na(estimates_alg_FDiv_benthos[,3]) !=T)),]), 
         p.mat = res1$p, sig.level = .05,
         insig = "pch")  # or blank to none
dev.off()

## average values of indexes

rbind (apply(estimates_alg_benthos,2,mean,na.rm=T),
       apply(estimates_alg_FRic_benthos,2,mean,na.rm=T),
       apply(estimates_alg_FEve_benthos,2,mean,na.rm=T),
       apply(estimates_alg_FDiv_benthos,2,mean,na.rm=T)
)

## sd
rbind (apply(estimates_alg_benthos,2,sd,na.rm=T),
   apply(estimates_alg_FRic_benthos,2,sd,na.rm=T),
   apply(estimates_alg_FEve_benthos,2,sd,na.rm=T),
   apply(estimates_alg_FDiv_benthos,2,sd,na.rm=T)
)

# ----------------------------------------- #
#             covariates
# ------------------------------------------#

colnames(res_lomolino)[1:3] <- c("nsamples","obs_ss", "EST.Rich")
res_lomolino <- data.frame(res_lomolino,Site = sites_fish_complete)
colnames(res_lomolino_bentos) <- c("nsamples","obs_ss", "EST.Rich")
res_lomolino_bentos <- data.frame(res_lomolino_bentos,Site =sites_bentos_complete)

## subsetting lomolino's data

res_lomolino <- res_lomolino[which(res_lomolino[,1] < res_lomolino[,2]),]
res_lomolino_bentos <- res_lomolino_bentos[which(res_lomolino_bentos[,1] < res_lomolino_bentos[,2]),]

# and sites
sites_bentos_complete <- sites_bentos_complete [which(res_lomolino_bentos[,1] < res_lomolino_bentos[,2])]
sites_fish_complete <- sites_fish_complete [which(res_lomolino[,1] < res_lomolino[,2])]

# ---------------------------------------- #
# organizing data of fishes
# list of SR results
list_data_fish_SR <- list (res_table_samples, res_sp_accum_fish_asymptote,
                           res_lomolino)
# object do indicate NAs
list_data_fish_SR<-lapply (list_data_fish_SR[1:2], function (i)
   
   cbind (i,
          maintain = ifelse (is.na(i$EST.Rich)==T, "No","Yes")
   )
)
## object to indicate imprecise sample size and richness estimate of lomolino's model
res_lomolino$maintain <- ifelse (res_lomolino$nsamples <  res_lomolino$obs_ss, "Yes","No")

list_data_fish_SR <- c(list_data_fish_SR,
                       list(res_lomolino))

# list of FD results
list_data_fish_FD <- list (av_FRic_fish_MSS, av_FRic_fish_PSS,
                           av_FRic_fish_lomolino)

# list of FEve results
list_data_fish_FEve <- list (av_FEve_fish_MSS, av_FEve_fish_PSS,
                             av_FEve_fish_lomolino)

# list of FDiv results
list_data_fish_FDiv <- list (av_FDiv_fish_MSS, av_FDiv_fish_PSS,
                             av_FDiv_fish_lomolino)


# method of rarefaction
method <- c("MSS","PSSi","Lomolino")

# organizing data across lists

org_data_fish <- lapply (seq(1,length(list_data_fish_FD)), function (i) {      
   # covariates_site$sea_data[[1]] is for fishes, [[2]] is for benthos
   cov_fish <- covariates_site$sea_data[[1]][which(list_data_fish_SR[[i]]$maintain == "Yes"),]
   # bind richness estimate
   cov_fish <- data.frame(cov_fish, 
                          EST.rich = list_data_fish_SR[[i]]$EST.Rich[which(list_data_fish_SR[[i]]$maintain == "Yes")])
   # bind FD estimate
   cov_fish <- cbind(cov_fish, 
                     FD = list_data_fish_FD[[i]])
   
   # bind Feve
   cov_fish <- cbind(cov_fish, 
                     FEve = list_data_fish_FEve[[i]])
   # bind FDiv
   cov_fish <- cbind(cov_fish, 
                     FDiv = list_data_fish_FDiv[[i]])
   
   # bind coordinates
   cov_fish <- cbind(cov_fish,
                     covariates_site$coord$coord_peixes [match(rownames(cov_fish),covariates_site$coord$coord_peixes$Group.1),c("Lon","Lat","distance")])
   
   # bind organism
   cov_fish <- cbind(cov_fish,
                     Organism="Fishes")
   
   # bind method
   cov_fish <- cbind(cov_fish,
                     Method=method[i])
   
   # ------------------------
   # organize and standardize covariates
   # inverse of temperature, standardize covariates, insert region and locality
   # ------------------------
   
   boltzmann_factor <- 8.62e-5
   #a<-20:25 + 273.15
   #plot(a, (1/boltzmann_factor * (1/mean (a) - 1/(a))))
   
   # FISH
   ## MMS
   temp_kelvin<- cov_fish$BO2_tempmean_ss + 273.15
   inv_temp <- 1 / boltzmann_factor * (1 / mean(temp_kelvin) - 1 / temp_kelvin)
   cov_fish <- cbind(cov_fish,
                     inv_temp=inv_temp)
   
   ## standardize other variables
   cov_fish$BO2_tempmean_ss_std <- (cov_fish$BO2_tempmean_ss-mean(cov_fish$BO2_tempmean_ss))/sd(cov_fish$BO2_tempmean_ss)
   cov_fish$BO2_ppmean_ss_std <- (cov_fish$BO2_ppmean_ss-mean(cov_fish$BO2_ppmean_ss))/sd(cov_fish$BO2_ppmean_ss)
   cov_fish$BO2_salinitymean_ss_std <- (cov_fish$BO2_salinitymean_ss-mean(cov_fish$BO2_salinitymean_ss))/sd(cov_fish$BO2_salinitymean_ss)
   cov_fish$BO_damean_std <- (cov_fish$BO_damean-mean(cov_fish$BO_damean))/sd(cov_fish$BO_damean)
   cov_fish$inv_temp_std <- (cov_fish$inv_temp-mean(cov_fish$inv_temp))/sd(cov_fish$inv_temp)
   cov_fish$distanceLog <- log(cov_fish$distance)
   cov_fish$distance_std <- (cov_fish$distanceLog-mean(cov_fish$distanceLog))/sd(cov_fish$distanceLog)
   
   ## region
   cov_fish <- cbind(cov_fish,
                     Region=covariates_site$region [match(rownames(cov_fish),rownames(covariates_site$region))])
   cov_fish$Region <- as.factor (cov_fish$Region)
   
   ## locality
   locality <- unlist(lapply (strsplit (covariates_site$site_names,"\\."), function (i) i[1]))
   cov_fish <- cbind(cov_fish,
                     Locality=locality [match(rownames(cov_fish),rownames(covariates_site$region))])
   cov_fish$Locality <- as.factor (cov_fish$Locality)
   
   ; # return
   
   cov_fish
}
)

# ---------------------------------------- #
# organizing data of benthos
# list of SR results
list_data_benthos_SR <- list (res_table_samples_bentos, 
                              res_sp_accum_bentos_asymptote,
                              res_lomolino_bentos)
# object do indicate NAs
list_data_benthos_SR<-lapply (list_data_benthos_SR[1:2], function (i)
   
   cbind (i,
          maintain = ifelse (is.na(i$EST.Rich)==T, "No","Yes")
   )
)
## object to indicate imprecise sample size and richness estimate of lomolino's model
res_lomolino_bentos$maintain <- ifelse (res_lomolino_bentos$nsamples <  res_lomolino_bentos$obs_ss, "Yes","No")

list_data_benthos_SR <- c(list_data_benthos_SR,
                          list(res_lomolino_bentos))

# list of FD results
list_data_benthos_FD <- list (av_FRic_benthos_MSS, 
                              av_FRic_benthos_PSS,
                              av_FRic_benthos_lomolino)

# list of FEve results
list_data_benthos_FEve <- list (av_FEve_benthos_MSS, 
                                av_FEve_benthos_PSS,
                                av_FEve_benthos_lomolino)

# list of FDiv results
list_data_benthos_FDiv <- list (av_FDiv_benthos_MSS, 
                                av_FDiv_benthos_PSS,
                                av_FDiv_benthos_lomolino)

# organizing data across lists

org_data_benthos <- lapply (seq(1,length(list_data_benthos_FD)), function (i) {      
   # covariates_site$sea_data[[1]] is for fishes, [[2]] is for benthos
   cov_benthos <- covariates_site$sea_data[[2]][which(list_data_benthos_SR[[i]]$maintain == "Yes"),]
   # bind richness estimate
   cov_benthos <- data.frame(cov_benthos, 
                             EST.rich = list_data_benthos_SR[[i]]$EST.Rich[which(list_data_benthos_SR[[i]]$maintain == "Yes")])
   # bind FD estimate
   cov_benthos <- cbind(cov_benthos, 
                        FD = list_data_benthos_FD[[i]])
   
   # bind Feve
   cov_benthos <- cbind(cov_benthos, 
                        FEve = list_data_benthos_FEve[[i]])
   
   # bind FDiv
   cov_benthos <- cbind(cov_benthos, 
                        FDiv = list_data_benthos_FDiv[[i]])
   
   # bind coordinates
   cov_benthos <- cbind(cov_benthos,
                        covariates_site$coord$coord_bentos [match(rownames(cov_benthos),covariates_site$coord$coord_bentos$Group.1),c("Lon","Lat","distance")])
   
   # bind organism
   cov_benthos <- cbind(cov_benthos,
                        Organism="Benthos")
   
   # bind method
   cov_benthos <- cbind(cov_benthos,
                        Method=method[i])
   
   # ------------------------
   # organize and standardize covariates
   # inverse of temperature, standardize covariates, insert region and locality
   # ------------------------
   
   boltzmann_factor <- 8.62e-5
   #a<-20:25 + 273.15
   #plot(a, (1/boltzmann_factor * (1/mean (a) - 1/(a))))
   
   # FISH
   ## MMS
   temp_kelvin<- cov_benthos$BO2_tempmean_ss + 273.15
   inv_temp <- 1 / boltzmann_factor * (1 / mean(temp_kelvin) - 1 / temp_kelvin)
   cov_benthos <- cbind(cov_benthos,
                        inv_temp=inv_temp)
   
   ## standardize other variables
   cov_benthos$BO2_tempmean_ss_std <- (cov_benthos$BO2_tempmean_ss-mean(cov_benthos$BO2_tempmean_ss))/sd(cov_benthos$BO2_tempmean_ss)
   cov_benthos$BO2_ppmean_ss_std <- (cov_benthos$BO2_ppmean_ss-mean(cov_benthos$BO2_ppmean_ss))/sd(cov_benthos$BO2_ppmean_ss)
   cov_benthos$BO2_salinitymean_ss_std <- (cov_benthos$BO2_salinitymean_ss-mean(cov_benthos$BO2_salinitymean_ss))/sd(cov_benthos$BO2_salinitymean_ss)
   cov_benthos$BO_damean_std <- (cov_benthos$BO_damean-mean(cov_benthos$BO_damean))/sd(cov_benthos$BO_damean)
   cov_benthos$inv_temp_std <- (cov_benthos$inv_temp-mean(cov_benthos$inv_temp))/sd(cov_benthos$inv_temp)
   cov_benthos$distanceLog <- log(cov_benthos$distance)
   cov_benthos$distance_std <- (cov_benthos$distanceLog-mean(cov_benthos$distanceLog))/sd(cov_benthos$distanceLog)
   
   ## region
   cov_benthos <- cbind(cov_benthos,
                        Region=covariates_site$region [match(rownames(cov_benthos),rownames(covariates_site$region))])
   cov_benthos$Region <- as.factor (cov_benthos$Region)
   
   ## locality
   locality <- unlist(lapply (strsplit (covariates_site$site_names,"\\."), function (i) i[1]))
   cov_benthos <- cbind(cov_benthos,
                        Locality=locality [match(rownames(cov_benthos),rownames(covariates_site$region))])
   cov_benthos$Locality <- as.factor (cov_benthos$Locality)
   
   ## coordinates
   
   
   ; # return
   
   cov_benthos
}
)

## aggregate list of results across organisms by rarefaction method

complete_results_for_fig2 <- lapply (seq (1,length(org_data_fish)), function (i)
   
   rbind(org_data_fish [[i]], org_data_benthos[[i]])
   
)

# -----------------------------------------
## correlation between variables

lapply (org_data_fish, function (i)
   cor (i [,c("BO2_tempmean_ss_std","BO2_ppmean_ss_std","BO2_salinitymean_ss_std",
              "BO_damean_std","inv_temp_std", "distance_std")])
)

# 
lapply (org_data_benthos, function (i)
   cor (i [,c("BO2_tempmean_ss_std","BO2_ppmean_ss_std","BO2_salinitymean_ss_std",
              "BO_damean_std","inv_temp_std", "distance_std")])
   )

### go with (correlation < 0.8) for both datasets

# inv_temp_std, BO2_ppmean_ss_std, BO2_salinitymean_ss_std, BO_damean_std, distance_std

## save results for Figure 2
save (complete_results_for_fig2, 
      file=here("output", "complete_results_for_fig2.RData"))

############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FD RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################

load (here ("output", "complete_results_for_fig2.RData"))

# run model (ancova)
model.ancova <-lapply (complete_results_for_fig2, function (i)
   
   brm (FD ~ EST.rich*Organism,
        data=i,
        family = gaussian (link="identity"),
        chains=3,
        iter = 30000,
        warmup = 20000,
        thin=20)
)
# summary of results
lapply (model.ancova,summary)

# plotting
lapply (seq(1,3), function (i) {

   plot(conditional_effects(model.ancova[[i]],
                         method="fitted",
                         re_formula=NA,
                         robust=T,
                         effects = "EST.rich:Organism",
                         points=T,
                         prob = 0.95),
     
     theme = theme_classic() +
        theme (axis.title = element_text(size=15),
               axis.text = element_text(size=12),
               legend.position = "none"),
     points=T)

ggsave (here("output",paste (i,"ancovaplot.pdf")),height = 5,width=5)

})

# compare slopes
m.lst <- lapply (model.ancova, emtrends, "Organism", var="EST.rich")

m.lst_tab <- lapply(m.lst,summary, point.est = mean)
do.call(rbind,m.lst_tab)

save (model.ancova,m.lst,
      file=here("output", "ancova.RData"))

############################################################################
# -------------------------------------------------------------------------
#          Linear models
#
#        Environmental drivers
#
#           RICHNESS
#
# -------------------------------------------------------------------------
############################################################################

# ----------------
#      BENTHOS SPECIES RICHNESS
#-----------------

# MCMC settings
ni <- 30000
nb <- 20000
nt <- 20
nc <- 3

# set formula (the same for benthos and fishes)
formula <- brms::bf(log(EST.rich) ~ inv_temp +
                       
                       BO2_ppmean_ss_std +
                      
                       BO2_salinitymean_ss_std +
                       
                       BO_damean_std,# +
                       
                       #distance_std,
                    
                    nl = F) # nonlinear

# priors for benthos
priors_benthos <- lapply (org_data_benthos, function (i) 
   
   set_prior("normal(0,5)", class = "b") #+ 
     
     #set_prior("normal(-2,5)", coef = "inv_temp")
      
            
            )

# run MCMC chains for benthos DS
MCMC_runs_benthos_SR <- lapply ( seq (1,length (org_data_benthos)), function (i)
   
   
   brms::brm(formula,
          
          data = org_data_benthos[[i]], 
          
          family = gaussian(),
          
          prior = priors_benthos[[i]], 
          
          chains = nc,
          
          iter = ni,
          
          warmup = nb,
          
          thin=nt
          
          )
   )

## save
save (MCMC_runs_benthos_SR, file=here ("output", "MCMC_runs_benthos_SR.Rdata"))

# ----------------
#      FISH SPECIES RICHNESS
#-----------------

# priors for fish
priors_fish <- lapply (org_data_fish, function (i) 
   
   set_prior("normal(0,5)", class = "b")# + 
     
     #set_prior("normal(-2,5)", coef = "inv_temp")
   
)

# run MCMC chains for fish DS
MCMC_runs_fishes_SR <- lapply ( seq (1,length (org_data_fish)), function (i)
   
   
   brms::brm(formula,
             
             data = org_data_fish[[i]], 
             
             family = gaussian(),
             
             prior = priors_fish[[i]], 
             
             chains = nc,
             
             iter = ni,
             
             warmup = nb,
             
             thin=nt
             
   )
)
# save
save (MCMC_runs_fishes_SR, file=here ("output", "MCMC_runs_fishes_SR.Rdata"))

# ----------------
#      BENTHOS FRIC
#-----------------

# set formula (the same for benthos and fishes)
formula_FRic <- brms::bf((FD) ~ inv_temp +
                       
                       BO2_ppmean_ss_std +
                       
                       BO2_salinitymean_ss_std +
                       
                       BO_damean_std,# +
                       
                       #distance_std,
                    
                    nl = F) # nonlinear

# priors for benthos
priors_benthos <- lapply (org_data_benthos, function (i) 
   
   set_prior("normal(0,5)", class = "b")# + 
      
      #set_prior("normal(-2,5)", coef = "inv_temp")
   
)

# run MCMC chains for benthos DS
MCMC_runs_benthos_FRic <- lapply ( seq (1,length (org_data_benthos)), function (i)
   
   
   brms::brm(formula_FRic,
             
             data = org_data_benthos[[i]], 
             
             family = gaussian(),
             
             prior = priors_benthos[[i]], 
             
             chains = nc,
             
             iter = ni,
             
             warmup = nb,
             
             thin=nt
             
   )
)

## save
save (MCMC_runs_benthos_FRic, file=here ("output", "MCMC_runs_benthos_FRic.Rdata"))


# ----------------
#      FISH FRIC
#-----------------

# priors for fish
priors_fish <- lapply (org_data_fish, function (i) 
   
   set_prior("normal(0,5)", class = "b") #+ 
      
      #set_prior("normal(-2,5)", coef = "inv_temp")
   
)

# run MCMC chains for fish DS
MCMC_runs_fishes_FRic <- lapply ( seq (1,length (org_data_fish)), function (i)
   
   
   brms::brm(formula_FRic,
             
             data = org_data_fish[[i]], 
             
             family = gaussian(),
             
             prior = priors_fish[[i]], 
             
             chains = nc,
             
             iter = ni,
             
             warmup = nb,
             
             thin=nt
             
   )
)
# save
save (MCMC_runs_fishes_FRic, file=here ("output", "MCMC_runs_fishes_FRic.Rdata"))

# ----------------
#      BENTHOS FEVE
#-----------------

# set formula (the same for benthos and fishes)
formula_FEve <- brms::bf((FEve) ~ inv_temp +
                            
                            BO2_ppmean_ss_std +
                            
                            BO2_salinitymean_ss_std +
                            
                            BO_damean_std,# +
                            
                            #distance_std,
                         
                         nl = F) # nonlinear

# priors for benthos
priors_benthos <- lapply (org_data_benthos, function (i) 
   
   set_prior("normal(0,5)", class = "b")# + 
      
     # set_prior("normal(-2,5)", coef = "inv_temp")
   
)

# run MCMC chains for benthos DS
MCMC_runs_benthos_FEve <- lapply ( seq (1,length (org_data_benthos)), function (i)
   
   
   brms::brm(formula_FEve,
             
             data = org_data_benthos[[i]], 
             
             family = gaussian(),
             
             prior = priors_benthos[[i]], 
             
             chains = nc,
             
             iter = ni,
             
             warmup = nb,
             
             thin=nt
             
   )
)

## save
save (MCMC_runs_benthos_FEve, file=here ("output", "MCMC_runs_benthos_FEve.Rdata"))


# ----------------
#      FISH FEVE
#-----------------

# priors for fish
priors_fish <- lapply (org_data_fish, function (i) 
   
   set_prior("normal(0,5)", class = "b") #+ 
      
     # set_prior("normal(-2,5)", coef = "inv_temp")
   
)

# run MCMC chains for fish DS
MCMC_runs_fishes_FEve <- lapply ( seq (1,length (org_data_fish)), function (i)
   
   
   brms::brm(formula_FEve,
             
             data = org_data_fish[[i]], 
             
             family = gaussian(),
             
             prior = priors_fish[[i]], 
             
             chains = nc,
             
             iter = ni,
             
             warmup = nb,
             
             thin=nt
             
   )
)
# save
save (MCMC_runs_fishes_FEve, file=here ("output", "MCMC_runs_fishes_FEve.Rdata"))


# ----------------
#      BENTHOS FDiv
#-----------------

# set formula (the same for benthos and fishes)
formula_FDiv <- brms::bf((FDiv) ~ inv_temp +
                            
                            BO2_ppmean_ss_std +
                            
                            BO2_salinitymean_ss_std +
                            
                            BO_damean_std,# +
                         
                         #distance_std,
                         
                         nl = F) # nonlinear

# priors for benthos
priors_benthos <- lapply (org_data_benthos, function (i) 
   
   set_prior("normal(0,5)", class = "b")# + 
   
   # set_prior("normal(-2,5)", coef = "inv_temp")
   
)

# run MCMC chains for benthos DS
MCMC_runs_benthos_FDiv <- lapply ( seq (1,length (org_data_benthos)), function (i)
   
   
   brms::brm(formula_FDiv,
             
             data = org_data_benthos[[i]], 
             
             family = gaussian(),
             
             prior = priors_benthos[[i]], 
             
             chains = nc,
             
             iter = ni,
             
             warmup = nb,
             
             thin=nt
             
   )
)

## save
save (MCMC_runs_benthos_FDiv, file=here ("output", "MCMC_runs_benthos_FDiv.Rdata"))


# ----------------
#      FISH FDIV
#-----------------

# priors for fish
priors_fish <- lapply (org_data_fish, function (i) 
   
   set_prior("normal(0,5)", class = "b") #+ 
   
   # set_prior("normal(-2,5)", coef = "inv_temp")
   
)

# run MCMC chains for fish DS
MCMC_runs_fishes_FDiv <- lapply ( seq (1,length (org_data_fish)), function (i)
   
   
   brms::brm(formula_FDiv,
             
             data = org_data_fish[[i]], 
             
             family = gaussian(),
             
             prior = priors_fish[[i]], 
             
             chains = nc,
             
             iter = ni,
             
             warmup = nb,
             
             thin=nt
             
   )
)
# save
save (MCMC_runs_fishes_FDiv, file=here ("output", "MCMC_runs_fishes_FDiv.Rdata"))


# -----------------------------------------------

## List of results for figure 3

## fish results
load (here ("output","MCMC_runs_fishes_SR.RData"))
load (here ("output","MCMC_runs_fishes_FRic.RData"))
load (here ("output","MCMC_runs_fishes_FEve.RData"))
load (here ("output","MCMC_runs_fishes_FDiv.RData"))

# benthos
load (here ("output","MCMC_runs_benthos_SR.RData"))
load (here ("output","MCMC_runs_benthos_FRic.RData"))
load (here ("output","MCMC_runs_benthos_FEve.RData"))
load (here ("output","MCMC_runs_benthos_FDiv.RData"))

# FISH
## coefficients
list_of_results_fish <- list(MCMC_runs_fishes_SR,
                             MCMC_runs_fishes_FRic, 
                             MCMC_runs_fishes_FEve,
                             MCMC_runs_fishes_FDiv
)


# set names in indexes
names(list_of_results_fish) <- c("SR", "FRic", "FEve", "FDiv")
# set names of algorithms within indexes 
list_of_results_fish <- lapply (list_of_results_fish, function (i) {
   names (i) <- c("MSS","PSS", "Lomolino");
   i
})

## organizing results (coeficients)

org_results <- lapply (seq(1,length (list_of_results_fish)), function (i)
   lapply (seq (length(list_of_results_fish[[1]])), function (k)
      
      data.frame (
         fixef(list_of_results_fish[[i]][[k]],
               summary = TRUE,
               robust = FALSE,
               probs = c(0.025, 0.975)
         ), 
         "Algorithm" = names(list_of_results_fish[[i]])[k],
         "Index" = names(list_of_results_fish)[i])
   ))

# melt list 
org_results <- do.call(rbind, 
                       
                       lapply (org_results, function (i)
                          
                          do.call(rbind , i)
                          
                       )
)

# bind parameter
org_results$Parameter <- c("Intercept",
                           "Temperature",
                           "Productivity",
                           "Salinity",
                           "Turbidity")
# bind organism
org_results$Organism <- "Fishes"

# -----------------------------------------------
# BENTHOS

## coefficients
list_of_results_benthos <- list(MCMC_runs_benthos_SR,
                                MCMC_runs_benthos_FRic, 
                                MCMC_runs_benthos_FEve,
                                MCMC_runs_benthos_FDiv
)


# set names in indexes
names(list_of_results_benthos) <- c("SR", "FRic", "FEve", "FDiv")
# set names of algorithms within indexes 
list_of_results_benthos <- lapply (list_of_results_benthos, function (i) {
   names (i) <- c("MSS","PSS", "Lomolino");
   i
})

## organizing results (coeficients)

org_results_benthos <- lapply (seq(1,length (list_of_results_benthos)), function (i)
   lapply (seq (length(list_of_results_benthos[[1]])), function (k)
      
      data.frame (
         fixef(list_of_results_benthos[[i]][[k]],
               summary = TRUE,
               robust = FALSE,
               probs = c(0.025, 0.975)
         ), 
         "Algorithm" = names(list_of_results_benthos[[i]])[k],
         "Index" = names(list_of_results_benthos)[i])
   ))

# melt list 
org_results_benthos <- do.call(rbind, 
                               
                               lapply (org_results_benthos, function (i)
                                  
                                  do.call(rbind , i)
                                  
                               )
)

# bind parameter
org_results_benthos$Parameter <- c("Intercept",
                                   "Temperature",
                                   "Productivity",
                                   "Salinity",
                                   "Turbidity")
# bind organism
org_results_benthos$Organism <- "Benthos"

## bind fish and benthos data

complete_results <- rbind (org_results,
                           org_results_benthos)

# adjusting
colnames (complete_results)[3:4] <- c("lower","upper")

complete_results$Index <- factor (complete_results$Index,
                                  levels = c("SR", "FRic","FEve", "FDiv"))

## save results
save (complete_results, file=here("output", "complete_results_for_fig3.RData"))


### check of goodness of fit

## fit we want to check
# fishes
# MSS
loo_test_fish_MSS <- lapply (list(MCMC_runs_fishes_SR[[1]],
             MCMC_runs_fishes_FRic[[1]],
             MCMC_runs_fishes_FEve[[1]],
             MCMC_runs_fishes_FDiv[[1]]),
             
             loo)


do.call(rbind, 
        lapply (loo_test_fish_MSS, function (i) i$estimates [2,]))

# all pareto's k lower than 0.5

lapply (loo_test_fish_MSS, pareto_k_table)

# PSSi
loo_test_fish_PSSi <- lapply (list(MCMC_runs_fishes_SR[[2]],
             MCMC_runs_fishes_FRic[[2]],
             MCMC_runs_fishes_FEve[[2]],
             MCMC_runs_fishes_FDiv[[2]]),
             
             loo)

do.call(rbind, 
        lapply (loo_test_fish_PSSi, function (i) i$estimates [2,]))

# all pareto's k lower than 0.5

lapply (loo_test_fish_PSSi, pareto_k_table)

## Lomolino's 

# PSSi
loo_test_fish_lomolino <- lapply (list(MCMC_runs_fishes_SR[[3]],
                                   MCMC_runs_fishes_FRic[[3]],
                                   MCMC_runs_fishes_FEve[[3]],
                                   MCMC_runs_fishes_FDiv[[3]]),
                              
                              loo)

do.call(rbind, 
        lapply (loo_test_fish_lomolino, function (i) i$estimates [2,]))

# all pareto's k lower than 0.5

lapply (loo_test_fish_lomolino, pareto_k_table)

## benthos

# MSS
loo_test_benthos_MSS <- lapply (list(MCMC_runs_benthos_SR[[1]],
                                  MCMC_runs_benthos_FRic[[1]],
                                  MCMC_runs_benthos_FEve[[1]],
                                  MCMC_runs_benthos_FDiv[[1]]),
                                
                                loo)


do.call(rbind, 
        lapply (loo_test_benthos_MSS, function (i) i$estimates [2,]))

# all pareto's k lower than 0.5

lapply (loo_test_benthos_MSS, pareto_k_table)

# PSSi
loo_test_benthos_PSSi <- lapply (list(MCMC_runs_benthos_SR[[2]],
                                   MCMC_runs_benthos_FRic[[2]],
                                   MCMC_runs_benthos_FEve[[2]],
                                   MCMC_runs_benthos_FDiv[[2]]),
                                 
                                 loo)

do.call(rbind, 
        lapply (loo_test_benthos_PSSi, function (i) i$estimates [2,]))

# all pareto's k lower than 0.5

lapply (loo_test_benthos_PSSi, pareto_k_table)

## lomolino

loo_test_benthos_lomolino <- lapply (list(MCMC_runs_benthos_SR[[3]],
                                      MCMC_runs_benthos_FRic[[3]],
                                      MCMC_runs_benthos_FEve[[3]],
                                      MCMC_runs_benthos_FDiv[[3]]),
                                 
                                 loo)

do.call(rbind, 
        lapply (loo_test_benthos_lomolino, function (i) i$estimates [2,]))

# all pareto's k lower than 0.5

lapply (loo_test_benthos_lomolino, pareto_k_table)
