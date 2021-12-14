
## ------------------------------------------------------------ ##
##                 Analysis of functional diversity             ##
## ------------------------------------------------------------ ##
## laoding packages
source("R/packages.R")
source("R/functions.R")

## function to test space quality (from Maire et al. 2015)
source("R/quality_funct_space_fromdist2.R")

# ------------------------------------------ #
# Load rarefied data of fishes and benthos
# sites x transect / video x species
# ------------------------------------------ #

load (here ("output","random_composition_bentos.RData"))
load (here ("output","random_composition_fish.RData"))

# ---------------------------------- #
# Load fish trait data
# ---------------------------------- #

fish_traits <- read.csv(here("data","traits","Atributos_especies_Atlantico_&_Pacifico_Oriental_2020_04_28.csv"), 
                        header=TRUE, sep=";")
### Adjust spp names in trait dataset
fish_traits$Name <- tolower(gsub(" ", ".", fish_traits$Name)) 

### transforming characters into numeric
fish_traits$Body_size <- as.numeric (gsub (",",".",fish_traits$Body_size))
fish_traits$TempPref_min <- as.numeric (gsub (",",".",fish_traits$TempPref_min))
fish_traits$TemPref_mean <- as.numeric (gsub (",",".",fish_traits$TemPref_mean))
fish_traits$TempPref_max <- as.numeric (gsub (",",".",fish_traits$TempPref_max))
fish_traits$Aspect_ratio <- as.numeric (gsub (",",".",fish_traits$Aspect_ratio))
fish_traits$Trophic_level <- as.numeric (gsub (",",".",fish_traits$Trophic_level))
fish_traits$PLD <- as.numeric (gsub (",",".",fish_traits$PLD))
fish_traits$Vulnerability <- as.numeric (gsub (",",".",fish_traits$Vulnerability))
fish_traits$Depth_min <- as.numeric (gsub (",",".",fish_traits$Depth_min))
fish_traits$Depth_max <- as.numeric (gsub (",",".",fish_traits$Depth_max))
fish_traits$Depth_range <- as.numeric (gsub (",",".",fish_traits$Depth_range))

### editing diet (removing NA and empty cells)
fish_traits <- fish_traits [which(is.na(fish_traits$Diet) ==F),]
fish_traits<-fish_traits [-which(fish_traits$Diet == ""),]

### Select fish species from Atlantic database to match UVCs
UVS_spp <- unique(unlist(lapply (rdm_composition_complete, colnames)))
fish_traits <- fish_traits[which(fish_traits$Name %in% UVS_spp),]

### df with sp names
fishes <- data.frame (Name=UVS_spp)
### subsetting traits based on defined groups
fish_traits_subset <- fish_traits
### matching
fish_traits2 <- merge(fishes, fish_traits_subset, by="Name")

### Select traits of interest and remove NAs
fish_traits2 <- na.omit(fish_traits2[,c("Name","Home_range","Diel_activity","Size_group","Body_size",
                                        "Level_water","Body_shape","Diet","Caudal_fin","Mouth_position")])
# dim(fish_traits2)

# Adjust traits into ordered categories
### body size is continuous
rownames(fish_traits2) <- fish_traits2$Name
fish_traits2 <- fish_traits2[,-1]
fish_traits2 <- fish_traits2[order(rownames(fish_traits2)),]
fish_traits2$Body_size <- as.numeric(fish_traits2$Body_size)

### transform body size
size <- (fish_traits2$Body_size) # log the size 
# discreticize
size  <- sapply(fish_traits2$Body_size, function(x) {if (x<=50) {1} 
   else if (x>=50&x<100) {2} 
   else if (x>=100&x<150) {3}
   else if (x>=150) {4}}
)

### here transform each ordered category in 1, 2, 3.... 
mobility  <- sapply(fish_traits2$Home_range, function(x) {if (x=="sed") {1} 
   else if (x=="mob") {2} 
   else if (x=="vmob") {3}}
)

schooling <- sapply(fish_traits2$Size_group, function(x) {if (x=="sol") {1} 
   else if (x=="pair") {2} 
   else if (x=="smallg") {3} 
   else if (x=="medg") {4} 
   else if (x=="largeg") {5}}
)

level  <- sapply(fish_traits2$Level_water, function(x) {if (x=="bottom") {1} 
   else if (x=="low") {2} 
   else if (x=="high") {3}}
)

### ordering trait values
mobility  <-ordered (mobility)
level  <-ordered (level) 
schooling <-ordered (schooling) 

### all traits into a dataframe
fish_traits_ord <-data.frame (Size=ordered(size), 
                              Mobility=ordered (mobility), 
                              #Activity=as.factor (fish_traits2$Diel_activity), 
                              Schooling=ordered (schooling), 
                              #Level=ordered (level), 
                              #Diet=as.factor(fish_traits2$Diet), 
                              Body_shape=as.factor(fish_traits2$Body_shape),
                              Caudal_fin=as.factor(fish_traits2$Caudal_fin), 
                              Mouth_position=as.factor(fish_traits2$Mouth_position))

rownames(fish_traits_ord) <- rownames(fish_traits2)

#----------------------------------#
#      Load benthos trait dataset
# ---------------------------------#

bent_traits <- read.csv(here("data","traits","Database_benthos.csv"), 
                        header=TRUE, sep=";")

# Remove benthos from dataset as below
to_remove <- c("areia.e.cascalho","desconhecido","estrela","ourico1","ourico2",
               "outra.ascidia","outro.anthozoa", "outro.echinoderma",
               "outro.hydrozoa","quadrado","outro.crustaceo","sombra")

bent2 <- bent_traits[-which(bent_traits$groups %in% to_remove),]

# list of benthic spp
benthos_spp <- data.frame(groups=bent2$groups)

## matching these data   
bent_traits2 <- merge(bent2, benthos_spp, by="groups") # merge traits with groups in occurrence data
rownames(bent_traits2) <- bent_traits2$groups
bent_traits2 <- bent_traits2[2:8]

# Transform benthos body size in ordered category

bent_traits2 <- bent_traits2[order(rownames(bent_traits2)),]
benthos_body_size <- sapply(bent_traits2$body_size, function(x) {
   if (x=="S") {1} 
   else if (x=="M") {2} 
   else if (x=="L") {3} 
   else if (x=="XL") {4}}
)

benthos_body_size  <-ordered (benthos_body_size) #order traits

# Set data for benthos traits with ordered categories
bent_traits_ord <-data.frame (body_size=ordered (benthos_body_size), 
                              growth_form=as.factor(bent_traits2$growth_form),
                              modularity=as.factor(bent_traits2$modularity), 
                              mobility=as.factor(bent_traits2$mobility),
                              reproductive_mode=as.factor(bent_traits2$reproductive_mode), 
                              carbonate=as.factor(bent_traits2$carbonate.accretion))

rownames(bent_traits_ord) <- rownames(bent_traits2)

# --------------------------------------------------------------------- #
#    Functional diversity per site, data set, and sample size def
# --------------------------------------------------------------------- #

# ------
# Fish
#-------

### tests using composition obtained by the minimum sample size (MSS)
nc <- 3
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores
clusterEvalQ(cl, library(FD))
clusterEvalQ(cl, library(cluster))
clusterEvalQ(cl, library(clue))

# export your data and function
clusterExport(cl, c("rdm_composition_complete", 
                    "fish_traits_ord",
                    "function_FD_fish_abundW",
                    "quality_funct_space_fromdist"))

FD_fish_MSS_abundW <- parLapply (cl, 
                          rdm_composition_complete, function (i) {
   
                           tryCatch ( 
                              function_FD_fish_abundW (i, fish_traits_ord),
                                      error = function (e)
                                        return (e)
                           )

                          })

stopCluster (cl)

# save it
save (FD_fish_MSS_abundW, file=here("output","FD_fish_MSS_abundW.RData"))

# ------
# Benthos
#-------

### tests using composition obtained by the minimum sample size (MSS)
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores
clusterEvalQ(cl, library(FD))
clusterEvalQ(cl, library(cluster))
clusterEvalQ(cl, library(clue))

# export your data and function
clusterExport(cl, c("rdm_composition_complete_bentos", 
                    "bent_traits_ord",
                    "function_FD_benthos_abundW",
                    "quality_funct_space_fromdist"))

FD_benthos_MSS_abundW <- parLapply (cl, 
                           rdm_composition_complete_bentos, function (i){
                                 
                                 tryCatch (function_FD_benthos_abundW (i, bent_traits_ord),
                                            error = function (e) 
                                              return (e))
                                 }
                           
)

stopCluster (cl)

save (FD_benthos_MSS_abundW, file=here("output","FD_benthos_MSS_abundW.RData"))

