
## ------------------------------------------------------------ ##
##                 Analysis of functional diversity             ##
## ------------------------------------------------------------ ##

## chamar os pacotes
source("R/packages.R")
source("R/functions.R")
## function to test space quality
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
UVS_spp <- unique(unlist(lapply (rdm_composition_asymptote_complete, colnames)))
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
size <-log (fish_traits2$Body_size) # log the size 

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
fish_traits_ord <-data.frame (Size=as.numeric(size), 
                              Mobility=ordered (mobility), 
                              Activity=as.factor (fish_traits2$Diel_activity), 
                              Schooling=ordered (schooling), 
                              Level=ordered (level), 
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

# function

function_FD <- function (site.data, spp.trait.data) {
      
      #---------------------------------------------------#
      # Adjust fish abundance data to relative abundance
      #-------------------------------------------------- #
      
      abund_matrix <- site.data[,which(colnames(site.data) %in% rownames(spp.trait.data))] 
      rel_abund <- abund_matrix/rowSums(abund_matrix) # relative abundance/cover
      quais_sitios_manter <-which(is.na (rowSums (rel_abund))!= T) 
      rel_abund <- rel_abund [quais_sitios_manter,]
      rel_abund <- rel_abund [,which(colSums (rel_abund) > 0)]
      
      # and match fish spp in both datasets
      fish_traits_ord_match <- spp.trait.data[which(rownames(spp.trait.data) %in% colnames (rel_abund)),]
      
      ## ---------------------------------- ##
      ## Building the gower distance matrix ##
      ## ---------------------------------- ##
      
      gower_matrix <- daisy (fish_traits_ord_match, metric=c("gower")) 
      
      ## ---------------------------------------------- ##
      ## Building the functional space based on a PCOA 
      ## quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
      ## ---------------------------------------------- ##
      
      pco<-dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10) 
      
      ### barplot of eigenvalues for each axis 
      barplot(pco$eig) 
      
      ### calculate and show percentage of inertia explained by the two first axes
      (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) 
      
      ## ----------------------------------------------------- ##
      # Testing the Quality of the functional space based on the method 
      # of Maire et al. (2015) in GEB, i.e. how many axes do we need to keep to 
      # faithfully represent the original Gower's distances.
      # The function will test the quality of the space from 2 to n axes using 
      # dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
      ## ------------------------------------------------------ ##
      
      quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
                                             plot="quality_funct_space_I") 
      
      ### the minimal value corresponds to the best space to use (min SD)
      axes_to_choose <- which(quality$meanSD == min(quality$meanSD)) 
      ### calculate and show the percentage of inertia explained by the choosen axes
      (Inertia7<- (sum(pco$eig[1:axes_to_choose])) /(sum(pco$eig))) 
      
      ### Method with dbFD() 
      fd <- dbFD (gower_matrix, rel_abund, 
                  corr ="none", 
                  w.abun = T,
                  m = 5, 
                  stand.FRic = TRUE,
                  print.pco = TRUE,
                  calc.FGR = F,
                  clust.type = "kmeans")
      
      ### in the case of samples with too few species
      fd_indexes  <- data.frame (matrix (NA, ncol=8,nrow=nrow(rel_abund),
                                         dimnames = list (NULL,
                                                          c("nbsp",
                                                            "sing",
                                                            "FRic",
                                                            "qual.FRic",
                                                            "FEve",
                                                            "FDiv",
                                                            "FDis",
                                                            "RaoQ"))))
      
      fd_indexes [quais_sitios_manter,1] <- fd$nbsp
      fd_indexes[quais_sitios_manter,2] <- fd$sing.sp
      fd_indexes[quais_sitios_manter,3] <- fd$FRic
      fd_indexes[quais_sitios_manter,4] <- fd$qual.FRic
      fd_indexes[quais_sitios_manter,5] <- fd$FEve
      fd_indexes[quais_sitios_manter,6] <- fd$FDiv
      fd_indexes[quais_sitios_manter,7] <- fd$FDis
      fd_indexes[quais_sitios_manter,8] <- fd$RaoQ
      
      ### list with results
      results <- list (Fdindexes = fd_indexes,
                       chosenAxes = axes_to_choose,
                       InertiaPCO=Inertia2,
                       InertiaQuality=Inertia7,
                       explanAxes=pco$eig/sum(pco$eig),
                       convexHullVolume=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$vol,
                       convexHullArea=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$area,
                       axesPCO=fd$x.axes)

      return (results)
      
}

# ------
# Fish
#-------

### tests using composition obtained by the minimum sample size (MSS)
nc <- 3
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores
clusterEvalQ(cl, library(FD))
clusterEvalQ(cl, library(cluster))

# export your data and function
clusterExport(cl, c("rdm_composition_complete", 
                    "fish_traits_ord",
                    "function_FD",
                    "quality_funct_space_fromdist"))

FD_fish_MSS <- parLapply (cl, 
                          rdm_composition_complete, function (i)
   
                             function_FD (i, fish_traits_ord)

                          )

stopCluster (cl)

save (FD_fish_MSS, file=here("output","FD_fish_MSS.RData"))

### tests using composition obtained by the asymptotic sample size (ASSi)
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores
clusterEvalQ(cl, library(FD))
clusterEvalQ(cl, library(cluster))

# export your data and function
clusterExport(cl, c("rdm_composition_asymptote_complete", 
                    "fish_traits_ord",
                    "function_FD",
                    "quality_funct_space_fromdist"))

FD_fish_ASSi <- parLapply (cl, 
                           rdm_composition_asymptote_complete, function (i)
                             
                             function_FD (i, fish_traits_ord)
                          
)

stopCluster (cl)

save (FD_fish_ASSi, file=here("output","FD_fish_ASSi.RData"))

# ------
# Benthos
#-------

### tests using composition obtained by the minimum sample size (MSS)
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores
clusterEvalQ(cl, library(FD))
clusterEvalQ(cl, library(cluster))

# export your data and function
clusterExport(cl, c("rdm_composition_complete_bentos", 
                    "bent_traits_ord",
                    "function_FD",
                    "quality_funct_space_fromdist"))

FD_benthos_MSS <- parLapply (cl, 
                           rdm_composition_complete_bentos, function (i)
                              
                              function_FD (i, bent_traits_ord)
                           
)

stopCluster (cl)

save (FD_benthos_MSS, file=here("output","FD_benthos_MSS.RData"))

### tests using composition obtained by the asymptotic sample size (ASSi)

cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores
clusterEvalQ(cl, library(FD))
clusterEvalQ(cl, library(cluster))

# export your data and function
clusterExport(cl, c("rdm_composition_asymptote_complete_bentos", 
                    "bent_traits_ord",
                    "function_FD",
                    "quality_funct_space_fromdist"))

FD_benthos_ASSi <- parLapply (cl, 
                             rdm_composition_asymptote_complete_bentos, function (i)
                                
                                function_FD (i, bent_traits_ord)
                             
)

stopCluster (cl)

save (FD_benthos_ASSi, file=here("output","FD_benthos_ASSi.RData"))

