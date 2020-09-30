##### Analysis of functional diversity patterns of Brazilian reefs ###

## chamar os pacotes
source("R/packages.R")
source("R/functions.R")
## function to test space quality
source("R/quality_funct_space_fromdist2.R")

#### Load fish and benthos occurrence data
load (here ("output","data_drivers_analysis.RData"))

#### Run script on "codigo_organizacao_dados.R" until 
# dados_peixes_bentos <- list(peixes = comp_peixes, bentos = comp_bentos)

fish <- dados_peixes_bentos[[1]]
bent <- dados_peixes_bentos[[2]]

#### Load trait data for fish and benthos
## Fish traits

fish_traits <- read.csv(here("data","traits","Atributos_especies_Atlantico_&_Pacifico_Oriental_2020_04_28.csv"), 
                        header=TRUE, sep=";")
fish_traits$Name <- tolower(gsub(" ", ".", fish_traits$Name)) #Adjust spp names in trait dataset

## transforming characters into numeric
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

# Select fish species from Atlantic database to match UVCs

fishes <- data.frame(Name=colnames(fish))
# matching considering all species, and only herbivorous and sessile/mob inv

list_fish_diet <- list( all = c("hd","hm","om","pk","fc","im","is"),
                        herbivorous = c("hd","hm"),
                        sessileInv = c("is","im"),
                        nonHS = c("om","pk","fc"))

## doing that for all three sets of species at once 
## two functions
### 1 - test the overall functional space belonging to each group
### 2 - test the contribution of each group to the overall functional space

# function 1
FD_results_f1 <- lapply (list_fish_diet, function (i) {
   
   # subsetting traits
   fish_traits_subset <- fish_traits [which (fish_traits$Diet %in% i),]
   # matching
   fish_traits2 <- merge(fishes, fish_traits_subset, by="Name")
   # dim(fish_traits2)
   
   # Select traits of interest and remove NAs
   fish_traits2 <- na.omit(fish_traits2[,c("Name","Home_range","Diel_activity","Size_group","Body_size",
                                           "Level_water","Body_shape","Diet","Caudal_fin","Mouth_position")])
   # dim(fish_traits2)
   
   # Adjust traits into ordered categories
   # body size is continuous
   rownames(fish_traits2) <- fish_traits2$Name
   fish_traits2 <- fish_traits2[,-1]
   fish_traits2 <- fish_traits2[order(rownames(fish_traits2)),]
   fish_traits2$Body_size <- as.numeric(fish_traits2$Body_size)
   ## transform body size
   size <-log (fish_traits2$Body_size) # log the size 
   
   # here transform each ordered category in 1, 2, 3.... 
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
   
   ## ordering trait values
   mobility  <-ordered (mobility)
   level  <-ordered (level) 
   schooling <-ordered (schooling) 
   
   ## all traits into a dataframe
   fish_traits_ord <-data.frame (Size=as.numeric(size), 
                                 Mobility=ordered (mobility), 
                                 Activity=as.factor (fish_traits2$Diel_activity), 
                                 Schooling=ordered (schooling), 
                                 Level=ordered (level), 
                                 Diet=as.factor(fish_traits2$Diet), 
                                 Body_shape=as.factor(fish_traits2$Body_shape),
                                 Caudal_fin=as.factor(fish_traits2$Caudal_fin), 
                                 Mouth_position=as.factor(fish_traits2$Mouth_position))
   
   rownames(fish_traits_ord) <- rownames(fish_traits2)
   # dim(fish_traits_ord)
   
   ######## Building the gower distance matrix ########
   
   # first calculate gower distance on traits
   gower_matrix <- daisy (fish_traits_ord, metric=c("gower")) 
   
   # Building the functional space based on a PCOA 
   pco<-dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
   barplot(pco$eig) # barplot of eigenvalues for each axis 
   (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) # percentage of inertia explained by the two first axes
   
   # Testing the Quality of the functional space based on the method of Maire et al. (2015) in GEB, i.e. how many axes do we need to keep to faithfully represent the original Gower's distances
   quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
                                          plot="quality_funct_space_I") # The function will test the quality of the space from 2 to 10 axes using dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
   axes_to_choose <- which(quality$meanSD == min(quality$meanSD)) # the minimal value corresponds to the best space to use, here 7 axes 
   (Inertia7<- (sum(pco$eig[1:axes_to_choose])) /(sum(pco$eig))) # percentage of inertia explained by the 6 first axes = 76%
   
   # Functional diversity indices for fish #######
   # Adjust fish abundance data to relative abundance
   fish_abund <- fish[,which(colnames(fish) %in% rownames(fish_traits_ord))] # match fish spp in both datasets
   rownames(fish_abund) <- fish$locality_site
   fish_rel <- fish_abund/rowSums(fish_abund) # relative abundance for fish
   
   # Method with dbFD() 
   fd <- dbFD (gower_matrix, fish_rel, 
               corr ="none", 
               m = axes_to_choose, 
               stand.FRic = TRUE,
               print.pco = TRUE)
   
   # building a DF with results
   fish_fd <- data.frame (nbsp=fd$nbsp, sing.sp=fd$sing.sp, FRic=fd$FRic, qual.FRic=fd$qual.FRic,
                          FEve=fd$FEve, FDiv=fd$FDiv, FDis=fd$FDis, RaoQ=fd$RaoQ)
   # and then a list
   results <- list (Fdindexes = fish_fd,
                    chosenAxes = axes_to_choose,
                    InertiaPCO=Inertia2,
                    InertiaQuality=Inertia7,
                    convexHullVolum=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$vol,
                    convexHullArea=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$area,
                    axesPCO=pco$li)
   ; ## return
   
   results
   
   
})

# function 2
FD_results_f2 <- lapply (list_fish_diet, function (i) {
   
   # matching
   fish_traits2 <- merge(fishes, fish_traits, by="Name")
   # dim(fish_traits2)
   
   # Select traits of interest and remove NAs
   fish_traits2 <- na.omit(fish_traits2[,c("Name","Home_range","Diel_activity","Size_group","Body_size",
                                           "Level_water","Body_shape","Diet","Caudal_fin","Mouth_position")])
   # dim(fish_traits2)
   
   # Adjust traits into ordered categories
   # body size is continuous
   rownames(fish_traits2) <- fish_traits2$Name
   fish_traits2 <- fish_traits2[,-1]
   fish_traits2 <- fish_traits2[order(rownames(fish_traits2)),]
   fish_traits2$Body_size <- as.numeric(fish_traits2$Body_size)
   ## transform body size
   size <-log (fish_traits2$Body_size) # log the size 
   
   # here transform each ordered category in 1, 2, 3.... 
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
   
   ## ordering trait values
   mobility  <-ordered (mobility)
   level  <-ordered (level) 
   schooling <-ordered (schooling) 
   
   ## all traits into a dataframe
   fish_traits_ord <-data.frame (Size=as.numeric(size), 
                                 Mobility=ordered (mobility), 
                                 Activity=as.factor (fish_traits2$Diel_activity), 
                                 Schooling=ordered (schooling), 
                                 Level=ordered (level), 
                                 Diet=as.factor(fish_traits2$Diet), 
                                 Body_shape=as.factor(fish_traits2$Body_shape),
                                 Caudal_fin=as.factor(fish_traits2$Caudal_fin), 
                                 Mouth_position=as.factor(fish_traits2$Mouth_position))
   
   rownames(fish_traits_ord) <- rownames(fish_traits2)
   # dim(fish_traits_ord)
   
   ######## Building the gower distance matrix ########
   
   # first calculate gower distance on traits
   gower_matrix <- daisy (fish_traits_ord, metric=c("gower")) 
   
   ## subsetting the gower matrix 
   fish_traits_subset <- fish_traits [which (fish_traits$Diet %in% i),"Name"]
   ## 
   gower_matrix_subset <- as.dist(as.matrix(gower_matrix) [which(attr (gower_matrix,"Labels") %in% fish_traits_subset),
                            which(attr (gower_matrix,"Labels") %in% fish_traits_subset)])
   
   # Building the functional space based on a PCOA 
   pco<-dudi.pco(quasieuclid(gower_matrix_subset), scannf=F, nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
   barplot(pco$eig) # barplot of eigenvalues for each axis 
   (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) # percentage of inertia explained by the two first axes
   
   # Testing the Quality of the functional space based on the method of Maire et al. (2015) in GEB, i.e. how many axes do we need to keep to faithfully represent the original Gower's distances
   quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
                                          plot="quality_funct_space_I") # The function will test the quality of the space from 2 to 10 axes using dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
   axes_to_choose <- which(quality$meanSD == min(quality$meanSD)) # the minimal value corresponds to the best space to use, here 7 axes 
   (Inertia7<- (sum(pco$eig[1:axes_to_choose])) /(sum(pco$eig))) # percentage of inertia explained by the 6 first axes = 76%
   
   # Functional diversity indices for fish #######
   # Adjust fish abundance data to relative abundance
   fish_abund <- fish[,which(colnames(fish) %in% rownames(pco$tab))] # match fish spp in both datasets
   rownames(fish_abund) <- fish$locality_site
   fish_rel <- fish_abund/rowSums(fish_abund) # relative abundance for fish
   
   # Method with dbFD() 
   fd <- dbFD (gower_matrix_subset, fish_rel, 
               corr ="none", 
               m = axes_to_choose, 
               stand.FRic = TRUE,
               print.pco = TRUE)
   
   # building a DF with results
   fish_fd <- data.frame (nbsp=fd$nbsp, 
                          sing.sp=fd$sing.sp, 
                          FRic=fd$FRic, 
                          qual.FRic=fd$qual.FRic,
                          FEve=fd$FEve, 
                          FDiv=fd$FDiv, 
                          FDis=fd$FDis, 
                          RaoQ=fd$RaoQ)
   # and then a list
   results <- list (Fdindexes = fish_fd,
                    chosenAxes = axes_to_choose,
                    InertiaPCO=Inertia2,
                    InertiaQuality=Inertia7,
                    convexHullVolum=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$vol,
                    convexHullArea=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$area,
                    axesPCO=pco$li)
   ; ## return
   
   results
   
})

## compare the volum of the two functions
## check area and volum (function1)
round(sapply (FD_results_f1, "[[","convexHullArea"),6)
round(sapply (FD_results_f1, "[[","convexHullVolum"),6)

## check area and volum (function2)
round(sapply (FD_results_f2, "[[","convexHullArea"),6)
round(sapply (FD_results_f2, "[[","convexHullVolum"),6)

## check explained variation (function1)
round(sapply (FD_results_f1, "[[","InertiaQuality"),6)
round(sapply (FD_results_f1, "[[","InertiaPCO"),6)

## check explained variation (function2)
round(sapply (FD_results_f2, "[[","InertiaQuality"),5)
round(sapply (FD_results_f2, "[[","InertiaPCO"),5)

# write.table (fish_fd, "FIndex_fish.csv", sep=";", dec=".")

#################################################################

# Analysis for benthos

# Remove benthos from dataset as below

bent2 <- bent[ , -which(names(bent) %in% c("Areia.e.Cascalho","Desconhecido","Estrela","ourico1","ourico2",
                                           "Outra.ascidia","Outro.anthozoa", "Outro.echinoderma",
                                           "Outro.hydrozoa","Quadrado","Outro.crustaceo","Sombra"))]

colnames(bent2) <- tolower(gsub(" ", ".", colnames(bent2))) #adjust colnames

benthos_spp <- as.data.frame(colnames(bent2[2:96]))
colnames(benthos_spp) <- paste("groups") 

# Work in benthos traits

# Load trait data
bent_traits <- read.csv(here("data","traits","Database_benthos.csv"), 
                        header=TRUE, sep=",")

bent_traits$groups <- tolower(gsub("_", ".", bent_traits$groups)) #adjust spp names in trait dataset

## creating a list of traits we are interested in
list_benthos_level <- list (all = unique(bent_traits$trophic_type),
                            autot = c("A","S"),
                            mixo = "AC")

## doing that for all three sets of species at once 
## two functions
### 1 - test the overall functional space belonging to each group
### 2 - test the contribution of each group to the overall functional space

# function 1
FD_results_f1_bentos <- lapply (list_benthos_level, function (i) {
   
   ## subsetting traits according to predefined groups
   bent_traits_subset <- bent_traits [which(bent_traits$trophic_type %in% i), ]
   
   ## matching these data   
   bent_traits2 <- merge(bent_traits_subset, benthos_spp, by="groups") # merge traits with groups in occurrence data
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
                                 #modularity=as.factor(bent_traits2$modularity), 
                                 #mobility=as.factor(bent_traits2$mobility),
                                 reproductive_mode=as.factor(bent_traits2$reproductive_mode), 
                                 carbonate=as.factor(bent_traits2$carbonate.accretion))
   
   rownames(bent_traits_ord) <- rownames(bent_traits2)
   dim(bent_traits_ord)
   
   ######## Building the gower distance matrix ########
   
   # first calculate gower distance on traits
   
   gower_benthos <- daisy (bent_traits_ord, metric=c("gower")) 
   
   # Building the functional space based on a PCOA 
   
   pco<-dudi.pco(quasieuclid(gower_benthos), scannf=F, nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
   barplot(pco$eig) # barplot of eigenvalues for each axis 
   (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) # percentage of inertia explained by the two first axes = 62%
   
   # Testing the Quality of the functional space
   
   quality<-quality_funct_space_fromdist(gower_benthos,  nbdim=10,   plot="quality_funct_space_I") # The function will test the quality of the space from 2 to 10 axes using dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
   axes_to_choose <- which(quality$meanSD == min(quality$meanSD)) # the minimal value corresponds to the best space to use, here 7 axes 
   (Inertia4 <- (sum(pco$eig[1:axes_to_choose])) /(sum(pco$eig))) # percentage of inertia explained by the 6 first axes = 59%
   
   # Functional diversity indices for benthos #######
   
   # Adjust benthos abundance data to relative abundance
   
   bent_abund <- bent2[,which(colnames(bent2) %in% rownames(bent_traits_ord))] # match fish spp in both datasets
   rownames(bent_abund) <- bent2$locality_site
   
   bent_rel <- bent_abund/rowSums(bent_abund) # relative abundance for fish
   bent_rel <- bent_rel [ which(is.na (rowSums (bent_rel))!= T),]
   # Method with dbFD() 
   
   fd_benthos <- dbFD (gower_benthos, 
                       bent_rel, 
                       corr ="none", 
                       m =  axes_to_choose, 
                       stand.FRic = TRUE,
                       print.pco = TRUE)
   
   bent_fd <- data.frame (nbsp.bent=fd_benthos$nbsp, 
                          sing.sp.bent=fd_benthos$sing.sp, 
                          FRic.bent=fd_benthos$FRic, 
                          qual.FRic.bent=fd_benthos$qual.FRic,
                          FEve.bent=fd_benthos$FEve, 
                          FDiv.bent=fd_benthos$FDiv, 
                          FDis.bent=fd_benthos$FDis, 
                          RaoQ.bent=fd_benthos$RaoQ)

   # and then a list
   results <- list (Fdindexes = bent_fd,
                    chosenAxes = axes_to_choose,
                    InertiaPCO=Inertia2,
                    InertiaQuality=Inertia4,
                    convexHullVolum=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$vol,
                    convexHullArea=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$area,
                    axesPCO=pco$li)
   ; ## return
   
   results
   
   })

## function 2
FD_results_f2_bentos <- lapply (list_benthos_level, function (i) {
   
   ## matching these data   
   bent_traits2 <- merge(bent_traits, benthos_spp, by="groups") # merge traits with groups in occurrence data
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
                                 #modularity=as.factor(bent_traits2$modularity), 
                                 #mobility=as.factor(bent_traits2$mobility),
                                 reproductive_mode=as.factor(bent_traits2$reproductive_mode), 
                                 carbonate=as.factor(bent_traits2$carbonate.accretion))
   
   rownames(bent_traits_ord) <- rownames(bent_traits2)
   dim(bent_traits_ord)
   
   ######## Building the gower distance matrix ########
   
   # first calculate gower distance on traits
   
   gower_benthos <- daisy (bent_traits_ord, metric=c("gower")) 
   
   ## subsetting the gower matrix 
   bent_traits_subset <- bent_traits [which (bent_traits$trophic_type %in% i),"groups"]
   ## 
   gower_matrix_subset <- as.dist(as.matrix(gower_benthos) [which(attr (gower_benthos,"Labels") %in% bent_traits_subset ),
                                                           which(attr (gower_benthos,"Labels") %in% bent_traits_subset)])
   
   
   # Building the functional space based on a PCOA 
   
   pco<-dudi.pco(quasieuclid(gower_matrix_subset), scannf=F, nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
   barplot(pco$eig) # barplot of eigenvalues for each axis 
   (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) # percentage of inertia explained by the two first axes = 62%
   
   # Testing the Quality of the functional space
   
   quality<-quality_funct_space_fromdist(gower_benthos,  nbdim=10,   plot="quality_funct_space_I") # The function will test the quality of the space from 2 to 10 axes using dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
   axes_to_choose <- which(quality$meanSD == min(quality$meanSD)) # the minimal value corresponds to the best space to use, here 7 axes 
   (Inertia4 <- (sum(pco$eig[1:axes_to_choose])) /(sum(pco$eig))) # percentage of inertia explained by the 6 first axes = 59%
   
   # Functional diversity indices for benthos #######
   
   # Adjust benthos abundance data to relative abundance
   
   bent_abund <- bent2[,which(colnames(bent2) %in% rownames(pco$tab))] # match fish spp in both datasets
   rownames(bent_abund) <- bent2$locality_site
   bent_rel <- bent_abund/rowSums(bent_abund) # relative abundance for fish
   quais_sitios_manter <-which(is.na (rowSums (bent_rel))!= T) 
   bent_rel <- bent_rel [quais_sitios_manter,]
   # Method with dbFD() 
   
   fd_benthos <- dbFD (gower_matrix_subset , 
                       bent_rel, 
                       corr ="none", 
                       m = axes_to_choose, 
                       stand.FRic = TRUE,
                       print.pco = TRUE)
   
   bent_fd <- data.frame (nbsp.bent=fd_benthos$nbsp, 
                          sing.sp.bent=fd_benthos$sing.sp, 
                          FRic.bent=fd_benthos$FRic, 
                          qual.FRic.bent=fd_benthos$qual.FRic,
                          FEve.bent=fd_benthos$FEve, 
                          FDiv.bent=fd_benthos$FDiv, 
                          FDis.bent=fd_benthos$FDis, 
                          RaoQ.bent=fd_benthos$RaoQ)
   
   # and then a list
   results <- list (Fdindexes = bent_fd,
                    chosenAxes = axes_to_choose,
                    InertiaPCO=Inertia2,
                    InertiaQuality=Inertia4,
                    convexHullVolum=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$vol,
                    convexHullArea=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$area,
                    axesPCO=pco$li,
                    quais_sitios_manter=quais_sitios_manter)
   ; ## return
   
   results
   
})

## compare the volum of the two functions
## check area and volum (function1)
round(sapply (FD_results_f1_bentos, "[[","convexHullArea"),6)
round(sapply (FD_results_f1_bentos, "[[","convexHullVolum"),6)

## check area and volum (function2)
round(sapply (FD_results_f2_bentos, "[[","convexHullArea"),6)
round(sapply (FD_results_f2_bentos, "[[","convexHullVolum"),6)

## check explained variation (function1)
round(sapply (FD_results_f1_bentos, "[[","InertiaQuality"),6)
round(sapply (FD_results_f1_bentos, "[[","InertiaPCO"),6)

## check explained variation (function2)
round(sapply (FD_results_f2_bentos, "[[","InertiaQuality"),6)
round(sapply (FD_results_f2_bentos, "[[","InertiaPCO"),6)

# write.table (bent_fd, "FIndex_bent.csv", sep=";", dec=".")

################ Plots

## relationship between richness and effort
# benthos
par(mfrow=c(1,2),mar = c (6,4,8,1))
plot(covariates_effort$effort$n_videos_bentos,FD_results_f1_bentos[[1]]$Fdindexes$nbsp.bent,
     ylab= "Number of species",
     xlab= "Number of videos per site",
     main = "Benthos",pch=19,col="coral")
abline(lm (FD_results_f1_bentos[[1]]$Fdindexes$nbsp.bent~covariates_effort$effort$n_videos_bentos),
       lwd=2,col="gray70")
# summary(lm (FD_results_f1_bentos$all$Fdindexes$nbsp.bent~covariates_effort$effort$n_videos_bentos))
text (x=17.5,y=17, labels=expression (paste("R"^2, "=-0.021")),cex=0.7)

# fishes
plot(covariates_effort$effort$n_transectos_peixes,FD_results_f1[[1]]$Fdindexes$nbsp,
     xlab= "Number of transects per site",
     ylab= "",
     main = "Fishes",pch=19,col="cyan3")
abline(lm (FD_results_f1[[1]]$Fdindexes$nbsp~covariates_effort$effort$n_transectos_peixes),
       lwd=2,col="gray70")
# summary(lm (FD_results_f1$all$Fdindexes$nbsp~covariates_effort$effort$n_transectos_peixes))
text (x=80,y=65, labels=expression (paste("R"^2, "=0.15**")),cex=0.7)

### rarefaction for richness of fishes and benthos
## 286 was sample size suggested by the function
# plot(FD_results_f2[[1]]$Fdindexes$nbsp,rarefy (fish,286))

### RELATIONSHIP BETWEEN FRIC OF BENTHOS AND FISHES
## FRIC

whole <- lapply (seq (1,4), function (k)

       (data.frame (FRicBenthos = FD_results_f2_bentos[[1]]$Fdindexes$FRic,
            FRicFishes = FD_results_f2[[k]]$Fdindexes$FRic,
            Group= names(FD_results_f2)[k],
            Benthos = "Whole"))
)

algae.filt <- lapply (seq (1,4), function (k)
      
       (data.frame (FRicBenthos = FD_results_f2_bentos[[2]]$Fdindexes$FRic,
                        FRicFishes = FD_results_f2[[k]]$Fdindexes$FRic,
                        Group= names(FD_results_f2)[k],
                        Benthos = "AlgaeFilterers"))
   )
   
corals <- lapply (seq (1,4), function (k)
      
       (data.frame (FRicBenthos = FD_results_f2_bentos[[3]]$Fdindexes$FRic,
                        FRicFishes = FD_results_f2[[k]]$Fdindexes$FRic [FD_results_f2_bentos[[3]]$quais_sitios_manter],
                        Group= names(FD_results_f2)[k],
                        Benthos="Corals"))
   )

## data for relationship between indexes
data_corr_index <- do.call (rbind, c( whole, algae.filt,corals))
## add comparison between benthic components, and fish diet
data_corr_index$Comp <- paste (data_corr_index$Group, data_corr_index$Benthos,sep=".")
   
## 
ggplot (data_corr_index, aes (x=FRicBenthos,y=FRicFishes)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Benthos+Group, scales = "free_y") + 
   theme_light()

## fRIC 

par (mfrow=c(2,2),mar=c(4,4,2,1))

plot(FD_results_f2_bentos$all$Fdindexes$nbsp.bent, 
     FD_results_f2_bentos$all$Fdindexes$FRic.bent, 
     ylim=c(0,0.83), xlim=c(0,20), pch=21, col="coral", bg="coral",
     ylab="Benthos FRic", xlab="Species Richness", cex=1.2)
plot(FD_results_f2[[1]]$Fdindexes$nbsp, 
     FD_results_f2[[1]]$Fdindexes$FRic, ylim=c(0,0.4), xlim=c(10,70), pch=21, col="cyan3", bg="cyan3",
     ylab="Fish FRic", xlab="Species Richness", cex=1.2)
plot(FD_results_f2_bentos$all$Fdindexes$FRic.bent, 
     FD_results_f2[[1]]$Fdindexes$FRic, ylim=c(0,0.3), xlim=c(0,0.83), pch=21, col="goldenrod", bg="goldenrod",
     ylab="Fish FRic", xlab="Benthos FRic", cex=1.2)

par (mfrow=c(2,2),mar=c(4,4,2,1))
plot(FD_results_f2_bentos$all$Fdindexes$FEve.bent~FD_results_f2_bentos$all$Fdindexes$nbsp.bent, ylim=c(0,0.8), xlim=c(0,30), pch=21, col="darkgray", bg="darkgray",
     ylab="Benthos FEveness", xlab="Benthos Richness", cex=1.2)
plot(FD_results_f2_bentos$all$Fdindexes$FDiv.bent~FD_results_f2_bentos$all$Fdindexes$nbsp.bent, ylim=c(0,1), xlim=c(0,30), pch=21, col="darkgray", bg="darkgray",
     ylab="Benthos FDivergence", xlab="Benthos Richness", cex=1.2)
plot(FD_results_f2[[1]]$Fdindexes$FEve~rarefy(fish,286), ylim=c(0,0.8), xlim=c(0,70), pch=21, col="darkgray", bg="darkgray",
     ylab="Fish FEveness", xlab="Fish Richness", cex=1.2)
plot(FD_results_f2[[1]]$Fdindexes$FDiv~rarefy(fish,286), ylim=c(0,1), xlim=c(0,70), pch=21, col="darkgray", bg="darkgray",
     ylab="Fish FDivergence", xlab="Fish Richness", cex=1.2)


### functional volume
par (mfrow=c(2,2),mar=c(4,4,2,1))
plot(FD_results_f2_bentos$all$axesPCO [,1:2],xlim =c(-0.6,1))
Plot_ConvexHull (FD_results_f2_bentos$all$axesPCO [,1],
                 FD_results_f2_bentos$all$axesPCO [,2],
                 lcolor="black")
Plot_ConvexHull (FD_results_f2_bentos$autot$axesPCO [,1],
                 FD_results_f2_bentos$autot$axesPCO [,2],
                 lcolor="green")
Plot_ConvexHull (FD_results_f2_bentos$mixo$axesPCO [,1],
                 FD_results_f2_bentos$mixo$axesPCO [,2],
                 lcolor="goldenrod")
## fishes
plot(FD_results_f1$all$axesPCO [,1:2],xlim =c(-0.6,1))
Plot_ConvexHull (FD_results_f2$all$axesPCO [,1],
                 FD_results_f2$all$axesPCO [,2],
                 lcolor=rgb(red = 0.5, green = 0.8, blue = 1, alpha = 1))
Plot_ConvexHull (FD_results_f2$herbivorous$axesPCO [,1],
                 FD_results_f2$herbivorous$axesPCO [,2],
                 lcolor=rgb(red = 0, green = 1, blue = 0, alpha = 1))
Plot_ConvexHull (FD_results_f2$sessileInv$axesPCO [,1],
                 FD_results_f2$sessileInv$axesPCO [,2],
                 lcolor= rgb(red = 0, green = 1, blue = 1, alpha = 0.5))
Plot_ConvexHull (FD_results_f2$nonHS$axesPCO [,1],
                 FD_results_f2$nonHS$axesPCO [,2],
                 lcolor=rgb(red = 0, green =0, blue = 1, alpha = 0.5))

## Load environmental variables

## Plots with environmental data
## building a DF with data

df_results <- data.frame (Nspec = FD_results_f2$all$Fdindexes$nbsp,
                           FRic = FD_results_f2$all$Fdindexes$FRic,
                          FEve = FD_results_f2$all$Fdindexes$FEve,
                          FDiv = FD_results_f2$all$Fdindexes$FDiv,
                          Group="Fishes")
# bind dataframe with benthos results
df_results <- rbind (df_results, 
                     data.frame (Nspec = FD_results_f2_bentos$all$Fdindexes$nbsp.bent,
                          FRic = FD_results_f2_bentos$all$Fdindexes$FRic.bent,
                          FEve = FD_results_f2_bentos$all$Fdindexes$FEve.bent,
                          FDiv = FD_results_f2_bentos$all$Fdindexes$FDiv.bent,
                          Group="Benthos"))
## melt these results
df_results <- melt(df_results)
## and bind covariates
df_results <- cbind(df_results,
         ReefType=covariates_site$biog_reef,
         Region=covariates_site$region,
         covariates_site$coord,
         covariates_site$sea_data)
# warning due to rownames of covariates_site$sea_data

### Latitude
ggplot (df_results, aes (x=Lat, y=value)) + 
   geom_point() + 
   geom_smooth() +
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Latitude") + 
   ylab("Index")

### Temperature
ggplot (df_results, aes (x=BO2_tempmean_ss, y=value)) + 
   geom_point() + 
   geom_smooth() +
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Average sea surface temperature (ºC)") + 
   ylab("Index")

### temperature range
ggplot (df_results, aes (x=BO2_temprange_ss, y=value)) + 
   geom_point() + 
   geom_smooth() +
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Range of sea surface temperature (ºC)") + 
   ylab("Index")

### primary productivity
ggplot (df_results, aes (x=BO2_ppmean_ss, y=value)) + 
   geom_point() + 
   geom_smooth() +
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Average primary productivity") + 
   ylab("Index")

### primary productivity range
ggplot (df_results, aes (x=BO2_pprange_ss, y=value)) + 
   geom_point() + 
   geom_smooth() +
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Primary productivity range") + 
   ylab("Index")

### Salinity
ggplot (df_results, aes (x=BO2_salinitymean_ss, y=value)) + 
   geom_point() + 
   geom_smooth() +
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Average salinity") + 
   ylab("Index")

# Salinity range
ggplot (df_results, aes (x=BO2_salinityrange_ss, y=value)) + 
geom_point() + 
   geom_smooth() +
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Salinity range") + 
   ylab("Index")

# ### Chlorophyl
ggplot (df_results, aes (x=BO2_chlomean_ss, y=value)) + 
   geom_point() + 
   geom_smooth() +
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Average chlorophyl") + 
   ylab("Index")

# Salinity range
ggplot (df_results, aes (x=BO2_chlorange_ss, y=value)) + 
   geom_point() + 
   geom_smooth() +
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Chlorophyl range") + 
   ylab("Index")

## region
ggplot (df_results, aes (x=Region, y=value)) + 
   geom_boxplot()+
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Region") + 
   ylab("Index")

## Reef type
ggplot (df_results, aes (x=ReefType, y=value)) + 
   geom_boxplot()+
   facet_wrap(~Group+variable,scales="free_y",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Reef type") + 
   ylab("Index")

