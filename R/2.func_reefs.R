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
                                 #Diet=as.factor(fish_traits2$Diet), 
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
                    axesPCO=pco$l1)
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
                                 #Diet=as.factor(fish_traits2$Diet), 
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
                    axesPCO=pco$l1)
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
list_benthos_level <- list (all =c("AC","S","C","G", "A"),
                            nonAut = c("AC","S","C","G"),
                            nonMixAut = c("S","C","G"),
                            corals = c("AC"))

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
                                 modularity=as.factor(bent_traits2$modularity), 
                                 #mobility=as.factor(bent_traits2$mobility),
                                 #reproductive_mode=as.factor(bent_traits2$reproductive_mode), 
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
   quais_sitios_manter <-which(is.na (rowSums (bent_rel))!= T) 
   bent_rel <- bent_rel [quais_sitios_manter,]
   # Method with dbFD() 
   
   fd_benthos <- dbFD (gower_benthos, 
                       bent_rel, 
                       corr ="none", 
                       m =  axes_to_choose, 
                       stand.FRic = TRUE,
                       print.pco = TRUE)
   
   ## some sites had no species
   
   bent_fd <- data.frame (matrix (NA, ncol=8,nrow=38,
           dimnames = list (NULL,
                            c("nbsp.bent",
                              "sing.sp.bent",
                              "FRic.bent",
                              "qual.FRic.bent",
                              "FEve.bent",
                              "FDiv.bent",
                              "FDis.bent",
                              "RaoQ.bent"))))
   
   bent_fd[quais_sitios_manter,1] <- fd_benthos$nbsp
   bent_fd[quais_sitios_manter,2] <- fd_benthos$sing.sp
   bent_fd[quais_sitios_manter,3] <- fd_benthos$FRic
   bent_fd[quais_sitios_manter,4] <- fd_benthos$qual.FRic
   bent_fd[quais_sitios_manter,5] <- fd_benthos$FEve
   bent_fd[quais_sitios_manter,6] <- fd_benthos$FDiv
   bent_fd[quais_sitios_manter,7] <- fd_benthos$FDis
   bent_fd[quais_sitios_manter,8] <- fd_benthos$RaoQ
   
   # and then a list
   results <- list (Fdindexes = bent_fd,
                    chosenAxes = axes_to_choose,
                    InertiaPCO=Inertia2,
                    InertiaQuality=Inertia4,
                    convexHullVolum=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$vol,
                    convexHullArea=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$area,
                    axesPCO=pco$l1)
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
                                 modularity=as.factor(bent_traits2$modularity), 
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
   
   ## some sites had no species
   bent_fd <- data.frame (matrix (NA, ncol=8,nrow=38,
                                  dimnames = list (NULL,
                                                   c("nbsp.bent",
                                                     "sing.sp.bent",
                                                     "FRic.bent",
                                                     "qual.FRic.bent",
                                                     "FEve.bent",
                                                     "FDiv.bent",
                                                     "FDis.bent",
                                                     "RaoQ.bent"))))
   bent_fd[quais_sitios_manter,1] <- fd_benthos$nbsp
   bent_fd[quais_sitios_manter,2] <- fd_benthos$sing.sp
   bent_fd[quais_sitios_manter,3] <- fd_benthos$FRic
   bent_fd[quais_sitios_manter,4] <- fd_benthos$qual.FRic
   bent_fd[quais_sitios_manter,5] <- fd_benthos$FEve
   bent_fd[quais_sitios_manter,6] <- fd_benthos$FDiv
   bent_fd[quais_sitios_manter,7] <- fd_benthos$FDis
   bent_fd[quais_sitios_manter,8] <- fd_benthos$RaoQ
   
   # and then a list
   results <- list (Fdindexes =  bent_fd,
                    chosenAxes = axes_to_choose,
                    InertiaPCO=Inertia2,
                    InertiaQuality=Inertia4,
                    convexHullVolum=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$vol,
                    convexHullArea=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$area,
                    axesPCO=pco$l1)
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
par(mfrow=c(2,2),mar = c (4,4,2,1))
plot(covariates_effort$effort$n_videos_bentos,FD_results_f1_bentos[[1]]$Fdindexes$nbsp.bent,
     ylab= "Number of species",
     xlab= "Number of videos per site",
     main = "Benthos",pch=19,col="coral",
     cex.axis=0.6,cex.lab=0.8)
abline(lm (FD_results_f1_bentos[[1]]$Fdindexes$nbsp.bent~covariates_effort$effort$n_videos_bentos),
       lwd=2,col="gray70")
# summary(lm (FD_results_f1_bentos$all$Fdindexes$nbsp.bent~covariates_effort$effort$n_videos_bentos))
text (x=17.5,y=17, labels=expression (paste("R"^2, "=-0.021")),cex=0.7)

# fishes
plot(covariates_effort$effort$n_transectos_peixes,FD_results_f1[[1]]$Fdindexes$nbsp,
     xlab= "Number of transects per site",
     ylab= "",
     main = "Fishes",pch=19,col="cyan3",
     cex.axis=0.6,cex.lab=0.8)
abline(lm (FD_results_f1[[1]]$Fdindexes$nbsp~covariates_effort$effort$n_transectos_peixes),
       lwd=2,col="gray70")
# summary(lm (FD_results_f1$all$Fdindexes$nbsp~covariates_effort$effort$n_transectos_peixes))
text (x=80,y=65, labels=expression (paste("R"^2, "=0.15**")),cex=0.7)

### rarefaction for richness of fishes and benthos
## 286 was sample size suggested by the function; lower than the minimum row max across sites (max abundance)
rarefied_richness <- rarefy (fish,100)
plot(covariates_effort$effort$n_transectos_peixes,rarefied_richness,
     ylab= "Rarefied fish richness",
     xlab = "Number of transects per site",
     cex.axis=0.6,cex.lab=0.8,
     pch=19, col = "blue")
abline (lm(rarefied_richness ~ covariates_effort$effort$n_transectos_peixes),
        lwd=2,col="gray70")
## summary (lm(rarefy (fish,286) ~ covariates_effort$effort$n_transectos_peixes))
text (x=80,y=28, labels=expression (paste("R"^2, "=-0.006")),cex=0.7)

### RELATIONSHIP BETWEEN FRIC OF BENTHOS AND FISHES
## FRIC

whole_data <- lapply (seq (1,length(FD_results_f2)), function (k)
   do.call (rbind, lapply (seq (1,length(FD_results_f2_bentos)), function (i)

       data.frame (
            RicBenthos = FD_results_f2_bentos[[i]]$Fdindexes$nbsp.bent,
            RicFishes = FD_results_f2[[k]]$Fdindexes$nbsp,
            RareRicFishes = rarefied_richness,
            FRicBenthos = FD_results_f2_bentos[[i]]$Fdindexes$FRic,
            FRicFishes = FD_results_f2[[k]]$Fdindexes$FRic,
            FEveBenthos = FD_results_f2_bentos[[i]]$Fdindexes$FEve.bent,
            FEveFishes = FD_results_f2[[k]]$Fdindexes$FEve,
            FDivBenthos = FD_results_f2_bentos[[i]]$Fdindexes$FDiv.bent,
            FDivFishes = FD_results_f2[[k]]$Fdindexes$FDiv,
            Group= names(FD_results_f2)[k],
            Benthos = names(list_benthos_level)[i])
)))

## data for relationship between indexes
data_corr_index <- do.call (rbind, whole_data)

## Relationship between FRic of Benthos and Fishes
ggplot (data_corr_index, aes (x=FRicBenthos,y=FRicFishes)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Benthos+Group, scales = "free") + 
   theme_light()

## Relationship between FRic and richness of Benthos and Fishes
# Rich benthos, FRic Fishes
ggplot (data_corr_index, aes (x=RicBenthos,y=FRicFishes)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Benthos+Group, scales = "free") + 
   theme_light()
# Rich benthos, rich fishes
ggplot (data_corr_index, aes (x=RicBenthos,y=RicFishes)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Benthos+Group, scales = "free") + 
   theme_light()
# FRic Fishes,  FRic benthos
ggplot (data_corr_index, aes (x=FRicBenthos,y=FRicFishes)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Benthos+Group, scales = "free") + 
   theme_light()
# FEve benthos,  Feve fishes
ggplot (data_corr_index, aes (x=FEveBenthos,y=FEveFishes)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Benthos+Group, scales = "free") + 
   theme_light()
# FDiv benthos,  FRic fishes
ggplot (data_corr_index, aes (x=FDivBenthos,y=FRicFishes)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Benthos+Group, scales = "free") + 
   theme_light()
# FRic benthos,  Feve benthos
ggplot (data_corr_index, aes (x=FRicBenthos,y=FEveBenthos)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Benthos, scales = "free") + 
   theme_light()
# FRic benthos,  FDiv benthos
ggplot (data_corr_index, aes (x=FRicBenthos,y=FDivBenthos)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Benthos, scales = "free") + 
   theme_light()
# FRic fishes,  Feve fishes
ggplot (data_corr_index, aes (x=FRicFishes,y=FEveFishes)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Group, scales = "free") + 
   theme_light()
# FRic fishes,  FDiv fishes
ggplot (data_corr_index, aes (x=FRicFishes,y=FDivFishes)) + 
   geom_point() + 
   geom_smooth()+
   facet_wrap(~ Group, scales = "free") + 
   theme_light()

####################################################################
### functional volume
par (mfrow=c(2,2),mar=c(4,4,2,1))
plot(FD_results_f2_bentos$all$axesPCO [,1:2],xlim =c(-5,5),ylim =c(-5,5))
Plot_ConvexHull (FD_results_f2_bentos$all$axesPCO [,1],
                 FD_results_f2_bentos$all$axesPCO [,2],
                 lcolor="black")
Plot_ConvexHull (FD_results_f2_bentos[[2]]$axesPCO [,1],
                 FD_results_f2_bentos[[2]]$axesPCO [,2],
                 lcolor="green")
Plot_ConvexHull (FD_results_f2_bentos[[3]]$axesPCO [,1],
                 FD_results_f2_bentos[[3]]$axesPCO [,2],
                 lcolor="goldenrod")
Plot_ConvexHull (FD_results_f2_bentos[[4]]$axesPCO [,1],
                 FD_results_f2_bentos[[4]]$axesPCO [,2],
                 lcolor="blue")

## fishes
plot(FD_results_f2$all$axesPCO [,1:2],xlim =c(-5,5),ylim =c(-5,5))
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

df_results <- data.frame (Nspec = rarefied_richness,
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
   geom_smooth(method='gam',formula = y ~ s(x, bs = "cs")) +
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Latitude") + 
   ylab("Index")

### Temperature
ggplot (df_results, aes (x=BO2_tempmean_ss, y=value)) + 
   geom_point() + 
   geom_smooth(method='gam',formula = y ~ s(x, bs = "cs")) +
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Average sea surface temperature (ºC)") + 
   ylab("Index")

### temperature range
ggplot (df_results, aes (x=BO2_temprange_ss, y=value)) + 
   geom_point() + 
   geom_smooth(method='gam',formula = y ~ s(x, bs = "cs")) +
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Range of sea surface temperature (ºC)") + 
   ylab("Index")

### primary productivity
ggplot (df_results, aes (x=BO2_ppmean_ss, y=value)) + 
   geom_point() + 
   geom_smooth(method='gam',formula = y ~ s(x, bs = "cs")) +
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Average primary productivity") + 
   ylab("Index")

### primary productivity range
ggplot (df_results, aes (x=BO2_pprange_ss, y=value)) + 
   geom_point() + 
   geom_smooth(method='gam',formula = y ~ s(x, bs = "cs")) +
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Primary productivity range") + 
   ylab("Index")

### Salinity
ggplot (df_results, aes (x=BO2_salinitymean_ss, y=value)) + 
   geom_point() + 
   geom_smooth(method='gam',formula = y ~ s(x, bs = "cs")) +
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Average salinity") + 
   ylab("Index")

# Salinity range
ggplot (df_results, aes (x=BO2_salinityrange_ss, y=value)) + 
geom_point() + 
   geom_smooth(method='gam',formula = y ~ s(x, bs = "cs")) +
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Salinity range") + 
   ylab("Index")

# ### Chlorophyl
ggplot (df_results, aes (x=BO2_chlomean_ss, y=value)) + 
   geom_point() + 
   geom_smooth(method='gam',formula = y ~ s(x, bs = "cs")) +
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Average chlorophyl") + 
   ylab("Index")

# Salinity range
ggplot (df_results, aes (x=BO2_chlorange_ss, y=value)) + 
   geom_point() + 
   geom_smooth(method='gam',formula = y ~ s(x, bs = "cs")) +
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Chlorophyl range") + 
   ylab("Index")

## region
ggplot (df_results, aes (x=Region, y=value)) + 
   geom_boxplot()+
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Region") + 
   ylab("Index")

## Reef type
ggplot (df_results, aes (x=ReefType, y=value)) + 
   geom_boxplot()+
   facet_wrap(~Group+variable,scales="free",
              ncol=4,nrow=2) + 
   theme_classic() + 
   xlab("Reef type") + 
   ylab("Index")

############################ 
## MAPPING

## creating a dataframe for each result
# fish
fish_richness <- do.call (cbind, 
   lapply (seq(1,length(FD_results_f2)), function (i)

      FD_results_f2[[i]]$Fdindexes$nbsp
      
))
colnames (fish_richness) <- names(FD_results_f2)

# benthos
benthos_richness <- do.call (cbind, 
                             
    lapply (seq(1,length(FD_results_f2_bentos)), function (i)
            
            FD_results_f2_bentos[[i]]$Fdindexes$nbsp.bent
            
         ))
colnames (benthos_richness) <- names(FD_results_f2_bentos)
benthos_richness <- data.frame(benthos_richness)
benthos_richness$aut <- benthos_richness$all - benthos_richness$nonAut 
benthos_richness[is.na(benthos_richness)] <- 0

## binding coordinates (fishes for now)
fish_richness <- cbind (fish_richness, covariates_site$coord)
## jittering these coordinates
fish_richness$LonJitter <- jitter (fish_richness$Lon,factor=400)
fish_richness$LatJitter <- jitter (fish_richness$Lat,factor=600)

## just moving charts of benthic communities to the left
benthos_richness <- cbind(benthos_richness,LatJitter = fish_richness$LatJitter,
                          LonJitter=fish_richness$LonJitter+6.5)

## maps
# mapa mundi
world <- ne_countries(scale = "medium", returnclass = "sf")

# cortar o mapa para ver a america do Sul e parte da central
wm <- ggplot() + 
   geom_sf (data=world, size = 0.1, 
            fill= "gray90",colour="gray90") +
   coord_sf (xlim = c(-50, -20),  ylim = c(-30, 2), expand = FALSE) +
   theme_bw() + #xlab ("Longitude")  + ylab ("Latitude") +
   theme(panel.border = element_blank(), 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "lightcyan",#darkslategray1
                                         colour = "lightcyan"),
         axis.text.x = element_text(size=6),
         axis.ticks.x=element_line(size=1),
         axis.text.y = element_text(size=6),
         axis.ticks.y=element_line(size=1),
         axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         title = element_blank()) + ggtitle("Number of fish and benthic species")+
   xlab("Longitude") + ylab("Latitude")


## advise to jitter : https://stackoverflow.com/questions/52806580/pie-charts-in-geom-scatterpie-overlapping
## pie: http://www.spectdata.com/index.php/2018/10/25/how-to-use-ggplot-to-plot-pie-charts-on-a-map/

wm_pie <- wm + geom_scatterpie(aes(x=LonJitter, y=LatJitter, r=all/max(all)),alpha=0.5,
                               data = fish_richness,
                               cols = c(
                                        "herbivorous",
                                        "sessileInv",
                                        "nonHS"),
                               #pie_scale = 0.1,
                               size=0,
                               sorted_by_radius = F,
                               legend_name = "Groups") 

### adding benthos richness

wm_pie_col_benthos <- wm_pie + geom_scatterpie(aes(x=LonJitter, y=LatJitter,r=all/max(all,na.rm=T)),alpha=0.5,
                data = benthos_richness,
                cols = c("aut",
                         "nonAut",
                         "nonMixAut",
                         "corals"
                   ),
                pie_scale = 1,
                size=0,
                sorted_by_radius = F,
                legend_name = "Benthos") + 
   
   theme (legend.title = element_text(size=7),
          legend.text = element_text(size=7),
          legend.position = c(0.37, 0.9),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6,6,6,6),
          legend.background = element_blank(),
          title=element_text(size=7)) 


wm_pie_col_benthos <- wm_pie_col_benthos + scale_fill_manual(
        
                                       breaks=c("aut", "nonAut", "nonMixAut","corals",
                                                "nonHS","sessileInv","herbivorous"),
                                       labels=c("Autotrophs", "Non-autotrophs", 
                                                "Non-autotrophs and non-mixotrophs",
                                                "Mixotrophs",
                                                "Carnivores, planktivores, omnivores",
                                                "Invertivores",
                                                "Herbivores"),
                                       values= c("aut" = "#003300",
                                                 "nonAut" = "#009900",
                                                 "nonMixAut" = "#00CC66",
                                                 "corals" = "#B2FF66",
                                                 "nonHS" = "#990000",
                                                 "sessileInv" = "#CC0000",
                                                 "herbivorous" = "#FF6666"))

wm_pie_col_benthos <- wm_pie_col_benthos + geom_scatterpie_legend(fish_richness$all/max(fish_richness$all), 
                                                  x=-35, y=0.5, n=3, 
                                                  labeller=function(x) (x*max(fish_richness$all)))
wm_pie_col_benthos <- wm_pie_col_benthos + geom_scatterpie_legend(benthos_richness$all/max(benthos_richness$all,na.rm=T), 
                                                  x=-28, y=0.5, n=3, 
                                                  labeller=function(x) (x*max(benthos_richness$all,na.rm=T)))


### Latitude and FRIC
lat_fric <- ggplot (df_results [which (df_results$variable == "FRic"),], aes (x=Lat, y=value,
                                                                  col=Group)) + 
   geom_point() + 
   geom_smooth() +
   scale_color_manual(values=c('#009900','#ca0020'))+
   theme_classic() + 
   ylab(NULL) + 
   coord_flip() + 
   theme_classic()+
   ggtitle ("FRic")+
   
   theme (legend.title = element_text(size=6),
          legend.text = element_text(size=5),
          legend.position = c(1, 0.8),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6,6,6,6),
          legend.background = element_blank(),
          axis.title.x = element_text(size=8),
          axis.text.x = element_text(size=6),
          axis.title.y =  element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          title=element_text(size=7),
          panel.background = element_rect(fill = 'lightcyan', 
                                          colour = 'lightcyan'),
          plot.margin = unit(c(0.2,0.2,0.3,-0.5), "cm")) 

## LATITUDE AND FEVE
lat_feve <- ggplot (df_results [which (df_results$variable == "FEve"),], aes (x=Lat, y=value,
                                                                              col=Group)) + 
   geom_point() + 
   geom_smooth() +
   scale_color_manual(values=c('#009900','#ca0020'))+
   theme_classic() + 
   ylab(NULL) + 
   coord_flip() + 
   theme_classic()+
   ggtitle ("FEve")+
   
   theme (legend.position = "none",
          axis.title.x = element_text(size=6),
          axis.text.x = element_text(size=5),
          axis.title.y =  element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          title=element_text(size=7),
          panel.background = element_rect(fill = 'lightcyan', 
                                          colour = 'lightcyan'),
          plot.margin = unit(c(0.2,0.2,0.3,-0.1), "cm")) 


## LATITUDE AND FDIV
lat_fdiv <- ggplot (df_results [which (df_results$variable == "FDiv"),], aes (x=Lat, y=value,
                                                                              col=Group)) + 
   geom_point() + 
   geom_smooth() +
   scale_color_manual(values=c('#009900','#ca0020'))+
   theme_classic() + 
   ylab(NULL) + 
   coord_flip() + 
   theme_classic()+
   ggtitle ("FDiv")+
   
   theme (legend.position = "none",
          axis.title.x = element_text(size=6),
          axis.text.x = element_text(size=5),
          axis.title.y =  element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          title=element_text(size=7),
          panel.background = element_rect(fill = 'lightcyan', 
                                          colour = 'lightcyan'),
          plot.margin = unit(c(0.2,0.2,0.3,-0.1), "cm")) 


## panel with plot and map
pdf(here ("output", "MapFD.pdf"),width = 10,heigh=7)
grid.arrange(wm_pie_col_benthos, 
             lat_fric,
             lat_feve,
             lat_fdiv,
             ncol=9,nrow=6,
             layout_matrix = rbind (c(1,1,1,1,1,1,2,3,4),
                                    c(1,1,1,1,1,1,2,3,4),
                                    c(1,1,1,1,1,1,2,3,4),
                                    c(1,1,1,1,1,1,2,3,4),
                                    c(1,1,1,1,1,1,2,3,4),
                                    c(1,1,1,1,1,1,2,3,4)))
dev.off()

