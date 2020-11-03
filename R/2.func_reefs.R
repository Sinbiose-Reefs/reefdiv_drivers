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
# editing diet (removing NA and empty cells)
fish_traits <- fish_traits [which(is.na(fish_traits$Diet) ==F),]
fish_traits<-fish_traits [-which(fish_traits$Diet == ""),]

# Select fish species from Atlantic database to match UVCs

fishes <- data.frame(Name=colnames(fish))
# matching considering all species, and only herbivorous and sessile/mob inv

list_fish_diet <- list( all = c("hd","hm","om","pk","fc","im","is"),
                        herbivorous = c("hd","hm"),
                        sessileInv = c("is","im"),
                        om = c("om"),
                        fc = "fc",
                        pk = "pk"
                        )

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
   
   ################################################################
   ## or activate this if you want a subset of gower matrix
   # first calculate gower distance on traits
   #gower_matrix <- daisy (fish_traits_ord, metric=c("gower")) 

   ## subsetting the gower matrix 
   #fish_traits_subset <- fish_traits [which (fish_traits$Diet %in% i),"Name"]
   ## 
   #gower_matrix <- as.dist(as.matrix(gower_matrix) [which(attr (gower_matrix,"Labels") %in% fish_traits_subset),
   #                                                           which(attr (gower_matrix,"Labels") %in% fish_traits_subset)])
   ###################################################################
   
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
   quais_sitios_manter <-which(is.na (rowSums (fish_rel))!= T) 
   fish_rel <- fish_rel [quais_sitios_manter,]
   # Method with dbFD() 
   fd <- dbFD (gower_matrix, fish_rel, 
               corr ="none", 
               m = axes_to_choose, 
               stand.FRic = TRUE,
               print.pco = TRUE,
               calc.FGR = F,
               clust.type = "kmeans")
   
   fish_fd  <- data.frame (matrix (NA, ncol=8,nrow=nrow(fish_abund),
                                  dimnames = list (NULL,
                                                   c("nbsp",
                                                     "sing",
                                                     "FRic",
                                                     "qual.FRic",
                                                     "FEve",
                                                     "FDiv",
                                                     "FDis",
                                                     "RaoQ"))))
   
   fish_fd [quais_sitios_manter,1] <- fd$nbsp
   fish_fd[quais_sitios_manter,2] <- fd$sing.sp
   fish_fd[quais_sitios_manter,3] <- fd$FRic
   fish_fd[quais_sitios_manter,4] <- fd$qual.FRic
   fish_fd[quais_sitios_manter,5] <- fd$FEve
   fish_fd[quais_sitios_manter,6] <- fd$FDiv
   fish_fd[quais_sitios_manter,7] <- fd$FDis
   fish_fd[quais_sitios_manter,8] <- fd$RaoQ
   
   ## FD total
   # Method with dbFD() 
   total <- aggregate (fish_abund,
                       by = list(covariates_site$region[,'Region']),
                       FUN=sum)[,-1]
   total <- total/rowSums(total)
   
   ## calculating total FD           
   fd_total <- dbFD (gower_matrix, total, 
               corr ="none", 
               m = axes_to_choose, 
               stand.FRic = TRUE,
               print.pco = TRUE,
               calc.FGR = F,
               clust.type = "kmeans")
   
   # and then a list
   results <- list (Fdindexes = fish_fd,
                    fd_total = fd_total,
                    chosenAxes = axes_to_choose,
                    InertiaPCO=Inertia2,
                    InertiaQuality=Inertia7,
                    convexHullVolum=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$vol,
                    convexHullArea=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$area,
                    axesPCO=fd$x.axes)
   ; ## return
   
   results
   
   
})

## compare the volum of the two functions
## check area and volum (function1)
round(sapply (FD_results_f1, "[[","convexHullArea"),6)
round(sapply (FD_results_f1, "[[","convexHullVolum"),6)

## check explained variation (function1)
round(sapply (FD_results_f1, "[[","InertiaQuality"),6)
round(sapply (FD_results_f1, "[[","InertiaPCO"),6)

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
## mante coral e alga

list_benthos_level <- list (all =c("AC","S","C","G","A"),
                            aut = "A", # produtores
                            nonAut = c("S","C","G"), ## consumidores
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
   
   #################################################################
   ## or activate this if you want a subset of the distance matrix
   #gower_benthos <- daisy (bent_traits_ord, metric=c("gower")) 
   ## subsetting the gower matrix 
   #bent_traits_subset <- bent_traits [which (bent_traits$trophic_type %in% i),"groups"]
   ## 
   #gower_matrix_ <- as.dist(as.matrix(gower_benthos) [which(attr (gower_benthos,"Labels") %in% bent_traits_subset ),
   #                                                         which(attr (gower_benthos,"Labels") %in% bent_traits_subset)])
   ###################################################################
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
                       print.pco = TRUE,
                       calc.FGR = F,
                       clust.type = "kmeans")
   
   ## some sites had no species
   
   bent_fd <- data.frame (matrix (NA, ncol=8,nrow=nrow(bent_abund),
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
   
   ## FD total
   # Method with dbFD() 
   # Method with dbFD() 
   total <- aggregate (bent_abund,
                       by = list(covariates_site$region[,'Region'][which(covariates_effort$effort$n_videos_bentos>0)]),
                       FUN=sum)[,-1]
   total <- total/rowSums(total)
   
   ## calculating total FD           
   fd_total <- dbFD (gower_benthos, total, 
                     corr ="none", 
                     m = axes_to_choose, 
                     stand.FRic = TRUE,
                     print.pco = TRUE,
                     calc.FGR = F,
                     clust.type = "kmeans")
   
   
   # and then a list
   results <- list (Fdindexes = bent_fd,
                    fd_total = fd_total,
                    chosenAxes = axes_to_choose,
                    InertiaPCO=Inertia2,
                    InertiaQuality=Inertia4,
                    convexHullVolum=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$vol,
                    convexHullArea=convhulln(pco$tab[,1:axes_to_choose],output.options = T)$area,
                    axesPCO=fd_benthos$x.axes)
   ; ## return
   
   results
   
   })

## compare the volum of the two functions
## check area and volum (function1)
round(sapply (FD_results_f1_bentos, "[[","convexHullArea"),6)
round(sapply (FD_results_f1_bentos, "[[","convexHullVolum"),6)

## check explained variation (function1)
round(sapply (FD_results_f1_bentos, "[[","InertiaQuality"),6)
round(sapply (FD_results_f1_bentos, "[[","InertiaPCO"),6)

## finding the number of functional entities
bent_traits_subset <- bent_traits [which(bent_traits$groups %in% colnames (bent2)[-1]),]

total_funct_ent <- unique(
   unlist (
      lapply (seq (1,nrow (bent_traits_subset)), function (i) {
   
   FE <- paste(unname(bent_traits_subset[i,-1] [1]),
               unname(bent_traits_subset[i,-1] [2]),
               unname(bent_traits_subset[i,-1] [3]),
               unname(bent_traits_subset[i,-1] [4]),
               unname(bent_traits_subset[i,-1] [5]),
               unname(bent_traits_subset[i,-1] [6]),
               sep=".")
    }  
  )
 )
)

# site entities
# site composition
list_sites_bentos <- lapply (seq (1,nrow (bent2)), function (i){
   
   sub1 <- bent2[i,-1]
   sub2 <- sub1[,which(sub1>0)] 
   
})

## trait composition per site
traits_site <- lapply (list_sites_bentos, function (i)
   
   bent_traits [which(bent_traits$groups %in% colnames (i)),]
)

### finding site functional entities
site_funct_ent <- lapply (traits_site, function (k)
   unique(
      unlist (lapply (seq (1,nrow (k)), function (i) {
         
         FE <- paste(unname(bent_traits_subset[i,-1] [1]),
                     unname(bent_traits_subset[i,-1] [2]),
                     unname(bent_traits_subset[i,-1] [3]),
                     unname(bent_traits_subset[i,-1] [4]),
                     unname(bent_traits_subset[i,-1] [5]),
                     unname(bent_traits_subset[i,-1] [6]),
                     sep=".")
      }  
      )
   )
 )
)

###############################
### fishes functional entities
fish_traits2 <- merge(fishes, fish_traits, by="Name")

# Select traits of interest and remove NAs
fish_traits2 <- na.omit(fish_traits2[,c("Name","Home_range","Diel_activity","Size_group","Body_size",
                                        "Level_water","Body_shape","Diet","Caudal_fin","Mouth_position")])

# here transform each ordered category in 1, 2, 3.... 
fish_traits2$Body_size  <- sapply(fish_traits2$Body_size, function(x) {
   if (x<=10) {"S"} 
   else if (x>10&x<=50) {"M"} 
   else if (x>50&x<=400) {"L"}
   else if (x>400) {"XL"}
   }
)

## all traits into a dataframe
fish_traits_ord <-data.frame (Size=as.character(fish_traits2$Body_size), 
                              Mobility=as.character(fish_traits2$Home_range), 
                              Activity= as.character(fish_traits2$Diel_activity), 
                              Schooling=as.character(fish_traits2$Size_group), 
                              Level=as.character(fish_traits2$Level_water), 
                              #Diet=as.factor(fish_traits2$Diet), 
                              Body_shape=as.character(fish_traits2$Body_shape),
                              Caudal_fin=as.character(fish_traits2$Caudal_fin), 
                              Mouth_position=as.character(fish_traits2$Mouth_position))

rownames(fish_traits_ord) <- (fish_traits2$Name)

## finding total FE
total_funct_ent_fish <- unique(
   unlist ( lapply (seq(1,nrow(fish_traits_ord)), function (i) {

      FE <- paste(unname(fish_traits_ord[i,] [1]),
            unname(fish_traits_ord[i,] [2]),
            unname(fish_traits_ord[i,] [3]),
            unname(fish_traits_ord[i,] [4]),
            unname(fish_traits_ord[i,] [5]),
            unname(fish_traits_ord[i,] [6]),
            unname(fish_traits_ord[i,] [7]),
            unname(fish_traits_ord[i,] [8]),
            sep=".")
}
)))

# site entities
# site composition
list_sites_fishes <- lapply (seq (1,nrow (fish)), function (i){
   
   sub1 <- fish[i,-1]
   sub2 <- sub1[,which(sub1>0)] 
   
})

## trait composition per site
traits_site_fishes <- lapply (list_sites_fishes, function (i)
   
   fish_traits_ord [which(rownames(fish_traits_ord) %in% colnames (i)),]
)

### finding site functional entities
site_FE_fishes <- lapply (traits_site_fishes, function (k)
   unique(
      unlist (lapply (seq (1,nrow (k)), function (i) {
         
         FE <- paste(unname(fish_traits_ord[i,] [1]),
                     unname(fish_traits_ord[i,] [2]),
                     unname(fish_traits_ord[i,] [3]),
                     unname(fish_traits_ord[i,] [4]),
                     unname(fish_traits_ord[i,] [5]),
                     unname(fish_traits_ord[i,] [6]),
                     unname(fish_traits_ord[i,] [7]),
                     unname(fish_traits_ord[i,] [8]),
                     sep=".")
      }  
      )
      )
   )
)

###
fishes_FE <- unlist(lapply (site_FE_fishes, length))
total_fishes_FE <- length(total_funct_ent_fish )
###
bentos_FE <- unlist(lapply (site_funct_ent , length))
total_bentos_FE <- length(total_funct_ent )

## saving these results
save (bent_traits, fish_traits,
      bent2, fish,
      list_benthos_level,list_fish_diet,
      covariates_site,
      covariates_effort,
      FD_results_f1_bentos,FD_results_f1,
      fishes_FE,total_fishes_FE,
      bentos_FE,total_bentos_FE,
      file = here("output", "results_FD_analyses.RData"))
