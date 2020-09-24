##### Analysis of functional diversity patterns of Brazilian reefs ###

#### Load fish and benthos occurrence data
#### Run script on "codigo_organizacao_dados.R" until 

# dados_peixes_bentos <- list(peixes = comp_peixes, bentos = comp_bentos)

fish <- dados_peixes_bentos[[1]]
bent <- dados_peixes_bentos[[2]]

#### Load trait data for fish and benthos

## Fish traits

fish_traits <- read.csv(here("data","traits","Atributos_especies_Atlantico_&_Pacifico_Oriental_2020_04_28.csv"), 
                        header=TRUE, sep=";")
fish_traits$Name <- tolower(gsub(" ", ".", fish_traits$Name)) #Adjust spp names in trait dataset

# Select fish species from Atlantic database to match UVCs

fishes <- as.data.frame(colnames(fish))
colnames(fishes) <- paste("Name")

fish_traits2 <- merge(fishes, fish_traits, by="Name")
dim(fish_traits2)

# Select traits of interest and remove NAs

fish_traits2 <- na.omit(fish_traits2[,c("Name","Home_range","Diel_activity","Size_group","Body_size",
                                "Level_water","Body_shape","Diet","Caudal_fin","Mouth_position")])
dim(fish_traits2)

# Adjust traits into ordered categories

rownames(fish_traits2) <- fish_traits2$Name
fish_traits2 <- fish_traits2[,-1]
fish_traits2 <- fish_traits2[order(rownames(fish_traits2)),]
fish_traits2$Body_size <- as.numeric(fish_traits2$Body_size)

# here transform each ordered category in 1, 2, 3.... 
mobility  <- sapply(fish_traits2$Home_range,     function(x) {if (x=="sed")         {1} else if (x=="mob")       {2} else if (x=="vmob")        {3}})
schooling <- sapply(fish_traits2$Size_group,      function(x) {if (x=="sol")         {1} else if (x=="pair")          {2} else if (x=="smallg")           {3} else if (x=="medg")           {4} else if (x=="largeg")            {5}})
level  <- sapply(fish_traits2$Level_water,    function(x) {if (x=="bottom")         {1} else if (x=="low")          {2} else if (x=="high")           {3}})

size      <-log (fish_traits2$Body_size) # log the size 
mobility  <-ordered (mobility) # Order
level  <-ordered (level) 
schooling <-ordered (schooling) 

fish_traits_ord <-data.frame (Size=size,  Mobility=mobility, Activity=fish_traits2$Diel_activity, 
                             Schooling=schooling, Level=level, Diet=fish_traits2$Diet, Body_shape=fish_traits2$Body_shape,
                             Caudal_fin=fish_traits2$Caudal_fin, Mouth_position=fish_traits2$Mouth_position)

rownames(fish_traits_ord) <- rownames(fish_traits2)
dim(fish_traits_ord)

######## Building the gower distance matrix ########

# first calculate gower distance on traits

gower_matrix <- daisy (fish_traits_ord, metric=c("gower")) 

# Building the functional space based on a PCOA 

pco<-dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
barplot(pco$eig) # barplot of eigenvalues for each axis 
Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig)) # percentage of inertia explained by the two first axes

# Testing the Quality of the functional space based on the method of Maire et al. (2015) in GEB, i.e. how many axes do we need to keep to faithfully represent the original Gower's distances
source("quality_funct_space_fromdist2.R")

quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   plot="quality_funct_space_I") # The function will test the quality of the space from 2 to 10 axes using dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
quality$meanSD # the minimal value corresponds to the best space to use, here 7 axes 
Inertia7<- (sum(pco$eig[1:7])) /(sum(pco$eig)) # percentage of inertia explained by the 6 first axes = 73%

# Functional diversity indices for fish #######

# Adjust fish abundance data to relative abundance

fish_abund <- fish[,which(colnames(fish) %in% rownames(fish_traits_ord))] # match fish spp in both datasets
rownames(fish_abund) <- fish$locality_site

fish_rel <- fish_abund/rowSums(fish_abund) # relative abundance for fish


# Method with dbFD() 

fd <- dbFD (gower_matrix, fish_rel, corr ="none", m = 7, stand.FRic = TRUE,
            print.pco = TRUE)

fish_fd <- data.frame (nbsp=fd$nbsp, sing.sp=fd$sing.sp, FRic=fd$FRic, qual.FRic=fd$qual.FRic,
                      FEve=fd$FEve, FDiv=fd$FDiv, FDis=fd$FDis, RaoQ=fd$RaoQ)

write.table (fish_fd, "FIndex_fish.csv", sep=";", dec=".")

#################################################################

# Analysis for benthos

# Remove benthos from dataset as below

bent2 <- bent[ , -which(names(bent) %in% c("Areia.e.Cascalho","Desconhecido","Estrela","ourico1","ourico2","Outra.ascidia","Outro.anthozoa",
                                  "Outro.echinoderma","Outro.hydrozoa","Quadrado","Outro.crustaceo","Sombra"))]

colnames(bent2) <- tolower(gsub(" ", ".", colnames(bent2))) #adjust colnames

benthos_spp <- as.data.frame(colnames(bent2[2:96]))
colnames(benthos_spp) <- paste("groups") 

# Work in benthos traits

# Load trait data
bent_traits <- read.csv(here("data","traits","Database_benthos.csv"), 
                        header=TRUE, sep=",")

bent_traits$groups <- tolower(gsub("_", ".", bent_traits$groups)) #adjust spp names in trait dataset
bent_traits2 <- merge(bent_traits, benthos_spp, by="groups") # merge traits with groups in occurrence data
rownames(bent_traits2) <- bent_traits2$groups
bent_traits2 <- bent_traits2[2:8]


# Transform benthos body size in ordered category

bent_traits2 <- bent_traits2[order(rownames(bent_traits2)),]
benthos_body_size <- sapply(bent_traits2$body_size,     function(x) {if (x=="S")         {1} else if (x=="M")       {2} else if (x=="L")        {3} else if (x=="XL")            {4}})

benthos_body_size  <-ordered (benthos_body_size) #order traits


# Set data for benthos traits with ordered categories
bent_traits_ord <-data.frame (body_size=benthos_body_size, growth_form=bent_traits2$growth_form,
                              modularity=bent_traits2$modularity, mobility=bent_traits2$mobility,
                              reproductive_mode=bent_traits2$reproductive_mode, carbonate=bent_traits2$carbonate.accretion)
                
rownames(bent_traits_ord) <- rownames(bent_traits2)
dim(bent_traits_ord)

######## Building the gower distance matrix ########

# first calculate gower distance on traits

gower_benthos <- daisy (bent_traits_ord, metric=c("gower")) 

# Building the functional space based on a PCOA 

pco<-dudi.pco(quasieuclid(gower_benthos), scannf=F, nf=10) # quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
barplot(pco$eig) # barplot of eigenvalues for each axis 
Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig)) # percentage of inertia explained by the two first axes

# Testing the Quality of the functional space

quality<-quality_funct_space_fromdist(gower_benthos,  nbdim=10,   plot="quality_funct_space_I") # The function will test the quality of the space from 2 to 10 axes using dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
quality$meanSD # the minimal value corresponds to the best space to use, here 4 axes 
Inertia4<- (sum(pco$eig[1:4])) /(sum(pco$eig)) # percentage of inertia explained by the 6 first axes = 73%

# Functional diversity indices for benthos #######

# Adjust benthos abundance data to relative abundance

bent_abund <- bent2[,which(colnames(bent2) %in% rownames(bent_traits_ord))] # match fish spp in both datasets
rownames(bent_abund) <- bent2$locality_site

bent_rel <- bent_abund/rowSums(bent_abund) # relative abundance for fish


# Method with dbFD() 

fd_benthos <- dbFD (gower_benthos, bent_rel, corr ="none", m = 4, stand.FRic = TRUE,
            print.pco = TRUE)

bent_fd <- data.frame (nbsp.bent=fd_benthos$nbsp, sing.sp.bent=fd_benthos$sing.sp, FRic.bent=fd_benthos$FRic, qual.FRic.bent=fd_benthos$qual.FRic,
                       FEve.bent=fd_benthos$FEve, FDiv.bent=fd_benthos$FDiv, FDis.bent=fd_benthos$FDis, RaoQ.bent=fd_benthos$RaoQ)

write.table (bent_fd, "FIndex_bent.csv", sep=";", dec=".")


################ Plots

fd_reefs <- cbind(bent_fd, fish_fd)

par(mfrow=c(1,3))
plot(fd_reefs$nbsp.bent, fd_reefs$FRic.bent, ylim=c(0,0.8), xlim=c(0,20), pch=21, col="coral", bg="coral",
     ylab="Benthos FRic", xlab="Species Richness", cex=1.2)
plot(fd_reefs$nbsp, fd_reefs$FRic, ylim=c(0,0.4), xlim=c(10,70), pch=21, col="cyan3", bg="cyan3",
     ylab="Fish FRic", xlab="Species Richness", cex=1.2)
plot(fd_reefs$FRic.ben, fd_reefs$FRic, ylim=c(0,0.3), xlim=c(0,0.8), pch=21, col="goldenrod", bg="goldenrod",
     ylab="Fish FRic", xlab="Benthos FRic", cex=1.2)

par(mfrow=c(2,2))
plot(fd_reefs$FEve.bent~fd_reefs$nbsp.bent, ylim=c(0,0.8), xlim=c(0,30), pch=21, col="darkgray", bg="darkgray",
     ylab="Benthos FEveness", xlab="Benthos Richness", cex=1.2)
plot(fd_reefs$FDiv.bent~fd_reefs$nbsp.bent, ylim=c(0,1), xlim=c(0,30), pch=21, col="darkgray", bg="darkgray",
     ylab="Benthos FDivergence", xlab="Benthos Richness", cex=1.2)
plot(fd_reefs$FEve~fd_reefs$nbsp, ylim=c(0,0.8), xlim=c(0,70), pch=21, col="darkgray", bg="darkgray",
     ylab="Fish FEveness", xlab="Fish Richness", cex=1.2)
plot(fd_reefs$FDiv~fd_reefs$nbsp, ylim=c(0,1), xlim=c(0,70), pch=21, col="darkgray", bg="darkgray",
     ylab="Fish FDivergence", xlab="Fish Richness", cex=1.2)


## Load environmental variables

env <- load(here("output","data_drivers_analysis.RData"))
env <- cbind(covariates_site$coord, covariates_site$sea_data, covariates_site$biog_reef)

data_fd <- cbind(fd_reefs, env)

save (data_fd, ### functional indices and environmental variables
      file=here("output","functional_indices.RData"))

## Plots with environmental data (latitude)
#fish
par(mfrow=c(2,3))
plot(data_fd$FRic~data_fd$Lat, ylim=c(0,0.3), xlim=c(0,-30), pch=21, col="cyan2", bg="cyan2",
     ylab="Fish FRic", xlab="", xaxt="n", cex=1.2)
plot(data_fd$FEve~data_fd$Lat, ylim=c(0,1), xlim=c(0,-30), pch=21, col="cyan2", bg="cyan2",
     ylab="Fish FEveness", xlab="", xaxt="n", cex=1.2)
plot(data_fd$FDiv~data_fd$Lat, ylim=c(0,1), xlim=c(0,-30), pch=21, col="cyan2", bg="cyan2",
     ylab="Fish FDivergence", xlab="", xaxt="n", cex=1.2)

#benthos
plot(data_fd$FRic.bent~data_fd$Lat, ylim=c(0,1), xlim=c(0,-30), pch=21, col="plum3", bg="plum3",
     ylab="Benthos FRic", xlab="Latitude", cex=1.2)
plot(data_fd$FEve.bent~data_fd$Lat, ylim=c(0,1), xlim=c(0,-30), pch=21, col="plum3", bg="plum3",
     ylab="Benthos FEveness", xlab="Latitude", cex=1.2)
plot(data_fd$FDiv.bent~data_fd$Lat, ylim=c(0,1), xlim=c(0,-30), pch=21, col="plum3", bg="plum3",
     ylab="Benthos FDivergence", xlab="Latitude", cex=1.2)

## Plots with environmental data (temperature)
#fish
par(mfrow=c(2,3))
plot(data_fd$FRic~data_fd$BO2_tempmean_ss, ylim=c(0,0.3), xlim=c(20,30), pch=21, col="cyan2", bg="cyan2",
     ylab="Fish FRic", xlab="", xaxt="n", cex=1.2)
plot(data_fd$FEve~data_fd$BO2_tempmean_ss, ylim=c(0,1), xlim=c(20,30), pch=21, col="cyan2", bg="cyan2",
     ylab="Fish FEveness", xlab="", xaxt="n", cex=1.2)
plot(data_fd$FDiv~data_fd$BO2_tempmean_ss, ylim=c(0,1), xlim=c(20,30), pch=21, col="cyan2", bg="cyan2",
     ylab="Fish FDivergence", xlab="", xaxt="n", cex=1.2)

#benthos
plot(data_fd$FRic.bent~data_fd$BO2_tempmean_ss, ylim=c(0,1), xlim=c(0,30), pch=21, col="plum3", bg="plum3",
     ylab="Benthos FRic", xlab="Mean SST", cex=1.2)
plot(data_fd$FEve.bent~data_fd$BO2_tempmean_ss, ylim=c(0,1), xlim=c(0,30), pch=21, col="plum3", bg="plum3",
     ylab="Benthos FEveness", xlab="Mean SST", cex=1.2)
plot(data_fd$FDiv.bent~data_fd$BO2_tempmean_ss, ylim=c(0,1), xlim=c(0,30), pch=21, col="plum3", bg="plum3",
     ylab="Benthos FDivergence", xlab="MEan SST", cex=1.2)
