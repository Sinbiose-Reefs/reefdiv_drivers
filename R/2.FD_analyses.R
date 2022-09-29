
## ------------------------------------------------------------ ##
##                 Analysis of functional diversity             ##
## ------------------------------------------------------------ ##
## loading packages
source("R/packages.R")
source("R/functions.R")

## function to test space quality (from Maire et al. 2015)
source("R/quality_funct_space_fromdist2.R")

# ------------------------------------------ #
# Load rarefied data of fishes and benthos
# sites x transect / video x species
# ------------------------------------------ #

load (here ("data","modeling_data.RData"))


# ---------------------------------- #
# Load fish trait data
# ---------------------------------- #

fish_traits <- read.csv(here("data",
                             "traits",
                             "Atributos_especies_Atlantico_&_Pacifico_Oriental_2020_04_28.csv"), 
                        header=TRUE, sep=";")


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
UVS_spp <- colnames(comp_fish)
fish_traits <- fish_traits[which(fish_traits$Name %in% UVS_spp),]

### df with sp names
fishes <- data.frame (Name=UVS_spp)
### subsetting traits based on defined groups
fish_traits_subset <- fish_traits
### matching
fish_traits2 <- merge(fishes, 
                      fish_traits_subset, 
                      by="Name")

### Select traits of interest and remove NAs
fish_traits2 <- na.omit(fish_traits2[,c("Name",
                                        "TemPref_mean",
                                        "Depth_range",
                                        "Size_group",
                                        "Body_size",
                                        "Level_water",
                                        "Body_shape",
                                        "Home_range",
                                        "Mouth_position",
                                        "Caudal_fin")])
#
# Adjust traits into ordered categories
### body size is continuous
rownames(fish_traits2) <- fish_traits2$Name
#fish_traits2 <- fish_traits2[,-1]
fish_traits2 <- fish_traits2[order(rownames(fish_traits2)),]
#fish_traits2$Body_size <- as.numeric(fish_traits2$Body_size)


### transform body size
size <- (fish_traits2$Body_size) # log the size 

### here transform each ordered category in 1, 2, 3.... 
mobility  <- sapply(fish_traits2$Home_range, function(x) {if (
            x=="sed") {1} 
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
fish_traits_ord <-data.frame (Size=(size), 
                              Temp = fish_traits2$TemPref_mean,
                              Depth = fish_traits2$Depth_range,
                              Mobility=ordered (mobility), 
                              Schooling=ordered (schooling), 
                              #Caudal_fin=as.factor(fish_traits2$Caudal_fin), 
                              #Mouth_position=as.factor(fish_traits2$Mouth_position),
                              Body_shape=as.factor(fish_traits2$Body_shape)
                              )

rownames(fish_traits_ord) <- rownames(fish_traits2)

# filtering out spp not in trait dataset
comp_fish <- comp_fish[,which(colnames(comp_fish) %in% rownames(fish_traits_ord) )]# rm spp not in trait dataset



#----------------------------------#
#      Load benthos trait dataset
# ---------------------------------#



bent_traits <- read.csv(here("data","traits","Database_benthos.csv"), 
                        header=TRUE, sep=";")


# adjust names

bent_traits$taxon <- gsub ( "\\.", " ", bent_traits$taxon)


# Transform benthos body size in ordered category

benthos_body_size <- sapply(bent_traits$body_size, function(x) {
   if (x=="S") {1} 
   else if (x=="M") {2} 
   else if (x=="L") {3} 
   else if (x=="XL") {4}}
)

benthos_body_size  <-ordered (benthos_body_size) #order traits



# Set data for benthos traits with ordered categories
bent_traits_ord <-data.frame (body_size=ordered (benthos_body_size), 
                              growth_form=as.factor(bent_traits$growth_form),
                              modularity=as.factor(bent_traits$modularity), 
                              mobility=as.factor(bent_traits$mobility),
                              reproductive_mode=as.factor(bent_traits$reproductive_mode), 
                              carbonate=as.factor(bent_traits$carbonate.accretion))

rownames(bent_traits_ord) <- bent_traits$groups


# adjusting colnames and filtering species not in the trait dataset
colnames (comp_benthos) <- tolower (colnames(comp_benthos))
comp_benthos <- comp_benthos[,which(colnames(comp_benthos) %in% rownames(bent_traits_ord) )]# rm spp not in trait dataset
benthos_traits <- bent_traits_ord[which(rownames (bent_traits_ord) %in% colnames(comp_benthos)),]
# order
comp_benthos<- comp_benthos [,order(colnames(comp_benthos))]
benthos_traits<- benthos_traits [order(rownames(benthos_traits)),]


# separate corals and algae
corals <- c("agaricia.fragilis",
            "agaricia.humilis",
            "agaricia.sp", 
            "favia.gravida",
            "madracis.decactis",
            "meandrina.brasiliensis",
            "millepora.alcicornis",
            "millepora.incrusting",
            "millepora.nitida",
            "millepora.sp",
            "montastraea.cavernosa", 
            "mussismilia.braziliensis",
            "mussismilia.harttii",
            "mussismilia.hispida", 
            "mussismilia.leptophylla",
            "mussismilia.spp", 
            "porites.astreoides",
            "porites.branneri",
            "porites.sp",
            "siderastrea.spp")



# composition corals
comp_corals <- comp_benthos[,which(colnames(comp_benthos) %in% corals)]
comp_corals<- comp_corals [,order(colnames(comp_corals))]

# corals traits from Bleuel et al. (unpublished data)
require(openxlsx)
coral_traits <- read.xlsx (here ("data", 
                                "traits", 
                                "Planilha_functional_traits_CORAL_TRAITS.xlsx"),
                           sheet=1)

# removing rows without data
coral_traits <- coral_traits[is.na(coral_traits$spp_modified) != T,]
coral_traits$Depth.lower.limit <- as.numeric(coral_traits$Depth.lower.limit)

# modified spp names
coral_traits$spp_modified <- tolower (gsub (" ", ".", coral_traits$spp_modified))
# match names
coral_traits$spp_modified[which(coral_traits$spp_modified == "meandrina.braziliensis")] <- "meandrina.brasiliensis"
coral_traits$spp_modified[which(coral_traits$spp_modified == "mussismilia.hartii")] <- "mussismilia.harttii"
# genera
coral_traits$genera <- sapply (strsplit(coral_traits$spp_modified, "\\."), "[[",1)

# coral traits in coral community
coral_traits_subset <- coral_traits[which (coral_traits$spp_modified %in% corals),
                                    c("spp_modified", "Sexual_system", "Growth_rate", "growth_form", "reprod_mode", "Depth.lower.limit")]

# genera level
require(dplyr)
coral_traits_to_bind<-coral_traits %>% 
  group_by(genera) %>%
  summarise (Sexual_system = getmode(Sexual_system),
             Growth_rate = getmode(Growth_rate),
             growth_form = getmode(growth_form),
             reprod_mode = getmode(reprod_mode),
             Depth.lower.limit = mean(Depth.lower.limit,na.rm=T))


# what's missing
missing_corals <- corals [which(corals %in% coral_traits$spp_modified ==F)]
missing_corals_genera <- sapply (strsplit(missing_corals, "\\."), "[[",1)
# match with trait data to bind
coral_traits_to_bind <- coral_traits_to_bind [match (missing_corals_genera, coral_traits_to_bind$genera),]
# change colnames to match
coral_traits_to_bind <-coral_traits_to_bind %>% 
  rename (spp_modified = genera)
coral_traits_to_bind$spp_modified <- missing_corals

# bind all data
coral_traits_complete <-bind_rows(coral_traits_subset,coral_traits_to_bind)

# order
coral_traits_complete<- coral_traits_complete [order((coral_traits_complete$spp_modified)),]

# names
rownames(coral_traits_complete) <- coral_traits_complete$spp_modified
coral_traits_complete<-coral_traits_complete[,-1]
# change levels
coral_traits_complete$growth_form[which(coral_traits_complete$growth_form == "massive/phaceloid")] <- "massive"
coral_traits_complete$growth_form[which(coral_traits_complete$growth_form == "branching, other")] <- "branching"
coral_traits_complete$growth_form[which(coral_traits_complete$growth_form == "massive, other")] <- "massive"
coral_traits_complete$growth_form[which(coral_traits_complete$growth_form == "plate/other")] <- "plate"
# NA depth
coral_traits_complete$Depth.lower.limit [is.na(coral_traits_complete$Depth.lower.limit)] <- coral_traits_complete$Depth.lower.limit[which(rownames(coral_traits_complete) == "agaricia.sp")]
# binary traits as numeric
coral_traits_complete$Sexual_system <- as.factor(coral_traits_complete$Sexual_system)
coral_traits_complete$Growth_rate <- as.factor(coral_traits_complete$Growth_rate)
coral_traits_complete$reprod_mode <- as.factor(coral_traits_complete$reprod_mode)



# algae


# composition algae
algae <- c("bryopsis.pennata",
           "calcareous.articulate.algae",
           "calcareous.turf",
           "caulerpa.racemosa",
           "caulerpa.sp",
           "caulerpa.verticillata",
           "chaetomorpha.sp",
           "champia.parvula",
           "codium.intertextum",
           "codium.spp",
           "colpomenia.sinuosa",
           "corticated.algae",
           "crostose.coralline.algae",
           "dictyopteris",
           "dictyopteris.plagiogramma",
           "dictyota.sp",
           "foliaceous.algae",
           "galaxaura.sp",
           "gelidiella.acerosa",
           "gelidiopsis",
           "gelidium.floridanum",
           "green.filamentous.algae",
           "halimeda",
           "hypnea.musciformis",
           "jania.amphiroa",
           "laurencia.sp",
          "leathery.algae",
          "padina.sp",
         "sargassum.sp",
         "stypopodium",
         "tricleocarpa.cylindrica",
         "udotea",
         "ulvophyceae",
         "ventricaria.ventricosa",
         "wrangelia")



# composition matrix
comp_algae <- comp_benthos[,which(colnames(comp_benthos) %in% algae)]
algae_traits <- bent_traits_ord[which(rownames(bent_traits_ord) %in% colnames(comp_algae)),]

# order
comp_algae<- comp_algae [,order(colnames(comp_algae))]
algae_traits<- algae_traits [order(rownames(algae_traits)),
                             c("body_size", "growth_form", "carbonate")]





# --------------------------------------------------------------------- #
#    Functional diversity per site, data set, and sample size def
# --------------------------------------------------------------------- #

# ------
# Fish
#-------
require(ape)
tab_ord <- fish_traits_ord
tab_ord$Body_shape <- as.numeric(tab_ord$Body_shape)
tab_ord$Size<- as.numeric(tab_ord$Size)
tab_ord$Temp<- as.numeric(tab_ord$Temp)
tab_ord$Depth<- as.numeric(tab_ord$Depth)
tab_ord$Mobility<- as.numeric(tab_ord$Mobility)
tab_ord$Schooling<- as.numeric(tab_ord$Schooling)
tab_ord$Body_shape<- as.numeric(tab_ord$Body_shape)


# explanation
pcoa(vegdist (tab_ord, "gower"), correction = "cailliez")$values$Rel_corr_eig[1:3]

# composition fish
FD_fish <- dbFD ((fish_traits_ord),
      comp_fish,
      w.abun = T,
      corr = "cailliez",
      calc.FRic = T,
      m = "max",
      scale.RaoQ = F,
      stand.FRic = T,
      calc.CWM = F,
      calc.FDiv = F,
      print.pco = T)



# --------
# algae
#---------


# explanation
tab_ord <- algae_traits
tab_ord$body_size <- as.numeric(tab_ord$body_size)
tab_ord$growth_form <- as.numeric(tab_ord$growth_form)
tab_ord$carbonate <- as.numeric(tab_ord$carbonate)

pcoa(vegdist (tab_ord, "gower"), correction = "cailliez")$values$Rel_corr_eig[1:3]

# composition

FD_algae <- dbFD ((algae_traits),
                    data.matrix(comp_algae),
                  w.abun = T,
                  corr = "cailliez",
                  calc.FRic = T,
                  m = "max",
                  scale.RaoQ = F,
                  stand.FRic = T,
                  calc.CWM = F,
                  calc.FDiv = F,
                  print.pco = T)



# ---------------------------------------


# corals

comm_with_corals <- which(rowSums(comp_corals>0) > 0)
comp_corals_sub <- comp_corals[comm_with_corals,] # rm comm with no sp


# explanation
tab_ord <- coral_traits_complete
tab_ord$Sexual_system <- as.numeric(tab_ord$Sexual_system)
tab_ord$Growth_rate <- as.numeric(tab_ord$Growth_rate)
tab_ord$growth_form <- as.numeric(as.factor(tab_ord$growth_form))
tab_ord$reprod_mode <- as.numeric(tab_ord$reprod_mode)
tab_ord$Depth.lower.limit <- as.numeric(tab_ord$Depth.lower.limit)


pcoa(vegdist (tab_ord, "gower"), correction = "cailliez")$values$Rel_corr_eig[1:3]


# run
FD_corals <- dbFD ((coral_traits_complete),
                   data.matrix(comp_corals_sub),
                   w.abun = T,
                   corr = "cailliez",
                   calc.FRic = T,
                   m = "max",
                   scale.RaoQ = F,
                   stand.FRic = T,
                   calc.CWM = F,
                   calc.FDiv = F,
                   print.pco = T)

# missing data to zero
FD_corals$FRic[is.na(FD_corals$FRic)]<-0

# bind zeros in the communities with no coral
df_corals <-data.frame (
  site = seq (1,nrow (comp_corals)),  
  SR = rowSums (comp_corals>0),
  FRic= 0,
  RaoQ =0
)

df_corals[comm_with_corals,"FRic"] <- FD_corals$FRic
df_corals[comm_with_corals,"RaoQ"] <- FD_corals$RaoQ


# save 
save (comp_fish,
      FD_fish,
      FD_corals,
      comm_with_corals,
      comp_corals,
      df_corals,
      FD_algae,
      comp_algae,
      file=here ("output", "FD_results.RData"))


