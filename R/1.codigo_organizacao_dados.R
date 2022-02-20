
# ------------------------------------------------------------------- #
#         Organizing data for the MS "Diversity patterns and drivers" 
#-------------------------------------------------------------------- #

## load packages and functions
source("R/packages.R")
source("R/functions.R")
#source("R/function_lomolino_richness.R")

# Load data

# ------------------------------ #
#        Benthic data
# ------------------------------ #

bentos <- read.xlsx(here("data","detection","Updated_compiled_quadrats_allsites.xlsx"),
                    sheet = 1, colNames = TRUE,detectDates=F)

## adjust data - bug of openxlsx
bentos$eventDate <-convertToDate(bentos$eventDate)

# Remove benthos from dataset as below
rm_sp <- c("Areia.e.Cascalho","Desconhecido","Estrela","ourico1","ourico2",
           "Outra.ascidia","Outro.anthozoa", "Outro.echinoderma",
           "Outro.hydrozoa","Quadrado","Outro.crustaceo","Sombra")

bentos <- bentos[which(bentos$Taxon %in% rm_sp == F),]
# to lower case  
bentos$Taxon <- tolower (bentos$Taxon)

# proportion of not identified taxa (family was considered identified)
1-(74/95)

# -----------------------------
#       Fish data
# -----------------------------

peixes <- read.xlsx(here("data","detection","UpdatedData_RMorais_et_al_2017.xlsx"),
                    sheet = 1, colNames = TRUE,detectDates=F)

## adjust data
peixes$eventDate <-convertToDate(peixes$eventDate)

# obtain an event ID (it is the ID of each site) 
bentos$eventID_MOD <- substr(bentos$eventID, 1,nchar(as.character(bentos$eventID))-5) 
peixes$eventID_MOD  <- substr(peixes$eventID, 1,nchar(as.character(peixes$eventID))-5) 

## sensitivity analyzes could be done for sites present in both data sets
peixes_subset <- peixes [which(peixes$eventID_MOD %in% bentos$eventID_MOD),]
#peixes_subset <- peixes # but now we work with all data per taxa

## the same for benthos
bentos_subset <- bentos [which(bentos$eventID_MOD %in% peixes$eventID_MOD),]
#bentos_subset <- bentos # complete dataset
##  Observer ID (number)
peixes_subset$ID.observer <- as.numeric(as.factor(peixes_subset$Observer))

## ID for sites
peixes_subset$locality_site <- paste(peixes_subset$Locality,  # fishes
                                     peixes_subset$Site,
                                     peixes_subset$eventDepth,
                                     sep=".")
#
bentos_subset$locality_site <- paste(bentos_subset$Locality, # benthos
                                     bentos_subset$Site,
                                     bentos_subset$eventDepth,
                                     sep=".")

# REMOVE ISLAND SITES
peixes_subset <- peixes_subset[which(peixes_subset$Region != "oc_isl"),]
bentos_subset <- bentos_subset[which(bentos_subset$Region != "oc_isl"),]

# total number of belt transects
length(unique(peixes_subset$Transect_id))
# total number of videos
sum(
  rowSums(
    table(bentos_subset$locality_site, 
          bentos_subset$Video_number)>0)
  )
  
#---------------------------------------------------#
#             composition data                      #
#---------------------------------------------------#

#########################################
## fish composition (site x spp)
#########################################

comp_peixes <-  cast(peixes_subset,
                     formula = locality_site  ~ ScientificName,
                     #formula = Locality  ~ ScientificName,
                     value= "IndCounting",
                     fun.aggregate = sum)

comp_peixes <- comp_peixes[order (comp_peixes$locality_site,decreasing=F),]
#comp_peixes <- comp_peixes[order (comp_peixes$Locality,decreasing=F),]

###########################################################################
## data for SAC (species accumulation curve, site  x transect x spp)
###########################################################################

# unuique sites
sites_fish_complete <- unique (peixes_subset$locality_site)
sites_fish_complete <- sites_fish_complete[order(sites_fish_complete,decreasing=F)]

# tables per site
rar_peixes <- lapply (sites_fish_complete, function (i)
    
    cast (peixes_subset [which (peixes_subset$locality_site == i),],
          formula = Transect_id ~ ScientificName,
          value="IndCounting",
          fun.aggregate = sum)[,-1]
    )

# which sites have more than 2 transects
enough_transects <- which(unlist(lapply(rar_peixes,nrow))>2)

# rm sites with too few samples
sites_fish_complete <- sites_fish_complete[enough_transects]

# also from rarefaction data
rar_peixes <- rar_peixes[enough_transects]

# SAC based on the number of samples,
# nrandom samples  (for all analysis)
niter <- niter
rarefied_richness_fish <- lapply (rar_peixes, 
                                  specaccum,
                                  method="random", 
                                  permutations=niter)

### plotting 
par (mfrow=c(1,1))
plot(NA,
     xlim=c(0,150),
     ylim=c(-5,100),
     xlab = "Number of samples",
     ylab = "Fish richness")

lapply (rarefied_richness_fish,plot,add=T)
abline(v=5,lwd=2,col="gray50",lty=2)

##---------------------------------- ## 
# richness estimate based on the minimum number of samples
## ---------------------------------- ##

min_samples <- sapply (rarefied_richness_fish,"[[","sites")
# finding the minimum number of samples across sites
min_samples<- min(unlist(sapply (min_samples,max,simplify=F)))

# finding estimated richness
est_rich <- sapply (rarefied_richness_fish,"[[","richness")
# finding the richness for min_samples
est_rich <- lapply (est_rich, function (i) i [min_samples])

# finding sd
# finding estimated richness
sd_rich <- sapply (rarefied_richness_fish,"[[","sd")
# finding the richness for min_samples
sd_rich <- lapply (sd_rich, function (i) i [min_samples])

# toatl number of samples per site
total_samples <- sapply (rarefied_richness_fish,"[[","sites")
total_samples <- unlist (lapply (total_samples, max))

# res table

res_table_samples <- data.frame (Site = sites_fish_complete,
                                 EST.Rich=unlist(est_rich),
                                 SD.Rich=unlist(sd_rich),
                                 n_samples = min_samples,
                                 n_total= total_samples)

# finding site composition based on min_samples 
# a random sample of five transects

# -------------------------------#
# 1) random sample based on min samples
# -------------------------------#

# obtain random compositions

ncores <- ncores
cl <- makeCluster(ncores) ## number of cores = generally ncores -1

# exportar pacote para os cores

# export your data and function
clusterExport(cl, c("niter", 
                    "fc_random_composition",
                    "rar_peixes",
                    "min_samples"))

rdm_composition <- parLapply (cl, seq(1,niter), function (k)
                           lapply(rar_peixes, function (i)
  
                fc_random_composition(data = i,
                                      nsamples = min_samples,
                                      replace= F)
                

                     )
                  )
                

stopCluster (cl)

# list of species to adjust cols
list_spp <- unique(unlist(lapply (rar_peixes, function (k)
   colnames(k))))

# matrix to bind (with missing names)

rdm_composition_complete <-lapply (rdm_composition, function (i)
  do.call(rbind, lapply (i, function (k) {

    bind_to_rdm_composition <- matrix(0, 
                                      ncol = length(which (list_spp %in% colnames (k)  == F)),
                                 dimnames = list (NULL,
                                                  list_spp[which (list_spp %in% colnames (k)  == F)]
                                 ))
    # cbind
    bind_rdm_composition <- cbind(k,bind_to_rdm_composition)
    # ordering by sp name to avoid problems when combining samples from diff sites
    bind_rdm_composition <- bind_rdm_composition [,order (colnames(bind_rdm_composition))];
    bind_rdm_composition

})))

# saving

save (comp_peixes,# composition per site (not rarefied)
      rar_peixes,# observed data
      rdm_composition_complete, # random composition based on min samples
      res_table_samples, # results of estimates based on min samples
      sites_fish_complete, # analyzed sites
      file = here("output", "random_composition_fish.RData"))


#########################################
## benthos composition (site  x spp)
#########################################

comp_bentos <-  cast(bentos_subset,
                     formula = locality_site  ~ Taxon,
                     #formula = Locality  ~ Taxon,
                     value= "Cover",
                     fun.aggregate = mean)

comp_bentos<- comp_bentos[order(comp_bentos$locality_site,decreasing=F),]
#comp_bentos<- comp_bentos[order(comp_bentos$Locality,decreasing=F),]

# save benthos and fish composition (non-rarefied)
#save (comp_peixes, comp_bentos,
#      file=here ("output","composition_non_rarefied.RData"))

###################################################################
## data for SAC (species accumulation curve, site  x video x spp)
###################################################################
sites_bentos_complete <- unique (bentos_subset$locality_site)
sites_bentos_complete <- sites_bentos_complete [order(sites_bentos_complete,decreasing = F)]

# create one table of samples per speices x site 
rar_bentos <- lapply (sites_bentos_complete, function (i)
  
  cast (bentos_subset [which (bentos_subset$locality_site == i),],
        formula = Video_number ~ Taxon,
        value="Cover",
        fun.aggregate = max)[,-1])

# which sites have more than 2 videos
enough_videos <- which(unlist(lapply(rar_bentos,nrow))>2)

# rm sites with too few samples
sites_bentos_complete <- sites_bentos_complete[enough_videos]

# rm also from rarefaction data
rar_bentos <- rar_bentos[enough_videos]

## SAC based on the number of samples

rarefied_richness_bentos<- lapply (rar_bentos, 
                                  specaccum,
                                  method="random", ## method="rarefaction" - based on individuals
                                  permutations=niter)

### plotting 
par (mfrow=c(1,1))
plot(NA,
     xlim=c(0,25),
     ylim=c(-2,40),
     xlab = "Number of samples",
     ylab = "Benthos richness")

lapply (rarefied_richness_bentos,plot,add=T)
abline(v=3,lwd=2,col="gray50",lty=2)

##---------------------------------- ## 
# based on minimum number of samples
## ---------------------------------- ##
min_samples_bentos <- sapply (rarefied_richness_bentos,"[[","sites")
# finding the minimum number of samples across sites
min_samples_bentos<- min(unlist(sapply (min_samples_bentos,max,simplify=F)))

# finding estimated richness
est_rich <- sapply (rarefied_richness_bentos,"[[","richness")
# finding the richness for min_samples
est_rich <- lapply (est_rich, function (i) i [min_samples_bentos])

# finding sd
# finding estimated richness
sd_rich <- sapply (rarefied_richness_bentos,"[[","sd")
# finding the richness for min_samples
sd_rich <- lapply (sd_rich, function (i) i [min_samples_bentos])

# toatl number of samples per site
total_samples <- sapply (rarefied_richness_bentos,"[[","sites")
total_samples <- unlist (lapply (total_samples, max))

# res table
res_table_samples_bentos <- data.frame (Site = sites_bentos_complete,
                                 EST.Rich=unlist(est_rich),
                                 SD.Rich=unlist(sd_rich),
                                 n_samples = min_samples_bentos,
                                 n_total= total_samples)

# finding site composition based on min_samples 
# a random sample of five transects

# ----------------------------------------#
# 1) random sample based on min samples
# ----------------------------------------#

# obtain random compositions

cl <- makeCluster(ncores) ## number of cores = generally ncores -1

# exportar pacote para os cores

# export your data and function
clusterExport(cl, c("niter", 
                    "fc_random_composition",
                    "rar_bentos",
                    "min_samples_bentos"))

rdm_composition_bentos <- parLapply (cl, seq(1,niter), function (k)
  lapply(rar_bentos, function (i)
    fc_random_composition(data = i,
                          nsamples = min_samples_bentos,
                          replace= F)
  ))

stopCluster (cl)

# list of species to adjust cols
list_spp <- unique(unlist(lapply (rar_bentos, function (k)
  colnames(k))))

# matrix to bind (with missing names)

rdm_composition_complete_bentos <-lapply (rdm_composition_bentos, function (i)
  do.call(rbind, lapply (i, function (k) {
    
	if (length(which (list_spp %in% colnames (k)  == F))==0){k}
	
	else {
	
	    bind_to_rdm_composition <- matrix(0, 
                                      ncol = length(which (list_spp %in% colnames (k)  == F)),
                                      dimnames = list (NULL,
                                                       list_spp[which (list_spp %in% colnames (k)  == F)]                                      ))
    
    # cbind
    bind_rdm_composition <- cbind(k,bind_to_rdm_composition)
    # ordering by sp name to avoid problems when combining samples from diff sites
    bind_rdm_composition <- bind_rdm_composition [,order (colnames(bind_rdm_composition))];
    bind_rdm_composition
	
	}
 ;k
    
  })))

# rm first column
# save

save (comp_bentos,
      rar_bentos,
      rdm_composition_complete_bentos, # random composition based on min samples
      res_table_samples_bentos, # results of estimates based on min samples
      sites_bentos_complete, # analyzed sites
      file = here("output", "random_composition_bentos.RData"))


#######################################################################################
################         ENVIRONMENTAL PREDICTORS               #######################
#######################################################################################

## REGION
# northeastern
nord <- cast(peixes_subset,
             formula = locality_site ~ Region,
             value= "IndCounting",
             fun.aggregate = sum)
nord<- nord[which(nord$locality_site %in% sites_fish_complete),]
regiao <- ifelse (nord$ne_reefs>0,"northeastern","southern")

# depth
prof <- cast(peixes_subset,
                     formula = locality_site ~ eventDepth,
                     value= "IndCounting",
                     fun.aggregate = sum)
prof<- prof[which(prof$locality_site %in% sites_fish_complete),]
prof <- ifelse (prof$fundo>0,"fundo","raso")

# ---------
# effort

# unique transection IDs, checking length == n trans
n_transeccoes <- lapply (sites_fish_complete, function (i)
  length (unique (peixes_subset [which(peixes_subset$locality_site == i),"Transect_id"]))
)

# unique videos per site
# fish site names as comparison
n_videos <- lapply (sites_bentos_complete, function (i)
  length (unique (bentos_subset [which(bentos_subset$locality_site == i),"Video_number"]))
)

# effort table
tabela_esforco <- list(fish=data.frame (sites_fish_complete = sites_fish_complete,
                              n_transectos_peixes = unlist(n_transeccoes)),
                       benthos = data.frame (sites_bentos_complete = sites_bentos_complete,
                              n_videos_bentos = unlist(n_videos)))
## list of data
# complete
comp_peixes <- comp_peixes[which(comp_peixes$locality_site %in% sites_fish_complete),]
comp_bentos <- comp_bentos [which(comp_bentos$locality_site %in% sites_bentos_complete),]
dados_peixes_bentos <- list(peixes = comp_peixes,
                            riq_peixes = rowSums (comp_peixes[,-1]>0),
                            bentos = comp_bentos,
                            riq_bentos=rowSums (comp_bentos[,-1]>0))

## geo coordinates
# fish
coordenadas_peixes <- aggregate(peixes_subset, 
                                by= list (peixes_subset$locality_site), 
                                FUN=mean)[c("Group.1","Lon","Lat")]
coordenadas_peixes<- coordenadas_peixes[which(coordenadas_peixes$Group.1 %in% sites_fish_complete),]

# benthos
coordenadas_bentos <- aggregate(bentos_subset, 
                         by= list (bentos_subset$locality_site), 
                         FUN=mean)[c("Group.1","Lon","Lat")]
coordenadas_bentos<- coordenadas_bentos[which(coordenadas_bentos$Group.1 %in% sites_bentos_complete),]

# -----------------------
# BiO Oracle - extracting covariate data
# Explore datasets in the package
# devtools::install_github("lifewatch/sdmpredictors")
layers <- list_layers()
#View (layers [grep ("Bio-ORACLE",layers$dataset_code),])

# Download specific layers to the current directory
# set prefered folder (to download data)
options(sdmpredictors_datadir=here ("data","environment"))

## chlorophil has different extent - loading and extracting in two steps         
layers_oracle <- load_layers(c("BO2_tempmean_ss",
                               #"BO2_temprange_ss",
                               "BO2_ppmean_ss", 
                               #"BO2_pprange_ss",
                               "BO2_salinitymean_ss", 
                               #"BO2_salinityrange_ss",
                               #"BO_damax",
                               "BO_damean"
                               #,"BO_damin"
                               ))

## these data have different extent
#layers_oracle_Chl <- load_layers (c("BO2_chlomean_ss","BO2_chlorange_ss"))

#f <- system.file("data/environment/BO2_tempmean_ss_lonlat.tif", 
#                 package="raster")

## coordinates to spatial points
# adjusting one coordinate
coordenadas_peixes [grep("perua",coordenadas_peixes$Group.1),2] <- as.numeric(-35.082658) ## perua preta

# list of points
sp_points <- list(coordenadas_peixes,
                  coordenadas_bentos)
# coord to sppoints df
spdf <- lapply (sp_points, function (i) 
  
  SpatialPointsDataFrame(coords = i[,2:3], data = i,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")))

## extracting data

extracted_sea_data <- lapply (spdf, function (i) 
  extract (layers_oracle, i,method='simple', fun=mean))
extracted_sea_data<-lapply (seq(1,length(extracted_sea_data)), function (i) {
  rownames(extracted_sea_data[[i]]) <- sp_points[[i]]$Group.1;
  extracted_sea_data[[i]]
})

#extracted_sea_data_Chl <- extract (layers_oracle_Chl, spdf,method='simple', fun=mean)
#rownames(extracted_sea_data_Chl ) <- sp_points$Group.1

## binding these dfs
#extracted_sea_data <- cbind(extracted_sea_data,

# ------------------------------
#  distance offshore
# BR coastline, download from here https://mapcruzin.com/free-brazil-arcgis-maps-shapefiles.htm

BR <- readOGR(dsn=here("data", "environment","brazil-coastline"), "brazil_coastline")
crs(BR) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
BR <- spTransform(BR, CRS("+init=epsg:4326"))

# use dist2Line from geosphere - only works for WGS84 
# data benthos
sp_data <- SpatialPoints(coordenadas_bentos[,2:3])
crs(sp_data) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
sp_data <- spTransform(sp_data, CRS("+init=epsg:4326"))

# measuring the distance
dist_bentos <- geosphere::dist2Line(p = sp_data, 
                             line = (BR))
# binding coords
coordenadas_bentos <- cbind (coordenadas_bentos, dist_bentos)
# plot
ggplot(data = world) +
  geom_sf() +
  geom_point(data=coordenadas_bentos, 
             aes(x=Lon, y=Lat, col =distance))# +
  #coord_sf(xlim = c(-55, -20), ylim = c(-33,0 ), expand = FALSE)

## fish
sp_data <- SpatialPoints(coordenadas_peixes[,2:3])
crs(sp_data) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
sp_data <- spTransform(sp_data, CRS("+init=epsg:4326"))

# measuring distance
dist_peixes <- geosphere::dist2Line(p = sp_data, 
                                    line = (BR))

# bind coords and distance
coordenadas_peixes <- cbind (coordenadas_peixes, dist_peixes)
#plot to check
ggplot(data = world) +
  geom_sf() +
  geom_point(data=coordenadas_peixes, 
             aes(x=Lon, y=Lat, col =distance))# +
  #coord_sf(xlim = c(-55, -20), ylim = c(-33,0 ), expand = FALSE)

# ====================== 
# magris data (reef area)

shapefiles<-list.files(here ("data","environment","magris_reef_map"),pattern=".shp")
# load all at once
shapes <- lapply (shapefiles, function (shp) 
  
      readOGR (here ("data","environment","magris_reef_map"),
                   gsub(".shp","",shp)) 
)

# subset of habitats
list_habitats <- list ("amazon" = "AO11",# amazon
                      "eastern"= c("EC11",# eastern
                         "EC12",
                         "EC13",
                         "EC14",
                         "EC15",
                         "EC16"),
                      "northeastern" = c("NC11",# northeastern
                                         "NC12",
                                         "NC13"),
                       "riogrande"="RC11", # rio grande
                      "southeastern"="SC11" #southeastern
                       )

# reef location
BR_reefs <- lapply (seq (1,length(shapes)), function (shp)
  # extract codes  
  shapes[[shp]][which(shapes[[shp]]$habitat %in% list_habitats[[shp]]),]
  
)


# bind 
BR_reefs <- bind(BR_reefs[[1]],
                 BR_reefs[[2]],
                 BR_reefs[[3]],
                 BR_reefs[[4]],
                 BR_reefs[[5]])

# create a grid for extracting data
# based on the extent of extracted data
grd_df <- expand.grid(x = seq(from = extent (BR_reefs)[1]-1,
                              to = extent (BR_reefs)[2]+1, 
                              by = .5),
                      y = seq(from = extent (BR_reefs)[3]-1,                                           
                              to = extent (BR_reefs)[4]+1, 
                              by = .5))  # expand points to grid

# Convert grd object to a matrix and then turn into a spatial
# points object
coordinates(grd_df) <- ~x + y

# Sp points into raster
grd_raster <- (raster(grd_df,resolution = .5))
crs(grd_raster) <-crs(BR_reefs)
values (grd_raster) <- runif(n=ncell(grd_raster))

# project
# reproject to laea # https://weiming-hu.github.io/projection-in-R/
BR_reefs <- spTransform(BR_reefs, 
                        "+proj=laea +lat_0=0 +lon_0=0 +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
to_raster <-  raster (extent(BR_reefs),
                      crs = "+proj=laea +lat_0=0 +lon_0=0 +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
                      res=55000)
grd_raster <- projectRaster(grd_raster, to_raster)

# southe maerica map
southAme<- readOGR(dsn= here("data","South_America"),encoding="latin1", 
                   layer="South_America")
southAme <- spTransform(southAme, 
                        "+proj=laea +lat_0=0 +lon_0=0 +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# extract
require("exactextractr")

# extract
prop_area <- coverage_fraction(grd_raster, 
                               st_combine(st_as_sf(BR_reefs)))[[1]]
# points
pts <- spTransform(spdf[[1]], 
                   "+proj=laea +lat_0=0 +lon_0=0 +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# plot
plot1 <- gplot(prop_area) +  
  geom_tile(aes(fill=(value)),alpha=0.8) +
  scale_fill_gradient(low = "white", high = "darkred",
                     name="Reef area") +
   theme_classic() 

# 
plot2 <- plot1 + geom_polygon(data=southAme, 
                                    aes(x=long, y=lat, group=group),
                                    size = 0.1, fill="gray60", 
                                    colour="gray75",alpha=0.1) + 
  xlab("Longitude") + ylab("Latitude") +
  coord_fixed (xlim = c(-5000000,-2500000), 
               ylim = c(-4000000, 0), ratio = 1) +
  theme(legend.position = c(0.8,0.7),
        axis.text.x = element_text(angle=45))

pdf (here ("output_no_resampling","reef_area.pdf"),width=3,height = 4)
plot2
dev.off()

# values in site coords
extracted_area <- extract (prop_area, 
                           pts,
                           method='simple', fun=mean)

# spdf[[1]]$Lon == spdf[[2]]$Lon
## list with predictors
covariates_site <- list (site_names = sites_fish_complete,
                         regiao = regiao,
                         prof = prof,
                         reef_area = extracted_area,
                         coord = list(coord_bentos=coordenadas_bentos,
                                      coord_peixes = coordenadas_peixes),
                         sea_data = extracted_sea_data)
# effort
covariates_effort <- list(effort = tabela_esforco)


cor(data.frame (covariates_site$sea_data[[1]],
                area = covariates_site$reef_area,
                dist = covariates_site$coord$coord_bentos$distance))

# SAVE

save (covariates_site , ### dados de covariaveis das bases de dados
      covariates_effort, ## variaveis de esforco
      dados_peixes_bentos, ## dados de composicao 
      file=here ("output","env_data.RData"))


rm(list=ls())



