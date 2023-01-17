
# ------------------------------------------------------------------- #
#         Organizing data for the MS "Diversity patterns and drivers" 
#-------------------------------------------------------------------- #

## load packages and functions
source("R/packages.R")
source("R/functions.R")


# Load data

# ------------------------------ #
#        Benthic data
# ------------------------------ #





## ============================================
## benthic snapshots (Aued et al. 2018)






benthos_DF_eMOF_aued <- read.csv(here ("data",
                                       "detection",
                                       "AAued_spatialData",
                                       "DF_eMOF.csv"),sep=",",
                                 #encoding= "UTF-8",
                                 row.names=NULL)

# event core
benthos_event_core_aued <-  read.csv(here ("data",
                                           "detection",
                                           "AAued_spatialData",
                                           "event_core.csv"),sep=","#, 
                                     #encoding= "UTF-8"
                                     )




# matching event IDs to find site and locality (variables_we_want)
benthos_SN_data_aued <- benthos_event_core_aued [match (benthos_DF_eMOF_aued$eventID,
                                                        benthos_event_core_aued$eventID),]
# bind the occurrence data
benthos_SN_data_aued<- cbind (benthos_SN_data_aued,
                              benthos_DF_eMOF_aued)

# define classes of depth
benthos_SN_data_aued$cat_depth <- ifelse (benthos_SN_data_aued$maximumDepthinMeters<=7,
                                         "shallow",
                                         "deep")

# create sites by binding site, locality, and depth
benthos_SN_data_aued$site_analysis <- paste (benthos_SN_data_aued$site,
                                             benthos_SN_data_aued$locality,sep = "-")


## proportion of not identified taxa (family was considered identified)

rowSums(table (benthos_SN_data_aued$family, benthos_SN_data_aued$scientificName)>0)
table(is.na(benthos_SN_data_aued$family))

## proportion of not identified taxa (family was considered identified)

table (
    rowSums(
      
      table (benthos_SN_data_aued$verbatimIdentification, # family per each identified taxa
                  benthos_SN_data_aued$family)
      >0) # if identified at family level, the sum will be > 0
) / length(unique(benthos_SN_data_aued$verbatimIdentification)) # all taxa in the dataset


require(dplyr)

check<-benthos_SN_data_aued[,-1]  %>% 
  
    filter (site %in% c("ilhasc_norte", "ilhasc_sul") 
                                  ) %>%
  group_by(site,locality,verbatimIdentification) %>%
  summarize(siz=sum(measurementValue,na.rm=T)) %>%
  filter (siz>0) 



#-----------------------------------------------------------------------------
# Fish snapshots (Morais et al. 2017)





fish_DF_eMOF_morais <- read.csv(here ("data",
                                      "detection",
                                      "RMorais_spatialData",
                                      "DF_eMOF.csv"),sep=",",
                                encoding= "UTF-8",
                                row.names=NULL)

# event core
fish_event_core_morais <-  read.csv(here ("data",
                                          "detection",
                                          "RMorais_spatialData",
                                          "event_core.csv"),sep=",", 
                                    encoding= "UTF-8")



# matching event IDs to find site and locality (variables_we_want)
fish_SN_data_morais <- fish_event_core_morais [match (fish_DF_eMOF_morais$eventID,
                                                      fish_event_core_morais$eventID),]
# bind the occurrence data
fish_SN_data_morais<- cbind (fish_SN_data_morais,
                             fish_DF_eMOF_morais)

# filtering abundance data
fish_SN_data_morais<-fish_SN_data_morais[which(fish_SN_data_morais$measurementType == "abundance"),]

# define classes of depth
fish_SN_data_morais$cat_depth <- ifelse (fish_SN_data_morais$maximumDepthinMeters<=7,
                                         "shallow",
                                          "deep")

# create sites by binding site, locality, and depth
fish_SN_data_morais$site_analysis <- paste (fish_SN_data_morais$site,
                                            fish_SN_data_morais$locality,sep = "-")



# filtering sites with data in both datasets
dataset_fish <- fish_SN_data_morais[which(fish_SN_data_morais$site_analysis %in% 
                                            benthos_SN_data_aued$site_analysis),]
dataset_benthos <- benthos_SN_data_aued[which(benthos_SN_data_aued$site_analysis %in% 
                                                fish_SN_data_morais$site_analysis),]


# check if sites are in both datasets
(unique(dataset_benthos$site_analysis) %in% unique(dataset_fish$site_analysis))


# total number of belt transects
length(unique(dataset_fish$eventID))

# average
mean (apply(table (dataset_fish$eventID, 
                   dataset_fish$site_analysis)>0,
            2,sum))
# sd
sd (apply(table (dataset_fish$eventID, 
                 dataset_fish$site_analysis)>0,
          2,sum))

# total number of plots
length(unique(dataset_benthos$eventID))
# average
mean (apply(table (dataset_benthos$eventID, 
           dataset_benthos$site_analysis)>0,
      2,sum))
# sd
sd (apply(table (dataset_benthos$eventID, 
                   dataset_benthos$site_analysis)>0,
            2,sum))



#---------------------------------------------------#
#             composition data                      #
#---------------------------------------------------#
## fish & benthic composition (site x spp)

comp_fish <-  cast(dataset_fish, # fish
                     formula = site_analysis  ~ scientificName,
                     value= "measurementValue",
                     fun.aggregate = sum,
                     drop =F)
comp_fish <- comp_fish[order (comp_fish$site_analysis,decreasing=F),] # order

#
comp_benthos <-  cast(dataset_benthos, # benthos
                     formula = site_analysis  ~ verbatimIdentification,
                     value= "measurementValue",
                     fun.aggregate = mean,
                     drop =F)
comp_benthos <- comp_benthos[order (comp_benthos$site_analysis,decreasing=F),] # order


# define sites for all analyses
comp_fish$site_analysis ==  comp_benthos$site_analysis
sites <- comp_fish$site_analysis

# relative counts of fish
comp_fish<- comp_fish[,-1]
comp_fish <- comp_fish/rowSums (comp_fish)
# rm first col of benthic composition (site names)
comp_benthos<- comp_benthos[,-1]


# -----------------------------------
#   EFFORT PREDICTORS
# -----------------------------------
#fish
effort_fish <-  cast(dataset_fish, # fish
                   formula = site_analysis  ~ eventID,
                   value= "measurementValue",
                   fun.aggregate = sum,
                   drop =F)[,-1]
# benthos
effort_benthos <-  cast(dataset_benthos, # fish
                     formula = site_analysis  ~ eventID,
                     value= "measurementValue",
                     fun.aggregate = sum,
                     drop =F)[,-1]

# effort
effort_dataframe <- data.frame (sites = sites ,
                     fish_effort = rowSums(effort_fish>0),
                     benthos_effort = rowSums(effort_benthos>0))




# -----------------------------------
#   ENVIRONMENTAL PREDICTORS
# -----------------------------------


# average depth
site_depth <- tapply (dataset_fish$maximumDepthinMeters,
                      list(dataset_fish$site_analysis),
                      mean)

# region
site_region <- tapply (dataset_fish$higherGeography,
                      list(dataset_fish$site_analysis),
                      getmode)


site_region <- ifelse (sites %in% c(
  "abrolhos-chapeirao",
  "abrolhos-portinho_norte" ,
  "abrolhos-siriba" ,
  "btds_santos-farol_da_barra",
  "btds_santos-frades" ,                     
  "btds_santos-pedra_cardinal",                 
  "btds_santos-poste_quatro"  ,                 
  "costa_corais-barra_da_gale" ,             
  "costa_corais-gales"      ,                
  "costa_corais-taocas",
  "manuel_luis-ana_cristina",
  "rgnor_natal-batente_das_agulhas",            
  "rgnor_natal-pedra_do_silva",                 
  "rgnor_parrachos-maracajau",               
  "rgnor_parrachos-parrachos_de_rio_do_fogo"),
  "BrazilianCoastNorth",
  site_region)

# southeast and south sites
site_region <- ifelse (site_region == "BrazilianCoast",
                       "BrazilianCoastSouth", 
                       site_region)



## geo coordinates
coordinates_sites <- dataset_fish [,-which(colnames(dataset_fish) %in% c("eventID", "X"))] %>% 
  
  group_by(site_analysis) %>%

    summarise(decimalLatitude = mean(decimalLatitude),
            decimalLongitude = mean(decimalLongitude))

# adjusting some coordinates
coordinates_sites$decimalLatitude [which(coordinates_sites$site_analysis == "btds_santos-farol_da_barra")] <- -13.017295
coordinates_sites$decimalLongitude [which(coordinates_sites$site_analysis == "btds_santos-farol_da_barra")] <- -38.544638
coordinates_sites$decimalLatitude [grep("ilha_das_cabras", coordinates_sites$site_analysis)] <- -23.824199
coordinates_sites$decimalLongitude [grep("ilha_das_cabras", coordinates_sites$site_analysis)] <- -45.397819


# -----------------------
# NOAA Ocean Data




# parameters
SSTstartDate <- "2012-01-01"  ## define start date of your time series 

## set dataset source (monthly SST)
## the list of datasets is here
## https://coastwatch.pfeg.noaa.gov/erddap/search/index.html?page=1&itemsPerPage=1000&searchFor=griddap
# https://coastwatch.pfeg.noaa.gov/erddap/griddap/index.html?page=1&itemsPerPage=1000


# SST
SSTsource <- info("jplMURSST41mday")

# create a cluster of 'ncores', 
ncores <- 5
cl <- makeCluster (ncores)
# load data and functions in each core 
clusterExport(cl, c("SSTsource",
                    "SSTstartDate",
                    "coordinates_sites"))
# load packages in each core
clusterEvalQ(cl,library("rerddap"))

## Get sst 
SST <- parLapply (cl, seq (1,nrow (coordinates_sites)), function (i) {
  
  tryCatch(
    
    griddap(SSTsource, 
            time=c(SSTstartDate, "last"),
            longitude = c(coordinates_sites$decimalLongitude[i],
                          coordinates_sites$decimalLongitude[i]),
            latitude = c(coordinates_sites$decimalLatitude[i], 
                         coordinates_sites$decimalLatitude[i]), 
            fields = "sst", # to have the field you need to go to the graph option in the website (https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdVH2018r6718day.graph) and going to "color"
            fmt = "csv"),
    error = function(e) return ("NULL"))}#{print(e); print("retrying...")}
)
stopCluster (cl)




# turbidity
DEPTHsource <- info("erdMH1kd4901day")
cl <- makeCluster (ncores) #cluster
# load data and functions in each core 
clusterExport(cl, c("DEPTHsource", 
                    "SSTstartDate",
                    "coordinates_sites"))
# load packages in each core
clusterEvalQ(cl,library("rerddap"))

## Get depth
DEPTH <- parLapply (cl, seq (1,nrow (coordinates_sites)), function (i) {
  
  tryCatch(
    
    griddap(DEPTHsource, 
            time=c(SSTstartDate, "last"),
            longitude = c(coordinates_sites$decimalLongitude[i],
                          coordinates_sites$decimalLongitude[i]),
            latitude = c(coordinates_sites$decimalLatitude[i], 
                         coordinates_sites$decimalLatitude[i]), 
            fields = "k490", # to have the field you need to go to the graph option in the website (https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdVH2018r6718day.graph) and going to "color"
            fmt = "csv"),
    error = function(e) return ("NULL"))}#{print(e); print("retrying...")}
)
stopCluster (cl) # stop cluster





# productivity (Chlorophyll-a)
PRODsource <- info("erdMH1chla8day")
cl <- makeCluster (ncores) # cluster

# load data and functions in each core 
clusterExport(cl, c("PRODsource", 
                    "SSTstartDate",
                    "coordinates_sites"))
# load packages in each core
clusterEvalQ(cl,library("rerddap"))

## Get productivity
PRODUCTIVITY <- parLapply (cl, seq (1,nrow (coordinates_sites)), function (i) {
  
  tryCatch(
    
    griddap(PRODsource, 
            time=c(SSTstartDate, "last"),
            longitude = c(coordinates_sites$decimalLongitude[i],
                          coordinates_sites$decimalLongitude[i]),
            latitude = c(coordinates_sites$decimalLatitude[i], 
                         coordinates_sites$decimalLatitude[i]), 
            fields = "chlorophyll", # to have the field you need to go to the graph option in the website (https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdVH2018r6718day.graph) and going to "color"
            fmt = "csv"),
    error = function(e) return ("NULL"))}#{print(e); print("retrying...")}
)

stopCluster (cl)



# salinity source
SALINsource <- info("jplAquariusSSSDailyV5")

##
cl <- makeCluster (ncores)

## Get sst cl <- makeCluster (ncores)
# load data and functions in each core 
clusterExport(cl, c("SALINsource", 
                    "SSTstartDate",
                    "coordinates_sites"))
# load packages in each core
clusterEvalQ(cl,library("rerddap"))

# run
## Get productivity

SALINITY <- parLapply (cl, seq (1,nrow (coordinates_sites)), function (i) {
  
  tryCatch(
    
    griddap(SALINsource, 
            time=c(SSTstartDate, "last"),
            longitude = c(coordinates_sites$decimalLongitude[i],
                          coordinates_sites$decimalLongitude[i]),
            latitude = c(coordinates_sites$decimalLatitude[i], 
                         coordinates_sites$decimalLatitude[i]), 
            fields = "sss", # to have the field you need to go to the graph option in the website (https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdVH2018r6718day.graph) and going to "color"
            fmt = "csv"),
    error = function(e) return ("NULL"))}#{print(e); print("retrying...")}
)

stopCluster (cl)



## summarized  variables
# sst
SST_site <- unlist (lapply (SST, function (i) 
  
  i %>% summarise(sst = mean(sst,na.rm=T))
  
))
# turbidity
kd490_site <- unlist (lapply (DEPTH, function (i) 
  
  i %>% summarise(k490 = mean(k490,na.rm=T))
  
))
# PROD
prod_site <- unlist (lapply (PRODUCTIVITY, function (i) 
  
  i %>% summarise(chlorophyll = mean(chlorophyll,na.rm=T))
  
))
# salinity
salinity_site <- unlist (lapply (SALINITY, function (i) 
  
  i %>% summarise(sss = mean(sss,na.rm=T))
  
))


# site covs
site_covs <- data.frame (sites = sites,
            sst = SST_site,
            turbidity = kd490_site,
            productivity = prod_site,
            salinity = salinity_site,
            depth = site_depth)




# -----------------------

# try BioOracle too


# BiO Oracle - extracting covariate data
# Explore datasets in the package
# devtools::install_github("lifewatch/sdmpredictors")
layers <- list_layers()
#View (layers [grep ("Bio-ORACLE",layers$dataset_code),])

# Download specific layers to the current directory
# set prefered folder (to download data)
options(sdmpredictors_datadir=here ("data","environment"))

## chlorophil has different extent - loading and extracting in two steps         
#layers_oracle <- load_layers(c("BO2_tempmean_ss",
#                               #"BO2_temprange_ss",
#                               "BO2_ppmean_ss", 
#                               #"BO2_pprange_ss",
#                               "BO2_salinitymean_ss", 
#                               #"BO2_salinityrange_ss",
#                               #"BO_damax",
#                               "BO_damean"
#                               #,"BO_damin"
#))

biooracle_data <-  (list.files (here ("data","environment"),pattern = ".tif"))
biooracle_data <- lapply (biooracle_data, function (i) 
  
    raster (here ("data","environment",i))
)
biooracle_data<- stack (biooracle_data)#stack

# coord to sppoints df to extract
spdf <- SpatialPointsDataFrame(coords = coordinates_sites[,3:2], data = coordinates_sites,
                         proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

## extracting data

extracted_sea_data <- extract (biooracle_data, spdf,method='bilinear', 
                               fun=mean)
rownames (extracted_sea_data) <- sites

# correlation between datasets
cor (extracted_sea_data, 
     site_covs[,-1],
     use = "complete.obs")

# it's ok to use biooracle data (that is more complete)


# ------------------------------
#  distance offshore
# BR coastline, download from here https://mapcruzin.com/free-brazil-arcgis-maps-shapefiles.htm

BR <- readOGR(dsn=here("data", "environment","brazil-coastline"), "brazil_coastline")
crs(BR) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
BR <- spTransform(BR, CRS("+init=epsg:4326"))


# use dist2Line from geosphere - only works for WGS84 
sp_data <- spTransform(spdf, CRS("+init=epsg:4326"))


# measuring the distance
dist_bentos <- geosphere::dist2Line(p = sp_data, 
                             line = (BR))


# binding coords
coordinates_sites <- cbind (coordinates_sites, dist_bentos)


# plot
ggplot(data = world) +
  geom_sf() +
  geom_point(data=coordinates_sites, 
             aes(x=decimalLongitude  , y=decimalLatitude , col =distance))# +
  #coord_sf(xlim = c(-55, -20), ylim = c(-33,0 ), expand = FALSE)



# ====================== 
# magris data (reef area)

shapefiles<-list.files(here ("data","environment","magris_reef_map"),pattern=".shp")
shapefiles<-shapefiles[-grep ("cumulative",shapefiles)]
shapefiles<-shapefiles[-grep ("Threatened",shapefiles)]

# repeat noronha to have rocas atol
shapefiles <- c(shapefiles,"FN.shp")
shapefiles<-shapefiles[order(shapefiles)]

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
                      "noronha" = "FS13", # noronha
                      "atol" = "FS12", #rocas
                      "northeastern" = c("NC11",# northeastern
                                         "NC12",
                                         "NC13"),
                      "riogrande"="RC11", # rio grande
                      "southeastern"="SC11", #southeastern
                      "trindade" = "TS12" # trindade
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
                 BR_reefs[[5]],
                 BR_reefs[[6]],
                 BR_reefs[[7]],
                 BR_reefs[[8]])

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
BR_reefs_lambert <- spTransform(BR_reefs, 
                        "+proj=laea +lat_0=0 +lon_0=0 +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
to_raster <-  raster (extent(BR_reefs_lambert),
                      crs = "+proj=laea +lat_0=0 +lon_0=0 +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
                      res=55000)
grd_raster <- projectRaster(grd_raster, to_raster)

# southe maerica map
southAme<- readOGR(dsn= here("data","environment","South_America"),encoding="latin1", 
                   layer="South_America")
southAme <- spTransform(southAme, 
                        "+proj=laea +lat_0=0 +lon_0=0 +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# extract
require("exactextractr")

# extract
prop_area <- coverage_fraction(grd_raster, 
                               st_combine(st_as_sf(BR_reefs_lambert)))[[1]]
# points
pts <- spTransform(spdf, 
                   "+proj=laea +lat_0=0 +lon_0=0 +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# plot
plot1 <- gplot(prop_area) +  
  geom_tile(aes(fill=(value)),alpha=0.8) +
  scale_fill_viridis_c(begin = 0.1, 
                      end = 1,
                      direction=-1,
                      option = "magma",
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
plot2 + geom_point(data = data.frame (coordinates(pts)), 
                   aes (x=decimalLongitude,
                                           y=decimalLatitude),
                   shape=1,size =5,col = "gray")
dev.off()

# values in site coords
extracted_area <- extract (prop_area, 
                           pts,
                           method='simple', fun=mean)



# all covariate data into a dataframe

site_covs <- data.frame (sites = sites,
                         region = site_region,
                         sst = extracted_sea_data[,'BO2_tempmean_ss_lonlat'],
                         turbidity = extracted_sea_data[,'BO_damean_lonlat'],
                         productivity = extracted_sea_data[,'BO2_ppmean_ss_lonlat'],
                         salinity = extracted_sea_data[,'BO2_salinitymean_ss_lonlat'],
                         depth = site_depth,
                         reef_area = extracted_area,
                         offshore_distance = coordinates_sites$distance,
                         decimalLatitude = coordinates_sites$decimalLatitude,
                         decimalLongitude = coordinates_sites$decimalLongitude)


# SAVE

save (site_covs , ### site covariates
      effort_dataframe, ## effort
      comp_fish, ## fish composition
      comp_benthos, ## benthic composition
      
      file=here ("data","modeling_data.RData"))


rm(list=ls())

