
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
                                 encoding= "latin1",
                                 row.names=NULL)

# event core
benthos_event_core_aued <-  read.csv(here ("data",
                                           "detection",
                                           "AAued_spatialData",
                                           "event_core.csv"),sep=",", 
                                     encoding= "latin1"
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
#coordinates_sites$decimalLatitude [which(coordinates_sites$site_analysis == "btds_santos-farol_da_barra")] <- -13.017295
#coordinates_sites$decimalLongitude [which(coordinates_sites$site_analysis == "btds_santos-farol_da_barra")] <- -38.544638
#coordinates_sites$decimalLatitude [grep("ilha_das_cabras", coordinates_sites$site_analysis)] <- -23.824199
#coordinates_sites$decimalLongitude [grep("ilha_das_cabras", coordinates_sites$site_analysis)] <- -45.397819


# -----------------------

#  BioOracle 
# codes to download and save in the folder "environment"
# BiO Oracle - extracting covariate data
# Explore datasets in the package
# devtools::install_github("lifewatch/sdmpredictors")
# layers <- list_layers()
# View (layers [grep ("Bio-ORACLE",layers$dataset_code),])
# Download specific layers to the current directory
# set prefered folder (to download data)
# options(sdmpredictors_datadir=here ("data","environment"))
## chlorophil has different extent - loading and extracting in two steps         
# layers_oracle <- load_layers(c("BO2_tempmean_ss",
#                               "BO2_ppmean_ss", 
#                               "BO2_salinitymean_ss", 
#                               "BO_damean"
#                              ))

biooracle_data <-  (list.files (here ("data","environment"),pattern = ".tif"))
biooracle_data <- lapply (biooracle_data, function (i) 
  
    raster (here ("data","environment",i))
)
biooracle_data<- stack (biooracle_data)#stack

# coord to sppoints df to extract
spdf <- SpatialPointsDataFrame(coords = coordinates_sites[,3:2], data = coordinates_sites,
                         proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

## extracting data

extracted_sea_data <- raster::extract (biooracle_data, 
                                       spdf,method='bilinear', 
                               fun=mean)
rownames (extracted_sea_data) <- sites


# ------------------------------
#  distance offshore
# BR coastline, download from here https://mapcruzin.com/free-brazil-arcgis-maps-shapefiles.htm

#BR <- readOGR(dsn=here("data", "environment","brazil-coastline"), "brazil_coastline")
BR <- st_read(("data/environment/brazil-coastline/brazil_coastline.shp"))
#crs(BR) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
#BR <- spTransform(BR, CRS("+init=epsg:4326"))
st_crs(BR) <- st_crs(4326) # assign crs

# use dist2Line from geosphere - only works for WGS84 
sp_data <- st_as_sf(spdf, coords = c("decimalLongitude", "decimalLatitude"))
st_crs(sp_data) <- st_crs(4326) # assign crs

# measuring the distance
#dist_bentos <- st_distance(x = sp_data, 
#                             y = (BR))
# help : https://gis.stackexchange.com/questions/243994/how-to-calculate-distance-from-point-to-linestring-in-r-using-sf-library-and-g
dist_bentos <- geosphere::dist2Line(p = st_coordinates(sp_data), 
                     line = st_coordinates(BR)[,1:2])

# binding coords
coordinates_sites <- cbind (coordinates_sites, dist_bentos)


# plot
ggplot(data = world) +
  geom_sf() +
  geom_point(data=coordinates_sites, 
             aes(x=decimalLongitude  , y=decimalLatitude , col =distance))# +
  #coord_sf(xlim = c(-55, -20), ylim = c(-33,0 ), expand = FALSE)



# all covariate data into a dataframe

site_covs <- data.frame (sites = sites,
                         region = site_region,
                         sst = extracted_sea_data[,'BO2_tempmean_ss_lonlat'],
                         turbidity = extracted_sea_data[,'BO_damean_lonlat'],
                         productivity = extracted_sea_data[,'BO2_ppmean_ss_lonlat'],
                         salinity = extracted_sea_data[,'BO2_salinitymean_ss_lonlat'],
                         depth = site_depth,
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

