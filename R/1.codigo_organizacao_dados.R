## organizacaoo de dados, e analises para o manuscrito "Diversity patterns and drivers" 

## chamar os pacotes
source("R/packages.R")
source("R/functions.R")

## CODIGO PARA CARREGAR OS DADOS PARA ANALISES DE DIVERSIDADE

## dados dos bentos
bentos <- read.xlsx(here("data","detection","Updated_compiled_quadrats_allsites.xlsx"),
                    sheet = 1, colNames = TRUE,detectDates=F)
## converter data em data - bug do pacote openxlsx
bentos$eventDate <-convertToDate(bentos$eventDate)

## dados dos peixes
peixes <- read.xlsx(here("data","detection","UpdatedData_RMorais_et_al_2017.xlsx"),
                    sheet = 1, colNames = TRUE,detectDates=F)
## converter data em data - bug do pacote openxlsx
peixes$eventDate <-convertToDate(peixes$eventDate)

## filtrar os dados de peixes de acordo com o minimo de ano de bentos
# peixes <- peixes [which(peixes$eventYear >= min (bentos$eventYear) & peixes$eventYear <= max (bentos$eventYear)),]

## modificar o eventID removendo o ano, desde que escolhemos um lapso temporal que cobre 2010 a 2014, 
## de modo que possamos analisar os dados se peixes foi coletado em 2012 e bentos em 2014, por exemplo 

bentos$eventID_MOD <- substr(bentos$eventID, 1,nchar(as.character(bentos$eventID))-5) 
peixes$eventID_MOD  <- substr(peixes$eventID, 1,nchar(as.character(peixes$eventID))-5) 


## fazer um histograma pra saber o numero de eventIDS por ano

barplot_function(df1=bentos,#bentos
                 df2= peixes #peixes
                 )

## se quiser o subset de sitios que correspondem em ambos os data sets
## subset entre os datasets, desde que peixes ou bentos foram amostrados nos mesmos locais
# peixes_subset <- peixes [which(peixes$eventID_MOD %in% bentos$eventID_MOD),]

## da mesma forma, pegar o subset de ID de bentos que estao nos dados de peixes
# bentos_subset <- bentos [which(bentos$eventID_MOD %in% peixes$eventID_MOD),]

##  criar uma ID numerica para o observador
peixes_subset$ID.observer <- as.numeric(as.factor(peixes_subset$Observer))

## ID para o sitio (peixes)

peixes_subset$locality_site <- paste(peixes_subset$Locality, 
                                     peixes_subset$Site,sep=".")


## ID para o sitio (bentos)

bentos_subset$locality_site <- paste(bentos_subset$Locality, 
                                     bentos_subset$Site,sep=".")

### mapa para conferir a posicao dos sitios - com base no subset dos dados que correspondem
## tem um monte de pontos porque tem todos os eventIDs ai

initial_map_function (df1 = bentos_subset,# red
                      df2 = peixes_subset) # black

#######################################################################################
#######################################################################################
################         dados de composicao                    #######################
#######################################################################################
#######################################################################################

## composicao de peixes
comp_peixes <-  cast(peixes_subset,
                     formula = locality_site  ~ ScientificName,
                     value= "IndCounting",
                     fun.aggregate = sum)

## composicao de bentos
comp_bentos <-  cast(bentos_subset,
                     formula = locality_site  ~ Taxon,
                     value= "Cover",
                     fun.aggregate = mean)

#######################################################################################
#######################################################################################
################         covariaveis de sitio                   #######################
#######################################################################################
#######################################################################################

## de sitio (aquelas que variam por sitio)

## sitios do nordeste
nord <- cast(peixes_subset,
             formula = locality_site  ~ Region,
             value= "IndCounting",
             fun.aggregate = sum)

## regiao 
ilhas <-  ifelse(nord [,3] >0,"oc.isl",0)
nor <-  ifelse(nord [,2] >0,"nord",0)
sud <-  ifelse(nord [,4] >0,"sud",0)
regiao <-  c(ilhas[which(ilhas !=0)],
             nor[which(nor !=0)],
             sud[which(sud !=0)])

## tipo de recife
tipo_recife <- cast(bentos_subset,
                    formula = locality_site  ~ Reef,
                    value= "Cover",
                    fun.aggregate = max)

recife_biog <- as.factor (ifelse (tipo_recife  [,2] >0, "1","0") )

## numero de transeccoes e videos (esforço)

lista_sitios <- unique(comp_peixes$locality_site)
# encontrando o ID das transeccoes unicas, e vendo seu length == numero de transectos
n_transeccoes <- lapply (lista_sitios, function (i)
  length (unique (peixes_subset [which(peixes_subset$locality_site == i),"Transect_id"]))
)
# encontrando numero de videos unicos por sitio
n_videos <- lapply (lista_sitios, function (i)
  length (unique (bentos_subset [which(bentos_subset$locality_site == i),"Video_number"]))
)

# montando uma tabela
tabela_esforco <- data.frame (locality_site = lista_sitios,
                              n_transectos_peixes = unlist(n_transeccoes),
                              n_videos_bentos = unlist(n_videos))
## lista de dados

dados_peixes_bentos <- list(peixes = comp_peixes,
                            bentos = comp_bentos)

## coordenadas geográficas para os mapas e analises

coordenadas <- aggregate(bentos_subset, 
                         by= list (bentos_subset$locality_site), 
                         FUN=mean)[c("Group.1","Lon","Lat")]

### BiO Oracle covariates
# Explore datasets in the package
layers <- list_layers()
#View (layers [grep ("Bio-ORACLE",layers$dataset_code),])

# Download specific layers to the current directory
# set prefered folder (to download data)
options(sdmpredictors_datadir=here ("data","environment"))

## chlorophil has different extent - loading and extracting in two steps         
layers_oracle <- load_layers(c("BO2_tempmean_ss","BO2_temprange_ss",
                               "BO2_ppmean_ss", "BO2_pprange_ss",
                               "BO2_salinitymean_ss", "BO2_salinityrange_ss"))

## these data have different extent
layers_oracle_Chl <- load_layers (c("BO2_chlomean_ss","BO2_chlorange_ss"))

## coordinates to spatial points
sp_points <- coordenadas
spdf <- SpatialPointsDataFrame(coords = sp_points[,2:3], data = sp_points,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

## extracting data

extracted_sea_data <- extract (layers_oracle, spdf,method='simple', fun=mean)
rownames(extracted_sea_data) <- sp_points$Group.1
extracted_sea_data_Chl <- extract (layers_oracle_Chl, spdf,method='simple', fun=mean)
rownames(extracted_sea_data_Chl ) <- sp_points$Group.1

## binding these dfs
extracted_sea_data <- cbind(extracted_sea_data,
                            extracted_sea_data_Chl)

## lista de covariaveis

covariates_site <- list (biog_reef = recife_biog,
                         region = regiao,
                         site_names = lista_sitios,
                         coord = coordenadas,
                         sea_data = extracted_sea_data)

covariates_effort <- list(effort = tabela_esforco)

###############################################
# # # # # # # # #  SAVE # # # # # # # # # # # #
###############################################

### save data - para os modelos de coral

save (covariates_site , ### dados de covariaveis das bases de dados
      covariates_effort, ## variaveis de esforco
      dados_peixes_bentos, ## dados de composicao 
      file=here ("output","data_drivers_analysis.RData"))



