## pacotes necessários

require(here) ## para transitar entre pastas
require(reshape) ## para funcao tabela dinamica e melt (formato longo)
require(reshape2) ## para funcao tabela dinamica e melt (formato longo)

## pacotes para mapas
require(rnaturalearth)
require(rnaturalearthdata)
require(ggplot2)
require(gridExtra)
require(ggrepel)
require (scatterpie)
#require(sf)
require(rgeos)

## pacotes para carregar os dados
require(xlsx)
require(openxlsx)
require(dplyr)

## pacote para colar dados em arrays
require(abind)

### 
## pacotes para vizinhanca
require(rgdal)
require(raster)
require(rgeos)
require (spdep) 
require(maps)

## pacote para processamento paralelo
require(parallel)

# funcoes diversas (padronizacao)
require(vegan)
require(FD)
require(cluster)
require(clue)

# 
require("sdmpredictors")
require("leaflet")
library(mgcv)
require("gamm4")
require(nlme)
require(emmeans)

require(betareg)
require(MASS)
require(effects)

require(brms)
require(loo)
require (corrplot)

# plot
require(viridis)
require(rasterVis)

# maps and noaa data

# load packages
library(readr)
library(rerddap)
require(parallel)
library(lubridate)
library(dplyr)
library(flexdashboard)
library(reshape2)
library(leaflet)
library(ggplot2)
library(vegan)
library(xts)
library(dygraphs)
library(plotly)
library(mapdata)
library(RColorBrewer)
palette(brewer.pal(8, "Set2"))
