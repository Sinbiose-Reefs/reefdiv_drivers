
## --------------------
#       FIGURES
# --------------------

## call packages
source("R/packages.R")
source("R/functions.R")

### -------------------------
#  FIG - supporting information
# correlations between variables
load(here("output_no_resampling", "data_to_modeling_GLM.RData"))
corrplot (
  cor (cov_fish [,c("BO2_tempmean_ss_std","BO2_ppmean_ss_std","BO2_salinitymean_ss_std",
                    "BO_damean_std","distance_std","area","SamplingArea_std")]), 
  p.mat = cor.mtest(cor (cov_fish [,c("BO2_tempmean_ss_std","BO2_ppmean_ss_std","BO2_salinitymean_ss_std",
                                      "BO_damean_std","distance_std","area","SamplingArea_std")]))$p)


# metrics
cov_fish[,c("Group.1","SR","FRic","FEve","FDiv")]
cov_benthos[,c("Group.1","SR","FRic","FEve","FDiv")]

# -----------------------------------------------------------

# load estimates 
load(here("output_no_resampling", "FD_obs_data_fish.RData"))
load(here("output_no_resampling", "FD_obs_data_benthos.RData"))

# Fig. S1.1
# Load rarefied data of fishes and benthos, and other data
# sites x transect / video x species

load (here ("output","random_composition_bentos.RData"))
load (here ("output","random_composition_fish.RData"))

# ------------------------------------------ #
#
#     Load environment data

load (here ("output","env_data.RData"))


# ------------------------------------------- #
# map of distribution of samples and sampling effort
# fish
ObsRich <- rowSums (dados_peixes_bentos$peixes[,-1] >0)
eff_data <- data.frame (Sites=covariates_site$site_names,
                        Lat=covariates_site$coord$coord_peixes[,"Lat"],
                        Lon=covariates_site$coord$coord_peixes[,"Lon"],
                        Reg =covariates_site$regiao,
                        UVC=covariates_effort$effort$fish$n_transectos_peixes,
                        ObsRich = ObsRich)
# benthos
ObsRichBentos <- rowSums (dados_peixes_bentos$bentos[,-1] >0)
eff_data_bentos <- data.frame (Sites = rownames(covariates_site$sea_data[[2]]),
                               Lat=covariates_site$coord$coord_bentos[,"Lat"],
                               Lon=covariates_site$coord$coord_bentos[,"Lon"],
                               UVC=covariates_effort$effort$benthos$n_videos_bentos,
                               ObsRich = ObsRichBentos)

## plot fish variation in obs richness and effort of fish DS

# scaling factor for UVC data
scaleFactor <- max(eff_data$ObsRich) / max(eff_data$UVC)

# plot 
# fish

eff_obsrich_plot <- ggplot(eff_data, aes(x=Lat)) +
  geom_point(aes(y=ObsRich),
             alpha=0.2,col="#0e49b5") +
  geom_smooth(aes(y=ObsRich),fill="#0e49b5",alpha=0.1, 
              method = "lm", 
              formula = y ~ poly(x, 5), size = 1,col ="#0e49b5") +
  geom_point(aes(y=UVC*scaleFactor),
             alpha=0.2,col="black")+
  geom_smooth(aes(y=UVC*scaleFactor), 
              method = "lm", 
              formula = y ~ poly(x, 5), size = 1,col="black") +
  scale_y_continuous(name="Observed Richness", sec.axis=sec_axis(~./scaleFactor,name="Number of transects")) +
  theme_classic() +
  theme (axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         axis.text.x = element_text(size=8),
         axis.title.x = element_text(size=10),
         plot.margin = unit(c(0,-0.0,0,0.3), "cm")
  ) +
  coord_flip()

eff_obsrich_plot

#  benthos
scaleFactorB <- max(eff_data_bentos$ObsRich) / max(eff_data_bentos$UVC)

eff_obsrich_plot_benthos <- ggplot(eff_data_bentos, aes(x=Lat)) +
  geom_point(aes(y=ObsRich),
             alpha=0.2,col="#f2a154") +
  geom_smooth(aes(y=ObsRich),fill="#f2a154",alpha=0.1, 
              method = "lm", formula = y ~ poly(x, 5), size = 1,col ="#f2a154") +
  geom_point(aes(y=UVC*scaleFactorB),
             alpha=0.2,col="black")+
  geom_smooth(aes(y=UVC*scaleFactorB), 
              method = "lm", formula = y ~ poly(x, 5), size = 1,col="black") +
  scale_y_continuous(name="Observed Richness", sec.axis=sec_axis(~./scaleFactorB,name="Number of photo quadrats")) +
  theme_classic() +
  theme (axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         axis.text.x = element_text(size=8),
         axis.title.x = element_text(size=10),
         plot.margin = unit(c(0,0.3,0,0.3), "cm")
  ) +
  coord_flip()


eff_obsrich_plot_benthos

# temperature

data_lat_temp <- data.frame(Lat=covariates_site$coord$coord_bentos$Lat,
                            SST =covariates_site$sea_data[[1]][,"BO2_tempmean_ss"])
# plot
SST <- ggplot(data_lat_temp, aes(x=Lat,y=SST)) +
  
  geom_point(alpha=0.2,col="#8E0505") +
  
  geom_smooth(method = "glm", 
              formula = y ~ poly(x, 2), 
              size = 1,col ="#8E0505",
              fill="#8E0505",alpha=0.1) +
  
  scale_y_continuous(name="SST",sec.axis=sec_axis(~./scaleFactorB,name="Number of photo quadrats")) +
  theme_classic() +
  theme (axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         axis.text.x = element_text(size=8),
         axis.title.x = element_text(size=10),
         plot.margin = unit(c(0,0.3,0,0.3), "cm")
  ) +
  coord_flip()

SST

## maps
# mapa mundi
world <- ne_countries(scale = "medium", returnclass = "sf")

# cortar o mapa para ver a america do Sul e parte da central
wm <- ggplot() + 
  geom_sf (data=world, size = 0.1, 
           fill= "#aaaaaa",colour="#aaaaaa") +
  coord_sf (xlim = c(-50, -30),  ylim = c(-30, 4), expand = FALSE) +
  theme_bw() + #xlab ("Longitude")  + ylab ("Latitude") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#f4f9f9",#darkslategray1
                                        colour = "#f4f9f9"),
        axis.text.x = element_text(size=8),
        axis.ticks.x=element_line(size=1),
        axis.text.y = element_text(size=8),
        axis.ticks.y=element_line(size=1),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        title = element_blank(),
        plot.margin = unit(c(0,-0.8,0,0.3), "cm")) +
  xlab("Longitude") + ylab("Latitude")

map_peixes <- wm + geom_point(data=eff_data, aes(x=Lon, y=Lat),size=2,
                              col="#0e49b5")
map_peixes_bentos <- map_peixes + 
  geom_point (data=eff_data_bentos,aes(x=Lon+1.5, y=Lat+0.5),size=2,
              col="#f2a154")

# -------------------------------------------------------
# select of sites to plot (extremes in index values)
site_to_plot_alcatrazes <-c("alcatrazes.portinho_sudoeste.raso","alcatrazes.portinho_sudoeste.fundo")
site_to_plot_sc <- c("ilhasc_norte.deserta_norte.raso","ilhasc_norte.deserta_norte.fundo")
site_to_plot_ccorais <- c("btds_santos.farol_da_barra.raso","btds_santos.poste_quatro.fundo")
site_to_plot_es <- c("espirito_santo.escalvada.raso","espirito_santo.escalvada.fundo")

# now produce data for  polygons in the first plot  

df_plot<- lapply (list(site_to_plot_ccorais,
                       site_to_plot_alcatrazes,
                       site_to_plot_sc,
                       site_to_plot_es), function (t) { # over sites to choose
                         
                         
                           # select sites
                           site_raso <- which(sites_fish_complete == t[1]) # [2] always shallow
                           site_fundo <- which(sites_fish_complete == t[2]) # [1] always deep
                           
                           # the complete trait space per rarefaction dataset
                           data_to_space <- cbind(FD_obs_data_fish$axesPCO[,1:2],ext = F) # the complete space
                           complete_space <- data_to_space [chull(data_to_space[,1:2], y = NULL),] # convex hull on complete space
                           
                           # composition of shallow sites
                           s.raso<-colnames(comp_peixes[site_raso,][which(comp_peixes[site_raso,]>0)])
                           # composition of deep sites
                           s.fundo<-colnames(comp_peixes[site_fundo,][which(comp_peixes[site_fundo,]>0)])
                           
                           # subspace and convex hull of each rarefied data set
                           scores_subset_raso<-data_to_space [rownames(data_to_space)%in% (s.raso),] # subsetting species present into the sites
                           chull_site_raso <- (scores_subset_raso[chull(scores_subset_raso),])# Chull of this subset
                           scores_subset_fundo<-data_to_space [rownames(data_to_space)%in% (s.fundo),] # subsetting species present into the sites
                           chull_site_fundo <- (scores_subset_fundo[chull(scores_subset_fundo),])# Chull of this subset
                           
                           # create df to plot
                           df_plot <- rbind (data.frame (complete_space,data="complete",site= sites_fish_complete[site_fundo]),
                                             data.frame (chull_site_fundo,data="fundo",site= sites_fish_complete[site_fundo]),
                                             data.frame(chull_site_raso,data="raso",site= sites_fish_complete[site_fundo]))
                           
                           ;
                           df_plot
                           
                         })

# melt this df
df_plot<-do.call(rbind,df_plot)

# species from shallow and deep areas
comp_fish_deep <- comp_peixes[grep ("fundo", comp_peixes$locality_site),]
comp_fish_deep<-colnames(comp_fish_deep [,which(colSums(comp_fish_deep[,-1]) >0)])
comp_fish_shallow <- comp_peixes[grep ("raso", comp_peixes$locality_site),]
comp_fish_shallow<- colnames(comp_fish_shallow [,which(colSums(comp_fish_shallow[,-1]) >0)])
# find the species in both depths
fish_spp <- data.frame (species = colnames(comp_peixes)[-1],
                           depth = NA)
# depth
fish_spp[which(fish_spp$species %in% comp_fish_deep == T & fish_spp$species %in% comp_fish_shallow == F),"depth"] <- "fundo"
fish_spp[which(fish_spp$species %in% comp_fish_shallow == T & fish_spp$species %in% comp_fish_deep == F),"depth"] <- "raso"
fish_spp[which(fish_spp$species %in% comp_fish_shallow == T & fish_spp$species %in% comp_fish_deep == T),"depth"] <- "complete"
fish_spp$depth [is.na(fish_spp$depth)] <- "complete"
# ordering
FD_obs_data_fish$axesPCO <- FD_obs_data_fish$axesPCO [order(rownames(FD_obs_data_fish$axesPCO)),]
fish_spp<-fish_spp [which(fish_spp$species %in% rownames(FD_obs_data_fish$axesPCO)),]
fish_spp$species == rownames(FD_obs_data_fish$axesPCO)
# bind depth
FD_obs_data_fish$axesPCO <- cbind (FD_obs_data_fish$axesPCO,
                                      depth=fish_spp$depth)
FD_obs_data_fish$axesPCO$depth <- factor (FD_obs_data_fish$axesPCO$depth,
                                             levels = c("complete",
                                                        "fundo",
                                                        "raso"))
df_plot$data<-factor (df_plot$data,
                              levels = c("complete",
                                         "fundo",
                                         "raso"))

# index per depth

# results to annotate into the plot (average indexes across rarefaction samples)
indexes_per_sample_site_fish <- lapply (list(site_to_plot_ccorais,
                                             site_to_plot_alcatrazes,
                                             site_to_plot_sc,
                                             site_to_plot_es), function (t) { # over sites to choose
                                               # deep
                                               deep <- FD_obs_data_fish$Fdindexes[which(sites_fish_complete == t[1]),]
                                              shallow <- FD_obs_data_fish$Fdindexes[which(sites_fish_complete == t[2]),]
                                                 res <- (list(shallow=shallow,
                                                              deep=deep))
                                                 ;# return
                                                 res
                                               })
                                             
av_shallow_per_site_fish<-round(do.call(rbind,sapply (indexes_per_sample_site_fish, "[","shallow")),2)
# adjust names
rownames(av_shallow_per_site_fish)<-c(site_to_plot_ccorais[1],
                                      site_to_plot_alcatrazes[1],
                                      site_to_plot_sc[1],
                                      site_to_plot_es[1])
# SR

cbind(
  sites_fish_complete,
      rowSums(ifelse (comp_peixes > 0,1,0))
)

round(FD_obs_data_fish$explanAxes,3)
round(FD_obs_data_benthos$explanAxes,3)

# create a plot and set in an object (will be used bellow)
plotA_fish <- ggplot()+#data=FD_fish_MSS[[i]]$axesPCO, 
  #aes(A1, A2)) + # plot the complete space
  geom_polygon(data=df_plot, aes(x=A1,y=A2, group=paste(data,1,sep=""),
                                 fill=data), alpha=0.7) +
  geom_point(data=FD_obs_data_fish$axesPCO,
             aes (x=A1,y=A2,
                  colour = depth),
             alpha=0.7,
             size = 0.5,
             stroke = 1) +
  scale_fill_manual(values = c(complete = "#E6D5B8",
                               fundo = "#0E185F",
                               raso = "#56BBF1")) + # theme_bw() + # the points
  scale_colour_manual(values = c(complete = "#E6D5B8",
                                 fundo = "#56BBF1",
                                 raso = "#0E185F")) + # theme_bw() + # the points
  facet_wrap(~site,ncol= 1) + theme_classic()
plotA_fish
# benthos

df_plot_benthos<- lapply (list(site_to_plot_ccorais,
                       site_to_plot_alcatrazes,
                       site_to_plot_sc,
                       site_to_plot_es), function (t) { # over sites to choose
                         
                         
                         # select sites
                         site_raso <- which(sites_bentos_complete == t[1]) # [2] always shallow
                         site_fundo <- which(sites_bentos_complete == t[2]) # [1] always deep
                         
                         # the complete trait space per rarefaction dataset
                         data_to_space <- cbind(FD_obs_data_benthos$axesPCO[,1:2],ext = F) # the complete space
                         complete_space <- data_to_space [chull(data_to_space[,1:2], y = NULL),] # convex hull on complete space
                         
                         # composition of shallow sites
                         s.raso<-colnames(comp_bentos[site_raso,][which(comp_bentos[site_raso,]>0)])
                         # composition of deep sites
                         s.fundo<-colnames(comp_bentos[site_fundo,][which(comp_bentos[site_fundo,]>0)])
                         
                         # subspace and convex hull of each rarefied data set
                         scores_subset_raso<-data_to_space [rownames(data_to_space)%in% (s.raso),] # subsetting species present into the sites
                         chull_site_raso <- (scores_subset_raso[chull(scores_subset_raso),])# Chull of this subset
                         scores_subset_fundo<-data_to_space [rownames(data_to_space)%in% (s.fundo),] # subsetting species present into the sites
                         chull_site_fundo <- (scores_subset_fundo[chull(scores_subset_fundo),])# Chull of this subset
                         
                         # create df to plot
                         df_plot <- rbind (data.frame (complete_space,data="complete",site= sites_bentos_complete[site_fundo]),
                                           data.frame (chull_site_fundo,data="fundo",site= sites_bentos_complete[site_fundo]),
                                           data.frame(chull_site_raso,data="raso",site= sites_bentos_complete[site_fundo]))
                         
                         ;
                         df_plot
                         
                       })


# melt this df
df_plot_benthos<-do.call(rbind,df_plot_benthos)

# species from shallow and deep areas
comp_bentos_deep <- comp_bentos[grep ("fundo", comp_bentos$locality_site),]
comp_bentos_deep<-colnames(comp_bentos_deep [,which(colSums(comp_bentos_deep[,-1]) >0)])
comp_bentos_shallow <- comp_bentos[grep ("raso", comp_bentos$locality_site),]
comp_bentos_shallow<- colnames(comp_bentos_shallow [,which(colSums(comp_bentos_shallow[,-1]) >0)])
# find the species in both depths
benthos_spp <- data.frame (species = colnames(comp_bentos)[-1],
                           depth = NA)
# depth
benthos_spp[which(benthos_spp$species %in% comp_bentos_deep == T & benthos_spp$species %in% comp_bentos_shallow == F),"depth"] <- "fundo"
benthos_spp[which(benthos_spp$species %in% comp_bentos_shallow == T & benthos_spp$species %in% comp_bentos_deep == F),"depth"] <- "raso"
benthos_spp[which(benthos_spp$species %in% comp_bentos_shallow == T & benthos_spp$species %in% comp_bentos_deep == T),"depth"] <- "complete"
benthos_spp$depth [is.na(benthos_spp$depth)] <- "complete"
# ordering
FD_obs_data_benthos$axesPCO <- FD_obs_data_benthos$axesPCO [order(rownames(FD_obs_data_benthos$axesPCO)),]
benthos_spp<-benthos_spp [which(benthos_spp$species %in% rownames(FD_obs_data_benthos$axesPCO)),]
benthos_spp$species == rownames(FD_obs_data_benthos$axesPCO)
# bind depth
FD_obs_data_benthos$axesPCO <- cbind (FD_obs_data_benthos$axesPCO,
                                      depth=benthos_spp$depth)
FD_obs_data_benthos$axesPCO$depth <- factor (FD_obs_data_benthos$axesPCO$depth,
                                             levels = c("complete",
                                                        "fundo",
                                                        "raso"))
df_plot_benthos$data<-factor (df_plot_benthos$data,
                              levels = c("complete",
                                         "fundo",
                                         "raso"))

# create a plot and set in an object (will be used bellow)
plotA_benthos <- ggplot()+#data=FD_fish_MSS[[i]]$axesPCO, 
  #aes(A1, A2)) + # plot the complete space
  #geom_point() + theme_bw() + # the points
  geom_polygon(data=df_plot_benthos, aes(x=A1,y=A2, group=paste(data,1,sep=""),
                                 fill=data), 
               linetype=2, alpha=0.5) +
  #scale_fill_manual(values = c("gray70", "#FB9300", "#FAFCC2")) + 
  geom_point(data=FD_obs_data_benthos$axesPCO,
             aes (x=A1,y=A2,
                  colour = depth),
             alpha=0.7,
             size = 0.5,
             stroke = 1) +
  scale_fill_manual(values = c(complete = "#E6D5B8",
                                 fundo = "#FF5F00",
                               raso = "#F6F54D")) + # theme_bw() + # the points
  scale_colour_manual(values = c(complete = "#E6D5B8",
                                 fundo = "#F6F54D",
                                 raso = "#FF5F00")) + # theme_bw() + # the points
  facet_wrap(~site,ncol= 1) + theme_classic()
plotA_benthos

# save
pdf(file=here("output_no_resampling","Fig1_spaces_with_pts"),height=5,width=6.5)

grid.arrange(map_peixes_bentos, 
             plotA_benthos+theme(legend.position = "none"),
             plotA_fish+theme(legend.position = "none"),
             ncol=9,nrow=11,
             layout_matrix = rbind (c(1,1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,1,2,2,2,3,3,3)))
dev.off()

# ------------------------------------------------

# coefficients of drivers

load (here ("output_no_resampling", "MCMC_runs_multivariate_rarefied_no_aut.RData"))

# model selection

res$looic
res$param
res$best_model


# descripte statistics of fit
library("bayesplot")
library("ggplot2")
#library("rstanarm")   

# extract posterior probs
posterior <- as.array(res$best_model[[1]])
dimnames(posterior)

# plot
# here: https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html
pdf(here ("output_no_resampling","fig2_temp.pdf"),width=7,height=4)
color_scheme_set("blue")
mcmc_areas(posterior, pars = c("b_logSR_BO2_tempmean_ss_std",         
                               "b_logSRbenthos_BO2_tempmean_ss_std", 
                               "b_logFRic_BO2_tempmean_ss_std",
                               "b_logFRicbenthos_BO2_tempmean_ss_std", 
                                   "b_logFEve_BO2_tempmean_ss_std",
                               "b_logFEvebenthos_BO2_tempmean_ss_std",     
                               
                                   "b_logFDiv_BO2_tempmean_ss_std",
                                   
                               
                                   "b_logFDivbenthos_BO2_tempmean_ss_std"),
               prob = 0.5, # 85% intervals
               prob_outer = 0.95, # 95%
               point_est = "median",
           area_method = "equal height")

dev.off()

# area
pdf(here ("output_no_resampling","fig2_area.pdf"),width=7,height=4)
mcmc_areas(posterior, pars = c("b_logSR_area",         
                               "b_logSRbenthos_area", 
                               "b_logFRic_area",
                               "b_logFRicbenthos_area", 
                               "b_logFEve_area",
                               "b_logFEvebenthos_area",     
                               
                               "b_logFDiv_area",
                               
                               
                               "b_logFDivbenthos_area"),
           prob = 0.5, # 85% intervals
           prob_outer = 0.95, # 95%
           point_est = "mean",
           area_method = "equal height")
dev.off()
#sampling area
pdf(here ("output_no_resampling","fig2_samplingarea.pdf"),width=7,height=4)
mcmc_areas(posterior, pars = c("b_logSR_SamplingArea_std",         
                             "b_logSRbenthos_SamplingArea_std", 
                             "b_logFRic_SamplingArea_std",
                             "b_logFRicbenthos_SamplingArea_std", 
                             "b_logFEve_SamplingArea_std",
                             "b_logFEvebenthos_SamplingArea_std",     
                             
                               "b_logFDiv_SamplingArea_std",
                             
                               
                   "b_logFDivbenthos_SamplingArea_std"),
            prob = 0.5, # 85% intervals
            prob_outer = 0.95, # 95%
            point_est = "mean",
            area_method = "equal height")
dev.off()
# ====================================================== #

library(terra)
library(flexmix)
library(modeltools)
library(tidyverse)
library(brms)
library(ggdist)
require(here)
require(tidybayes)

# load results
#load(here ("output_sensitivity","MCMC_runs_multivariate_rarefied_no_aut.rdata"))

# recompile model with more iterations and chains
nd <- res$best_model[[1]]$data %>%
  dplyr::select(-BO2_tempmean_ss_std,
                -area,
                -SamplingArea_std) %>%
  dplyr::mutate(dplyr::across(tidyselect::everything(), log)) %>%
  dplyr::rename_with(~paste0("log", .x)) %>%
  dplyr::mutate(sst = res$best_model[[1]]$data$BO2_tempmean_ss_std,
                area = res$best_model[[1]]$data$area,
                samp.area = res$best_model[[1]]$data$SamplingArea_std)
model <- brms::brm(
  mvbind(logSR, logFRic, logFEve, logFDiv, logSR_benthos, logFRic_benthos, logFEve_benthos, logFDiv_benthos) ~ sst+area+samp.area,
  data = nd, cores = nc, chains = nc, iter = ni, warmup = nb
)

# bayesian R2
br2<- bayes_R2(model) # mean for each reponse
round(br2,2)
round(mean(br2 [,'Estimate']),2) # average across responses

# predicted values
all_preds <- exp(brms::posterior_predict(model))
cor_out <- vector(mode = "list", length = nrow(all_preds))
for (i in seq_len(nrow(all_preds))) {
  cor_i <- cor(as.matrix(all_preds[i, , ]))
  ind <- which(upper.tri(cor_i, diag = TRUE), arr.ind = TRUE)
  nn <- dimnames(cor_i)
  cor_out[[i]] <- data.frame(var_a = nn[[1]][ind[, 1]],
                             var_b = nn[[2]][ind[, 2]],
                             correlation_r = cor_i[ind]) %>%
    dplyr::filter(var_a != var_b) %>%
    dplyr::mutate(.draw = i)
}

cor_out <- do.call("rbind.data.frame", cor_out) %>%
  dplyr::mutate(dplyr::across(tidyselect::starts_with("var_"),
                              ~gsub("^log", "", .x))) %>%
  dplyr::mutate(var_a = case_when(!grepl("benthos", var_a) ~ paste0(var_a, "_fish"), TRUE ~ var_a),
                var_b = case_when(!grepl("benthos", var_b) ~ paste0(var_b, "_fish"), TRUE ~ var_b),
                cor_pair = paste0(var_a, " ~ ", var_b))

# select correlation in the diagonal (i.e., cor of the same metric between fish and benthos)
sel_cors <- c("SR_fish ~ SRbenthos",
              "FRic_fish ~ FRicbenthos",
              "FEve_fish ~ FEvebenthos",
              "FDiv_fish ~ FDivbenthos")

# selecting interesting correlations
cor_out_sel <- cor_out[which(cor_out$cor_pair %in% sel_cors),]
(cor_out_sel %>%
    dplyr::group_by(cor_pair) %>%
    dplyr::summarise(ggdist::median_hdci(correlation_r),
                     prob_pos = sum(correlation_r > 0) / n(),
                     prob_neg = sum(correlation_r < 0) / n())
)

# aggregating to get averages and plot a vertical line into the histogram
#aggregate (cor_out_sel,by=list(cor_out_sel$cor_pair),
#           FUN=mean)

vline.data <- cor_out %>%
  group_by(cor_pair) %>%
  summarize(z = mean(correlation_r ))

# plot posterior distribution of each pair as histogram
post_cors_plot <- ggplot(data = cor_out) +
  geom_histogram(mapping = aes(x = correlation_r), fill = "dodgerblue2",
                 colour = "grey30",
                 bins=60) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
  labs(x = "Pearson correlation", y = "Frequency") +
  facet_wrap(~cor_pair) +
  geom_vline(aes(xintercept = z), vline.data, 
             colour = "red",size=1.5)+ # certical line with average
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))
post_cors_plot

ggsave("post_cors_plot.pdf", post_cors_plot, width = 11.88, height = 8.81)

# summary table with median +/- 95% credible intervals
# probability of correlation being negative or positive.
(cor_out %>%
    dplyr::group_by(cor_pair) %>%
    dplyr::summarise(ggdist::median_hdci(correlation_r),
                     prob_pos = sum(correlation_r > 0) / n(),
                     prob_neg = sum(correlation_r < 0) / n())
)

# plot observed vs. mean predicted correlation
post_cor_means <- cor_out %>%
  dplyr::group_by(cor_pair) %>%
  dplyr::summarise(mean_post_r = mean(correlation_r))


cor_orig <- nd %>%
  dplyr::select(-sst) %>%
  dplyr::rename_with(~gsub("^log", "", .x)) %>%
  exp %>%
  cor

ind <- which(upper.tri(cor_orig, diag = TRUE), arr.ind = TRUE)
nn <- dimnames(cor_orig)
cor_orig <- data.frame(var_a = nn[[1]][ind[, 1]], var_b = nn[[2]][ind[, 2]],
                       obs_cor_r = cor_orig[ind]) %>%
  dplyr::filter(var_a != var_b) %>%
  dplyr::mutate(var_a = case_when(!grepl("benthos", var_a) ~ paste0(var_a, "_fish"), TRUE ~ gsub("_", "", var_a)),
                var_b = case_when(!grepl("benthos", var_b) ~ paste0(var_b, "_fish"), TRUE ~ gsub("_", "", var_b)),
                cor_pair = paste0(var_a, " ~ ", var_b))
merged_cor <- merge(post_cor_means, cor_orig[, -c(1:2)])
# factors to represent indexes
#merged_cor$cor_pair <- gsub ("EstRich", "SR",merged_cor$cor_pair )
merged_cor$cor_pair <- gsub ("benthos", "Benthos",merged_cor$cor_pair )
merged_cor$cor_pair <- gsub ("fish", "Fish",merged_cor$cor_pair )
merged_cor$cor_pair <- gsub ("_", "",merged_cor$cor_pair )

# bind residual correlations
require(reshape)
corr_res <- VarCorr(model)$residual__$cor # residals coor
corr_res <- melt(corr_res[,1,]) # melt
corr_res <- corr_res[which(corr_res$value != 1),]# rm corr == 1
# edit labels
# var1
#corr_res$Var1<-gsub ("_","",corr_res$Var1)
corr_res[,1]<-gsub ("log","",corr_res[,1])
corr_res[,1] <- gsub ("EstRich", "SR",corr_res[,1])
corr_res[,1] <- gsub ("benthos", "Benthos",corr_res[,1])
corr_res[,1] <- gsub ("fish", "Fish",corr_res[,1] )
corr_res[,1] <- gsub ("_", "",corr_res[,1] )
corr_res[,1][which(corr_res[,1] == "SR")] <- "SRFish"
corr_res[,1][which(corr_res[,1] == "FRic")] <- "FRicFish"
corr_res[,1][which(corr_res[,1] == "FEve")] <- "FEveFish"
corr_res[,1][which(corr_res[,1] == "FDiv")] <- "FDivFish"
# var2
corr_res[,2]<-gsub ("_","",corr_res[,2])
corr_res[,2]<-gsub ("log","",corr_res[,2])
corr_res[,2] <- gsub ("EstRich", "SR",corr_res[,2])
corr_res[,2] <- gsub ("benthos", "Benthos",corr_res[,2])
corr_res[,2] <- gsub ("fish", "Fish",corr_res[,2])
corr_res[,2] <- gsub ("_", "",corr_res[,2])
corr_res[,2] [which(corr_res[,2] == "SR")] <- "SRFish"
corr_res[,2][which(corr_res[,2] == "FRic")] <- "FRicFish"
corr_res[,2][which(corr_res[,2] == "FEve")] <- "FEveFish"
corr_res[,2][which(corr_res[,2] == "FDiv")] <- "FDivFish"
# bind
corr_res$var3 <- paste (corr_res[,1], corr_res[,2], sep= " ~ ")
corr_res <- corr_res[which(corr_res$var3 %in% merged_cor$cor_pair),]

# bind resid corr
merged_cor$corr_res <- corr_res$value[match (merged_cor$cor_pair,corr_res$var3)] #== merged_cor$cor_pair

# plot
scatter_cor <- ggplot(data = merged_cor) +
  geom_point(mapping = aes(y = obs_cor_r, x = mean_post_r), shape = 21,
             size = 3, fill = "grey90") +
  theme_classic() +
  geom_abline(slope = 1, linetype = 3) +
  labs(y = "Observed Pearson correlation",
       x = "Posterior mean of predicted Pearson's correlation") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) +
  xlim(-0.6, 1) + ylim(-0.6, 1)

require (ggrepel)
scatter_cor + geom_text_repel(data = merged_cor, aes (y = obs_cor_r, x = mean_post_r,
                                                      label = cor_pair),
                              box.padding = 0.6, 
                              max.overlaps = Inf,
                              size=3)

ggsave(here ("output_sensitivity","scatter_cor.pdf"), 
       scatter_cor, width = 6.5, height = 6)


# resid corr
# plot
scatter_cor_res <- ggplot(data = merged_cor) +
  geom_point(mapping = aes(y = obs_cor_r, x = corr_res), 
             shape = 21,
             size = 3, fill = "grey90") +
  theme_classic() +
  geom_abline(slope = 1, linetype = 3) +
  labs(y = "Observed Pearson correlation",
       x = "Posterior mean of Pearson's residual correlation") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) +
  xlim(-0.6, 1) + ylim(-0.6, 1)

require (ggrepel)
scatter_cor_res + geom_text_repel(data = merged_cor, aes (y = obs_cor_r, x = corr_res,
                                                          label = cor_pair),
                                  box.padding = 0.6, 
                                  max.overlaps = Inf,
                                  size=3)

# predicted vs residual
scatter_cor_res <- ggplot(data = merged_cor) +
  geom_point(mapping = aes(y = mean_post_r, x = corr_res), 
             shape = 21,
             size = 3, fill = "grey90") +
  theme_classic() +
  geom_abline(slope = 1, linetype = 3) +
  labs(y = "Posterior mean of predicted Pearson's correlation",
       x = "Posterior mean of Pearson's residual correlation") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) +
  xlim(-0.6, 1) + ylim(-0.6, 1)

require (ggrepel)
scatter_cor_res + geom_text_repel(data = merged_cor, aes (y = mean_post_r, 
                                                          x = corr_res,
                                                          label = cor_pair),
                                  box.padding = 0.6, 
                                  max.overlaps = Inf,
                                  size=3)

# absolute diff

merged_cor<- cbind(merged_cor, 
                   diff = (abs(merged_cor$mean_post_r) - abs(merged_cor$corr_res)))
merged_cor<- cbind(merged_cor, 
                   diff2 = (abs(merged_cor$mean_post_r) - abs(merged_cor$obs_cor_r)))

merged_cor<-melt(merged_cor)
colnames(merged_cor)[which(colnames(merged_cor) == "variable")] <- "Correlation"
levels(merged_cor$Correlation)[which(levels(merged_cor$Correlation) == "obs_cor_r")] <- "Observed"
levels(merged_cor$Correlation)[which(levels(merged_cor$Correlation) == "corr_res")] <- "Residual"
levels(merged_cor$Correlation)[which(levels(merged_cor$Correlation) == "mean_post_r")] <- "Predicted"

# plot diff

a<- ggplot (merged_cor[which(merged_cor$Correlation %in% c("Observed", "Residual", "Predicted")),], 
            aes(x= value,
                y= reorder(cor_pair, - value),
                group = cor_pair,
                color = Correlation
            )) + 
  geom_point(size=2) + 
  theme_classic()  +
  xlab("Pearson's correlation") + 
  ylab ("Pairs of diversity metrics")+
  theme(legend.position = c(0.8,0.8)) + 
  geom_vline(aes(xintercept =0),alpha =0.5,size=1,col = "gray") +
  xlim(c(-1,1))


a

# ==============================================================
# param probability

vars_ext<-get_variables(model)[grep("sst",get_variables(model))]

# sst SR fish
model %>%
  spread_draws(b_logSR_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSR_sst),
                   prob_pos = sum(b_logSR_sst > 0) / n(),
                   prob_neg = sum(b_logSR_sst < 0) / n())

# sst SR benthos
model %>%
  spread_draws(b_logSRbenthos_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSRbenthos_sst),
                   prob_pos = sum(b_logSRbenthos_sst > 0) / n(),
                   prob_neg = sum(b_logSRbenthos_sst < 0) / n())

# sst Fric fish
model %>%
  spread_draws(b_logFRic_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRic_sst),
                   prob_pos = sum(b_logFRic_sst > 0) / n(),
                   prob_neg = sum(b_logFRic_sst < 0) / n())

# sst fric benthos
model %>%
  spread_draws(b_logFRicbenthos_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicbenthos_sst),
                   prob_pos = sum(b_logFRicbenthos_sst > 0) / n(),
                   prob_neg = sum(b_logFRicbenthos_sst < 0) / n())

# sst FEve fish
model %>%
  spread_draws(b_logFEve_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEve_sst),
                   prob_pos = sum(b_logFEve_sst > 0) / n(),
                   prob_neg = sum(b_logFEve_sst < 0) / n())

# sst FEve benthos
model %>%
  spread_draws(b_logFEvebenthos_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEvebenthos_sst),
                   prob_pos = sum(b_logFEvebenthos_sst > 0) / n(),
                   prob_neg = sum(b_logFEvebenthos_sst < 0) / n())

# sst Fdiv fish
model %>%
  spread_draws(b_logFDiv_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDiv_sst),
                   prob_pos = sum(b_logFDiv_sst > 0) / n(),
                   prob_neg = sum(b_logFDiv_sst < 0) / n())

# sst fdiv benthos
model %>%
  spread_draws(b_logFDivbenthos_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDivbenthos_sst),
                   prob_pos = sum(b_logFDivbenthos_sst > 0) / n(),
                   prob_neg = sum(b_logFDivbenthos_sst < 0) / n())


# reef area
vars_ext<-get_variables(model)[grep("area",get_variables(model))]

# sst SR fish
model %>%
  spread_draws(b_logSR_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSR_area),
                   prob_pos = sum(b_logSR_area > 0) / n(),
                   prob_neg = sum(b_logSR_area < 0) / n())

# sst SR benthos
model %>%
  spread_draws(b_logSRbenthos_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSRbenthos_area),
                   prob_pos = sum(b_logSRbenthos_area > 0) / n(),
                   prob_neg = sum(b_logSRbenthos_area < 0) / n())

# sst Fric fish
model %>%
  spread_draws(b_logFRic_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRic_area),
                   prob_pos = sum(b_logFRic_area > 0) / n(),
                   prob_neg = sum(b_logFRic_area < 0) / n())

# sst fric benthos
model %>%
  spread_draws(b_logFRicbenthos_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicbenthos_area),
                   prob_pos = sum(b_logFRicbenthos_area > 0) / n(),
                   prob_neg = sum(b_logFRicbenthos_area < 0) / n())

# sst FEve fish
model %>%
  spread_draws(b_logFEve_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEve_area),
                   prob_pos = sum(b_logFEve_area > 0) / n(),
                   prob_neg = sum(b_logFEve_area < 0) / n())

# sst FEve benthos
model %>%
  spread_draws(b_logFEvebenthos_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEvebenthos_area),
                   prob_pos = sum(b_logFEvebenthos_area > 0) / n(),
                   prob_neg = sum(b_logFEvebenthos_area < 0) / n())


# sst Fdiv fish
model %>%
  spread_draws(b_logFDiv_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDiv_area),
                   prob_pos = sum(b_logFDiv_area > 0) / n(),
                   prob_neg = sum(b_logFDiv_area < 0) / n())

# sst fdiv benthos
model %>%
  spread_draws(b_logFDivbenthos_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDivbenthos_area),
                   prob_pos = sum(b_logFDivbenthos_area > 0) / n(),
                   prob_neg = sum(b_logFDivbenthos_area < 0) / n())


# sampliung area

vars_ext<-get_variables(model)[grep("samp.area",get_variables(model))]

# sst SR fish
model %>%
  spread_draws(b_logSR_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSR_samp.area),
                   prob_pos = sum(b_logSR_samp.area > 0) / n(),
                   prob_neg = sum(b_logSR_samp.area < 0) / n())

# sst SR benthos
model %>%
  spread_draws(b_logSRbenthos_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSRbenthos_samp.area),
                   prob_pos = sum(b_logSRbenthos_samp.area > 0) / n(),
                   prob_neg = sum(b_logSRbenthos_samp.area < 0) / n())

# sst Fric fish
model %>%
  spread_draws(b_logFRic_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRic_samp.area),
                   prob_pos = sum(b_logFRic_samp.area > 0) / n(),
                   prob_neg = sum(b_logFRic_samp.area < 0) / n())

# sst fric benthos
model %>%
  spread_draws(b_logFRicbenthos_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicbenthos_samp.area),
                   prob_pos = sum(b_logFRicbenthos_samp.area > 0) / n(),
                   prob_neg = sum(b_logFRicbenthos_samp.area < 0) / n())

# sst FEve fish
model %>%
  spread_draws(b_logFEve_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEve_samp.area),
                   prob_pos = sum(b_logFEve_samp.area > 0) / n(),
                   prob_neg = sum(b_logFEve_samp.area < 0) / n())

# sst FEve benthos
model %>%
  spread_draws(b_logFEvebenthos_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEvebenthos_samp.area),
                   prob_pos = sum(b_logFEvebenthos_samp.area > 0) / n(),
                   prob_neg = sum(b_logFEvebenthos_samp.area < 0) / n())


# sst Fdiv fish
model %>%
  spread_draws(b_logFDiv_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDiv_samp.area),
                   prob_pos = sum(b_logFDiv_samp.area > 0) / n(),
                   prob_neg = sum(b_logFDiv_samp.area < 0) / n())

# sst fdiv benthos
model %>%
  spread_draws(b_logFDivbenthos_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDivbenthos_samp.area),
                   prob_pos = sum(b_logFDivbenthos_samp.area > 0) / n(),
                   prob_neg = sum(b_logFDivbenthos_samp.area < 0) / n())




