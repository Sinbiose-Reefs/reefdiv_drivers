
## --------------------
#       FIGURES
# --------------------

## call packages
source("R/packages.R")
source("R/functions.R")

# -----------------------------------------
## correlation between variables

load(here("output", "data_to_modeling_GLM.RData"))

# fish
pdf(here("output","vectorized", "correlation_variables.pdf"))
par(mfrow=c(1,2))
corrplot (
  cor (cov_fish [,c("BO2_tempmean_ss_std","BO2_ppmean_ss_std","BO2_salinitymean_ss_std",
                    "BO_damean_std","distance_std")]), 
  p.mat = cor.mtest(cor (cov_fish [,c("BO2_tempmean_ss_std","BO2_ppmean_ss_std","BO2_salinitymean_ss_std",
                                      "BO_damean_std","distance_std")]))$p)

# benthos 
corrplot (
  cor (cov_benthos [,c("BO2_tempmean_ss_std","BO2_ppmean_ss_std","BO2_salinitymean_ss_std",
                       "BO_damean_std","distance_std")]), 
  p.mat = cor.mtest(cor (cov_benthos [,c("BO2_tempmean_ss_std","BO2_ppmean_ss_std","BO2_salinitymean_ss_std",
                                         "BO_damean_std","distance_std")]))$p)

dev.off()


# ------------------------------------------ #
# Load rarefied data of fishes and benthos
# sites x transect / video x species

load (here ("output","random_composition_bentos.RData"))
load (here ("output","random_composition_fish.RData"))

# ------------------------------------------ #
#
#     Load environment data

load (here ("output","env_data.RData"))

# ------------------------------------------- #
# Fig S1
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

pdf(file=here("output","vectorized","FigS1.pdf"),height=5,width=9)

grid.arrange(map_peixes_bentos, 
             eff_obsrich_plot,
             eff_obsrich_plot_benthos,
             SST,
             ncol=10,nrow=11,
             layout_matrix = rbind (c(1,1,1,1,2,2,3,3,4,4),
                                    c(1,1,1,1,2,2,3,3,4,4),
                                    c(1,1,1,1,2,2,3,3,4,4),
                                    c(1,1,1,1,2,2,3,3,4,4),
                                    c(1,1,1,1,2,2,3,3,4,4),
                                    c(1,1,1,1,2,2,3,3,4,4),
                                    c(1,1,1,1,2,2,3,3,4,4),
                                    c(1,1,1,1,2,2,3,3,4,4),
                                    c(1,1,1,1,2,2,3,3,4,4)))

dev.off()


# ------------------------------------------ #
#
#     Load results of the rarefaction

# -----------------------------------------------

# MSS algorithm 

# ------------------------------------------ #

load (here ("output","FD_fish_MSS.RData"))
load (here ("output","FD_benthos_MSS.RData"))

# remove the ~15  random samples with problems (those with length==2)
FD_fish_MSS_mod<-FD_fish_MSS [-which(lapply (FD_fish_MSS,length) ==2)] # 2 is the length of problematic samples 
FD_benthos_MSS_mod<-FD_benthos_MSS [-which(lapply (FD_benthos_MSS,length) ==2)]

# also rm from composition
rdm_composition_complete_mod<-rdm_composition_complete [-which(lapply (FD_fish_MSS,length) ==2)]
rdm_composition_complete_bentos_mod<-rdm_composition_complete_bentos [-which(lapply (FD_benthos_MSS,length) ==2)]

# select of sites to plot (extremes in index values)
site_to_plot_alcatrazes <-c("alcatrazes.portinho_sudoeste.raso","alcatrazes.portinho_sudoeste.fundo")
site_to_plot_sc <- c("ilhasc_norte.deserta_norte.raso","ilhasc_norte.deserta_norte.fundo")
site_to_plot_ccorais <- c("btds_santos.farol_da_barra.raso","btds_santos.poste_quatro.fundo")
site_to_plot_es <- c("espirito_santo.escalvada.raso","espirito_santo.escalvada.fundo")

# now produce data for  polygons in the first plot  

df_plot<- lapply (list(site_to_plot_ccorais,
                       site_to_plot_alcatrazes,
                       site_to_plot_sc,
                       site_to_plot_es), function (t) # over sites to choose
  
  lapply (seq(1,length(FD_fish_MSS_mod)), function (k) { # over rarefied data
    
    # select sites
    site_raso <- which(sites_fish_complete == t[1]) # [2] always shallow
    site_fundo <- which(sites_fish_complete == t[2]) # [1] always deep
    
    # the complete trait space per rarefaction dataset
    data_to_space <- cbind(FD_fish_MSS_mod[[k]]$axesPCO[,1:2],ext = F) # the complete space
    complete_space <- data_to_space [chull(data_to_space[,1:2], y = NULL),] # convex hull on complete space
    
    # composition of shallow sites
    s.raso<-colnames(rdm_composition_complete_mod[[k]][site_raso,][which(rdm_composition_complete_mod[[k]][site_raso,]>0)])
    # composition of deep sites
    s.fundo<-colnames(rdm_composition_complete_mod[[k]][site_fundo,][which(rdm_composition_complete_mod[[k]][site_fundo,]>0)])
    
    # subspace and convex hull of each rarefied data set
    scores_subset_raso<-data_to_space [rownames(data_to_space)%in% (s.raso),] # subsetting species present into the sites
    chull_site_raso <- (scores_subset_raso[chull(scores_subset_raso),])# Chull of this subset
    scores_subset_fundo<-data_to_space [rownames(data_to_space)%in% (s.fundo),] # subsetting species present into the sites
    chull_site_fundo <- (scores_subset_fundo[chull(scores_subset_fundo),])# Chull of this subset
    
    # create df to plot
    df_plot <- rbind (data.frame (complete_space,data="complete",raref.set=k,site= sites_fish_complete[site_fundo]),
                      data.frame (chull_site_fundo,data="fundo",raref.set=k,site= sites_fish_complete[site_fundo]),
                      data.frame(chull_site_raso,data="raso",raref.set=k,site= sites_fish_complete[site_fundo]))
    
    ;
    df_plot
    
  }))

# melt this df
df_plot <- lapply (df_plot, function (i) do.call(rbind, i))
df_plot<-do.call(rbind,df_plot)

# create a plot and set in an object (will be used bellow)
plotA_fish <- ggplot()+#data=FD_fish_MSS[[i]]$axesPCO, 
                #aes(A1, A2)) + # plot the complete space
  #geom_point() + theme_bw() + # the points
  geom_polygon(data=df_plot, aes(x=A1,y=A2, group=paste(data,raref.set,sep=""),
                                 fill=data), alpha=0.02) +
  scale_fill_manual(values = c("gray70", "#046582", "#B5EAEA")) + 
  facet_wrap(~site) + theme_classic()

# results to annotate into the plot (average indexes across rarefaction samples)
indexes_per_sample_site_fish <- lapply (list(site_to_plot_ccorais,
                                        site_to_plot_alcatrazes,
                                        site_to_plot_sc,
                                        site_to_plot_es), function (t) { # over sites to choose
                                                        
                                         lapply (FD_fish_MSS_mod, function (i) {# for each rarefied dataset
                                              deep <- i$Fdindexes[which(sites_fish_complete == t[1]),]
                                              shallow <- i$Fdindexes[which(sites_fish_complete == t[2]),]
                                              res <- (list(shallow=shallow,
                                                                deep=deep))
                                              ;# return
                                              res
                                         })
                                    }
                                )

# get the index values to average and show
# shallow
av_shallow_per_site_fish <- lapply (indexes_per_sample_site_fish,function (i) # each site
  round(apply (do.call(rbind,sapply (i, "[","shallow")),2,mean,na.rm=T),2)
)
# adjust names
names(av_shallow_per_site_fish)<-c(site_to_plot_ccorais[1],
           site_to_plot_alcatrazes[1],
           site_to_plot_sc[1],
           site_to_plot_es[1])
# deep
av_deep_per_site_fish <- lapply (indexes_per_sample_site_fish,function (i) # each site
  round(apply (do.call(rbind,sapply (i, "[","deep")),2,mean,na.rm=T),2)
)
# adjust names
names(av_deep_per_site_fish)<-c(site_to_plot_ccorais[2],
                              site_to_plot_alcatrazes[2],
                              site_to_plot_sc[2],
                              site_to_plot_es[2])


## --------------------------------------------------------------
# Benthos
## --------------------------------------------------------------

# considering the list of selected sites
# now produce data for  polygons in the first plot  
df_plot<- lapply (list(site_to_plot_ccorais,
                       site_to_plot_alcatrazes,
                       site_to_plot_sc,
                       site_to_plot_es), function (t) # over sites to choose
                         
                         lapply (seq(1,length(FD_benthos_MSS_mod)), function (k) { # over rarefied data
                           
                           # select sites
                           site_fundo <- which(sites_bentos_complete == t[1]) # [1] always deep
                           site_raso <- which(sites_bentos_complete == t[2]) # [2] always shallow
                           
                           # the complete trait space per rarefaction dataset
                           data_to_space <- cbind(FD_benthos_MSS_mod[[k]]$axesPCO[,1:2],ext = F) # the complete space
                           complete_space <- data_to_space [chull(data_to_space[,1:2], y = NULL),] # convex hull on complete space
                           
                           # composition of deep sites
                           s.fundo<-colnames(rdm_composition_complete_bentos_mod[[k]][site_fundo,][which(rdm_composition_complete_bentos_mod[[k]][site_fundo,]>0)])
                           # composition of shallow sites
                           s.raso<-colnames(rdm_composition_complete_bentos_mod[[k]][site_raso,][which(rdm_composition_complete_bentos_mod[[k]][site_raso,]>0)])
                           
                           # subspace and convex hull of each rarefied data set
                           scores_subset_fundo<-data_to_space [rownames(data_to_space)%in% (s.fundo),] # subsetting species present into the sites
                           chull_site_fundo <- (scores_subset_fundo[chull(scores_subset_fundo),])# Chull of this subset
                           scores_subset_raso<-data_to_space [rownames(data_to_space)%in% (s.raso),] # subsetting species present into the sites
                           chull_site_raso <- (scores_subset_raso[chull(scores_subset_raso),])# Chull of this subset
                           
                           # create df to plot
                           df_plot <- rbind (data.frame (complete_space,data="complete",raref.set=k,site= sites_fish_complete[site_fundo]),
                                             data.frame (chull_site_fundo,data="fundo",raref.set=k,site= sites_fish_complete[site_fundo]),
                                             data.frame(chull_site_raso,data="raso",raref.set=k,site= sites_fish_complete[site_fundo]))
                           
                           ;
                           df_plot
                           
                         }))

# melt this df
df_plot <- lapply (df_plot, function (i) do.call(rbind, i))
df_plot<-do.call(rbind,df_plot)

# create a plot and set in an object (will be used bellow)
plotA_benthos <- ggplot()+#data=FD_fish_MSS[[i]]$axesPCO, 
  #aes(A1, A2)) + # plot the complete space
  #geom_point() + theme_bw() + # the points
  geom_polygon(data=df_plot, aes(x=A1,y=A2, group=paste(data,raref.set,sep=""),
                                 fill=data), alpha=0.02) +
  scale_fill_manual(values = c("gray70", "#FB9300", "#FAFCC2")) + 
  facet_wrap(~site) + theme_classic()

# results to annotate into the plot (average indexes across rarefaction samples)
indexes_per_sample_site <- lapply (list(site_to_plot_ccorais,
                                        site_to_plot_alcatrazes,
                                        site_to_plot_sc,
                                        site_to_plot_es), function (t) { # over sites to choose
                                          
                                          lapply (FD_benthos_MSS_mod, function (i) {# for each rarefied dataset
                                            deep <- i$Fdindexes[which(sites_bentos_complete == t[1]),]
                                            shallow <- i$Fdindexes[which(sites_bentos_complete == t[2]),]
                                            res <- (list(shallow=shallow,
                                                         deep=deep))
                                            ;# return
                                            res
                                          })
                                        }
)

# get the index values to average and show
# shallow
av_shallow_per_site <- lapply (indexes_per_sample_site,function (i) # each site
  round(apply (do.call(rbind,sapply (i, "[","shallow")),2,mean,na.rm=T),2)
)
# adjust names
names(av_shallow_per_site)<-c(site_to_plot_ccorais[2],
                              site_to_plot_alcatrazes[2],
                              site_to_plot_sc[2],
                              site_to_plot_es[2])
# deep
av_deep_per_site <- lapply (indexes_per_sample_site,function (i) # each site
  round(apply (do.call(rbind,sapply (i, "[","deep")),2,mean,na.rm=T),2)
)
# adjust names
names(av_deep_per_site)<-c(site_to_plot_ccorais[1],
                           site_to_plot_alcatrazes[1],
                           site_to_plot_sc[1],
                           site_to_plot_es[1])

# ---------------------------------------------------------------
# match functional spaces with maps of sites

pdf(file=here("output","vectorized","Fig1_MSS.pdf"),height=5,width=6.5)

grid.arrange(map_peixes_bentos, 
             plotA_benthos,
             plotA_fish,
             ncol=9,nrow=11,
             layout_matrix = rbind (c(1,1,1,1,2,2,2,2,2),
                                    c(1,1,1,1,2,2,2,2,2),
                                    c(1,1,1,1,2,2,2,2,2),
                                    c(1,1,1,1,2,2,2,2,2),
                                    c(1,1,1,1,3,3,3,3,3),
                                    c(1,1,1,1,3,3,3,3,3),
                                    c(1,1,1,1,3,3,3,3,3),
                                    c(1,1,1,1,3,3,3,3,3)))
dev.off()

# data to annotate into the plot
av_shallow_per_site_fish # fish shallow
av_deep_per_site_fish # fish deep
av_shallow_per_site # benthos shallow
av_deep_per_site # benthos deep

# average intertia per axis
#benthos
axis_inertia_benthos <-sapply(FD_benthos_MSS,"[[", "explanAxes")# extract
colMeans(do.call(rbind,sapply (axis_inertia_benthos, "[",1:2)))# average of first two axes across rarefied samples
# fish
axis_inertia_fish <-sapply(FD_fish_MSS,"[[", "explanAxes")# extract
colMeans(do.call(rbind,sapply (axis_inertia_fish, "[",1:2)))# average of first two axes across rarefied samples


### -------------------------
#  FIG 2
# coefficients of drivers

load (here ("output", "MCMC_runs_multivariate_rarefied_no_aut.RData"))

# model selection

res$looic
res$param
res$best_model

# extract model effects
# 95% credible interval
effects_model<- fixef(res$best_model[[1]],
      summary = TRUE,
      robust = FALSE,
      probs = c(0.025, 0.975))

#90% credible interval
effects_model90<-fixef(res$best_model[[1]],
        summary = TRUE,
        robust = FALSE,
        probs = c(0.15, 0.85))


# bind algorithm name
effects_model <- data.frame (effects_model,
                             effects_model90[,3:4])
  

# melt the list to have a df with all results
colnames(effects_model)[c(3:6)] <- c("lower95","upper95","lower90","upper90") # change colnames
# adjust these results to plot
complete_results<- data.frame (effects_model,
                               Parameter = NA,
                               Organism = NA,
                               Index = NA)
# adjust params
complete_results [grep ("Intercept",rownames (complete_results)),"Parameter"] <- "Intercept"
complete_results [grep ("BO2_tempmean_ss_std",rownames (complete_results)),"Parameter"] <- "Temperature"
# adjust organisms
complete_results [grep ("benthos",rownames (complete_results)),"Organism"] <- "Benthos"
complete_results [is.na(complete_results$Organism),"Organism"] <- "Fish"
# adjust indexes
complete_results [grep ("logEstRich",rownames (complete_results)),"Index"] <- "SR"
complete_results [grep ("logFRic",rownames (complete_results)),"Index"] <- "FRic"
complete_results [grep ("logFEve",rownames (complete_results)),"Index"] <- "FEve"
complete_results [grep ("logFDiv",rownames (complete_results)),"Index"] <- "FDiv"

# plot
dodge <- c(0.4,0.4)
pd <- position_dodge(dodge)
pdf_pt <- position_dodge(dodge)

a <- ggplot (complete_results[which (complete_results$Parameter == "Temperature"),], 
             
             aes  (y=Index, x=Estimate, 
                          colour=Organism, fill=Organism)) + 
  
  geom_errorbar(aes(xmin=lower95,xmax=upper95),width = 0.2,size=0.5,
                position=pd)  + 
  geom_errorbar(aes(xmin=lower90,xmax=upper90),width = 0.2,size=1,
                position=pd)  + 
  
  theme_classic() + 
  
  geom_point(position=(pdf_pt), 
             size=1.5)+ 
  
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray90", size=0.5)+
  
  #facet_wrap(~Algorithm+Index,scale="free",ncol=4) + 

  scale_color_manual(values=c("#f2a154", "#0e49b5")) + 
  
  xlab("Standardized effect size") + 
  
  ylab ("Metric") + 
  
  #xlim(-0.5,0.5) +
  
  theme(axis.text.x = element_text(angle = 45,size=7))

a

ggsave (file=here("output","vectorized","Fig2.pdf"),width=4,height=3)

# ---------------------------------------------
# correlation between metrics and CWM

# extract estimate of each trait level
# all levels
levels_cwm<-colnames(FD_fish_MSS[[1]]$cwm)

# remove problematic cwm estimates
prob<-(lapply (sapply(FD_fish_MSS,"[[","cwm"),nrow)) # null is prob
FD_fish_MSS <- FD_fish_MSS[which(unlist(lapply(prob,is.null)) == F)]

# extract cwm
cwm_fish <- lapply (levels_cwm, function (k) # for eeach trait level
  
  do.call(cbind, # cbind estimate of cwm for each trait level
          lapply (FD_fish_MSS, function (i) # in each rarefaction dataset
    
            i$cwm[,which(colnames(i$cwm) %in% k)] # extract estimate of CWM
      )
    )
  )

# obtain the average cwm across rarefaction datasets and trait levels
mean_cwm_fish <- do.call (cbind, 
                          lapply(cwm_fish,rowMeans)#  average across rarefied datasets
  )
colnames(mean_cwm_fish) <- levels_cwm

# the correlation between  functional metrics and cwm
correlation_metrics_cwm <- cor (cbind (cov_fish[,c("FRic","FEve", "FDiv")],
                mean_cwm_fish)
)

# the three first rows matters here
# also eliminating the first three cols (the metrics)
correlation_metrics_cwm <- correlation_metrics_cwm[1:3,-c(1:3)]
# ordering to get the highest values of correlation,if you want positive, and selecting the ntraits required

ntraits <- 2

#
corr_FRic <- correlation_metrics_cwm[1,][order(correlation_metrics_cwm[1,],decreasing=T)] [1:ntraits]
corr_FRic_neg <- correlation_metrics_cwm[1,][order(correlation_metrics_cwm[1,],decreasing=F)] [1:ntraits]
corr_FEve <- correlation_metrics_cwm[2,][order(correlation_metrics_cwm[2,],decreasing=T)] [1:ntraits]
corr_FEve_neg <- correlation_metrics_cwm[2,][order(correlation_metrics_cwm[2,],decreasing=F)] [1:ntraits]
corr_FDiv <- correlation_metrics_cwm[3,][order(correlation_metrics_cwm[3,],decreasing=T)] [1:ntraits]
corr_FDiv_neg <- correlation_metrics_cwm[3,][order(correlation_metrics_cwm[3,],decreasing=F)] [1:ntraits]

# now obtaining a df with the variation in cwm across the richness gradient (as done for each functional metric)
# first cbind richness and cwm correlated to each functional metric
# correlated to FRic
df_FRic_cwm <- data.frame (EstRich = cov_fish$BO2_tempmean_ss,
                          mean_cwm_fish[,which(colnames(mean_cwm_fish) %in% c(names(corr_FRic),names(corr_FRic_neg)))]
       )
df_FRic_cwm <- cbind(melt(df_FRic_cwm,"EstRich"),
                     Group="Fish",
                     Metric = "FRic")
                     
# correlated to FEve
df_FEve_cwm <- data.frame (EstRich = cov_fish$BO2_tempmean_ss,
                           mean_cwm_fish[,which(colnames(mean_cwm_fish) %in% c(names(corr_FEve),names(corr_FEve_neg)))]
)
df_FEve_cwm <- cbind(melt(df_FEve_cwm,"EstRich"),
                     Group="Fish",
                     Metric = "FEve")

# correlated to FDiv
df_FDiv_cwm <- data.frame (EstRich = cov_fish$BO2_tempmean_ss,
                           mean_cwm_fish[,which(colnames(mean_cwm_fish) %in% c(names(corr_FDiv),names(corr_FDiv_neg)))]
)
df_FDiv_cwm <- cbind(melt(df_FDiv_cwm,"EstRich"),
                     Group="Fish",
                     Metric = "FDiv")

# -----------------------------------------------                     
# do the same for benthos
# extract estimate of each trait level
# all levels
levels_cwm_benthos <- colnames(FD_benthos_MSS[[1]]$cwm)

# remove problematic cwm estimates
prob<-(lapply (sapply(FD_benthos_MSS,"[[","cwm"),nrow)) # null is prob
FD_benthos_MSS <- FD_benthos_MSS[which(unlist(lapply(prob,is.null)) == F)]

# extract cwm
cwm_benthos <- lapply (levels_cwm_benthos, function (k) # for eeach trait level
  
  do.call(cbind, # cbind estimate of cwm for each trait level
          lapply (FD_benthos_MSS, function (i) # in each rarefaction dataset
            
            i$cwm[,which(colnames(i$cwm) %in% k)] # extract estimate of CWM
          )
  )
)

# obtain the average cwm across rarefaction datasets and trait levels
mean_cwm_benthos <- do.call (cbind, 
                          lapply(cwm_benthos,rowMeans)#  average across rarefied datasets
)
colnames(mean_cwm_benthos) <- levels_cwm_benthos

# the correlation between  functional metrics and cwm
correlation_metrics_cwm_benthos <- cor (cbind (cov_benthos[,c("FRic","FEve", "FDiv")],
                                        mean_cwm_benthos)
)

# the three first rows matters here
# also eliminating the first three cols (the metrics)
correlation_metrics_cwm_benthos <- correlation_metrics_cwm_benthos[1:3,-c(1:3)]
# ordering to get the highest values
corr_FRic_benthos <- correlation_metrics_cwm_benthos[1,][order(correlation_metrics_cwm_benthos[1,],decreasing=T)] [1:ntraits]
corr_FRic_benthos_neg <- correlation_metrics_cwm_benthos[1,][order(correlation_metrics_cwm_benthos[1,],decreasing=F)] [1:ntraits]
corr_FEve_benthos <- correlation_metrics_cwm_benthos[2,][order(correlation_metrics_cwm_benthos[2,],decreasing=T)] [1:ntraits]
corr_FEve_benthos_neg <- correlation_metrics_cwm_benthos[2,][order(correlation_metrics_cwm_benthos[2,],decreasing=F)] [1:ntraits]
corr_FDiv_benthos <- correlation_metrics_cwm_benthos[3,][order(correlation_metrics_cwm_benthos[3,],decreasing=T)] [1:ntraits]
corr_FDiv_benthos_neg <- correlation_metrics_cwm_benthos[3,][order(correlation_metrics_cwm_benthos[3,],decreasing=F)] [1:ntraits]

# now plotting the variation in cwm across the richness gradient (as done for each functional metric)
# first cbind richness and cwm correlated to each functional metric
# correlated to FRic
df_FRic_cwm_benthos <- data.frame (EstRich = cov_benthos$BO2_tempmean_ss,
                                  mean_cwm_benthos[,which(colnames(mean_cwm_benthos) %in% c(names(corr_FRic_benthos),names(corr_FRic_benthos_neg)))]
)
df_FRic_cwm_benthos <- cbind(melt(df_FRic_cwm_benthos,"EstRich"),
                             Group="Benthos",
                             Metric = "FRic")
                             
# correlated to FEve
df_FEve_cwm_benthos <- data.frame (EstRich = cov_benthos$BO2_tempmean_ss,
                                  mean_cwm_benthos[,which(colnames(mean_cwm_benthos) %in% c(names(corr_FEve_benthos),names(corr_FEve_benthos_neg)))]
)
df_FEve_cwm_benthos <- cbind(melt(df_FEve_cwm_benthos,"EstRich"),
                             Group="Benthos",
                             Metric = "FEve")
                             
# correlated to FDiv
df_FDiv_cwm_benthos <- data.frame (EstRich = cov_benthos$BO2_tempmean_ss,
                                 mean_cwm_benthos[,which(colnames(mean_cwm_benthos) %in% c(names(corr_FDiv_benthos),names(corr_FDiv_benthos_neg)))]
)
df_FDiv_cwm_benthos <- cbind(melt(df_FDiv_cwm_benthos,"EstRich"),
                             Group="Benthos",
                             Metric = "FDiv")

# now rbind all these data
df_plot_cwm <- rbind(df_FRic_cwm,
      df_FRic_cwm_benthos,
      df_FEve_cwm,
      df_FEve_cwm_benthos,
      df_FDiv_cwm,
      df_FDiv_cwm_benthos)


# adjusting levels of the factor metric
df_plot_cwm$Metric<-factor(df_plot_cwm$Metric,
                           levels=c("FRic",
                                    "FEve",
                                    "FDiv"))

# pasting a signal  to differentiate positive from negative correlations
df_plot_cwm$variable <- ifelse (df_plot_cwm$variable %in% c(names(corr_FRic_neg),
                                    names(corr_FEve_neg),
                                    names(corr_FDiv_neg),
                                    names(corr_FRic_benthos_neg),
                                    names(corr_FEve_benthos_neg),
                                    names(corr_FDiv_benthos_neg)),
        paste(df_plot_cwm$variable,"(-)",sep=""), 
        paste(df_plot_cwm$variable,"(+)",sep=""))

# changing the names of trait levels
df_plot_cwm$variable<-gsub ("body_size_","Size.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Size_","Size.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Mouth_position_","Mouth.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Mobility_","Mobility.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("mobility_M","Mobile",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("mobility_S","Sessile",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Body_shape_","BShape.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("reproductive_mode_S","Sexual",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("growth_form_M","Massive",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("growth_form_B","Branching",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("growth_form_EF","Encr/filam",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("growth_form_E","Encrusting",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("carbonate_C","Add",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("carbonate_N","Not add",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Caudal_fin_","Fin.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Schooling_","Schooling.",(df_plot_cwm$variable))

# plot

# split data for each "facet"
leg_facet <- split(df_plot_cwm,f = df_plot_cwm$Group)

# plot
p1<-ggplot(df_plot_cwm, aes(x=EstRich,y=value,
                            group = variable,
                            fill=Group,
                            color=variable))+
  #geom_smooth(alpha=0.1,method="lm",se=T) +  #"loess"
  theme_classic()+
  theme(legend.text=element_text(size=8),
        legend.key.width= unit(2, 'cm'),
        axis.text = element_text(size=8),
        axis.title = element_text(size=11),
        strip.text = element_text(size=14),
        strip.background = element_rect(
          color="white", fill="white", size=1.5, linetype="solid"
        ))+
  facet_wrap(~Metric)+
  xlab("Temperature") + 
  ylab ("Trait frequency (CWM)") +
  scale_size(range=c(1, 1), guide=FALSE)


# do it for each "facet"
# benthos
p2 <- p1 %+% leg_facet$Benthos + 
  
  binomial_smooth(
    formula = y ~ splines::ns(x, 2),
    se=T,
    #size=1,
    fill = "orange",
    alpha=0.02,
    aes(linetype=variable)) + 
  
  scale_color_viridis_d(option = "inferno") + 
  
  theme(axis.title.x = element_blank())


# fish
p3 <- p1 %+% leg_facet$Fish +
  
  binomial_smooth(
    formula = y ~ splines::ns(x, 2),
    se=T,
    size=1,
    fill = "blue",
    alpha=0.02,
    aes(linetype=variable)) + 
  
  scale_color_viridis_d(option = "cividis")

# here the plot
pdf(here ("output","vectorized","cwm_richness.pdf"),height=7,width=10)

grid.arrange(p2,p3)

dev.off()
