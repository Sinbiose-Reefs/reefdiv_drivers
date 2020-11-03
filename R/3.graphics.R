####################################################################
############ GRAPHIC REPRESENTATION OF THE RESULTS

## load packages
source("R/packages.R")
source("R/functions.R")

## load results
load (here("output", "results_FD_analyses.RData"))

## Load environmental variables
## Plots with environmental data
## building a DF with data

## encontrar a riqueza por regiao
# peixes
reg_riq <- rowSums(aggregate (fish[,-1], by=list(covariates_site$region[,'Region']),
           FUN=sum)>0)

names(reg_riq) <- c("Northeastern","Oceanic Islands","Southeastern")
## match site and reg richness
reg_riq <- reg_riq [match  (covariates_site$region[,"Region"],names (reg_riq))]

# bentos
reg_riq_bentos <- rowSums(aggregate (bent2[,-1], by=list(covariates_site$region[which(covariates_effort$effort$n_videos_bentos>0),'Region']),
                              FUN=sum)>0)

names(reg_riq_bentos) <- c("Northeastern","Oceanic Islands","Southeastern")
## match site and reg richness
reg_riq_bentos <- reg_riq_bentos [match  (covariates_site$region[which(covariates_effort$effort$n_videos_bentos>0),"Region"],names (reg_riq_bentos))]

## regional functional divesrity
# fish
reg_FD <- FD_results_f1$all$fd_total$FRic
names(reg_FD) <- c("Northeastern","Oceanic Islands","Southeastern")
reg_FD <- reg_FD [match  (covariates_site$region,names (reg_FD))]
# benthos
reg_FD_benthos <- FD_results_f1_bentos$all$fd_total$FRic
names(reg_FD_benthos) <- c("Northeastern","Oceanic Islands","Southeastern")
reg_FD_benthos <- reg_FD_benthos [match  (covariates_site$region[which(covariates_effort$effort$n_videos_bentos>0),"Region"],names (reg_FD_benthos))]

## data
df_results <- data.frame (Nspec = FD_results_f1$all$Fdindexes$nbsp/reg_riq,
                          FRic = FD_results_f1$all$Fdindexes$FRic/reg_FD,
                          #FEve = FD_results_f1$all$Fdindexes$FEve,
                          #FDiv = FD_results_f1$all$Fdindexes$FDiv,
                          #FE = fishes_FE/total_fishes_FE,
                          Group="Fishes")

# bind dataframe with benthos results
df_results_bentos <-data.frame (Nspec = FD_results_f1_bentos$all$Fdindexes$nbsp.bent/reg_riq_bentos,
                                 FRic = FD_results_f1_bentos$all$Fdindexes$FRic.bent/reg_FD_benthos,
                                 #FEve = FD_results_f1_bentos$all$Fdindexes$FEve.bent,
                                 #FDiv = FD_results_f1_bentos$all$Fdindexes$FDiv.bent,
                                #FE = bentos_FE/total_bentos_FE,
                                 Group="Benthos")
df_results_bentos$FRic[df_results_bentos$FRic>1]<-1
## melt these results
df_results <- melt(df_results)
df_results_bentos<- melt(df_results_bentos)

## and bind covariates
df_results <- cbind(df_results,
                    Region=covariates_site$region,
                    covariates_site$coord$coord_peixes,
                    decostand (covariates_site$sea_data[[1]], "standardize"))
df_results_bentos <- cbind(df_results_bentos,
                    Region=covariates_site$region[which(covariates_effort$effort$n_videos_bentos>0)],
                    covariates_site$coord$coord_bentos,
                    decostand (covariates_site$sea_data[[2]],"standardize"))

# warning due to rownames of covariates_site$sea_data
### testing correlation between variables

cor (
  df_results [which (df_results$variable == "Nspec"),c("Lat","BO2_tempmean_ss","BO2_temprange_ss","BO2_ppmean_ss",
                                                       "BO2_pprange_ss","BO2_salinitymean_ss","BO2_salinityrange_ss")])

## finding the main drivers through models

indexes <- c("Nspec", "FRic")

model_results <- lapply (indexes, function (i) {
  mod_lat <- gam(value ~ s(Lat*-1, bs="cr"), data=df_results[which(df_results$variable==i),])
  mod_temp <- gam(value ~ s(BO2_tempmean_ss, bs="cr"), data=df_results[which(df_results$variable==i),])
  mod_prod <- gam(value ~ s(BO2_ppmean_ss, bs="cr"), data=df_results[which(df_results$variable==i),])
  mod_sal <- gam(value ~ s(BO2_salinitymean_ss, bs="cr"), data=df_results[which(df_results$variable==i),])

  res <- list (summary (mod_lat),
               summary(mod_temp),
               summary(mod_prod),
               summary(mod_sal))
  
})

## peixes
### Latitude
pdf(file=here("output","fish.pdf"),height=4,width=6)
ggplot (df_results, aes (x=Lat*-1, y=value)) + 
  geom_point(aes (col=Region)) + 
  scale_color_manual(values=c("#FFC93C","#07689F","#A2D5F2"))+
  geom_smooth() +
  facet_wrap(~variable,scales="free",
             ncol=2,nrow=2) + 
  theme_bw() + 
  xlab("Latitude") + 
  ylab("Index") + 
  theme(legend.position = "top")

### Temperature
ggplot (df_results, aes (x=(BO2_tempmean_ss), y=value)) + 
  geom_point(aes (col=Region)) + 
  scale_color_manual(values=c("#FFC93C","#07689F","#A2D5F2"))+
  geom_smooth() +
  facet_wrap(~variable,scales="free",
             ncol=2,nrow=2) + 
  theme_bw() + 
  xlab("Average sea surface temperature (ºC)") + 
  ylab("Index")+ 
  theme(legend.position = "top")

### primary productivity
ggplot (df_results, aes (x=BO2_ppmean_ss, y=value)) + 
  geom_point(aes (col=Region)) + 
  scale_color_manual(values=c("#FFC93C","#07689F","#A2D5F2"))+
  geom_smooth() +
  facet_wrap(~variable,scales="free",
             ncol=2,nrow=2) + 
  theme_bw() + 
  xlab("Average primary productivity") + 
  ylab("Index")+ 
  theme(legend.position = "top")

### salinity productivity
ggplot (df_results, aes (x=BO2_salinitymean_ss, y=value)) + 
  geom_point(aes (col=Region)) + 
  scale_color_manual(values=c("#FFC93C","#07689F","#A2D5F2"))+
  geom_smooth() +
  facet_wrap(~variable,scales="free",
             ncol=2,nrow=2) + 
  theme_bw() + 
  xlab("Average salinity") + 
  ylab("Index")+ 
  theme(legend.position = "top")


dev.off()

##### benthos

model_results_bentos <- lapply (indexes, function (i) {
  mod_lat <- gam(value ~ s(Lat*-1, bs="cr"), data=df_results_bentos[which(df_results_bentos$variable==i),])
  mod_temp <- gam(value ~ s(BO2_tempmean_ss, bs="cr"), data=df_results_bentos[which(df_results_bentos$variable==i),])
  mod_prod <- gam(value ~ s(BO2_ppmean_ss, bs="cr"), data=df_results_bentos[which(df_results_bentos$variable==i),])
  mod_sal <- gam(value ~ s(BO2_salinitymean_ss, bs="cr"), data=df_results_bentos[which(df_results_bentos$variable==i),])
  res <- list (summary (mod_lat),
               summary(mod_temp),
               summary(mod_prod),
               summary(mod_sal))
  
})

#
ggplot (df_results_bentos, aes (x=Lat*-1, y=value)) + 
  geom_point(aes (col=Region)) + 
  scale_color_manual(values=c("#FFC93C","#07689F","#A2D5F2"))+
  geom_smooth() +
  facet_wrap(~variable,scales="free",
             ncol=2,nrow=2) + 
  theme_bw() + 
  xlab("Latitude") + 
  ylab("Index")+ 
  theme(legend.position = "top")


### Temperature
pdf(file=here("output","bent_temp.pdf"),height=4,width=6)
ggplot (df_results_bentos, aes (x=BO2_tempmean_ss, y=value)) + 
  geom_point(aes (col=Region)) + 
  scale_color_manual(values=c("#FFC93C","#07689F","#A2D5F2"))+
  geom_smooth() +
  facet_wrap(~variable,scales="free",
             ncol=2,nrow=2) + 
  theme_bw() + 
  xlab("Average sea surface temperature (ºC)") + 
  ylab("Index")+ 
  theme(legend.position = "top")
dev.off()

### primary productivity
ggplot (df_results_bentos, aes (x=BO2_ppmean_ss, y=value)) + 
  geom_point(aes (col=Region)) + 
  scale_color_manual(values=c("#FFC93C","#07689F","#A2D5F2"))+
  geom_smooth() +
  facet_wrap(~variable,scales="free",
             ncol=2,nrow=2) + 
  theme_bw() + 
  xlab("Average primary productivity") + 
  ylab("Index")+ 
  theme(legend.position = "top")

### salinity productivity
pdf(file=here("output","bent_salt.pdf"),height=4,width=6)
ggplot (df_results_bentos, aes (x=BO2_salinitymean_ss, y=value)) + 
  geom_point(aes (col=Region)) + 
  scale_color_manual(values=c("#FFC93C","#07689F","#A2D5F2"))+
  geom_smooth() +
  facet_wrap(~variable,scales="free",
             ncol=2,nrow=2) + 
  theme_bw() + 
  xlab("Average salinity") + 
  ylab("Index")+ 
  theme(legend.position = "top")


dev.off()

################################
### ric of each group and FRic relative to total
percent_rich <- df_results[which(df_results$variable == "Nspec"),c(1:3,5)]
percent_rich_FD <- cbind (percent_rich, 
                          FRic = df_results[which(df_results$variable == "FRic"),3])

percent_rich_bentos <- df_results_bentos[which(df_results_bentos$variable == "Nspec"),c(1:3,5)]
percent_rich_bentos <- cbind (percent_rich_bentos, 
                          FRic = df_results_bentos[which(df_results_bentos$variable == "FRic"),3])
## bind these data
percent_rich_FD <- rbind(percent_rich_FD,
  percent_rich_bentos)

percent_rich_FD$Region <- covariates_site$region [match (percent_rich_FD$Group.1, rownames (covariates_site$region)),"Region"]

pdf(here ("output", "scatterSR_FD.pdf"),width = 6,heigh=4)
## div tax vs functional
ggplot (percent_rich_FD, aes (x=value, y=FRic, group= Group,col=Group)) + 
  geom_point(aes (size=Region),alpha=0.3) + 
  geom_smooth(method="gam") +
  theme_classic() + 
  xlab("Species richness (%)") + 
  ylab("Functional diversity (%)") + 
  labs(
    colour = "Group",
    shape = "Region"
  ) + scale_colour_manual(
    
    values= c("Benthos" = "#00CC66",
              "Fishes" = "#990000")
  ) + ylim (0,1) + xlim (0,0.5) 

dev.off()

############################ 
## MAPPING

## creating a dataframe for each result
# fish
fish_richness <- do.call (cbind, 
                          lapply (seq(1,length(FD_results_f1)), function (i)
                            
                            FD_results_f1[[i]]$Fdindexes$nbsp
                            
                          ))
colnames (fish_richness) <- names(FD_results_f1)

# benthos
benthos_richness <- do.call (cbind, 
                             
                             lapply (seq(1,length(FD_results_f1_bentos)), function (i)
                               
                               FD_results_f1_bentos[[i]]$Fdindexes$nbsp.bent
                               
                             ))
colnames (benthos_richness) <- names(FD_results_f1_bentos)
benthos_richness <- data.frame(benthos_richness)
benthos_richness$aut <- benthos_richness$all - benthos_richness$nonAut 
benthos_richness[is.na(benthos_richness)] <- 0

## binding coordinates (fishes for now)
fish_richness <- cbind (fish_richness, covariates_site$coord$coord_peixes)
fish_richness[is.na(fish_richness)] <- 0
## jittering these coordinates
fish_richness$LonJitter <- jitter (fish_richness$Lon,factor=400)
fish_richness$LatJitter <- jitter (fish_richness$Lat,factor=600)

## just moving charts of benthic communities to the left
benthos_richness <- cbind(benthos_richness,
                          LatJitter = fish_richness$LatJitter[which(covariates_effort$effort$n_videos_bentos>0)],
                          LonJitter= fish_richness$LonJitter[which(covariates_effort$effort$n_videos_bentos>0)]+5)


#### aggregating sites
# fish
site_fish <- unlist(lapply (strsplit(covariates_site$site_names,"\\."), "[",1))
agg_fish_rich <- aggregate(fish_richness[,-which(colnames(fish_richness)=="Group.1")], 
          by = list(site_fish),
          FUN=mean)

# benthos
site_bentos <- unlist(lapply (strsplit(covariates_site$site_names,"\\."), "[",1))[which(covariates_effort$effort$n_videos_bentos>0)]
agg_benthos_rich <- aggregate(benthos_richness, 
                           by = list(site_bentos),
                           FUN=mean)

## maps
# mapa mundi
world <- ne_countries(scale = "medium", returnclass = "sf")

# cortar o mapa para ver a america do Sul e parte da central
wm <- ggplot() + 
  geom_sf (data=world, size = 0.1, 
           fill= "gray90",colour="gray90") +
  coord_sf (xlim = c(-50, -20),  ylim = c(-30, 2), expand = FALSE) +
  theme_bw() + #xlab ("Longitude")  + ylab ("Latitude") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "aliceblue",#darkslategray1
                                        colour = "aliceblue"),
        axis.text.x = element_text(size=6),
        axis.ticks.x=element_line(size=1),
        axis.text.y = element_text(size=6),
        axis.ticks.y=element_line(size=1),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        title = element_blank()) +
  xlab("Longitude") + ylab("Latitude")


## advise to jitter : https://stackoverflow.com/questions/52806580/pie-charts-in-geom-scatterpie-overlapping
## pie: http://www.spectdata.com/index.php/2018/10/25/how-to-use-ggplot-to-plot-pie-charts-on-a-map/

wm_pie <- wm + geom_scatterpie(aes(x=LonJitter, y=LatJitter, r=all/max(all)),alpha=0.5,
                               data = agg_fish_rich,
                               cols = c(
                                 "herbivorous",
                                 "sessileInv",
                                 "om",
                                 "fc",
                                 "pk"),
                               #pie_scale = 0.1,
                               size=0.5,
                               sorted_by_radius = F,
                               legend_name = "Groups") + 
  ggtitle("Fishes and benthos richness")

### adding benthos richness

wm_pie_col_benthos <- wm_pie + geom_scatterpie(aes(x=LonJitter, y=LatJitter,r=all/max(all,na.rm=T)),alpha=0.5,
                                               data = agg_benthos_rich,
                                               cols = c("aut",
                                                        "nonAut",
                                                        "corals"
                                               ),
                                               pie_scale = 1,
                                               size=0.5,
                                               sorted_by_radius = F,
                                               legend_name = "Benthos") + 
  
  theme (legend.title = element_text(size=7),
         legend.text = element_text(size=7),
         legend.position = c(0.36, 0.9),
         legend.justification = c("right", "top"),
         legend.box.just = "right",
         legend.margin = margin(6,6,6,6),
         legend.background = element_blank(),
         title=element_text(size=10)) 


wm_pie_col_benthos <- wm_pie_col_benthos + scale_fill_manual(
  
  breaks=c("aut", 
           "nonAut",
           "corals",
           "herbivorous",
           "om",
           "pk",
           "sessileInv",
           "fc"),
  labels=c("Algae", 
           "Filter feeders, carnivores, grazers",
           "Corals",
           "Herbivores",
           "Omnivores",
           "Planktivores",
           "Invertivores",
           "Piscivores"),
  values= c("aut" = "#003300",
            "nonAut" = "#00CC66",
            "corals" = "#B2FF66",
            "herbivorous" = "#132743",
            "om" = "#F8EFD4",
            "pk" = "#FE7171",
            "sessileInv" = "#F0A500",
            "fc" = "#931A25")
)

wm_pie_col_benthos <- wm_pie_col_benthos + geom_scatterpie_legend(fish_richness$all/max(fish_richness$all), 
                                                                  x=-30, y=-28, n=3, 
                                                                  labeller=function(x) round(x*max(fish_richness$all),2))
wm_pie_col_benthos <- wm_pie_col_benthos + geom_scatterpie_legend(benthos_richness$all/max(benthos_richness$all,na.rm=T), 
                                                                  x=-25, y=-28, n=3, 
                                                                  labeller=function(x) round(x*max(benthos_richness$all,na.rm=T),2))

pdf(here ("output", "MapSR.pdf"),width = 10,heigh=7)
wm_pie_col_benthos
dev.off()

### FRic

## creating a dataframe for each result
# fish
fish_fric <- data.frame (FRic = FD_results_f1[[1]]$Fdindexes$FRic,
                         FEve = FD_results_f1[[1]]$Fdindexes$FEve,
                         FDiv = FD_results_f1[[1]]$Fdindexes$FDiv)

# benthos
benthos_fric <- data.frame (FRic = FD_results_f1_bentos[[1]]$Fdindexes$FRic.bent,
                            FEve = FD_results_f1_bentos[[1]]$Fdindexes$FEve.bent,
                            FDiv = FD_results_f1_bentos[[1]]$Fdindexes$FDiv.bent)
benthos_fric[is.na(benthos_fric)] <- 0

## aggregate data 
agg_fish_FRic <- aggregate(fish_fric, 
                              by = list(site_fish),
                              FUN=mean)
agg_benthos_FRic <- aggregate(benthos_fric, 
                              by = list(site_bentos),
                              FUN=mean)

## binding already jiterred coordinates (fishes)
agg_fish_FRic <- cbind (agg_fish_FRic, 
                    #richness = fish_richness$all,
                    LonJitter=agg_fish_rich$LonJitter,
                    LatJitter = agg_fish_rich$LatJitter)
## binding already jiterred coordinates (benthos)
agg_benthos_FRic <- cbind (agg_benthos_FRic, 
                        #richness = fish_richness$all,
                        LonJitter=agg_benthos_rich$LonJitter,
                        LatJitter = agg_benthos_rich$LatJitter)

## plot 

wm_pie <- wm + geom_point (data = agg_fish_FRic, aes(x=LonJitter, y = LatJitter, size=FRic),
                           col = "#990000", alpha= 0.5) + 
                scale_size(range=c(2,8))

wm_pie_FRIC <- wm_pie + geom_point (data = agg_benthos_FRic, aes(x=LonJitter, y = LatJitter, 
                                                                 size=FRic),
                               col="#00CC66",alpha=0.5) + 
            ggtitle ("Functional richness (FRic)") + 
      
    theme (legend.title = element_text(size=8),
           legend.text = element_text(size=8),
           legend.position = c(0.25, 0.9),
           legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6,6,6,6),
          legend.background = element_blank(),
          title=element_text(size=10)) + 
  scale_size(range=c(2,8))

wm_pie_FRIC

### 

## selecao baseada no bentos (numbero de especies de bentos em cada grupo)
sel_sites <- covariates_effort$effort$locality_site[which(covariates_effort$effort$n_videos_bentos>0)]#19,,
sel_sites <- sel_sites [c (32,30,12,3,16,6,21)]
# subset coordenadas
subset_coord <- covariates_site$coord$coord_peixes[which(covariates_site$coord$coord_peixes$Group.1 %in% sel_sites),]
subset_coord <- subset_coord [order (subset_coord$Lat,decreasing=T),]

## annotate in the plot
wm_pie_FRIC_space <- wm_pie_FRIC + 
  geom_segment(aes(x=subset_coord$Lon,y=subset_coord$Lat,
                xend=-20,yend=subset_coord$Lat+0.1),
                color = "gray", size=1,alpha=0.5,
               ) + 
 annotate("text", x = subset_coord$Lon[1]+11, y = subset_coord$Lat[1]+1.0, label = "Rocas Atoll",size=3) + 
  annotate("text", x = subset_coord$Lon[2]+11, y = subset_coord$Lat[2]-0.5, label = "Rio Grande do Norte",size=3)+
  annotate("text", x = subset_coord$Lon[3]+12.2, y = subset_coord$Lat[3]-0.5, label = "Coral Coast",size=3)+
  annotate("text", x = subset_coord$Lon[4]+16, y = subset_coord$Lat[4]+0.8, label = "Abrolhos",size=3)+
  annotate("text", x = subset_coord$Lon[5]+17, y = subset_coord$Lat[5]-0.5, label = "Espírito Santo",size=3)+
  annotate("text", x = subset_coord$Lon[6]+18.5, y = subset_coord$Lat[6]-0.6, label = "Arraial do Cabo",size=3)+
  annotate("text", x = subset_coord$Lon[7]+25.2, y = subset_coord$Lat[7]+1.1, label = "Santa Catarina",size=3)

wm_pie_FRIC_space

#################################################

### groups of fish species based on diet
fish_traits_subset <- lapply (list_fish_diet, function (i)
  fish_traits [which (fish_traits$Diet %in% i),])

# matching
subset_PCO <- lapply (fish_traits_subset, function (i) 
  (rownames(FD_results_f1$all$axesPCO) %in% i$Name))

## subsetting composition data
site_composition <- lapply (sel_sites , function (i) {
  
 subs1 <- fish [which(fish$locality_site == i),]
 subs2 <- subs1 [,which(subs1 >0)];
 subs2
  
})
  
## quais sp estao nos sitios selecionados
site_subset_composition <- lapply (site_composition, function (i)
  which(rownames(FD_results_f1$all$axesPCO) %in% colnames(i))
)

# subset dos eixos da ordenacao
axes_spp_composition <- lapply (site_subset_composition, function (i)
  FD_results_f1$all$axesPCO[i,1:2]
)

## convex hull data
hull.data <- lapply (seq (1,length(axes_spp_composition )), function (i) 
  cbind(axes_spp_composition[[i]] , 
        herb=subset_PCO$herbivorous[site_subset_composition[[i]]],
        inv = subset_PCO$sessileInv[site_subset_composition[[i]]],
        om=subset_PCO$om[site_subset_composition[[i]]],
        fc=subset_PCO$fc[site_subset_composition[[i]]],
        pk=subset_PCO$pk[site_subset_composition[[i]]])
)

## get plots 
vol_plot_fish <- lapply (seq (1,length(hull.data)), function (i) {
  
  hull.data <- fortify(hull.data[[i]])
  ## todas as spp
  a <- hull.data [chull(hull.data[,1:2], y = NULL),]
  ## herb
  herb <-hull.data [which(hull.data$herb==T),]
  b <- herb [chull(herb, y = NULL),]
  inv <-hull.data [which(hull.data$inv==T),]
  c <- inv [chull(inv, y = NULL),]
  om <-hull.data [which(hull.data$om==T),]
  d <- om [chull(om, y = NULL),]
  fc <-hull.data [which(hull.data$fc==T),]
  e <- fc [chull(fc, y = NULL),]
  pk <-hull.data [which(hull.data$pk==T),]
  f <- pk [chull(pk, y = NULL),]
  
  ## hull for complete data
  
  vol_plot <- ggplot(axes_spp_composition[[i]], aes(A1, A2)) + 
    geom_point() + theme_bw()+
    geom_polygon(data=a, aes (A1,A2),alpha=0.1,fill="#808080") + 
    geom_polygon(data=b, aes (A1,A2,group=herb, fill=herb),fill="#132743",alpha=0.5)+
    geom_polygon(data=c, aes (A1,A2,group=inv, fill=inv),fill="#F0A500",alpha=0.5)+
    geom_polygon(data=d, aes (A1,A2,group=om, fill=om),fill="#F8EFD4",alpha=0.5) +
    geom_polygon(data=e, aes (A1,A2,group=fc, fill=fc),fill="#931A25",alpha=0.5) +
    geom_polygon(data=f, aes (A1,A2,group=pk, fill=pk),fill="#FE7171",alpha=0.5) +
    
     theme (axis.title = element_blank(),
           axis.text = element_blank())
  ;
  vol_plot
})


##
### groups of benthos species based on diet
benthos_traits_subset <- lapply (list_benthos_level, function (i)
  bent_traits [which (bent_traits$trophic_type %in% i),])

# matching
subset_PCO_benthos <- lapply (benthos_traits_subset, function (i) 
  (rownames(FD_results_f1_bentos$all$axesPCO) %in% i$groups))

## subsetting composition data from previously selected sites
site_subset_composition_benthos <- lapply (sel_sites , function (i) {
  
  subs1 <- bent2 [which(bent2$locality_site == i),]
  subs2 <- subs1 [,which(subs1 >0)];
  subs2
  
})

## quais sp estao nos sitios selecionados
site_subset_composition_benthos <- lapply (site_subset_composition_benthos, function (i)
  which(rownames(FD_results_f1_bentos$all$axesPCO) %in% colnames(i))
)

# subset dos eixos da ordenacao
axes_spp_composition_benthos <- lapply (site_subset_composition_benthos, function (i)
  FD_results_f1_bentos$all$axesPCO[i,1:2]
)

## convex hull data
hull.data_benthos <- lapply (seq (1,length(axes_spp_composition_benthos)), function (i) 
  cbind(axes_spp_composition_benthos[[i]] , 
        aut = subset_PCO_benthos$aut[site_subset_composition_benthos[[i]]],
        nonAut = subset_PCO_benthos$nonAut[site_subset_composition_benthos[[i]]],
        corals = subset_PCO_benthos$corals[site_subset_composition_benthos[[i]]])
)

## get plots 
vol_plot_benthos <- lapply (seq (1,length(hull.data_benthos)), function (i) {
  
  hull.data <- fortify(hull.data_benthos[[i]])
  ## todas as spp
  a <- hull.data [chull(hull.data[,1:2], y = NULL),]
  ## herb
  aut <-hull.data [which(hull.data$aut==T),]
  b <- aut [chull(aut, y = NULL),]
  nonAut <-hull.data [which(hull.data$nonAut==T),]
  c <- nonAut [chull(nonAut, y = NULL),]
  corals <-hull.data [which(hull.data$corals==T),]
  d <- corals [chull(corals, y = NULL),]
  
  ## hull for complete data
  
  vol_plot <- ggplot(axes_spp_composition_benthos[[i]], aes(A1, A2)) + 
    geom_point() + theme_bw()+
    geom_polygon(data=a, aes (A1,A2),alpha=0.1, fill="#808080") + 
    geom_polygon(data=b, aes (A1,A2,group=aut, fill=aut),fill="#003300",alpha=0.5)+
  geom_polygon(data=c, aes (A1,A2,group=nonAut, fill=nonAut),fill="#00CC66",alpha=0.5)+
  geom_polygon(data=d, aes (A1,A2,group=corals, fill=corals),fill="#B2FF66",alpha=0.5) +
    
    theme (axis.title = element_blank(),
           axis.text = element_blank()) 
    
    
  ;
  vol_plot
})


## colocar titulo nos dois primeiros graficos
vol_plot_fish[[1]] <- vol_plot_fish[[1]] + ggtitle ("Fishes")+
  theme (title = element_text(size=8))

vol_plot_benthos[[1]] <- vol_plot_benthos[[1]]+ggtitle("Benthos") +
  theme (title = element_text(size=8))

## panel with plot and map
pdf(here ("output", "MapFD.pdf"),width = 10,heigh=7)
grid.arrange(wm_pie_FRIC_space, 
             vol_plot_fish[[1]],vol_plot_benthos[[1]],
             vol_plot_fish[[2]],vol_plot_benthos[[2]],
             vol_plot_fish[[3]],vol_plot_benthos[[3]],
             vol_plot_fish[[4]],vol_plot_benthos[[4]],
             vol_plot_fish[[5]],vol_plot_benthos[[5]],
             vol_plot_fish[[6]],vol_plot_benthos[[6]],
             vol_plot_fish[[7]],vol_plot_benthos[[7]],
             ncol=7,nrow=7,
             layout_matrix = rbind (c(1,1,1,1,1,2,3),
                                    c(1,1,1,1,1,4,5),
                                    c(1,1,1,1,1,6,7),
                                    c(1,1,1,1,1,8,9),
                                    c(1,1,1,1,1,10,11),
                                    c(1,1,1,1,1,12,13),
                                    c(1,1,1,1,1,14,15)))
dev.off()




# largura




### de molho


### RELATIONSHIP BETWEEN FRIC OF BENTHOS AND FISHES
## FRIC

whole_data <- lapply (seq (1,length(FD_results_f1)), function (k)
  do.call (rbind, lapply (seq (1,length(FD_results_f1_bentos)), function (i)
    
    data.frame (
      RicBenthos = FD_results_f1_bentos[[i]]$Fdindexes$nbsp.bent,
      RicFishes = FD_results_f1[[k]]$Fdindexes$nbsp,
      RareRicFishes = rarefied_richness[[k]],
      FRicBenthos = FD_results_f1_bentos[[i]]$Fdindexes$FRic,
      FRicFishes = FD_results_f1[[k]]$Fdindexes$FRic,
      FEveBenthos = FD_results_f1_bentos[[i]]$Fdindexes$FEve.bent,
      FEveFishes = FD_results_f1[[k]]$Fdindexes$FEve,
      FDivBenthos = FD_results_f1_bentos[[i]]$Fdindexes$FDiv.bent,
      FDivFishes = FD_results_f1[[k]]$Fdindexes$FDiv,
      Group= names(FD_results_f1)[k],
      Benthos = names(list_benthos_level)[i])
  )))

## data for relationship between indexes
data_corr_index <- do.call (rbind, whole_data)

## Relationship between FRic of Benthos and Fishes
ggplot (data_corr_index, aes (x=FRicBenthos,y=FRicFishes)) + 
  geom_point() + 
  geom_smooth()+
  facet_wrap(~ Benthos+Group, scales = "free") + 
  theme_light()

## Relationship between FRic and richness of Benthos and Fishes
# Rich benthos, FRic Fishes
ggplot (data_corr_index, aes (x=RicBenthos,y=FRicFishes)) + 
  geom_point() + 
  geom_smooth()+
  facet_wrap(~ Benthos+Group, scales = "free") + 
  theme_light()

# Rich benthos, rich fishes
ggplot (data_corr_index, aes (x=RicBenthos,y=RareRicFishes)) + 
  geom_point() + 
  geom_smooth()+
  facet_wrap(~ Benthos+Group, scales = "free") + 
  theme_light()

# FEve benthos,  Feve fishes
ggplot (data_corr_index, aes (x=FEveBenthos,y=FEveFishes)) + 
  geom_point() + 
  geom_smooth()+
  facet_wrap(~ Benthos+Group, scales = "free") + 
  theme_light()

# FDiv benthos,  FRic fishes
ggplot (data_corr_index, aes (x=FDivBenthos,y=FRicFishes)) + 
  geom_point() + 
  geom_smooth()+
  facet_wrap(~ Benthos+Group, scales = "free") + 
  theme_light()

# FRic benthos,  Feve benthos
ggplot (data_corr_index, aes (x=FRicBenthos,y=FEveBenthos)) + 
  geom_point() + 
  geom_smooth()+
  facet_wrap(~ Benthos, scales = "free") + 
  theme_light()
# FRic benthos,  FDiv benthos
ggplot (data_corr_index, aes (x=FRicBenthos,y=FDivBenthos)) + 
  geom_point() + 
  geom_smooth()+
  facet_wrap(~ Benthos, scales = "free") + 
  theme_light()
# FRic fishes,  Feve fishes
ggplot (data_corr_index, aes (x=FRicFishes,y=FEveFishes)) + 
  geom_point() + 
  geom_smooth()+
  facet_wrap(~ Group, scales = "free") + 
  theme_light()
# FRic fishes,  FDiv fishes
ggplot (data_corr_index, aes (x=FRicFishes,y=FDivFishes)) + 
  geom_point() + 
  geom_smooth()+
  facet_wrap(~ Group, scales = "free") + 
  theme_light()





### Latitude and FRIC
lat_fric <- ggplot (df_results [which (df_results$variable == "FRic"),], aes (x=Lat, y=value,
                                                                              col=Group)) + 
  geom_point() + 
  geom_smooth() +
  scale_color_manual(values=c('#009900','#ca0020'))+
  theme_classic() + 
  ylab(NULL) + 
  coord_flip() + 
  theme_classic()+
  ggtitle ("FRic")+
  
  theme (legend.title = element_text(size=6),
         legend.text = element_text(size=5),
         legend.position = c(1, 0.8),
         legend.justification = c("right", "top"),
         legend.box.just = "right",
         legend.margin = margin(6,6,6,6),
         legend.background = element_blank(),
         axis.title.x = element_text(size=8),
         axis.text.x = element_text(size=6),
         axis.title.y =  element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         title=element_text(size=7),
         panel.background = element_rect(fill = 'aliceblue', 
                                         colour = 'aliceblue'),
         plot.margin = unit(c(0,-0.8,0,0.3), "cm")) 