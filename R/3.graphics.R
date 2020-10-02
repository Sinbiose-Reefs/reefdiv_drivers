####################################################################
############ GRAPHIC REPRESENTATION OF THE RESULTS

## load packages
source("R/packages.R")
source("R/functions.R")

## load results
load (here("output", "results_FD_analyses.RData"))

## relationship between richness and effort
# benthos
par(mfrow=c(2,2),mar = c (4,4,2,1))
plot(covariates_effort$effort$n_videos_bentos,FD_results_f1_bentos[[1]]$Fdindexes$nbsp.bent,
     ylab= "Number of species",
     xlab= "Number of videos per site",
     main = "Benthos",pch=19,col="coral",
     cex.axis=0.6,cex.lab=0.8)
abline(lm (FD_results_f1_bentos[[1]]$Fdindexes$nbsp.bent~covariates_effort$effort$n_videos_bentos),
       lwd=2,col="gray70")
# summary(lm (FD_results_f1_bentos$all$Fdindexes$nbsp.bent~covariates_effort$effort$n_videos_bentos))
text (x=17.5,y=17, labels=expression (paste("R"^2, "=-0.021")),cex=0.7)

# fishes
plot(covariates_effort$effort$n_transectos_peixes,FD_results_f1[[1]]$Fdindexes$nbsp,
     xlab= "Number of transects per site",
     ylab= "",
     main = "Fishes",pch=19,col="cyan3",
     cex.axis=0.6,cex.lab=0.8)
abline(lm (FD_results_f1[[1]]$Fdindexes$nbsp~covariates_effort$effort$n_transectos_peixes),
       lwd=2,col="gray70")
# summary(lm (FD_results_f1$all$Fdindexes$nbsp~covariates_effort$effort$n_transectos_peixes))
text (x=80,y=65, labels=expression (paste("R"^2, "=0.15**")),cex=0.7)

### rarefaction for richness of fishes and benthos
## 286 was sample size suggested by the function; lower than the minimum row max across sites (max abundance)
rarefied_richness <- lapply (list_fish_diet, function (i)  {
  
  fish_traits_subset <- fish_traits [which (fish_traits$Diet %in% i),]
  # matching
  fish_data <- fish [, which(colnames(fish) %in% fish_traits_subset$Name)]
  ## rarefy
  rarefied_richness <- rarefy (fish_data,100)
  ;
  rarefied_richness
}
)

###
# sample-based rarefaction
# transects needed to sample local richness

#
plot(covariates_effort$effort$n_transectos_peixes,rarefied_richness[[1]],
     ylab= "Rarefied fish richness",
     xlab = "Number of transects per site",
     cex.axis=0.6,cex.lab=0.8,
     pch=19, col = "blue")
abline (lm(rarefied_richness[[1]] ~ covariates_effort$effort$n_transectos_peixes),
        lwd=2,col="gray70")
## summary (lm(rarefied_richness[[1]] ~ covariates_effort$effort$n_transectos_peixes))
text (x=80,y=27, labels=expression (paste("R"^2, "=-0.01")),cex=0.7)

## Load environmental variables
## Plots with environmental data
## building a DF with data

df_results <- data.frame (Nspec = FD_results_f1$all$Fdindexes$nbsp/max(FD_results_f1$all$Fdindexes$nbsp,na.rm=T),
                          FRic = FD_results_f1$all$Fdindexes$FRic/max(FD_results_f1$all$Fdindexes$FRic, na.rm=T),
                          FEve = FD_results_f1$all$Fdindexes$FEve,
                          FDiv = FD_results_f1$all$Fdindexes$FDiv,
                          Group="Fishes")
# bind dataframe with benthos results
df_results <- rbind (df_results, 
                     data.frame (Nspec = FD_results_f1_bentos$all$Fdindexes$nbsp.bent/max(FD_results_f1_bentos$all$Fdindexes$nbsp.bent,na.rm=T),
                                 FRic = FD_results_f1_bentos$all$Fdindexes$FRic.bent/max(FD_results_f1_bentos$all$Fdindexes$FRic.bent,na.rm=T),
                                 FEve = FD_results_f1_bentos$all$Fdindexes$FEve.bent,
                                 FDiv = FD_results_f1_bentos$all$Fdindexes$FDiv.bent,
                                 Group="Benthos"))
## melt these results
df_results <- melt(df_results)
## and bind covariates
df_results <- cbind(df_results,
                    ReefType=covariates_site$biog_reef,
                    Region=covariates_site$region,
                    covariates_site$coord,
                    covariates_site$sea_data)
df_results$Region[which(df_results$Region == "sud")] <- "Southeastern"
df_results$Region[which(df_results$Region == "nord")] <- "Northeastern"
df_results$Region[which(df_results$Region == "oc.isl")] <- "Oceanic Islands"

# warning due to rownames of covariates_site$sea_data

### testing correlation between variables

cor (
  df_results [which (df_results$variable == "Nspec"),c("BO2_tempmean_ss","BO2_temprange_ss","BO2_ppmean_ss",
                                                       "BO2_pprange_ss","BO2_salinitymean_ss","BO2_salinityrange_ss",
                                                       "BO2_chlomean_ss", "BO2_chlorange_ss")]
)

### Latitude
ggplot (df_results, aes (x=Lat, y=value)) + 
  geom_point() + 
  geom_smooth() +
  facet_wrap(~Group+variable,scales="free",
             ncol=4,nrow=2) + 
  theme_bw() + 
  xlab("Latitude") + 
  ylab("Index")

### Temperature
ggplot (df_results, aes (x=BO2_tempmean_ss, y=value)) + 
  geom_point() + 
  geom_smooth() +
  facet_wrap(~Group+variable,scales="free",
             ncol=4,nrow=2) + 
  theme_bw() + 
  xlab("Average sea surface temperature (ºC)") + 
  ylab("Index")

### primary productivity
ggplot (df_results, aes (x=BO2_ppmean_ss, y=value)) + 
  geom_point() + 
  geom_smooth() +
  facet_wrap(~Group+variable,scales="free",
             ncol=4,nrow=2) + 
  theme_bw() + 
  xlab("Average primary productivity") + 
  ylab("Index")

### ric of each group and FRic
percent_rich <- df_results[which(df_results$variable == "Nspec"),c(1:3,5)]
percent_rich_FD <- cbind (percent_rich, 
                          FRic = df_results[which(df_results$variable == "FRic"),3])

pdf(here ("output", "scatterSR_FD.pdf"),width = 7,heigh=5)
## div tax vs functional
ggplot (percent_rich_FD, aes (x=value, y=FRic, group= Group,col=Group)) + 
  geom_point(aes (size=Region)) + 
  geom_smooth(method="glm") +
  theme_classic() + 
  xlab("Species richness (%)") + 
  ylab("Functional diversity (%)") + 
  labs(
    colour = "Group",
    shape = "Region"
  ) + scale_colour_manual(
    
    values= c("Benthos" = "#00CC66",
              "Fishes" = "#990000")
  )

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
fish_richness <- cbind (fish_richness, covariates_site$coord)
## jittering these coordinates
fish_richness$LonJitter <- jitter (fish_richness$Lon,factor=400)
fish_richness$LatJitter <- jitter (fish_richness$Lat,factor=600)

## just moving charts of benthic communities to the left
benthos_richness <- cbind(benthos_richness,
                          LatJitter = fish_richness$LatJitter,
                          LonJitter=fish_richness$LonJitter+6.5)

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
                               data = fish_richness,
                               cols = c(
                                 "herbivorous",
                                 "sessileInv",
                                 "pred"),
                               #pie_scale = 0.1,
                               size=0.5,
                               sorted_by_radius = F,
                               legend_name = "Groups") + 
  ggtitle("Fishes and benthos richness")

### adding benthos richness

wm_pie_col_benthos <- wm_pie + geom_scatterpie(aes(x=LonJitter, y=LatJitter,r=all/max(all,na.rm=T)),alpha=0.5,
                                               data = benthos_richness,
                                               cols = c("aut",
                                                        "nonMixAut",
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
  
  breaks=c("aut", "nonMixAut","corals",
           "pred","sessileInv","herbivorous"),
  labels=c("Algae", 
           "Filter feeders, carnivores, grazers",
           "Corals",
           "Carnivores, planktivores, omnivores",
           "Invertivores",
           "Herbivores"),
  values= c("aut" = "#003300",
            "nonMixAut" = "#00CC66",
            "corals" = "#B2FF66",
            "pred" = "#990000",
            "sessileInv" = "#CC0000",
            "herbivorous" = "#FF6666")
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

## binding already jiterred coordinates (fishes)
fish_fric <- cbind (fish_fric, 
                    richness = fish_richness$all,
                    LonJitter=fish_richness$LonJitter,
                    LatJitter = fish_richness$LatJitter)
## binding already jiterred coordinates (benthos)
benthos_fric <- cbind (benthos_fric, 
                       richness = benthos_richness$all,
                       LonJitter=benthos_richness$LonJitter,
                       LatJitter = benthos_richness$LatJitter)

## plot 

wm_pie <- wm + geom_scatterpie(aes(x=LonJitter, y=LatJitter, r=richness/max(richness)),
                               alpha=0.5,
                               col="#990000",
                               data = fish_fric,
                               cols = c(
                                 "FRic",
                                 "FEve",
                                 "FDiv"),
                               #pie_scale = 0.1,
                               size=0.5,
                               sorted_by_radius = F,
                               legend_name = "Index") 

### adding benthos richness

wm_pie_col_benthos <- wm_pie + geom_scatterpie(aes(x=LonJitter, y=LatJitter,r=richness/max(richness,na.rm=T)),
                                               alpha=0.5,
                                               col="#006633",
                                               data = benthos_fric,
                                               cols = c("FRic",
                                                        "FEve",
                                                        "FDiv"
                                               ),
                                               #pie_scale = 1,
                                               size=0.5,
                                               sorted_by_radius = F) + 
  
  theme (legend.title = element_text(size=7),
         legend.text = element_text(size=7),
         legend.position = c(0.35, 0.8),
         legend.justification = c("right", "top"),
         legend.box.just = "right",
         legend.margin = margin(6,6,6,6),
         legend.background = element_blank(),
         title=element_text(size=10)) 


wm_pie_col_benthos <- wm_pie_col_benthos + scale_fill_manual(
  
  values= c("FRic" = "#000000",
            "FEve" = "#A0A0A0",
            "FDiv" = "#FFFFFF"
  )
) + 
  ggtitle ("Functional diversity of fish and benthic communities")

wm_pie_col_benthos <- wm_pie_col_benthos + geom_scatterpie_legend(fish_fric$richness/max(fish_fric$richness), 
                                                                  x=-30, y=-28, n=3, 
                                                                  labeller=function(x) round(x*max(fish_fric$richness),2))
wm_pie_col_benthos <- wm_pie_col_benthos + geom_scatterpie_legend(benthos_fric$richness/max(benthos_fric$richness,na.rm=T), 
                                                                  x=-25, y=-28, n=3, 
                                                                  labeller=function(x) round(x*max(benthos_fric$richness,na.rm=T),2))

wm_pie_col_benthos 

### 

## selecao baseada no bentos (numbero de especies de bentos em cada grupo)
sel_sites <- c (32,28,12,3,16,6,21)#19,,
# subset coordenadas
subset_coord <- covariates_site$coord[sel_sites,]

## annotate in the plot
wm_pie_col_benthos <- wm_pie_col_benthos + 
  geom_segment(aes(x=subset_coord$Lon+7,y=subset_coord$Lat+0.1,
                xend=-20,yend=subset_coord$Lat+0.1),
                color = "gray", size=1,alpha=0.5,
               ) + 
 annotate("text", x = subset_coord$Lon[1]+11, y = subset_coord$Lat[1]-0.5, label = "Rocas Atoll",size=3) + 
  annotate("text", x = subset_coord$Lon[2]+12.8, y = subset_coord$Lat[2]-0.5, label = "Parrachos",size=3)+
  annotate("text", x = subset_coord$Lon[3]+12.2, y = subset_coord$Lat[3]-0.5, label = "Coral Coast",size=3)+
  annotate("text", x = subset_coord$Lon[4]+16.5, y = subset_coord$Lat[4]-0.5, label = "Abrolhos",size=3)+
  annotate("text", x = subset_coord$Lon[5]+17, y = subset_coord$Lat[5]-0.5, label = "Espírito Santo",size=3)+
  annotate("text", x = subset_coord$Lon[6]+18.4, y = subset_coord$Lat[6]-0.5, label = "Arraial do Cabo",size=3)+
  annotate("text", x = subset_coord$Lon[7]+24.5, y = subset_coord$Lat[7]+1.2, label = "Santa Catarina",size=3)

wm_pie_col_benthos

#################################################

### groups of fish species based on diet
fish_traits_subset <- lapply (list_fish_diet, function (i)
  fish_traits [which (fish_traits$Diet %in% i),])

# matching
subset_PCO <- lapply (fish_traits_subset, function (i) 
  (rownames(FD_results_f1$all$axesPCO) %in% i$Name))

## subsetting composition data
site_composition <- lapply (sel_sites , function (i)
  colnames (fish[i,][which(fish[i,]>0)])[-1]# menos 1 pq eh o nome do sitio
)
## quais sp estao nos sitios selecionados
site_subset_composition <- lapply (site_composition, function (i)
  which(rownames(FD_results_f1$all$axesPCO) %in% i)
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
        pred=subset_PCO$pred[site_subset_composition[[i]]])
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
  pred <-hull.data [which(hull.data$pred==T),]
  d <- pred [chull(pred, y = NULL),]
  
  ## hull for complete data
  
  vol_plot <- ggplot(axes_spp_composition[[i]], aes(A1, A2)) + 
    geom_point() + theme_bw()+
    geom_polygon(data=a, aes (A1,A2),alpha=0.1,fill="#808080") + 
    geom_polygon(data=b, aes (A1,A2,group=herb, fill=herb),fill="#FF6666",alpha=0.5)+
    geom_polygon(data=c, aes (A1,A2,group=inv, fill=inv),fill="#CC0000",alpha=0.5)+
    geom_polygon(data=d, aes (A1,A2,group=pred, fill=pred),fill="#990000",alpha=0.5) +
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
site_composition_benthos <- lapply (sel_sites , function (i)
  colnames (bent2[i,][which(bent2[i,]>0)])[-1]# menos 1 pq eh o nome do sitio
)
## quais sp estao nos sitios selecionados
site_subset_composition_benthos <- lapply (site_composition_benthos, function (i)
  which(rownames(FD_results_f1_bentos$all$axesPCO) %in% i)
)

# subset dos eixos da ordenacao
axes_spp_composition_benthos <- lapply (site_subset_composition_benthos, function (i)
  FD_results_f1_bentos$all$axesPCO[i,1:2]
)

## convex hull data
hull.data_benthos <- lapply (seq (1,length(axes_spp_composition_benthos)), function (i) 
  cbind(axes_spp_composition_benthos[[i]] , 
        nonaut = subset_PCO_benthos$nonAut[site_subset_composition_benthos[[i]]],
        nonmix = subset_PCO_benthos$nonMixAut[site_subset_composition_benthos[[i]]],
        mix = subset_PCO_benthos$corals[site_subset_composition_benthos[[i]]])
)

## get plots 
vol_plot_benthos <- lapply (seq (1,length(hull.data_benthos)), function (i) {
  
  hull.data <- fortify(hull.data_benthos[[i]])
  ## todas as spp
  a <- hull.data [chull(hull.data[,1:2], y = NULL),]
  ## herb
  nonaut <-hull.data [which(hull.data$nonaut==T),]
  b <- nonaut [chull(nonaut, y = NULL),]
  nonmix <-hull.data [which(hull.data$nonmix==T),]
  c <- nonmix [chull(nonmix, y = NULL),]
  mix <-hull.data [which(hull.data$mix==T),]
  d <- mix [chull(mix, y = NULL),]
  
  ## hull for complete data
  
  vol_plot <- ggplot(axes_spp_composition_benthos[[i]], aes(A1, A2)) + 
    geom_point() + theme_bw()+
    geom_polygon(data=a, aes (A1,A2),alpha=0.1, fill="#808080") + 
    geom_polygon(data=b, aes (A1,A2,group=nonaut, fill=nonaut),fill="#003300",alpha=0.5)+
  geom_polygon(data=c, aes (A1,A2,group=nonmix, fill=nonmix),fill="#00CC66",alpha=0.5)+
  geom_polygon(data=d, aes (A1,A2,group=mix, fill=mix),fill="#B2FF66",alpha=0.5) +
    
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
grid.arrange(wm_pie_col_benthos, 
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