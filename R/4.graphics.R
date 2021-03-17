
## --------------------
# FIGURES

## chamar os pacotes
source("R/packages.R")
source("R/functions.R")

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
# Fig 1
# map of distribution of samples, and sampling effort

ObsRich <- rowSums (dados_peixes_bentos$peixes[,-1] >0)

eff_data <- data.frame (Lat=covariates_site$coord$coord_peixes[,"Lat"],
                        Lon=covariates_site$coord$coord_peixes[,"Lon"],
                        UVC=covariates_effort$effort$n_transectos_peixes,
                        ObsRich = ObsRich)
## 
ObsRichBentos <- rowSums (dados_peixes_bentos$bentos[,-1] >0)

eff_data_bentos <- data.frame (Lat=covariates_site$coord$coord_bentos[,"Lat"],
                               Lon=covariates_site$coord$coord_bentos[,"Lon"],
                               UVC=covariates_effort$effort$n_videos_bentos[which(covariates_effort$effort$n_videos_bentos >0)],
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
         plot.margin = unit(c(0,-0.8,0,0.3), "cm")
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

## maps
# mapa mundi
world <- ne_countries(scale = "medium", returnclass = "sf")

# cortar o mapa para ver a america do Sul e parte da central
wm <- ggplot() + 
  geom_sf (data=world, size = 0.1, 
           fill= "#aaaaaa",colour="#aaaaaa") +
  coord_sf (xlim = c(-50, -25),  ylim = c(-30, 4), expand = FALSE) +
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

pdf(file=here("output","vectorized","Fig1.pdf"),height=5,width=9)

grid.arrange(map_peixes_bentos, 
             eff_obsrich_plot,
             eff_obsrich_plot_benthos,
             ncol=8,nrow=11,
             layout_matrix = rbind (c(1,1,1,1,2,2,3,3),
                                    c(1,1,1,1,2,2,3,3),
                                    c(1,1,1,1,2,2,3,3),
                                    c(1,1,1,1,2,2,3,3),
                                    c(1,1,1,1,2,2,3,3),
                                    c(1,1,1,1,2,2,3,3),
                                    c(1,1,1,1,2,2,3,3),
                                    c(1,1,1,1,2,2,3,3),
                                    c(1,1,1,1,2,2,3,3)))

dev.off()

#### -------------------------
# Fig 2
# accumulation of FD relative to SR

load (here ("output", "complete_results_for_fig2.RData"))

# accumulation of diversity

## plotting
## Minimum Sample Size (MSS)

p1 <- ggplot (complete_results_for_fig2 [[1]], aes (x=EST.rich,y=FD, fill=Organism, 
                                              colour=Organism)) + 
  geom_point() +
  scale_color_manual(values=c("#f2a154", "#0e49b5")) +
  geom_smooth(method="lm", formula = y ~x) + 
  theme_classic() + xlab ("") + ylab ("Functional diversity")+
  theme (legend.position = "none") + 
  ggtitle ("Minimum Sample Size (MSS)")

## Precision-based Sample Size (PSSi)

p2 <- ggplot (complete_results_for_fig2 [[2]], aes (x=EST.rich,y=FD, fill=Organism, 
                                              colour=Organism)) + 
  geom_point() +
  scale_color_manual(values=c("#f2a154", "#0e49b5")) +
  geom_smooth(method="lm", formula = y ~x) + 
  theme_classic() + xlab ("Estimated Richness") + ylab ("")+
  theme (legend.position = "none") + 
  ggtitle ("Precision-based Sample Size (PSSi)")

## Model-based Sample Size (Lomolino's model)

p3 <- ggplot (complete_results_for_fig2 [[3]], aes (x=EST.rich,y=FD, fill=Organism, 
                                              colour=Organism)) + 
  geom_point() +
  scale_color_manual(values=c("#f2a154", "#0e49b5")) +
  geom_smooth(method="lm", formula = y ~x) + 
  theme_classic() + xlab ("") + ylab ("")+
  theme (legend.position = "none") + 
  ggtitle ("Lomolino's model")

pdf(file=here("output","vectorized","Fig2.pdf"),height=4,width=11)

grid.arrange(p1,p2,p3,
             ncol=9,nrow=6,
             layout_matrix = rbind (c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3),
                                    c(1,1,1,2,2,2,3,3,3)))
dev.off()

### -------------------------
#  FIG 3
# coefficients of drivers

load (here ("output", "complete_results_for_fig3.RData"))

# plot
dodge <- c(0.4,0.4)
pd <- position_dodge(dodge)
pdf_pt <- position_dodge(dodge)

a <- ggplot (complete_results[which(complete_results$Parameter != "Intercept"),], 
             
             aes  (y=Parameter, x=Estimate, 
                          colour=Algorithm, fill=Algorithm)) + 
  
  geom_errorbar(aes(xmin=lower,xmax=upper),width = 0.2,size=0.5,
                position=pd) + theme_classic() + 
  
  geom_point(position=(pdf_pt), 
             size=1.5)+ 
  
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray50", size=0.1)+
  
  facet_wrap(~Organism+Index,scale="free",ncol=4) + 

  scale_color_manual(values=c("gray90", "gray60", "gray40")) + 
  
  xlab("Standardized effect size") + 
  
  ylab ("Parameter") + 
  
  #xlim(-0.5,0.5) +
  
  theme(axis.text.x = element_text(angle = 45,size=7))

a

ggsave (file=here("output","vectorized","Fig3.pdf"),width=11,height=6)


