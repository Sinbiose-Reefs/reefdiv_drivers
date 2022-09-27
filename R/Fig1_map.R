
## --------------------
#       FIGURES
# --------------------

## call packages
source("R/packages.R")
source("R/functions.R")

### -------------------------
#  FIG - supporting information
# correlations between variables
load(here("output", "modeling_data.RData"))



#correlation plot (Fig. S1)
# correlation between variables
corrplot (
  cor(site_covs[,c("sst", "turbidity","productivity", "salinity", "offshore_distance")]),
  p.mat = cor.mtest(cor(site_covs[,c("sst", "turbidity","productivity", "salinity", "offshore_distance")]))$p
  

  )



# -----------------------------------------------------------

# load results of functional analyses  
load(here("output", "FD_results.RData"))




# ------------------------------------------- #
# map of distribution of samples and sampling effort
# fish

# create a DF with all data

bind_fish_benthos<- cbind (site_covs,
                           effort_dataframe[,-which(colnames(effort_dataframe) == "sites")],
                           # benthos
                           SR_corals =  (df_corals$SR),
                           FRic_corals = df_corals$FRic,
                           Rao_corals = df_corals$RaoQ,
                           # fish
                           SR_fish =  (FD_fish$nbsp),
                           FRic_fish = FD_fish$FRic,
                           Rao_fish = FD_fish$RaoQ,
                           # algae
                           SR_algae =  (FD_algae$nbsp),
                           FRic_algae = FD_algae$FRic,
                           Rao_algae = FD_algae$RaoQ
                           
)

# melt this data to fit to one plot
bind_fish_benthos_SR_plot <- bind_fish_benthos[,c("decimalLatitude",
                                                  "SR_corals",
                                                  "SR_algae",
                                                  "SR_fish")]

# melt
bind_fish_benthos_SR_plot <- melt (bind_fish_benthos_SR_plot,id.vars = "decimalLatitude")

## plot fish variation in obs richness and effort of fish DS

# scaling factor for UVC data
# scaleFactor <- max(bind_fish_benthos$SR_fish) / max(bind_fish_benthos$fish_effort)


# raw correlation between taxa

cor (bind_fish_benthos[,c("SR_corals",
                          "SR_algae",
                          "SR_fish")])




# plot 
# fish

eff_obsrich_plot <- ggplot(bind_fish_benthos_SR_plot, aes(x=decimalLatitude, 
                                                          group = variable,
                                                          fill=variable,
                                                          colour=variable)) +
  geom_point(aes(y=log(value+1)),
             alpha=0.2) +
  geom_smooth(aes(y=log(value+1)),fill="gray",alpha=0.1, 
              method = "lm", 
              formula = y ~ poly(x, 1), size = 1) +  
  scale_fill_manual(values = c ("SR_corals" = "#FEB139",
                                "SR_algae" = "#42855B",
                                "SR_fish" = "#0096FF"))+
  scale_colour_manual(values = c ("SR_corals" = "#FEB139",
                                "SR_algae" = "#42855B",
                                "SR_fish" = "#0096FF"))+
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


# environmental covariates in space
## sea surface temperature 
SST <- ggplot(site_covs, aes(x=decimalLatitude,y=sst)) +
  
  geom_point(alpha=0.2,col="#8E0505") +
  
  geom_smooth(method = "glm", 
              formula = y ~ poly(x, 2), 
              size = 1,col ="#8E0505",
              fill="#8E0505",alpha=0.1) +
  
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

# turbidity
kd490 <- ggplot(site_covs, aes(x=decimalLatitude,y=turbidity)) +
  
  geom_point(alpha=0.2,col="#8E0505") +
  
  geom_smooth(method = "glm", 
              formula = y ~ poly(x, 2), 
              size = 1,col ="#8E0505",
              fill="#8E0505",alpha=0.1) +
  
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

kd490


# salinity
salinity <- ggplot(site_covs, aes(x=decimalLatitude,y=salinity)) +
  
  geom_point(alpha=0.2,col="#8E0505") +
  
  geom_smooth(method = "glm", 
              formula = y ~ poly(x, 2), 
              size = 1,col ="#8E0505",
              fill="#8E0505",alpha=0.1) +
  
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

salinity

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

# map
# site covs ordered by latitude
site_covs_lat <- site_covs[order(site_covs$decimalLatitude,decreasing=T),]
# define site number to plot in the map
site_covs_lat$site_number <- seq (1,nrow(site_covs_lat))
site_covs$site_number <- site_covs_lat$site_number[match (site_covs$sites, site_covs_lat$sites)]

# plot
map_peixes <- wm + geom_point(data=site_covs, 
                              
                              aes(x=decimalLongitude, 
                                  y=decimalLatitude,
                                  shape = region),
                              size=3,
                              col="#0e49b5",
                              alpha = 0.4) + 
  geom_text_repel (data = site_covs, aes(x=decimalLongitude, 
                            y=decimalLatitude,
                            label = site_number),
                   size = 4,
                   min.segment.length = 0,
                   box.padding = 0.3,
                   max.overlaps=100)


map_peixes




# functional trait spaces
# build it for each dataset


build_FS <- function (FD_output, 
                      composition, 
                      site_covs,
                      sel_coordinates = c(0,-5,-10,-15,-20,-25),
                      point_color,
                      space_color = "black",
                      complete_space_color = "gray30") {
        
        # find axes in the output
        axes_fish <- FD_output$x.axes
        
        # community
        comm_fish <- composition
        
        # closest to four coords
        seq_search <- sel_coordinates
        closest_fish <- lapply (as.list(seq_search), function (i) 
          closest (site_covs$decimalLatitude,i))
        
        # find the selected comms
        sel_comms_fish <- lapply (closest_fish, function (i) 
          which(site_covs$decimalLatitude == i[1])[1])
        # which ones
        # sites [unlist(sel_comms_fish)]
        # range rel abundance to plot in the selected comm
        range_plot_fish <- range(comm_fish [unlist(sel_comms_fish),]) 
        
        
        # find spp in each community and built the plot 
        fish_space <- lapply (sel_comms_fish, function (i) {
          
          # community subset
          sel_comm_data <- comm_fish[i,]
          pres_fish <- sel_comm_data[which(sel_comm_data>0)]
          
          # complete trait space
          all <- cbind (axes_fish[,1:2],ext = F)
          a <- all [chull(all[,1:2], y = NULL),]
          
          # space occupied by the community
          setB<-cbind(all, ext1=ifelse(rownames(all) %in% names(pres_fish),
                                       F,
                                       T))
          pk <-setB[which(setB$ext1==F),]
          f <- pk [chull(pk, y = NULL),]
          # abundance on pk
          pk$abund <- as.numeric(pres_fish[match(rownames(pk), names(pres_fish))])
          
          # plot space
          plotA <- ggplot(a, aes(A1, A2)) + 
            geom_point(colour=complete_space_color) + theme_bw()+
            geom_polygon(data=a, aes (A1,A2),alpha=1,
                         fill=complete_space_color) + 
            geom_polygon(data=f, aes (A1,A2,group=ext1, fill=ext1),alpha=1,
                         fill=space_color,size=3) +
            xlim(min (a$A1)-0.2,max (a$A1)+0.2) + 
            annotate("text",
                     x=0,
                     y=0,
                     size=5,
                     label=paste ("Site=", site_covs$site_number[i], 
                                  "\nFRic=", round(FD_output$FRic[i],6),
                                  "\nRao's Q=", round(FD_output$RaoQ[i],4))
            ) +
            geom_point(data=pk,aes (A1,A2,size=(abund)),
                       alpha=0.8,col=point_color,stroke=1) + 
            scale_size(name="Relative\nabundance",
                       limits=c(0,1),
                       breaks=seq(0,1,0.2))+ 
            theme(axis.text = element_text(size=6),
                  axis.title=element_text(size=8))
          ; # return
          plotA
          
        })

}

# run the function
# fish trait space
fish_space <- build_FS(FD_output = FD_fish,
                       composition = comp_fish,
                       site_covs = site_covs,
                       sel_coordinates = c(0,-4,-10,-17,-22,-24),
                       point_color = "#FFF5E4",
                       space_color = "#97D2EC",
                       complete_space_color = "#003865")

# coral trait space
coral_space <- build_FS(FD_output = FD_corals,
                       composition = comp_corals,
                       site_covs = site_covs,
                       sel_coordinates = c(0,-4,-10,-17,-22,-24),
                       point_color = "#FFF5E4",
                       space_color = "#FF9551",
                       complete_space_color = "#FF5B00")

# algae trait space
algae_space <- build_FS(FD_output = FD_algae,
                       composition = comp_algae,
                       site_covs = site_covs,
                       sel_coordinates = c(0,-4,-10,-17,-22,-24),
                       point_color = "#FFF5E4",
                       space_color = "#90B77D",
                       complete_space_color = "#1C6758")



# arrange
# tropical comm
array_fish1 <- grid.arrange(fish_space[[1]]+theme(legend.position="none",
                                                  legend.text = element_text(size=2),
                                                  legend.title = element_text(size=2),
                                                 legend.direction = "horizontal",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()),
                           fish_space[[2]]+theme(legend.position="none",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()),
                           fish_space[[3]]+theme(legend.position="none",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()),
                           nrow=1)

# subtropics
array_fish2 <- grid.arrange(fish_space[[4]]+theme(legend.position="none",
                                                  legend.direction = "horizontal",
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank()),
                            fish_space[[5]]+theme(legend.position="none",
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank()),
                            fish_space[[6]]+theme(legend.position="none",
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank()),
                            nrow=1)

# corals
# tropics
array_corals1 <- grid.arrange(coral_space [[1]]+theme(legend.position="none",
                                                 legend.direction = "horizontal",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()),
                             coral_space[[2]]+theme(legend.position="none",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()),
                             coral_space[[3]]+theme(legend.position="none",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()),
                             nrow=1)


# subtropics
array_corals2 <- grid.arrange(coral_space [[4]]+theme(legend.position="none",
                                                      legend.direction = "horizontal",
                                                      axis.title.x = element_blank(),
                                                      axis.text.x = element_blank()),
                              coral_space[[5]]+theme(legend.position="none",
                                                     axis.title.x = element_blank(),
                                                     axis.text.x = element_blank()),
                              coral_space[[6]]+theme(legend.position="none",
                                                     axis.title.x = element_blank(),
                                                     axis.text.x = element_blank()),
                              nrow=1)

# algae
# tropics
array_algae1 <- grid.arrange(algae_space [[1]]+theme(legend.position="none",
                                                     legend.direction = "horizontal",
                                                     axis.title.x = element_blank(),
                                                     axis.text.x = element_blank()),
                            algae_space[[2]]+theme(legend.position="none",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                            algae_space[[3]]+theme(legend.position="none",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                             nrow=1)

# subtropics
array_algae2 <- grid.arrange(algae_space [[4]]+theme(legend.position="none",
                                                     legend.direction = "horizontal",
                                                     axis.title.x = element_blank(),
                                                     axis.text.x = element_blank()),
                             algae_space[[5]]+theme(legend.position="none",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                             algae_space[[6]]+theme(legend.position="none",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                             nrow=1)


# map and covariates

# save
pdf(file=here("output",
              "figures", 
              "Fig1_map.pdf"),height=7,width=9)


map_arrange<- grid.arrange(map_peixes+theme(legend.position = "none"), 
                           eff_obsrich_plot+theme(legend.position = "none"),
                           SST,
                           kd490,
                           salinity,
                           ncol=14,nrow=11,
                           layout_matrix = rbind (c(1,1,1,1,1,1,1,1,2,2,2,3,3,3),
                                                  c(1,1,1,1,1,1,1,1,2,2,2,3,3,3),
                                                  c(1,1,1,1,1,1,1,1,2,2,2,3,3,3),
                                                  c(1,1,1,1,1,1,1,1,2,2,2,3,3,3),
                                                  c(1,1,1,1,1,1,1,1,2,2,2,3,3,3),
                                                  c(1,1,1,1,1,1,1,1,4,4,4,5,5,5),
                                                  c(1,1,1,1,1,1,1,1,4,4,4,5,5,5),
                                                  c(1,1,1,1,1,1,1,1,4,4,4,5,5,5),
                                                  c(1,1,1,1,1,1,1,1,4,4,4,5,5,5),
                                                  c(1,1,1,1,1,1,1,1,4,4,4,5,5,5),
                                                  c(1,1,1,1,1,1,1,1,4,4,4,5,5,5)))

dev.off()

# arrange all
ncol_plot <- 4
nrow_plot <- 6

# save
pdf(file=here("output",
              "figures", 
              "spaces_with_pts.pdf"),height=6,width=14)

grid.arrange(array_fish1,array_fish2,
             array_corals1,array_corals2,
             array_algae1,array_algae2,
             ncol = ncol_plot,
             nrow=nrow_plot,
             layout_matrix = rbind (c(rep(1,ncol_plot),
                                    rep(2,ncol_plot)),
                                    c(rep(1,ncol_plot),
                                    rep(2,ncol_plot)),
                                    c(rep(3,ncol_plot),
                                    rep(4,ncol_plot)),
                                    c(rep(3,ncol_plot),
                                    rep(4,ncol_plot)),
                                    c (rep(5,ncol_plot),
                                    rep(6,ncol_plot)),
                                    c (rep(5,ncol_plot),
                                    rep(6,ncol_plot))
             ))



dev.off()

# end