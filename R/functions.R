## funcoes necessarias

## exploracao dos dados

## fazer um histograma pra saber o numero de eventIDS por ano
barplot_function <- function (df1,df2) { 
  
  barplot(
    table(as.numeric(substr (unique(df1$eventID), 
                             nchar (as.character(unique(df1$eventID)))-3,
                             nchar(as.character(unique(df1$eventID)))))),
    xlab="Year",col="green",main="Number of eventIDs per year",
    space=1.2,ylim=c(0,60
    ))
  
  barplot(
    table(as.numeric(substr (as.character(unique(df2$eventID)), 
                             nchar (as.character(unique(df2$eventID)))-3,
                             nchar(as.character(unique(df2$eventID)))))),
    col="green4",main="Number of eventIDs per year",add=T,
    space=c(2,1.2,1.2,1.2,1.2),xaxt="n")
  
  legend ("topleft",legend= c("Benthos", "Fishes"),
          pch=15,bty="n",
          col=c("green","green4"),
          cex=1.5)
  
}

## mapa exploratorio dos sitios

initial_map_function <- function (df1, df2) { 
  
  # mapa mundi
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  require(ggplot2)
  library(ggrepel)
  
  # cortar o mapa para ver a america do Sul e parte da central
  wm <- ggplot() + 
    geom_sf (data=world, size = 0.1, 
             fill= "gray90",colour="gray90") +
    coord_sf (xlim = c(-50, -20),  ylim = c(-30, 0), expand = FALSE) +
    theme_bw() + xlab(NULL) + ylab(NULL) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x=element_blank(),
          axis.text.y = element_blank(),axis.ticks.y=element_blank(),
          title = element_text(size=8)) 
  
  jitter <- position_jitter(width = 0.2, height = 0.5)
  
  bentos_coord <- wm + geom_point(data=df1,aes (x= Lon, y=Lat),
                                  stroke=1,shape=1, size=1, 
                                  position = jitter,col="red") 
  
  peixes_coord <- bentos_coord + geom_point(data=df2,aes (x= Lon, y=Lat),
                                            shape=19, size=0.1, 
                                            position = jitter)
  
  peixes_coord ## mostre o mapa
}


### 

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]}

# OUTPUTS:
# convex hull around data points in a particular color (specified by lcolor)

# FUNCTION:
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  
# END OF FUNCTION

# -------------------------
# random sampling of compositions
#----------------------

fc_random_composition <- function (data,nsamples,replace) {
  rdm_sample <- sample (seq(1,nrow(data)),nsamples, replace = replace)
  rdm_comp <- data[rdm_sample,]
  # aggregate
  rdm_comp <- aggregate (rdm_comp, by = list (rep(1,nrow(rdm_comp))), FUN=sum) # rm two first cols
  rdm_comp <- rdm_comp [,-c(1,2)]
  return(rdm_comp)
}


# --------------------------------------------------------------------- #
#    Functional diversity per site, data set, and sample size def
# --------------------------------------------------------------------- #

# function

function_FD <- function (site.data, spp.trait.data) {
      
      #---------------------------------------------------#
      # Adjust fish abundance data to relative abundance
      #-------------------------------------------------- #
      
      abund_matrix <- site.data[,which(colnames(site.data) %in% rownames(spp.trait.data))] 
      rel_abund <- abund_matrix/rowSums(abund_matrix) # relative abundance/cover
      quais_sitios_manter <-which(is.na (rowSums (rel_abund))!= T) 
      rel_abund <- rel_abund [quais_sitios_manter,]
      rel_abund <- rel_abund [,which(colSums (rel_abund) > 0)]
      
      # and match fish spp in both datasets
      fish_traits_ord_match <- spp.trait.data[which(rownames(spp.trait.data) %in% colnames (rel_abund)),]
      
      ## ---------------------------------- ##
      ## Building the gower distance matrix ##
      ## ---------------------------------- ##
      
      gower_matrix <- daisy (fish_traits_ord_match, metric=c("gower")) 
      
      ## ---------------------------------------------- ##
      ## Building the functional space based on a PCOA 
      ## quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
      ## ---------------------------------------------- ##
      
      pco<-dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10) 
      
      ### barplot of eigenvalues for each axis 
      barplot(pco$eig) 
      
      ### calculate and show percentage of inertia explained by the two first axes
      (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) 
      
      ## ----------------------------------------------------- ##
      # Testing the Quality of the functional space based on the method 
      # of Maire et al. (2015) in GEB, i.e. how many axes do we need to keep to 
      # faithfully represent the original Gower's distances.
      # The function will test the quality of the space from 2 to n axes using 
      # dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
      ## ------------------------------------------------------ ##
      
      quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
                                             plot="quality_funct_space_I") 
      
      ### the minimal value corresponds to the best space to use (min SD)
      axes_to_choose <- which(quality$meanSD == min(quality$meanSD)) 
      ### calculate and show the percentage of inertia explained by the choosen axes
      (Inertia7<- (sum(pco$eig[1:axes_to_choose])) /(sum(pco$eig))) 
      
      ### Method with dbFD() 
      fd <- dbFD (gower_matrix, rel_abund, 
                  corr ="none", 
                  w.abun = T,
                  m = axes_to_choose, 
                  stand.FRic = TRUE,
                  print.pco = TRUE,
                  calc.FGR = F,
                  clust.type = "kmeans")
      
      ### in the case of samples with too few species
      fd_indexes  <- data.frame (matrix (NA, ncol=8,nrow=nrow(rel_abund),
                                         dimnames = list (NULL,
                                                          c("nbsp",
                                                            "sing",
                                                            "FRic",
                                                            "qual.FRic",
                                                            "FEve",
                                                            "FDiv",
                                                            "FDis",
                                                            "RaoQ"))))
      
      fd_indexes [quais_sitios_manter,1] <- fd$nbsp
      fd_indexes[quais_sitios_manter,2] <- fd$sing.sp
      fd_indexes[quais_sitios_manter,3] <- fd$FRic
      fd_indexes[quais_sitios_manter,4] <- fd$qual.FRic
      fd_indexes[quais_sitios_manter,5] <- fd$FEve
      fd_indexes[quais_sitios_manter,6] <- fd$FDiv
      fd_indexes[quais_sitios_manter,7] <- fd$FDis
      fd_indexes[quais_sitios_manter,8] <- fd$RaoQ
      
      ### list with results
      results <- list (Fdindexes = fd_indexes,
                       chosenAxes = axes_to_choose,
                       InertiaPCO=Inertia2,
                       InertiaQuality=Inertia7,
                       explanAxes=pco$eig/sum(pco$eig),
                       convexHullVolume=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$vol,
                       convexHullArea=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$area,
                       axesPCO=fd$x.axes)

      return (results)
      
}


