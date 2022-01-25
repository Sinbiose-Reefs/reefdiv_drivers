# closest function

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
  rdm_sample <- sample (seq(1,nrow(data)),nsamples, replace = replace) # list all samples deployed into a site
  rdm_comp <- data[rdm_sample,] # get the sample
  
  # aggregate
  rdm_comp <- aggregate (rdm_comp, by = list (rep(1,nrow(rdm_comp))), FUN=sum) #aggregate to obtain site composition
  rdm_comp <- rdm_comp  # rm the first column (site name 'group.1')
  return(rdm_comp) # return
}


# --------------------------------------------------------------------- #
#    Functional diversity per site, data set, and sample size def
# --------------------------------------------------------------------- #

# function for fishes (it calculates relative abundance)

function_FD_fish <- function (site.data, spp.trait.data) {
      
      #---------------------------------------------------#
      # Adjust fish abundance data to relative abundance
      #-------------------------------------------------- #
      
      abund_matrix <- site.data[,which(colnames(site.data) %in% rownames(spp.trait.data))] 
      rel_abund <- abund_matrix/rowSums(abund_matrix) # relative abundance/cover
      quais_sitios_manter <-which(is.na (rowSums (rel_abund))!= T) 
      rel_abund <- rel_abund [quais_sitios_manter,]
      rel_abund <- rel_abund [,which(colSums (rel_abund) > 0)]
      rel_abund [rel_abund>0]<-1
      
      # and match fish spp in both datasets
      fish_traits_ord_match <- spp.trait.data[which(rownames(spp.trait.data) %in% colnames (rel_abund)),]
      
      ## CWM
      cwm <-  tryCatch (functcomp(fish_traits_ord_match, as.matrix(rel_abund),
                                  CWM.type = "all", bin.num = NULL),error = function (e)
                                    return (e))
      
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
      #barplot(pco$eig) 
      
      ### calculate and show percentage of inertia explained by the two first axes
      (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) 
      
      ## ----------------------------------------------------- ##
      # Testing the Quality of the functional space based on the method 
      # of Maire et al. (2015) in GEB, i.e. how many axes do we need to keep to 
      # faithfully represent the original Gower's distances.
      # The function will test the quality of the space from 2 to n axes using 
      # dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
      ## ------------------------------------------------------ ##
      
      dev.set(dev.next())
      
      quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
                                             plot="quality_funct_space_I_fish") 
      
      ### the minimal value corresponds to the best space to use (min SD)
      axes_to_choose <- which(quality$meanSD == min(quality$meanSD)) 
      ### calculate and show the percentage of inertia explained by the choosen axes
      (Inertia7<- (sum(pco$eig[1:axes_to_choose])) /(sum(pco$eig))) 
      
      ### Method with dbFD() 
      fd <- dbFD (gower_matrix, rel_abund, 
                  corr ="none", 
                  w.abun = F,
                  m = axes_to_choose, 
                  stand.FRic = TRUE,
                  print.pco = TRUE,
                  calc.FGR = F,
                  clust.type = "kmeans",
                  calc.CWM = TRUE)
      
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
                       axesPCO=fd$x.axes,
                       cwm = cwm)

      return (results)
      
}


# function for benthos (it uses cover)

function_FD_benthos <- function (site.data, spp.trait.data) {
  
  #---------------------------------------------------#
  # Adjust benthic incidence
  #-------------------------------------------------- #
  
  abund_matrix <- site.data[,which(colnames(site.data) %in% rownames(spp.trait.data))] 
  rel_abund <- abund_matrix [,which(colSums (abund_matrix) > 0)]
  rel_abund [rel_abund>0]<-1
  rel_abund<-data.matrix(rel_abund)
  
  # and match fish spp in both datasets
  benthos_traits_ord_match <- spp.trait.data[which(rownames(spp.trait.data) %in% colnames (rel_abund)),]
  
  ## CWM
  cwm <-  tryCatch (functcomp(benthos_traits_ord_match,(rel_abund), # avoid and report errors
                    CWM.type = "all", bin.num = NULL),error = function (e)
                      return (e))
  
  ## ---------------------------------- ##
  ## Building the gower distance matrix ##
  ## ---------------------------------- ##
  
  gower_matrix <- daisy (benthos_traits_ord_match, metric=c("gower")) 
  
  ## ---------------------------------------------- ##
  ## Building the functional space based on a PCOA 
  ## quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
  ## ---------------------------------------------- ##
  
  pco<-dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10) 
  
  ### barplot of eigenvalues for each axis 
  #barplot(pco$eig) 
  
  ### calculate and show percentage of inertia explained by the two first axes
  (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) 
  
  ## ----------------------------------------------------- ##
  # Testing the Quality of the functional space based on the method 
  # of Maire et al. (2015) in GEB, i.e. how many axes do we need to keep to 
  # faithfully represent the original Gower's distances.
  # The function will test the quality of the space from 2 to n axes using 
  # dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
  ## ------------------------------------------------------ ##
  
  dev.set(dev.next())
  
  quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
                                         plot="quality_funct_space_I_benthos") 
  
  ### the minimal value corresponds to the best space to use (min SD)
  axes_to_choose <- which(quality$meanSD == min(quality$meanSD)) 
  ### calculate and show the percentage of inertia explained by the choosen axes
  (Inertia7<- (sum(pco$eig[1:axes_to_choose])) /(sum(pco$eig))) 
  
  ### Method with dbFD() 
  fd <- dbFD (gower_matrix, (rel_abund), 
              corr ="none", 
              w.abun = F,
              m = axes_to_choose, 
              stand.FRic = TRUE,
              print.pco = TRUE,
              calc.FGR = F,
              clust.type = "kmeans",
              calc.CWM = TRUE)
  
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
  
  fd_indexes [,1] <- fd$nbsp
  fd_indexes[,2] <- fd$sing.sp
  fd_indexes[,3] <- fd$FRic
  fd_indexes[,4] <- fd$qual.FRic
  fd_indexes[,5] <- fd$FEve
  fd_indexes[,6] <- fd$FDiv
  fd_indexes[,7] <- fd$FDis
  fd_indexes[,8] <- fd$RaoQ
  
  ### list with results
  results <- list (Fdindexes = fd_indexes,
                   chosenAxes = axes_to_choose,
                   InertiaPCO=Inertia2,
                   InertiaQuality=Inertia7,
                   explanAxes=pco$eig/sum(pco$eig),
                   convexHullVolume=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$vol,
                   convexHullArea=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$area,
                   axesPCO=fd$x.axes,
                   cwm = cwm)
  
  return (results)
  
}



# --------------------------------------------------------------------- #
#    Functional diversity per site, data set, and sample size def
#                       weighted by abundance
# --------------------------------------------------------------------- #

# function for fishes (it calculates relative abundance)

function_FD_fish_abundW <- function (site.data, spp.trait.data) {
  
  #---------------------------------------------------#
  # Adjust fish abundance data to relative abundance
  #-------------------------------------------------- #
  
  abund_matrix <- site.data[,which(colnames(site.data) %in% rownames(spp.trait.data))] 
  rel_abund <- abund_matrix#/rowSums(abund_matrix) # relative abundance/cover
  quais_sitios_manter <-which(is.na (rowSums (rel_abund))!= T) 
  rel_abund <- rel_abund [quais_sitios_manter,]
  rel_abund <- rel_abund [,which(colSums (rel_abund) > 0)]
  
  # and match fish spp in both datasets
  fish_traits_ord_match <- spp.trait.data[which(rownames(spp.trait.data) %in% colnames (rel_abund)),]
  
  ## CWM
  cwm <-  tryCatch (functcomp(fish_traits_ord_match, as.matrix(rel_abund),
                              CWM.type = "all", bin.num = NULL),error = function (e)
                                return (e))
  
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
  #barplot(pco$eig) 
  
  ### calculate and show percentage of inertia explained by the two first axes
  (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) 
  
  ## ----------------------------------------------------- ##
  # Testing the Quality of the functional space based on the method 
  # of Maire et al. (2015) in GEB, i.e. how many axes do we need to keep to 
  # faithfully represent the original Gower's distances.
  # The function will test the quality of the space from 2 to n axes using 
  # dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
  ## ------------------------------------------------------ ##
  
  #dev.set(dev.next())
  
  quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
                                         plot="quality_funct_space_I_fish_abund") 
  
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
              clust.type = "kmeans",
              calc.CWM = TRUE)
  
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
                   axesPCO=fd$x.axes,
                   cwm = cwm)
  
  return (results)
  
}


# function for benthos (it uses cover)

function_FD_benthos_abundW <- function (site.data, spp.trait.data) {
  
  #---------------------------------------------------#
  # Adjust benthic cover data 
  #-------------------------------------------------- #
  
  abund_matrix <- site.data[,which(colnames(site.data) %in% rownames(spp.trait.data))] 
  rel_abund <- abund_matrix [,which(colSums (abund_matrix) > 0)]
  rel_abund<-data.matrix(rel_abund)
  
  # and match fish spp in both datasets
  benthos_traits_ord_match <- spp.trait.data[which(rownames(spp.trait.data) %in% colnames (rel_abund)),]
  
  ## CWM
  cwm <-  tryCatch (functcomp(benthos_traits_ord_match,(rel_abund), # avoid and report errors
                              CWM.type = "all", bin.num = NULL),error = function (e)
                                return (e))
  
  ## ---------------------------------- ##
  ## Building the gower distance matrix ##
  ## ---------------------------------- ##
  
  gower_matrix <- daisy (benthos_traits_ord_match, metric=c("gower")) 
  
  ## ---------------------------------------------- ##
  ## Building the functional space based on a PCOA 
  ## quasieuclid() transformation to make the gower matrix as euclidean. nf= number of axis 
  ## ---------------------------------------------- ##
  
  pco<-dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10) 
  
  ### barplot of eigenvalues for each axis 
  #barplot(pco$eig) 
  
  ### calculate and show percentage of inertia explained by the two first axes
  (Inertia2<-(pco$eig[1]+pco$eig[2]) /(sum(pco$eig))) 
  
  ## ----------------------------------------------------- ##
  # Testing the Quality of the functional space based on the method 
  # of Maire et al. (2015) in GEB, i.e. how many axes do we need to keep to 
  # faithfully represent the original Gower's distances.
  # The function will test the quality of the space from 2 to n axes using 
  # dudi.pco(quasieuclid(gower_matrix), scannf=F, nf=10)
  ## ------------------------------------------------------ ##
  
  #dev.set(dev.next())
  
  quality<-quality_funct_space_fromdist( gower_matrix,  nbdim=10,   
                                         plot="quality_funct_space_I_benthos_abund") 
  
  ### the minimal value corresponds to the best space to use (min SD)
  axes_to_choose <- which(quality$meanSD == min(quality$meanSD)) 
  ### calculate and show the percentage of inertia explained by the choosen axes
  (Inertia7<- (sum(pco$eig[1:axes_to_choose])) /(sum(pco$eig))) 
  
  ### Method with dbFD() 
  fd <- dbFD (gower_matrix, (rel_abund), 
              corr ="none", 
              w.abun = T,
              m = axes_to_choose, 
              stand.FRic = TRUE,
              print.pco = TRUE,
              calc.FGR = F,
              clust.type = "kmeans",
              calc.CWM = TRUE)
  
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
  
  fd_indexes [,1] <- fd$nbsp
  fd_indexes[,2] <- fd$sing.sp
  fd_indexes[,3] <- fd$FRic
  fd_indexes[,4] <- fd$qual.FRic
  fd_indexes[,5] <- fd$FEve
  fd_indexes[,6] <- fd$FDiv
  fd_indexes[,7] <- fd$FDis
  fd_indexes[,8] <- fd$RaoQ
  
  ### list with results
  results <- list (Fdindexes = fd_indexes,
                   chosenAxes = axes_to_choose,
                   InertiaPCO=Inertia2,
                   InertiaQuality=Inertia7,
                   explanAxes=pco$eig/sum(pco$eig),
                   convexHullVolume=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$vol,
                   convexHullArea=convhulln(fd$x.axes[,1:axes_to_choose],output.options = T)$area,
                   axesPCO=fd$x.axes,
                   cwm = cwm)
  
  return (results)
  
}


# function of binomial smooth
# function from here
# https://ggplot2.tidyverse.org/reference/geom_smooth.html

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

