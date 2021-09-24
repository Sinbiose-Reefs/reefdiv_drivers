
# ----------------------------------------------------------------------------#
#    routine to modeling fish and benthos SR and FD  relative to environment
#                          using GLM
# 	PS: RUN IT IN R GUI
# ----------------------------------------------------------------------------#

## chamar os pacotes
source("R/packages.R")
source("R/functions.R")

# ------------------------------------------ #
# Load data to modeling

load (here ("output","data_to_modeling_GLM_abundW.RData"))


############################################################################
# -------------------------------------------------------------------------
#
#          Multivariate Linear models in BRMS
#
#        Effect of Environmental drivers
#
#
# -------------------------------------------------------------------------
############################################################################


# -------------------------------------
# prepare data to multivariate models
# -------------------------------------

# create a DF with all data

bind_fish_benthos<- cbind (cov_fish, 
                            EstRich_benthos = cov_benthos$EstRich,
                            FRic_benthos = cov_benthos$FRic,
                            FEve_benthos = cov_benthos$FEve,
                            FDiv_benthos = cov_benthos$FDiv)

# average values of metrics to present in the Results
round(apply(bind_fish_benthos[,c(1:4,22:25)],2,mean,na.rm=T),3)
round(apply(bind_fish_benthos[,c(1:4,22:25)],2,sd,na.rm=T),3)

## ----------------------------------------
# no spatial autocorrelation
## ----------------------------------------

# set formula (the same for benthos and fishes)
# complete model
formula <- brms::bf(mvbind (log(EstRich),
                            log(FRic),
                            log(FEve),
                            log(FDiv),
                            log(EstRich_benthos),
                            log(FRic_benthos),
                            log(FEve_benthos),
                            log(FDiv_benthos))~ BO2_tempmean_ss_std +
                                                distance_std +
                                                BO_damean_std + 
                                                Depth,
                    
                    nl = F)+
  set_rescor(TRUE) # nonlinear


# alternative 1
formulaA1 <- brms::bf(mvbind (log(EstRich),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(EstRich_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std +
                                                    distance_std +
                                                    BO_damean_std,
                      nl = F) +
  set_rescor(TRUE)# nonlinear

# alternative 2
formulaA2 <- brms::bf(mvbind (log(EstRich),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(EstRich_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std + 
                                                    distance_std,
                      nl = F) +
  set_rescor(TRUE)# nonlinear

# alternative 3
formulaA3 <- brms::bf(mvbind (log(EstRich),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(EstRich_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std,
                      nl = F) +
  set_rescor(TRUE)# nonlinear

#setting priors 
priors <-  c(set_prior("normal(0,5)",class = "b",coef = "",
              resp=c("logEstRich","logFDiv","logFEve","logFRic",
                     "logEstRichbenthos","logFDivbenthos","logFEvebenthos","logFRicbenthos")),
    set_prior("normal(0,5)",class = "Intercept",coef = "",
              resp=c("logEstRich","logFDiv","logFEve","logFRic",
                     "logEstRichbenthos","logFDivbenthos","logFEvebenthos","logFRicbenthos"))	
  )


# MCMC settings
ni <- 10000 
nb <- 8000
nt <- 4
nc <- 3

# run MCMC chains across different organisms and models
MCMC_runs <- lapply (list(formula,formulaA1,formulaA2, formulaA3), function(k) #ACROSS MODELS
    
    
    brms::brm(k, # for each model / formula
              
              data = bind_fish_benthos, 
              
              family = gaussian(),
              
              prior = priors, 
              
              chains = nc,
              
              cores = nc,

              iter = ni,
              
              warmup = nb,
              
              thin=nt,

              save_pars = save_pars(all = TRUE),
              
              control = list(adapt_delta = 0.99) # to avoid disagreement in chain mixing
              
    )
  )



# LOO model fit checking   
# run loo fit test
loo_test <- lapply (MCMC_runs,loo, moment_match=T,reloo=T)
# extract estimates
loo_sel <- lapply (loo_test, function (i)
         i$estimates[which(rownames(i$estimates) == "looic"),"Estimate"])
# extract looic (like AIC)
tab_mod_sel <- do.call(rbind,lapply (loo_test, function (i)
         i$estimates[which(rownames(i$estimates) == "looic"),]))
# extract estimated number of parameters (model adequacy)
tab_mod_fit <- do.call(rbind,lapply (loo_test, function (i)
         i$estimates[which(rownames(i$estimates) == "p_loo"),]))
# select the model with lowest looic
sel_model <- MCMC_runs[which(loo_sel == min(unlist(loo_sel)))]
   
# list of results 
res <- list (looic = tab_mod_sel,
             param= tab_mod_fit,
             best_model = sel_model
                   )

## save
save ( MCMC_runs,res, 
       file=here ("output", 
                  "MCMC_runs_multivariate_rarefied_no_aut_abundW.Rdata"))


############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FD RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################

# run model (ancova)
model.ancova <- brm (FRic ~ EstRich*Organism,
                data=data_to_modeling_GLM,
                family = gaussian (link="identity"),
                chains=nc,
                iter = ni,
                warmup = nb,
                thin=nt)


# summary of results
summary (model.ancova)

# plotting
p1<-plot(conditional_effects(model.ancova,
                         method="fitted",
                         re_formula=NA,
                         robust=T,
                         effects = "EstRich:Organism",
                         points=T,
                         prob = 0.95),
     
     theme = theme_classic() +
      
        theme (axis.title = element_text(size=15),
               axis.text = element_text(size=12),
               legend.position = "top") ,
     points=T) [[1]] + 
  
  scale_color_manual(values=c("#f2a154", "#0e49b5")) + 
  
  xlab("Species richness") + 
  
  ylab ("FRic")
    
  
# compare slopes
m.lst <- emtrends (model.ancova, "Organism", var="EstRich")
m.lst_tab <- summary(m.lst,point.est = mean)

save (model.ancova,m.lst,
      file=here("output", "ancovaFRic_abundW.RData"))

############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FEve RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################

# run model (ancova)
model.ancova.FEve <-  brm (FEve ~ EstRich*Organism,
       data=data_to_modeling_GLM,
       family = gaussian (link="identity"),
       chains=nc,
       iter = ni,
       warmup = nb,
       thin=nt)


# summary of results
summary (model.ancova.FEve)



# plotting
p2<-plot(conditional_effects(model.ancova.FEve,
                           method="fitted",
                           re_formula=NA,
                           robust=T,
                           effects = "EstRich:Organism",
                           points=T,
                           prob = 0.95),
       
       theme = theme_classic() +
         theme (axis.title = element_text(size=15),
                axis.text = element_text(size=12),
                legend.position = "top"),
       points=T)[[1]] + 
  
  scale_color_manual(values=c("#f2a154", "#0e49b5")) + 
  
  xlab("Species richness") + 
  
  ylab ("FEve")

# compare slopes
m.lst.FEve <- emtrends (model.ancova.FEve,  "Organism", var="EstRich")
m.lst_tab.FEve <- summary(m.lst.FEve,point.est = mean)

save (model.ancova.FEve,m.lst.FEve,
      file=here("output", "ancovaFEve_abundW.RData"))


############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FDiv RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################


# run model (ancova)
model.ancova.FDiv <- brm (FDiv ~ EstRich*Organism,
       data=data_to_modeling_GLM,
       family = gaussian (link="identity"),
       chains=nc,
       iter = ni,
       warmup = nb,
       thin=nt)

# summary of results
summary (model.ancova.FDiv)

# plotting
p3<-plot(conditional_effects(model.ancova.FDiv,
                           method="fitted",
                           re_formula=NA,
                           robust=T,
                           effects = "EstRich:Organism",
                           points=T,
                           prob = 0.95),
       
       theme = theme_classic() +
         theme (axis.title = element_text(size=15),
                axis.text = element_text(size=12),
                legend.position = "top"),
       points=T)[[1]] + 
  
  scale_color_manual(values=c("#f2a154", "#0e49b5")) + 
  
  xlab("Species richness") + 
  
  ylab ("FDiv")
  
# compare slopes
m.lst.FDiv <- emtrends (model.ancova.FDiv, "Organism", var="EstRich")
m.lst_tab.FDiv <- summary(m.lst.FDiv, point.est = mean)

# save
save (model.ancova.FDiv,m.lst.FDiv,
      file=here("output", "ancovaFDiv_abundW.RData"))

## arrange these plots into a panel
pdf (here ("output","vectorized","glm_slope_abundW.pdf"),width=7,height=5)
grid.arrange(p1,p2,p3,
             ncol=3,nrow=2)

dev.off()
# ------------------------------------------------------

## ----------------------------------------
# with spatial autocorrelation
## ----------------------------------------

# MCMC settings
ni <- 100000 
nb <- 80000
nt <- 4
nc <- 3


# help with car grouping : https://discourse.mc-stan.org/t/error-when-fitting-a-car-model-in-brms/16687
# distance and neighborhood
# benthos
Grid <- expand.grid(1:nrow(bind_fish_benthos[[1]]), 
			  1:nrow(bind_fish_benthos[[1]]))
K <- nrow(Grid)

require(vegan)
# set up distance and neighbourhood matrices
distance <- data.matrix(vegdist(bind_fish_benthos[[1]][,c("Lat","Lon")],"euclidean"))
W <- ifelse(distance<=5,1,0)
#W_benthos<- list(W=W_benthos)
 
W <- lapply (bind_fish_benthos, function (i){
  
  W <- W[which(rownames(W) %in%  rownames(i)), # need of subsetting to match dims
          which(colnames(W) %in%  rownames(i))]
  rownames(W)<- colnames(W) <- i$Region;
  W
}
)

# set formula (the same for benthos and fishes)
# complete model
formula <- brms::bf(mvbind (log(EstRich),
                            log(FRic),
                            log(FEve),
                            log(FDiv),
                            log(EstRich_benthos),
                            log(FRic_benthos),
                            log(FEve_benthos),
                            log(FDiv_benthos))~ BO2_tempmean_ss_std +distance_std +BO_damean_std + Depth+ 
                      car(W, gr=Region),
                    
                    nl = F)+
 			 set_rescor(TRUE) # nonlinear


# alternative 1
formulaA1 <- brms::bf(mvbind (log(EstRich),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(EstRich_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std +distance_std +BO_damean_std + 
                        car(W,gr=Region),
                  nl = F) +
		     set_rescor(TRUE)# nonlinear

# alternative 2
formulaA2 <- brms::bf(mvbind (log(EstRich),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(EstRich_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std + distance_std + 
                        car(W,gr=Region),
	               nl = F) +
  			   set_rescor(TRUE)# nonlinear

# alternative 3
formulaA3 <- brms::bf(mvbind (log(EstRich),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(EstRich_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std + 
                        car(W, gr=Region),
	                nl = F) +
			   set_rescor(TRUE)# nonlinear

# Benthos
# run MCMC
#setting priors 
priors <- lapply (bind_fish_benthos, function (i) 
      
            c(set_prior("normal(0,5)",class = "b", coef = "",
                        resp=c("logEstRich","logFDiv","logFEve","logFRic",
                               "logEstRichbenthos","logFDivbenthos","logFEvebenthos","logFRicbenthos")),
         	    set_prior("normal(0,5)",class = "Intercept", coef = "",
         	              resp=c("logEstRich","logFDiv","logFEve","logFRic",
         	                     "logEstRichbenthos","logFDivbenthos","logFEvebenthos","logFRicbenthos")),
              
              # CAR priors
              set_prior("gamma(0.01,0.01)", class = "car", coef = "",
                        resp=c("logEstRich","logFDiv","logFEve","logFRic",
                               "logEstRichbenthos","logFDivbenthos","logFEvebenthos","logFRicbenthos"))
	        )
      )

get_prior(formula, data=bind_fish_benthos[[1]],
	data2=list(W=W[[1]]))

# student_t(3, 1, 2.5) did not work
# student_t(3, 5, 2.5) did not work

# run MCMC chains across different organisms and models
   
MCMC_runs_autocorrelation <- #lapply (seq (1,length (bind_fish_benthos)), function (i) #ACROSS ALGORITHMS
      
      lapply (list(formula,formulaA1,formulaA2, formulaA3), function(k) #ACROSS MODELS
         
      
         brms::brm(k, # for each model / formula
             
             data = bind_fish_benthos[[1]], 
		
		         data2 = list(W=W[[1]]),
             
             family = gaussian(),
             
             prior = priors[[1]], 
             
             chains = nc,

		cores = nc,
             
             iter = ni,
             
             warmup = nb,
             
             thin=nt,

              save_pars = save_pars(all = TRUE),


	        control = list(adapt_delta = 0.99) # to avoid disagreement in chain mixing
             
             )
      )
#)
   
# ------------------------------- #
# model selection, across rarefaction algorithms
   
sel_fit_check_autocorrelation <- lapply (MCMC_runs_autocorrelation, function (i) {
      
      # run loo fit test
      loo_test <- loo (i,moment_match=T,reloo=T)
      # extract estimates
      loo_sel <- lapply (loo_test, function (i)
         i$estimates[which(rownames(i$estimates) == "looic"),"Estimate"])
      # extract looic (like AIC)
      tab_mod_sel <- do.call(rbind,lapply (loo_test, function (i)
         i$estimates[which(rownames(i$estimates) == "looic"),]))
      # extract estimated number of parameters (model adequacy)
      tab_mod_fit <- do.call(rbind,lapply (loo_test, function (i)
         i$estimates[which(rownames(i$estimates) == "p_loo"),]))
      # select the model with lowest looic
      sel_model <- i[which(loo_sel == min(unlist(loo_sel)))]
   
      # list of results 
      res <- list (looic = tab_mod_sel,
                   param= tab_mod_fit,
                   best_model = sel_model
                   )
      ;
      res
      

      } # Close model check and sel 
      )
   
## save
save ( MCMC_runs_autocorrelation,sel_fit_check_autocorrelation, 
	file=here ("output", "MCMC_runs_multivariate_rarefied_aut.Rdata"))


#




























# -----------------------------------------------

## List of results for figure 3 (regression coefficients)

load(here ("output","MCMC_runs_SR.RData"))
load(here ("output","MCMC_runs_FRic.RData"))
load(here ("output","MCMC_runs_FEve.RData"))
load(here ("output","MCMC_runs_FDiv.RData"))

# goodness of fit

gof_models <- list(lapply (analysis_SR, function (i) sapply (i, "[","looic")),
                   lapply (analysis_FRic, function (i) sapply (i, "[","looic")), 
                   lapply (analysis_FEve, function (i) sapply (i, "[","looic")),
                   lapply (analysis_FDiv, function (i) sapply (i, "[","looic"))
)


# FISH
## coefficients
list_of_results <- list(lapply (analysis_SR, function (i) sapply (i, "[[","best_model")),
                        lapply (analysis_FRic, function (i) sapply (i, "[[","best_model")), 
                        lapply (analysis_FEve, function (i) sapply (i, "[[","best_model")),
                        lapply (analysis_FDiv, function (i) sapply (i, "[[","best_model"))
)
# set names in indexes
names(list_of_results) <- c("SR", "FRic", "FEve","FDiv")

# calculate bayes R2
bayes_R2_models <-lapply (list_of_results, function (i) # index
  lapply (i, function (k) # organism
    do.call(rbind,lapply (k, function (j) # algorithm
  
       bayes_R2(j)
       
      ))
    )
  )

## goodness of fit
gof_best <- list(lapply (analysis_SR, function (i) sapply (i, "[","param")),
                        lapply (analysis_FRic, function (i) sapply (i, "[","param")), 
                        lapply (analysis_FEve, function (i) sapply (i, "[","param")),
                        lapply (analysis_FDiv, function (i) sapply (i, "[","param"))
)

# set names of algorithms within indexes 
#list_of_results <- lapply (list_of_results, function (i) {
#   names (i) <- c("MSS","PSS", "Lomolino");
#   i
#})

## organizing results (coeficients)
algorithms <- c("MSS","PSS", "Lomolino")
org_results <- lapply (seq(1,length (list_of_results)), function (i)# across index
   lapply (seq (length(list_of_results[[1]])), function (k) # across taxa
     lapply (seq (length(list_of_results[[1]][[1]])), function (j) # across algorithms
      
      data.frame (
         fixef(list_of_results[[i]][[k]][[j]],
               summary = TRUE,
               robust = FALSE,
               probs = c(0.025, 0.975)
         ), 
         "Algorithm" = algorithms[j],
         "Index" = names(list_of_results)[i],
         "Organism" = names(list_of_results[[i]])[k])
   )))

# melt list 
org_results <- lapply (org_results, function (k) 
  
              do.call(rbind, 
                       
                       lapply (k, function (i)
                          
                          do.call(rbind , i)
                          
                       )
  )
)

names(org_results) <- c("SR", "FRic", "FEve","FDiv")

## save results
save (list_of_results,
      org_results, 
      bayes_R2_models,
      file=here("output", "summarized_results_drivers.RData"))

# --------------------------------------------- #
# predict values to the whole BR coast

#layers <- list_layers()
#View (layers [grep ("Bio-ORACLE",layers$dataset_code),])

# Download specific layers to the current directory
# set prefered folder (to download data)
options(sdmpredictors_datadir=here ("data","environment"))

## chlorophil has different extent - loading and extracting in two steps         
layers_oracle <- load_layers(c("BO2_tempmean_ss",#"BO2_temprange_ss",
                               "BO2_ppmean_ss", #"BO2_pprange_ss",
                               #"BO_damax",
                               "BO_damean"#,"BO_damin"
                               ))
# ----------------------
# load batimetry shape
batimetry <- readOGR(here ("data","environment","batimetria_lito"), "BATIMETRIA_SRTM_30_POLIGONO_SIRGAS_2000")

# subset batimetry
batimetry<-batimetry [which(batimetry$PROFUNDIDA == "0 a -20 m"),]

# transform batimetry
batimetry<-spTransform(batimetry, crs(layers_oracle))

# subset biooracle layer
masked_layers <- mask(layers_oracle, batimetry)
#crop it to decrase extent
masked_layers<-crop(masked_layers,batimetry)
# rna <- reclassify(masked_layers, cbind(NA, -999))

# get the centroid of each cell
centroids <- coordinates (masked_layers)
centroids <- SpatialPoints(centroids)
crs(centroids) <- crs (masked_layers)

# get the centroids witihin batimetry
over_pts <- (over(centroids,batimetry))
centroids_batimetry<-centroids[which(is.na(over_pts$PROFUNDIDA) != T)]

# extract covariates in these points
masked_layers_batimetry <- (extract (masked_layers,centroids_batimetry,
                            method='simple', fun=mean))

# distance to coast
# BR coast, download from here https://mapcruzin.com/free-brazil-arcgis-maps-shapefiles.htm

BR <- readOGR(dsn=here("data", "environment","brazil-coastline"), "brazil_coastline")
crs(BR) <- crs(centroids_batimetry)
#BR <- spTransform(BR, CRS("+init=epsg:4326"))

# use dist2Line from geosphere - only works for WGS84 
# measure distance
require(geosphere)
dist_cells_batimetry <- lapply (seq (1,length(centroids_batimetry)), function (i) 
  
                    dist2Line(p = centroids_batimetry[i], 
                              line = (BR))
                    )
#save(dist_cells_batimetry,file=here("output","dist_cells_batimetry.RData"))
load(here("output","dist_cells_batimetry.RData"))
dist_cells_batimetry_df<-do.call(rbind,dist_cells_batimetry)

# cbind distance
masked_layers_batimetry <- cbind (masked_layers_batimetry, 
                                  dist_cells_batimetry_df,
                                  centroids_batimetry@coords)

# rm islands (too far from coast)
masked_layers_batimetry <- masked_layers_batimetry [which(masked_layers_batimetry[,"distance"]<= 150000),]

# remove too southern sites
masked_layers_batimetry<-masked_layers_batimetry[which(masked_layers_batimetry[,"lat"]>-28.5),]

# plot
ggplot(data = world) +
  geom_sf() +
  geom_point(data=data.frame(masked_layers_batimetry), 
             aes(x=x, y=y, col =distance   )) +
  coord_sf(xlim = c(-55, -20), ylim = c(-40,6), expand = FALSE)+
  theme_classic()

# standardize covariates
#E log distance
masked_layers_batimetry[,"distance"] <- log (masked_layers_batimetry[,"distance"])
# std
masked_layers_batimetry_std <- data.frame (
  
  BO2_tempmean_ss_std=(masked_layers_batimetry[,"BO2_tempmean_ss"]-mean(masked_layers_batimetry[,"BO2_tempmean_ss"],na.rm=T))/sd(masked_layers_batimetry[,"BO2_tempmean_ss"],na.rm=T),
  BO_damean_std=(masked_layers_batimetry[,"BO_damean"]-mean(masked_layers_batimetry[,"BO_damean"],na.rm=T))/sd(masked_layers_batimetry[,"BO_damean"],na.rm=T),
  distance_std=(masked_layers_batimetry[,"distance"]-mean(masked_layers_batimetry[,"distance"],na.rm=T))/sd(masked_layers_batimetry[,"distance"],na.rm=T),
  Depthraso="fundo",
  masked_layers_batimetry[,c("x","y")]
  )

# rm NAs
masked_layers_batimetry_std<- masked_layers_batimetry_std[is.na(masked_layers_batimetry_std$BO2_tempmean_ss_std)==F,]

## model to predict

plots_predictions <-lapply (list_of_results, function (i) # indexes
  lapply (i, function (k) # organism
    lapply (k, function (j) { # algorithm
      
        # sel model    
        model_to_pred <- j
        # name to plot legend
        legend_label <- colnames(model_to_pred$data)[1]
        # find parameters
        params_in_the_model <- (fixef(model_to_pred)[,"Estimate"])
        # subset data according to data found in the model
        tab_subset <- data.frame(masked_layers_batimetry_std[,
                                                which(colnames(masked_layers_batimetry_std)
                                                             %in%
                                                        names(params_in_the_model))])
        colnames(tab_subset)<- colnames(masked_layers_batimetry_std)[which(colnames(masked_layers_batimetry_std)
                                                %in%
                                                  names(params_in_the_model))]
        # make predictions (using brms functions)
        colnames(tab_subset)[which(colnames(tab_subset) == "Depthraso")] <- "Depth"
        pred_model <- predict(model_to_pred,
                            newdata=tab_subset)
        
        # bind predictions in the df
        tab_subset <- cbind(masked_layers_batimetry_std,
                            pred_model=pred_model)
        tab_subset<-tab_subset[,c("x","y","pred_model.Estimate")]
        
        # coordinates
        coordinates(tab_subset) <- ~x+y
        
        # coerce to SpatialPixelsDataFrame
        gridded(tab_subset) <- TRUE
        # coerce to raster
        rasterDF <- raster(tab_subset)
        # expoent for richness
        if (length(grep ("Rich",model_to_pred$formula))!=0) {
          rasterDF<-exp(rasterDF)
          } else {rasterDF}
        
        # aggregate data at lower resolution
        rasterDF <- aggregate(rasterDF, fact=5)
        rasterDF<-clamp(rasterDF,lower=0)
        
        # plot
        plot_raster <- gplot(rasterDF) +
          geom_tile(aes(x=x, y=y, fill=round(value,2)),alpha=1) +
            coord_sf(xlim = c(-55, -30), ylim = c(-30,6), expand = FALSE)+
          scale_fill_viridis(direction=-1, begin = 0, 
                             breaks = seq(range(values(rasterDF),na.rm=T)[1],
                                          range(values(rasterDF),na.rm=T)[2],
                                          (range(values(rasterDF),na.rm=T)[2]-range(values(rasterDF),na.rm=T)[1])/5),
                             limits=c(range(values(rasterDF),na.rm=T)[1],
                                      range(values(rasterDF),na.rm=T)[2]),
                             na.value=NA)+
          theme_classic() + 
          labs(fill=legend_label)+
          theme (legend.position = "right")
        
        # results to report
        res_list <- list (plot=plot_raster,
                          predictions=pred_model,
                          raster=rasterDF)
        
        ;
        
        res_list
        
        }
        )
    )
  )


## plot across algorithms of rarefaction

algorithm <- c("MSS","PSS","lomolino")

# --------------------------------------------

lapply (seq (1,length(algorithm)), function (i) {
    
    # SPECIES RICHNESS
    # arrangment of predictions
    a<-plots_predictions[[1]][[1]][[i]]$plot
    b<-plots_predictions[[1]][[2]][[i]]$plot
    
    df_cor <- data.frame (b=plots_predictions[[1]][[1]][[i]]$predictions[,"Estimate"],
                          p=plots_predictions[[1]][[2]][[i]]$predictions[,"Estimate"])
                          
    cor_plot1<-ggplot (df_cor,aes (x=exp(b),y=exp(p))) + 
      geom_point()+
      geom_smooth(method="lm") + 
      theme_classic()+
      xlab("Benthic richness") + 
      ylab("Fish richness") + 
      annotate("text",x = min(exp(df_cor$b)), y = min(exp(df_cor$p)), 
               label = paste ("Pearson's r2=",
                              round (cor(df_cor)[2,1],2)),
               size=3) + 
      xlim (min(exp(df_cor$b)),max(exp(df_cor$b)))+
      ylim (min(exp(df_cor$p)),max(exp(df_cor$p)))+
      theme (axis.title = element_text(size=9),
             axis.text = element_text(size=7)
             )
    
    # arranging plots
    pdf(here ("output", "vectorized", paste ("richness", algorithm[i],".pdf")))
    grid.arrange(a,b,cor_plot1,
                 ncol=6,nrow=8,
                 layout_matrix = rbind (c(NA,NA,3,3,NA,NA),
                                        c(NA,NA,3,3,NA,NA),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2)))
    
    dev.off()
    # --------------------------------------------------
    # FUNCTIONAL RICHNESS
    # arrangment of predictions
    a<-plots_predictions[[2]][[1]][[i]]$plot
    b<-plots_predictions[[2]][[2]][[i]]$plot
    
    df_cor <- data.frame (b=plots_predictions[[2]][[1]][[i]]$predictions[,"Estimate"],
                          p=plots_predictions[[2]][[2]][[i]]$predictions[,"Estimate"])
    
    cor_plot1<-ggplot (df_cor,aes (x=(b),y=(p))) + 
      geom_point()+
      geom_smooth(method="lm") + 
      theme_classic()+
      xlab("Benthic FRic") + 
      ylab("Fish FRic") + 
      annotate("text",x = min((df_cor$b)), y = min((df_cor$p)), 
               label = paste ("Pearson's r2=",
                              round (cor(df_cor)[2,1],2)),
               size=3) + 
      xlim (min((df_cor$b)),max((df_cor$b)))+
      ylim (min((df_cor$p)),max((df_cor$p)))+
      theme (axis.title = element_text(size=9),
             axis.text = element_text(size=7)
      )
    
    # arranging plots
    pdf(here ("output", "vectorized", paste ("fric", algorithm[i],".pdf")))
    grid.arrange(a,b,cor_plot1,
                 ncol=6,nrow=8,
                 layout_matrix = rbind (c(NA,NA,3,3,NA,NA),
                                        c(NA,NA,3,3,NA,NA),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2)))
    
    dev.off()
    
    # --------------------------------------------------
    # FUNCTIONAL EVENNESS
    # arrangment of predictions
    a<-plots_predictions[[3]][[1]][[i]]$plot
    b<-plots_predictions[[3]][[2]][[i]]$plot
    
    df_cor <- data.frame (b=plots_predictions[[3]][[1]][[i]]$predictions[,"Estimate"],
                          p=plots_predictions[[3]][[2]][[i]]$predictions[,"Estimate"])
    
    cor_plot1<-ggplot (df_cor,aes (x=(b),y=(p))) + 
      geom_point()+
      geom_smooth(method="lm") + 
      theme_classic()+
      xlab("Benthic FEve") + 
      ylab("Fish FEve") + 
      annotate("text",x =min((df_cor$b)), y = min((df_cor$p)), 
               label = paste ("Pearson's r2=",
                              round (cor(df_cor)[2,1],2)),
               size=3) + 
      xlim (min((df_cor$b)),max((df_cor$b)))+
      ylim (min((df_cor$p)),max((df_cor$p)))+
      theme (axis.title = element_text(size=9),
             axis.text = element_text(size=7)
      )
    
    # arranging plots
    pdf(here ("output", "vectorized", paste ("feve", algorithm[i],".pdf")))
    grid.arrange(a,b,cor_plot1,
                 ncol=6,nrow=8,
                 layout_matrix = rbind (c(NA,NA,3,3,NA,NA),
                                        c(NA,NA,3,3,NA,NA),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2)))
    dev.off()
    
    # --------------------------------------------------
    # FUNCTIONAL DIVERGENCE
    # arrangment of predictions
    a<-plots_predictions[[4]][[1]][[i]]$plot
    b<-plots_predictions[[4]][[2]][[i]]$plot
    
    df_cor <- data.frame (b=plots_predictions[[4]][[1]][[i]]$predictions[,"Estimate"],
                          p=plots_predictions[[4]][[2]][[i]]$predictions[,"Estimate"])
    
    cor_plot1<-ggplot (df_cor,aes (x=(b),y=(p))) + 
      geom_point()+
      geom_smooth(method="lm") + 
      theme_classic()+
      xlab("Benthic FDiv") + 
      ylab("Fish FDiv") + 
      annotate("text",x = min((df_cor$b)), y = min((df_cor$p)), 
               label = paste ("Pearson's r2=",
                              round (cor(df_cor)[2,1],2)),
               size=3) + 
      xlim (min((df_cor$b)),max((df_cor$b)))+
      ylim (min((df_cor$p)),max((df_cor$p)))+
      theme (axis.title = element_text(size=9),
             axis.text = element_text(size=7)
      )
    
    # arranging plots
    pdf(here ("output", "vectorized", paste ("fdiv", algorithm[i],".pdf")))
    grid.arrange(a,b,cor_plot1,
                 ncol=6,nrow=8,
                 layout_matrix = rbind (c(NA,NA,3,3,NA,NA),
                                        c(NA,NA,3,3,NA,NA),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2),
                                        c(1,1,1,2,2,2)))
    
    dev.off()
    
    }
)
#



# extract all predictions from MSS
MSS <- do.call (cbind ,lapply (plots_predictions, function (i)
  do.call(cbind,lapply (i, function (k)
    
    k[[1]]$predictions[,1]
    
  
  )
  )))

# change colnames
colnames(MSS) <- c("SR.Benthos","SR.Fish",
                   "FRic.Benthos","FRic.Fish",
                   "FEve.Benthos","FEve.Fish",
                   "FDiv.Benthos","FDiv.Fish")
MSS <- MSS[,order(colnames(MSS),decreasing=T)]
corrplot(cor(MSS),
         method = "circle",
         is.corr = T,
         type="lower",
         p.mat = cor.mtest(MSS)$p,
         diag=F,
         sig.level = 0.05)


# extract all predictions from MSS
PSS <- do.call (cbind ,lapply (plots_predictions, function (i)
  do.call(cbind,lapply (i, function (k)
    
    k[[2]]$predictions[,1]
    
    
  )
  )))

# change colnames
colnames(PSS) <- c("SR.Benthos","SR.Fish",
                   "FRic.Benthos","FRic.Fish",
                   "FEve.Benthos","FEve.Fish",
                   "FDiv.Benthos","FDiv.Fish")
PSS <- PSS[,order(colnames(PSS),decreasing=T)]
corrplot(cor(PSS),
         method = "circle",
         is.corr = T,
         type="lower",
         p.mat = cor.mtest(PSS)$p,
         diag=F,
         sig.level = 0.05)

# extract all predictions from S50
S50 <- do.call (cbind ,lapply (plots_predictions, function (i)
  do.call(cbind,lapply (i, function (k)
    
    k[[3]]$predictions[,1]
    
    
  )
  )))

# change colnames
colnames(S50) <- c("SR.Benthos","SR.Fish",
                   "FRic.Benthos","FRic.Fish",
                   "FEve.Benthos","FEve.Fish",
                   "FDiv.Benthos","FDiv.Fish")
S50 <- S50[,order(colnames(S50),decreasing=T)]
corrplot(cor(S50),
         method = "circle",
         is.corr = T,
         type="lower",
         p.mat = cor.mtest(S50)$p,
         diag=F,
         sig.level = 0.05)













############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FD RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################

# run model (ancova)
model.ancova <- lapply (data_to_modeling_GLM, function (i)
  
  brm (FRic ~ EstRich*Organism,
       data=i,
       family = gaussian (link="identity"),
       chains=3,
       iter = 30000,
       warmup = 20000,
       thin=20)
)

# summary of results
lapply (model.ancova,summary)

# plotting
lapply (seq(1,3), function (i) {
  
  plot(conditional_effects(model.ancova[[i]],
                           method="fitted",
                           re_formula=NA,
                           robust=T,
                           effects = "EstRich:Organism",
                           points=T,
                           prob = 0.95),
       
       theme = theme_classic() +
         
         theme (axis.title = element_text(size=15),
                axis.text = element_text(size=12),
                legend.position = "none"),
       points=T)
  
  
  
  ggsave (here("output",paste (i,"ancovaplot.pdf")),height = 5,width=5)
  
})

# compare slopes
m.lst <- lapply (model.ancova, emtrends, "Organism", var="EstRich")

m.lst_tab <- lapply(m.lst,summary, point.est = mean)
do.call(rbind,m.lst_tab)

save (model.ancova,m.lst,
      file=here("output", "ancova_FRic.RData"))

############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FEve RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################

# run model (ancova)
model.ancova.FEve <-lapply (data_to_modeling_GLM, function (i)
  
  brm (FEve ~ EstRich*Organism,
       data=i,
       family = gaussian (link="identity"),
       chains=3,
       iter = 30000,
       warmup = 20000,
       thin=20)
)

# summary of results
lapply (model.ancova.FEve,summary)

# plotting
lapply (seq(1,3), function (i) {
  
  plot(conditional_effects(model.ancova.FEve[[i]],
                           method="fitted",
                           re_formula=NA,
                           robust=T,
                           effects = "EstRich:Organism",
                           points=T,
                           prob = 0.95),
       
       theme = theme_classic() +
         theme (axis.title = element_text(size=15),
                axis.text = element_text(size=12),
                legend.position = "none"),
       points=T)
  
  ggsave (here("output",paste (i,"ancovaplot_Feve.pdf")),height = 5,width=5)
  
})

# compare slopes
m.lst.FEve <- lapply (model.ancova.FEve, emtrends, "Organism", var="EstRich")

m.lst_tab.FEve <- lapply(m.lst.FEve,summary, point.est = mean)
do.call(rbind,m.lst_tab.FEve)

save (model.ancova.FEve,m.lst.FEve,
      file=here("output", "ancovaFEve.RData"))


############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FDiv RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################


# run model (ancova)
model.ancova.FDiv <-lapply (data_to_modeling_GLM, function (i)
  
  brm (FDiv ~ EstRich*Organism,
       data=i,
       family = gaussian (link="identity"),
       chains=3,
       iter = 30000,
       warmup = 20000,
       thin=20)
)

# summary of results
lapply (model.ancova.FDiv,summary)

# plotting
lapply (seq(1,3), function (i) {
  
  plot(conditional_effects(model.ancova.FDiv[[i]],
                           method="fitted",
                           re_formula=NA,
                           robust=T,
                           effects = "EstRich:Organism",
                           points=T,
                           prob = 0.95),
       
       theme = theme_classic() +
         theme (axis.title = element_text(size=15),
                axis.text = element_text(size=12),
                legend.position = "none"),
       points=T)
  
  ggsave (here("output",paste (i,"ancovaplotFDiv.pdf")),height = 5,width=5)
  
})

# compare slopes
m.lst.FDiv <- lapply (model.ancova.FDiv, emtrends, "Organism", var="EstRich")

m.lst_tab.FDiv <- lapply(m.lst.FDiv,summary, point.est = mean)
do.call(rbind,m.lst_tab.FDiv)

save (model.ancova.FDiv,m.lst.FDiv,
      file=here("output", "ancovaFDiv.RData"))
