
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

load (here ("output_no_resampling","data_to_modeling_GLM.RData"))


ggplot(data_to_modeling_GLM,aes(x=,BO2_tempmean_ss_std ,y=log(SR))) + geom_point() + 
  geom_smooth(method="lm",formula = y ~ x)

ggplot(data_to_modeling_GLM,aes(x=,area ,y=log(SR))) + geom_point() + 
  geom_smooth(method="lm",formula = y ~ x)

ggplot(data_to_modeling_GLM,aes(x=,SamplingArea_std,y=log(SR))) + geom_point() + 
  geom_smooth(method="lm",formula = y ~ x)


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
                            SR_benthos = cov_benthos$SR,
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
formula <- brms::bf(mvbind (log(SR),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(SR_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std +
                        distance_std +
                        area +  
                        BO_damean_std +
                      Depth + 
                      SamplingArea_std,
                      nl = F) +
  set_rescor(TRUE)# nonlinear

# alternative 1
formulaA1 <- brms::bf(mvbind (log(SR),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(SR_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std +
                                                    distance_std +
                                                    area +  
                                                    BO_damean_std + 
                        SamplingArea_std ,
                      nl = F) +
  set_rescor(TRUE)# nonlinear

# alternative 2
formulaA2 <- brms::bf(mvbind (log(SR),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(SR_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std + 
                                                    distance_std +
                                                    area + 
                        SamplingArea_std,
                      nl = F) +
  set_rescor(TRUE)# nonlinear

# alternative 3
formulaA3 <- brms::bf(mvbind (log(SR),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(SR_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std + area + 
                        SamplingArea_std,
                      nl = F) +
  set_rescor(TRUE)# nonlinear


# alternative 4
formulaA4 <- brms::bf(mvbind (log(SR),
                              log(FRic),
                              log(FEve),
                              log(FDiv),
                              log(SR_benthos),
                              log(FRic_benthos),
                              log(FEve_benthos),
                              log(FDiv_benthos)) ~ BO2_tempmean_ss_std + 
                        SamplingArea_std,
                      nl = F) +
  set_rescor(TRUE)# nonlinear


#setting priors 
priors <-  c(set_prior("normal(0,5)",class = "b",coef = "",
                resp=c("logSR","logFDiv","logFEve","logFRic",
                     "logSRbenthos","logFDivbenthos","logFEvebenthos","logFRicbenthos")),
              set_prior("normal(0,5)",class = "Intercept",coef = "",
                resp=c("logSR","logFDiv","logFEve","logFRic",
                     "logSRbenthos","logFDivbenthos","logFEvebenthos","logFRicbenthos"))
  )


# MCMC settings
ni <- ni 
nb <- nb
nt <- nt
nc <- nc

# run MCMC chains across different organisms and models
MCMC_runs <- lapply (list(formula,
                          formulaA1,
                          formulaA2, 
                          formulaA3,
                          formulaA4), function(k) #ACROSS MODELS
    
    
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
       file=here ("output_no_resampling", 
                  "MCMC_runs_multivariate_rarefied_no_aut.Rdata"))


############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FD RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################

# run model (ancova)
model.ancova <- brm (FRic ~ SR*Organism,
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
                         effects = "SR:Organism",
                         points=T,
                         prob = 0.95),
     
     theme = theme_classic() +
      
        theme (axis.title = element_text(size=15),
               axis.text = element_text(size=12),
               legend.position = "top") ,
     points=T) [[1]] + 
  
  scale_color_manual(values=c("#f2a154", "#0e49b5")) + 
  
  xlab("Species richness") + 
  
  ylab ("FRic") +
  
  ylim(c(-0.05,1.0))
    
p1

# compare slopes
m.lst <- emtrends (model.ancova, "Organism", var="SR")
m.lst_tab <- summary(m.lst,point.est = mean)

save (model.ancova,m.lst,
      file=here("output_no_resampling", "ancovaFRic.RData"))

############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FEve RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################

# run model (ancova)
model.ancova.FEve <-  brm (FEve ~ SR*Organism,
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
                           effects = "SR:Organism",
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
m.lst.FEve <- emtrends (model.ancova.FEve,  "Organism", var="SR")
m.lst_tab.FEve <- summary(m.lst.FEve,point.est = mean)

save (model.ancova.FEve,m.lst.FEve,
      file=here("output_no_resampling", "ancovaFEve.RData"))


############################################################################
# -------------------------------------------------------------------------
#          ANCOVA
#
#        ACCUMULATION OF FDiv RELATIVE TO RICHNESS
#
# -------------------------------------------------------------------------
############################################################################


# run model (ancova)
model.ancova.FDiv <- brm (FDiv ~ SR*Organism,
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
                           effects = "SR:Organism",
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
m.lst.FDiv <- emtrends (model.ancova.FDiv, "Organism", var="SR")
m.lst_tab.FDiv <- summary(m.lst.FDiv, point.est = mean)

# save
save (model.ancova.FDiv,m.lst.FDiv,
      file=here("output_no_resampling", "ancovaFDiv.RData"))

## arrange these plots into a panel
pdf (here ("output_no_resampling","glm_slope.pdf"),width=7,height=5)
grid.arrange(p1,p2,p3,
             ncol=3,nrow=2)

dev.off()
# ------------------------------------------------------

## ----------------------------------------
# with spatial autocorrelation
## ----------------------------------------


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

rm(list=ls())