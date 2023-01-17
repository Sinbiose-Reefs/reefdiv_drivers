
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

load (here ("output_no_resampling_sensitivity","data_to_modeling_GLM.RData"))


ggplot(data_to_modeling_GLM,aes(x=BO2_tempmean_ss_std ,y=log(SR))) + geom_point() + 
  geom_smooth(method="lm",formula = y ~ x)
ggplot(data_to_modeling_GLM,aes(x=area ,y=log(SR))) + geom_point() + 
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
                        SamplingArea_std,
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
                        area+
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
MCMC_runs <- lapply (list(formulaA1,
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
       file=here ("output_no_resampling_sensitivity", 
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
  
  ylab ("FRic")
    
  
# compare slopes
m.lst <- emtrends (model.ancova, "Organism", var="SR")
m.lst_tab <- summary(m.lst,point.est = mean)

save (model.ancova,m.lst,
      file=here("output_no_resampling_sensitivity", "ancovaFRic.RData"))

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
      file=here("output_no_resampling_sensitivity", "ancovaFEve.RData"))


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
      file=here("output_no_resampling_sensitivity", "ancovaFDiv.RData"))

## arrange these plots into a panel
pdf (here ("output_no_resampling_sensitivity","glm_slope.pdf"),width=7,height=5)
grid.arrange(p1,p2,p3,
             ncol=3,nrow=2)

dev.off()
# ------------------------------------------------------
