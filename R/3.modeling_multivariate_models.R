
# ----------------------------------------------------------------------------#
#    routine to modeling fish and benthos SR and FD  relative to environment
#                          using GLM
# 	PS: RUN IT IN R GUI
# ----------------------------------------------------------------------------#

## laoding packages
source("R/packages.R")
source("R/functions.R")

## function to test space quality (from Maire et al. 2015)
source("R/quality_funct_space_fromdist2.R")

# ------------------------------------------ #
# Load rarefied data of fishes and benthos
# sites x transect / video x species
# ------------------------------------------ #

load (here ("output","modeling_data.RData"))
load (here ("output","FD_results.RData"))

# standardize covariate data
site_covs$sst_std <- (site_covs$sst - mean(site_covs$sst))/sd(site_covs$sst)
site_covs$turbidity_std <- (site_covs$turbidity - mean(site_covs$turbidity))/sd(site_covs$turbidity)
site_covs$productivity_std <- (site_covs$productivity - mean(site_covs$productivity))/sd(site_covs$productivity)
site_covs$salinity_std <- (site_covs$salinity - mean(site_covs$salinity))/sd(site_covs$salinity)
site_covs$offshore_distance_std <- (site_covs$offshore_distance - mean(site_covs$offshore_distance))/sd(site_covs$offshore_distance)

# effort data
effort_dataframe$fish_effort_std <- (effort_dataframe$fish_effort - mean(effort_dataframe$fish_effort))/sd(effort_dataframe$fish_effort)
effort_dataframe$benthos_effort_std <- (effort_dataframe$benthos_effort - mean(effort_dataframe$benthos_effort))/sd(effort_dataframe$benthos_effort)




# -------------------------------------
# prepare data to multivariate models
# -------------------------------------

# create a DF with all data

bind_fish_benthos<- cbind (site_covs,
                           effort_dataframe,
                           # benthos
                           SR_corals = scale (df_corals$SR),
                           FRic_corals = df_corals$FRic,
                           Rao_corals = df_corals$RaoQ,
                           # fish
                           SR_fish = scale (FD_fish$nbsp),
                           FRic_fish = FD_fish$FRic,
                           Rao_fish = FD_fish$RaoQ,
                           # algae
                           SR_algae = scale (FD_algae$nbsp),
                           FRic_algae = FD_algae$FRic,
                           Rao_algae = FD_algae$RaoQ
                           
)

# average values of metrics to present in the Results
round(apply(bind_fish_benthos[c(3:16,18:30)],2,mean,na.rm=T),3)
round(apply(bind_fish_benthos[c(3:16,18:30)],2,sd,na.rm=T),3)


# correlation between variables
cor(bind_fish_benthos[c(3:9)])





############################################################################
# -------------------------------------------------------------------------
#
#          Multivariate Linear models in BRMS
#
#        Effect of Environmental drivers (run directly on R GUI)
#
#
# -------------------------------------------------------------------------
############################################################################


# example here:
# https://cran.r-project.org/web/packages/brms/vignettes/brms_multivariate.html
# set formula (the same for benthos and fishes)
# complete model
formula1_fish <- brms::bf(mvbind (log(FRic_fish+1),
                                  log(Rao_fish+1)) ~ SR_fish+sst_std + # fixed effects
                                                turbidity_std +
                                                salinity_std  +
                            (1|region), # random effects (intercept)
                      nl = F)


# complete model
formula1_algae <- brms::bf(mvbind (log(FRic_algae+1),
                                    log(Rao_algae+1)) ~ SR_algae+sst_std + # fixed effects
                                      turbidity_std +
                                      salinity_std  + 
                               
                              (1|region), # random effect (intercept) 
                          nl = F) 
# complete model
formula1_corals <- brms::bf(mvbind (log(FRic_corals+1),
                                     log(Rao_corals+1)) ~ SR_corals+sst_std + # fixed effects
                               turbidity_std +
                               salinity_std  + 
                               
                               (1|region), # random effect (intercept) 
                             nl = F) 

#setting priors 
priors <-  c(set_prior("normal(0,4)",class = "b",coef = "",
                       resp=c("logFRicfish1",
                              "logRaofish1",
                              "logFRicalgae1",
                              "logRaoalgae1",
                              "logFRiccorals1",
                              "logRaocorals1")),
             set_prior("normal(0,4)",class = "Intercept",coef = "",
                       resp=c("logFRicfish1",
                              "logRaofish1",
                              "logFRicalgae1",
                              "logRaoalgae1",
                              "logFRiccorals1",
                              "logRaocorals1"))
)


# fit the model
fit_complete <- brms::brm(mvbf(formula1_fish, formula1_algae,formula1_corals) + # bind multivariate response

                    set_rescor(TRUE), # explicitly asking to estimate the residual correlation 
            data = bind_fish_benthos, 
            chains = 3, cores = 3, # N chains and cores for parallel processing
            iter = 20000,
            warmup = 18000,
            thin=1,
            prior = priors,
            seed=1234,
		save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.99,
				   max_treedepth = 15))

#

save(fit_complete,file=here ("output","fit_complete.RData"))



# ******************************
# alternative model without salinity

# set formula (the same for benthos and fishes)
# complete model
formula1_fish <- brms::bf(mvbind (log(FRic_fish+1),
                                  log(Rao_fish+1)) ~ SR_fish+sst_std + # fixed effects
                                                turbidity_std +
                            (1|region), # random effects (intercept)
                      nl = F)


# complete model
formula1_algae <- brms::bf(mvbind (log(FRic_algae+1),
                                    log(Rao_algae+1)) ~ SR_algae+sst_std + # fixed effects
                                      turbidity_std + 
                               
                              (1|region), # random effect (intercept) 
                          nl = F) 
# complete model
formula1_corals <- brms::bf(mvbind (log(FRic_corals+1),
                                     log(Rao_corals+1)) ~ SR_corals+sst_std + # fixed effects
                               turbidity_std +
                               
					(1|region), # random effect (intercept) 
                             nl = F) 

#setting priors 
priors <-  c(set_prior("normal(0,4)",class = "b",coef = "",
                       resp=c("logFRicfish1",
                              "logRaofish1",
                              "logFRicalgae1",
                              "logRaoalgae1",
                              "logFRiccorals1",
                              "logRaocorals1")),
             set_prior("normal(0,4)",class = "Intercept",coef = "",
                       resp=c("logFRicfish1",
                              "logRaofish1",
                              "logFRicalgae1",
                              "logRaoalgae1",
                              "logFRiccorals1",
                              "logRaocorals1"))
)


# fit the model
fit2 <- brms::brm(mvbf(formula1_fish, formula1_algae,formula1_corals) + 
                    set_rescor(TRUE), 
            data = bind_fish_benthos, 
            chains = 3, cores = 3,
            iter = 20000,
            warmup = 18000,
            thin=1,
            prior = priors,
            seed=1234,
		save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.99,
				   max_treedepth = 15))

#

save(fit2,file=here ("output","fit2.RData"))







# ******************************
# alternative model without salinity and turbidity


# set formula (the same for benthos and fishes)
# complete model
formula1_fish <- brms::bf(mvbind (log(FRic_fish+1),
                                  log(Rao_fish+1)) ~ SR_fish+sst_std + # fixed effects
                                                
                            (1|region), # random effects (intercept)
                      nl = F)


# complete model
formula1_algae <- brms::bf(mvbind (log(FRic_algae+1),
                                    log(Rao_algae+1)) ~ SR_algae+sst_std + # fixed effects
                                                                    (1|region), # random effect (intercept) 
                          nl = F) 
# complete model
formula1_corals <- brms::bf(mvbind (log(FRic_corals+1),
                                     log(Rao_corals+1)) ~ SR_corals+sst_std + # fixed effects
                                                              (1|region), # random effect (intercept) 
                             nl = F) 


# fit the model
fit3 <- brms::brm(mvbf(formula1_fish, formula1_algae,formula1_corals) + 
                    set_rescor(TRUE), 
            data = bind_fish_benthos, 
            chains = 3, cores = 3,
            iter = 20000,
            warmup = 18000,
            thin=1,
            prior = priors,
            seed=1234,
		save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.99,
				   max_treedepth = 15))

#

save(fit3,file=here ("output","fit3.RData"))


# ********************************
# without SR and other covariates; only SST



# set formula (the same for benthos and fishes)
# complete model
formula1_fish <- brms::bf(mvbind (log(FRic_fish+1),
                                  log(Rao_fish+1)) ~ sst_std + # fixed effects
                                                
                            (1|region), # random effects (intercept)
                      nl = F)


# complete model
formula1_algae <- brms::bf(mvbind (log(FRic_algae+1),
                                    log(Rao_algae+1)) ~ sst_std + # fixed effects
                                                                    (1|region), # random effect (intercept) 
                          nl = F) 
# complete model
formula1_corals <- brms::bf(mvbind (log(FRic_corals+1),
                                     log(Rao_corals+1)) ~ sst_std + # fixed effects
                                                              (1|region), # random effect (intercept) 
                             nl = F) 


# fit the model
fit_simple <- brms::brm(mvbf(formula1_fish, formula1_algae,formula1_corals) + 
                    set_rescor(TRUE), 
            data = bind_fish_benthos, 
            chains = 3, cores = 3,
            iter = 20000,
            warmup = 18000,
            thin=1,
            prior = priors,
            seed=1234,
		save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.99,
				   max_treedepth = 15))

#

save(fit_simple,file=here ("output","fit_simple.RData"))



# set formula (the same for benthos and fishes)
# complete model
formula1_fish <- brms::bf(mvbind (log(FRic_fish+1),
                                  log(Rao_fish+1)) ~ SR_fish + # fixed effects
                                                
                            (1|region), # random effects (intercept)
                      nl = F)


# complete model
formula1_algae <- brms::bf(mvbind (log(FRic_algae+1),
                                    log(Rao_algae+1)) ~ SR_algae + # fixed effects
                                                                    (1|region), # random effect (intercept) 
                          nl = F) 
# complete model
formula1_corals <- brms::bf(mvbind (log(FRic_corals+1),
                                     log(Rao_corals+1)) ~ SR_corals + # fixed effects
                                                              (1|region), # random effect (intercept) 
                             nl = F) 


# fit the model
fit_SR <- brms::brm(mvbf(formula1_fish, formula1_algae,formula1_corals) + 
                    set_rescor(TRUE), 
            data = bind_fish_benthos, 
            chains = 3, cores = 3,
            iter = 20000,
            warmup = 18000,
            thin=1,
            prior = priors,
            seed=1234,
		save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.99,
				   max_treedepth = 15))

#

save(fit_SR,file=here ("output","fit_SR.RData"))




# model selection analysis 
load(here ("output","fit_complete.RData"))
load(here ("output","fit2.RData"))
load(here ("output","fit3.RData"))
load(here ("output","fit_simple.RData"))
load(here ("output","fit_SR.RData"))





# LOO model fit checking   
# run loo fit test
loo_test <- lapply (list (fit_complete, fit2, fit3,fit_simple,fit_SR),
	loo, moment_match=T,reloo=T)


# save fit test
save(loo_test,file = here ("output","loo_test.RData"))



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
sel_model <- list (fit_complete, fit2, fit3,fit_simple,fit_SR)[which(loo_sel == min(unlist(loo_sel)))]
   
# list of results 
res <- list (looic = tab_mod_sel,
             param= tab_mod_fit,
             best_model = sel_model
             )

## save
save ( res, 
       file=here ("output", 
                  "MCMC_selected_model.Rdata"))





# end



