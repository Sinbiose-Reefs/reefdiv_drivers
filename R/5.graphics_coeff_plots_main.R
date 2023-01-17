
## --------------------
#       FIGURES
# --------------------

## call packages
source("R/packages.R")
source("R/functions.R")

### -------------------------
#  FIG - supporting information
# correlations between variables
load(here("data", "modeling_data.RData"))
load(here ("output", "FD_results.RData"))



# standardize covariate data
site_covs$sst_std <- (site_covs$sst - mean(site_covs$sst))/sd(site_covs$sst)
site_covs$turbidity_std <- (site_covs$turbidity - mean(site_covs$turbidity))/sd(site_covs$turbidity)
site_covs$productivity_std <- (site_covs$productivity - mean(site_covs$productivity))/sd(site_covs$productivity)
site_covs$salinity_std <- (site_covs$salinity - mean(site_covs$salinity))/sd(site_covs$salinity)
site_covs$offshore_distance_std <- (site_covs$offshore_distance - mean(site_covs$offshore_distance))/sd(site_covs$offshore_distance)

# effort data
effort_dataframe$fish_effort_std <- (effort_dataframe$fish_effort - mean(effort_dataframe$fish_effort))/sd(effort_dataframe$fish_effort)
effort_dataframe$benthos_effort_std <- (effort_dataframe$benthos_effort - mean(effort_dataframe$benthos_effort))/sd(effort_dataframe$benthos_effort)



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


#correlation plot (Fig. S1)


# save
pdf(file=here("output",
              "figures", 
              "corrplot"),height=5,width=5)

# correlation between variables
corrplot (
  cor(bind_fish_benthos[,c("sst", "turbidity","productivity", "salinity", "offshore_distance",
                           "SR_algae", "SR_fish", "SR_corals")]),
  p.mat = cor.mtest(cor(bind_fish_benthos[,c("sst", "turbidity","productivity", "salinity", "offshore_distance",
                                             "SR_algae", "SR_fish", "SR_corals")]))$p
  
  
)

dev.off()








# ------------------------------------------------

# bayesian R2

load (here ("output", 
            "MCMC_selected_model.RData"))


br2<- bayes_R2(res$best_model[[1]], robust=T) # mean for each reponse
round(br2,2)
round(mean(br2 [,'Estimate']),2) # average across responses




# coefficients of drivers

# model selection
tab_models <-data.frame (res$looic,
         res$param,
         model = c ("mv(FRic, Rao)~SR+SST+KD490+salinity+region",
                    "mv(FRic, Rao)~SR+SST+KD490+region",
                    "mv(FRic, Rao)~SR+SST+region",
                    "mv(FRic, Rao)~SST+region"))

# order 
tab_models<-tab_models [order (tab_models$Estimate,decreasing=F),]



# the best ranked model
# descriptive statistics of fit
library("bayesplot")
library("ggplot2")
#library("rstanarm")   

# extract posterior probs
posterior <- as.array(res$best_model[[1]])
dimnames(posterior)

color_scheme_set("blue")
p_SST <- mcmc_intervals(posterior,  pars = c( # fish
  "b_logFRicfish1_sst_std",         
  "b_logRaofish1_sst_std", 
  
  # algae
  
  "b_logFRicalgae1_sst_std",         
  "b_logRaoalgae1_sst_std", 
  
  # corals
  
  "b_logFRiccorals1_sst_std",         
  "b_logRaocorals1_sst_std" 
),
prob = 0.8, # 80% intervals
prob_outer = 0.95, # 99%
point_est = "median",
point_size = 2.5)

# richness
p_SR <- mcmc_intervals(posterior, pars = c( # fish
  "b_logFRicfish1_SR_fish", 
  "b_logRaofish1_SR_fish",
  
  # algae
  
  "b_logFRicalgae1_SR_algae", 
  "b_logRaoalgae1_SR_algae",
  
  # corals
  
  "b_logFRiccorals1_SR_corals", 
  "b_logRaocorals1_SR_corals"),
  
  prob = 0.8, # 50% intervals
  prob_outer = 0.95, # 95%
  point_est = "median",
  point_size = 2.5)


# region
p_region <- mcmc_intervals(posterior, pars = c( # fish
  # fric
  "b_logFRicfish1_Intercept",
  "b_logFRicfish1_regionBrazilianCoastSouth", 
  "b_logFRicfish1_regionBrazilianOceanicIslands",
  
  # rao
  "b_logRaofish1_Intercept",
  "b_logRaofish1_regionBrazilianCoastSouth", 
  "b_logRaofish1_regionBrazilianOceanicIslands",
  
  
  # algae
  "b_logFRicalgae1_Intercept",
  "b_logFRicalgae1_regionBrazilianCoastSouth",
  "b_logFRicalgae1_regionBrazilianOceanicIslands",
  
  "b_logRaoalgae1_Intercept",
  "b_logRaoalgae1_regionBrazilianCoastSouth",
  "b_logRaoalgae1_regionBrazilianOceanicIslands",
  
  # corals
  "b_logFRiccorals1_Intercept",
  "b_logFRiccorals1_regionBrazilianCoastSouth",
  "b_logFRiccorals1_regionBrazilianOceanicIslands",
  "b_logRaocorals1_Intercept",
  "b_logRaocorals1_regionBrazilianCoastSouth",
  "b_logRaocorals1_regionBrazilianOceanicIslands"


),
  
  prob = 0.8, # 50% intervals
  prob_outer = 0.95, # 95%
  point_est = "median",
  point_size = 2.5) 


pdf (here("output", "figures", "fig3"),width=15,height=5,family="sans")

grid.arrange(p_SST,
             p_SR,
             p_region, 
             ncol=7,nrow=1,
             layout_matrix = rbind(c (1,1,2,2,3,3,3)))

dev.off()




# rhat
# fitting statistics
# help here : https://biol609.github.io/lectures/23c_brms_prediction.html
png(here ("output","figures", "fitting_statistics_rhat"),
    units="cm",res=300, width=15,height=20)

rhat_vals <- rhat(res$best_model[[1]])
mcmc_rhat_data(rhat_vals)
mcmc_rhat(rhat_vals) + theme_bw() + yaxis_text(size = 4)

dev.off()

# n eff samples

png(here ("output","figures", "fitting_statistics_neff"),
    units="cm",res=300, width=15,height=20)

neff_vals <- neff_ratio(res$best_model[[1]])
mcmc_neff_data(neff_vals)
mcmc_neff(neff_vals)  + theme_bw() + yaxis_text(size = 4)


dev.off()



# ==============================================================




# full posterior exceedance probabilities


# param probability
require(tidybayes)
vars_ext<-get_variables(res$best_model[[1]])[grep("sst_std",get_variables(res$best_model[[1]]))]

# select model
model <- res$best_model[[1]]

# sst fish
# fric
model %>%
  spread_draws(b_logFRicfish1_sst_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicfish1_sst_std),
                   prob_pos = sum(b_logFRicfish1_sst_std > 0) / n(),
                   prob_neg = sum(b_logFRicfish1_sst_std < 0) / n())
# rao
model %>%
  spread_draws(b_logRaofish1_sst_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaofish1_sst_std),
                   prob_pos = sum(b_logRaofish1_sst_std > 0) / n(),
                   prob_neg = sum(b_logRaofish1_sst_std < 0) / n())

# sst algae
# fric
model %>%
  spread_draws(b_logFRicalgae1_sst_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicalgae1_sst_std),
                   prob_pos = sum(b_logFRicalgae1_sst_std > 0) / n(),
                   prob_neg = sum(b_logFRicalgae1_sst_std < 0) / n())

# Rao
model %>%
  spread_draws(b_logRaoalgae1_sst_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaoalgae1_sst_std),
                   prob_pos = sum(b_logRaoalgae1_sst_std > 0) / n(),
                   prob_neg = sum(b_logRaoalgae1_sst_std < 0) / n())



# sst corals
# fric
model %>%
  spread_draws(b_logFRiccorals1_sst_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRiccorals1_sst_std),
                   prob_pos = sum(b_logFRiccorals1_sst_std > 0) / n(),
                   prob_neg = sum(b_logFRiccorals1_sst_std < 0) / n())

# Rao
model %>%
  spread_draws(b_logRaocorals1_sst_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaocorals1_sst_std),
                   prob_pos = sum(b_logRaocorals1_sst_std > 0) / n(),
                   prob_neg = sum(b_logRaocorals1_sst_std < 0) / n())




# ------------------------------------

# species richness

# param probability

vars_ext<-get_variables(res$best_model[[1]])[grep("SR_",get_variables(res$best_model[[1]]))]


# SR fish
# fric
model %>%
  spread_draws(b_logFRicfish1_SR_fish) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicfish1_SR_fish),
                   prob_pos = sum(b_logFRicfish1_SR_fish > 0) / n(),
                   prob_neg = sum(b_logFRicfish1_SR_fish < 0) / n())
# rao
model %>%
  spread_draws(b_logRaofish1_SR_fish) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaofish1_SR_fish),
                   prob_pos = sum(b_logRaofish1_SR_fish > 0) / n(),
                   prob_neg = sum(b_logRaofish1_SR_fish < 0) / n())

# SR algae
# fric
model %>%
  spread_draws(b_logFRicalgae1_SR_algae) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicalgae1_SR_algae),
                   prob_pos = sum(b_logFRicalgae1_SR_algae > 0) / n(),
                   prob_neg = sum(b_logFRicalgae1_SR_algae < 0) / n())

# Rao
model %>%
  spread_draws(b_logRaoalgae1_SR_algae) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaoalgae1_SR_algae),
                   prob_pos = sum(b_logRaoalgae1_SR_algae > 0) / n(),
                   prob_neg = sum(b_logRaoalgae1_SR_algae < 0) / n())



# SR corals
# fric
model %>%
  spread_draws(b_logFRiccorals1_SR_corals) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRiccorals1_SR_corals),
                   prob_pos = sum(b_logFRiccorals1_SR_corals > 0) / n(),
                   prob_neg = sum(b_logFRiccorals1_SR_corals < 0) / n())

# Rao
model %>%
  spread_draws(b_logRaocorals1_SR_corals) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaocorals1_SR_corals),
                   prob_pos = sum(b_logRaocorals1_SR_corals > 0) / n(),
                   prob_neg = sum(b_logRaocorals1_SR_corals < 0) / n())






# ------------------------------------

# region

# param probability
vars_ext<-get_variables(res$best_model[[1]])[grep("region",get_variables(res$best_model[[1]]))]
intercepts <- fixef(res$best_model[[1]])[grep("Intercept", rownames(fixef(res$best_model[[1]]))),"Estimate"]

# fric
# south
model %>%
  spread_draws(b_logFRicfish1_regionBrazilianCoastSouth) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicfish1_regionBrazilianCoastSouth),
                   prob_pos_SE.S = sum(intercepts[1]+b_logFRicfish1_regionBrazilianCoastSouth > intercepts[1]) / n(),
                   prob_neg_SE.S = sum(intercepts[1]+b_logFRicfish1_regionBrazilianCoastSouth < intercepts[1]) / n()
                   
                   )
# islands
model %>%
  spread_draws(b_logFRicfish1_regionBrazilianOceanicIslands) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicfish1_regionBrazilianOceanicIslands),
                   prob_pos_OI = sum(intercepts[1]+b_logFRicfish1_regionBrazilianOceanicIslands > intercepts[1]) / n(),
                   prob_neg_OI = sum(intercepts[1]+b_logFRicfish1_regionBrazilianOceanicIslands < intercepts[1]) / n()
  )

# rao
# south
model %>%
  spread_draws(b_logRaofish1_regionBrazilianCoastSouth) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaofish1_regionBrazilianCoastSouth),
                   prob_pos_SE.S = sum(intercepts[2]+b_logRaofish1_regionBrazilianCoastSouth > intercepts[2]) / n(),
                   prob_neg_SE.S = sum(intercepts[2]+b_logRaofish1_regionBrazilianCoastSouth < intercepts[2]) / n()
                   
                   )

# islands
model %>%
  spread_draws(b_logRaofish1_regionBrazilianOceanicIslands) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaofish1_regionBrazilianOceanicIslands),
                   prob_pos_SE.S = sum(intercepts[2]+b_logRaofish1_regionBrazilianOceanicIslands > intercepts[2]) / n(),
                   prob_neg_SE.S = sum(intercepts[2]+b_logRaofish1_regionBrazilianOceanicIslands < intercepts[2]) / n()
                   
  )


# ----------------------------------------------------------
# algae

# fric
# south
model %>%
  spread_draws(b_logFRicalgae1_regionBrazilianCoastSouth) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicalgae1_regionBrazilianCoastSouth),
                   prob_pos_SE.S = sum(intercepts[3]+b_logFRicalgae1_regionBrazilianCoastSouth > intercepts[3]) / n(),
                   prob_neg_SE.S = sum(intercepts[3]+b_logFRicalgae1_regionBrazilianCoastSouth < intercepts[3]) / n()
                   
  )

# islands
model %>%
  spread_draws(b_logFRicalgae1_regionBrazilianOceanicIslands) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicalgae1_regionBrazilianOceanicIslands),
                   prob_pos_OI = sum(intercepts[3]+b_logFRicalgae1_regionBrazilianOceanicIslands > intercepts[3]) / n(),
                   prob_neg_OI = sum(intercepts[3]+b_logFRicalgae1_regionBrazilianOceanicIslands < intercepts[3]) / n()
                   
  )

# rao
# south
model %>%
  spread_draws(b_logRaoalgae1_regionBrazilianCoastSouth) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaoalgae1_regionBrazilianCoastSouth),
                   prob_pos_SE.S = sum(intercepts[4]+b_logRaoalgae1_regionBrazilianCoastSouth > intercepts[4]) / n(),
                   prob_neg_SE.S = sum(intercepts[4]+b_logRaoalgae1_regionBrazilianCoastSouth < intercepts[4]) / n()
                   
  )

# islands
model %>%
  spread_draws(b_logRaoalgae1_regionBrazilianOceanicIslands) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaoalgae1_regionBrazilianOceanicIslands),
                   prob_pos_SE.S = sum(intercepts[4]+b_logRaoalgae1_regionBrazilianOceanicIslands > intercepts[4]) / n(),
                   prob_neg_SE.S = sum(intercepts[4]+b_logRaoalgae1_regionBrazilianOceanicIslands < intercepts[4]) / n()
                   
  )



# ----------------------------------------------------------
# corals

# fric
# south
model %>%
  spread_draws(b_logFRiccorals1_regionBrazilianCoastSouth) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRiccorals1_regionBrazilianCoastSouth),
                   prob_pos_SE.S = sum(intercepts[5]+b_logFRiccorals1_regionBrazilianCoastSouth > intercepts[5]) / n(),
                   prob_neg_SE.S = sum(intercepts[5]+b_logFRiccorals1_regionBrazilianCoastSouth < intercepts[5]) / n()
                   
  )
# islands
model %>%
  spread_draws(b_logFRiccorals1_regionBrazilianOceanicIslands) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRiccorals1_regionBrazilianOceanicIslands),
                   prob_pos_OI = sum(intercepts[5]+b_logFRiccorals1_regionBrazilianOceanicIslands > intercepts[5]) / n(),
                   prob_neg_OI = sum(intercepts[5]+b_logFRiccorals1_regionBrazilianOceanicIslands < intercepts[5]) / n()
                   
  )

# rao
# south
model %>%
  spread_draws(b_logRaocorals1_regionBrazilianCoastSouth) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaocorals1_regionBrazilianCoastSouth),
                   prob_pos_SE.S = sum(intercepts[6]+b_logRaocorals1_regionBrazilianCoastSouth > intercepts[6]) / n(),
                   prob_neg_SE.S = sum(intercepts[6]+b_logRaocorals1_regionBrazilianCoastSouth < intercepts[6]) / n()
                   
  )

# islands
model %>%
  spread_draws(b_logRaocorals1_regionBrazilianOceanicIslands) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaocorals1_regionBrazilianOceanicIslands),
                   prob_pos_SE.S = sum(intercepts[6]+b_logRaocorals1_regionBrazilianOceanicIslands > intercepts[6]) / n(),
                   prob_neg_SE.S = sum(intercepts[6]+b_logRaocorals1_regionBrazilianOceanicIslands < intercepts[6]) / n()
                   
  )


