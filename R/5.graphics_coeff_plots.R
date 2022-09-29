
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


br2<- bayes_R2(res$best_model[[1]], robust=T) # mean for each reponse
round(br2,2)
round(mean(br2 [,'Estimate']),2) # average across responses




# coefficients of drivers

load (here ("output", 
            "MCMC_selected_model.RData"))

# model selection
tab_models <-data.frame (res$looic,
         res$param,
         model = c ("mv(FRic, Rao)~SR+SST+KD490+salinity+(1|region)",
                    "mv(FRic, Rao)~SR+SST+KD490+(1|region)",
                    "mv(FRic, Rao)~SR+SST+(1|region)",
                    "mv(FRic, Rao)~SST+(1|region)",
                    "mv(FRic, Rao)~SR+(1|region)"))

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
prob = 0.5, # 50% intervals
prob_outer = 0.95, # 99%
point_est = "median",
point_size = 2.5)

# turbidity
p_turbidity <- mcmc_intervals(posterior, pars = c( # fish
  "b_logFRicfish1_turbidity_std", 
  "b_logRaofish1_turbidity_std",
  
  # algae
  
  "b_logFRicalgae1_turbidity_std", 
  "b_logRaoalgae1_turbidity_std",
  
  # corals
  
  "b_logFRiccorals1_turbidity_std", 
  "b_logRaocorals1_turbidity_std"),
  
  prob = 0.5, # 50% intervals
  prob_outer = 0.95, # 95%
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
  
  prob = 0.5, # 50% intervals
  prob_outer = 0.95, # 95%
  point_est = "median",
  point_size = 2.5)


pdf (here("output", "figures", "fig2"),width=7,height=3)

grid.arrange(p_SST+theme (axis.text.y = element_blank()),
             p_turbidity+theme (axis.text.y = element_blank()),
             p_SR+theme (axis.text.y = element_blank()), 
             ncol=3)

dev.off()


# add 80% intervals, manually in inkscape

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

# turbidity
p_turbidity <- mcmc_intervals(posterior, pars = c( # fish
  "b_logFRicfish1_turbidity_std", 
  "b_logRaofish1_turbidity_std",
  
  # algae
  
  "b_logFRicalgae1_turbidity_std", 
  "b_logRaoalgae1_turbidity_std",
  
  # corals
  
  "b_logFRiccorals1_turbidity_std", 
  "b_logRaocorals1_turbidity_std"),
  
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, # 95%
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
  
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, # 95%
  point_est = "median",
  point_size = 2.5)


pdf (here("output", "figures", "fig2_80p"),width=7,height=3)

grid.arrange(p_SST+theme (axis.text.y = element_blank()),
             p_turbidity+theme (axis.text.y = element_blank()),
             p_SR+theme (axis.text.y = element_blank()), 
             ncol=3)

dev.off()



# rhat
# fitting statistics
# help here : https://biol609.github.io/lectures/23c_brms_prediction.html
pdf(here ("output","figures", "fitting_statistics"),onefile = T)

rhat_vals <- rhat(res$best_model[[1]])
mcmc_rhat_data(rhat_vals)
mcmc_rhat(rhat_vals) + theme_bw() + yaxis_text(size = 4)


# n eff samples

neff_vals <- neff_ratio(res$best_model[[1]])
mcmc_neff_data(neff_vals)
mcmc_neff(neff_vals)  + theme_bw() + yaxis_text(size = 4)


dev.off()



# ==============================================================




# full posterior excedance probabilities



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

# turbidity



# param probability

vars_ext<-get_variables(res$best_model[[1]])[grep("turbidity_std",get_variables(res$best_model[[1]]))]


# turbidity fish
# fric
model %>%
  spread_draws(b_logFRicfish1_turbidity_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicfish1_turbidity_std),
                   prob_pos = sum(b_logFRicfish1_turbidity_std > 0) / n(),
                   prob_neg = sum(b_logFRicfish1_turbidity_std < 0) / n())
# rao
model %>%
  spread_draws(b_logRaofish1_turbidity_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaofish1_turbidity_std),
                   prob_pos = sum(b_logRaofish1_turbidity_std > 0) / n(),
                   prob_neg = sum(b_logRaofish1_turbidity_std < 0) / n())

# turbidity algae
# fric
model %>%
  spread_draws(b_logFRicalgae1_turbidity_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicalgae1_turbidity_std),
                   prob_pos = sum(b_logFRicalgae1_turbidity_std > 0) / n(),
                   prob_neg = sum(b_logFRicalgae1_turbidity_std < 0) / n())

# Rao
model %>%
  spread_draws(b_logRaoalgae1_turbidity_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaoalgae1_turbidity_std),
                   prob_pos = sum(b_logRaoalgae1_turbidity_std > 0) / n(),
                   prob_neg = sum(b_logRaoalgae1_turbidity_std < 0) / n())



# turbidity corals
# fric
model %>%
  spread_draws(b_logFRiccorals1_turbidity_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRiccorals1_turbidity_std),
                   prob_pos = sum(b_logFRiccorals1_turbidity_std > 0) / n(),
                   prob_neg = sum(b_logFRiccorals1_turbidity_std < 0) / n())

# Rao
model %>%
  spread_draws(b_logRaocorals1_turbidity_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaocorals1_turbidity_std),
                   prob_pos = sum(b_logRaocorals1_turbidity_std > 0) / n(),
                   prob_neg = sum(b_logRaocorals1_turbidity_std < 0) / n())




# ------------------------------------

# turbidity



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



#  random effects

# explore fixed and random effects
# help here
# https://biol609.github.io/lectures/23c_brms_prediction.html


library(brms)
library(tidyverse)



nd <- data.frame(
  SR_fish = 0, sst_std = 0, turbidity_std = 0, SR_algae = 0, SR_corals = 0,
  region = sort(unique(model$data$region))
)

# The rows are HMC draws, the columns correspond to the rows in `nd`
#   re_formula = NULL for incorporation of random effects.
int_preds <- brms::posterior_predict(model, newdata = nd, re_formula = NULL)
# The output is a 3D array, with each element on the third dimension
#  corresponding to a different response variable from the multivariate normal
#  model.

head(int_preds)
# The columns correspond to sort(unique(mod$data$region)), so first columns is
#   "BrazilianCoastNorth", second is "BrazilianCoastSouth" and third is
#   "BrazilianOceanicIslands".
# So now for each response all you have to do is calculate among-column
#   differences row-wise, and that will give you the difference among regions.
# Transform array to list, makes it easier to deal with output
resp_names <- dimnames(int_preds)[[3]]
int_preds <- seq_len(dim(int_preds)[3]) %>%
  lapply(function(x, marray) marray[ , , x], marray = int_preds)
names(int_preds) <- resp_names

# NB: 3 levels for a random effect is not recommended at all. You can barely
#   estimate a mean with 3 levels, let alone a variance.
diff_tab <- purrr::map_dfr(int_preds, function(x) {
  y <- cbind(x[, 1] - x[, 2], x[, 1] - x[, 3], x[, 2] - x[, 3])
  # exceedance probabilities.
  #   1 - What is the probability of North being greater than south?
  #   2 - What is the probability of North being greater than Oc. islands?
  #   3 - What is the probability of South being greater than Oc. islands?
  ex_probs <- c(sum(y[, 1] > 0), sum(y[, 2] > 0), sum(y[, 3] > 0)) / nrow(y)
  y %>%
    ggdist::median_hdci() %>%
    dplyr::mutate(diff = c("bn-bs", "bn-oi", "bs-oi"), ex_prob = ex_probs)
}, .id = "Response")

# Just one case where evidence for difference is > 0.9
# functional richness of algae is substantially higher in the south when
#   compared to the oceanic islands.
diff_tab %>%
  dplyr::filter(ex_prob > 0.9)










# aggregating to get averages and plot a vertical line into the histogram
#aggregate (cor_out_sel,by=list(cor_out_sel$cor_pair),
#           FUN=mean)

vline.data <- cor_out %>%
  group_by(cor_pair) %>%
  summarize(z = mean(correlation_r ))

# plot posterior distribution of each pair as histogram
post_cors_plot <- ggplot(data = cor_out) +
  geom_histogram(mapping = aes(x = correlation_r), fill = "dodgerblue2",
                 colour = "grey30",
                 bins=60) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
  labs(x = "Pearson correlation", y = "Frequency") +
  facet_wrap(~cor_pair) +
  geom_vline(aes(xintercept = z), vline.data, 
             colour = "red",size=1.5)+ # certical line with average
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))
post_cors_plot

ggsave("post_cors_plot.pdf", post_cors_plot, width = 11.88, height = 8.81)



# plot observed vs. mean predicted correlation
post_cor_means <- cor_out %>%
  dplyr::group_by(cor_pair) %>%
  dplyr::summarise(mean_post_r = ggdist::median_hdci(correlation_r))















summary(model)$rescor_pars

# check original vs residual correlation
merged_cor <- merge(post_cor_means, cor_orig[, -c(1:2)], by = "cor_pair")

# bind residual correlations
require(reshape)
corr_res <- VarCorr(res$best_model[[1]] )$residual__$cor # residuals coor
corr_res <- melt(corr_res[,1,]) # melt
corr_res <- corr_res[which(corr_res$value != 1),]# rm corr == 1

# edit labels
# var1
#corr_res$Var1<-gsub ("_","",corr_res$Var1)
corr_res[,1]<-gsub ("log","",corr_res[,1])
corr_res[,2]<-gsub ("log","",corr_res[,2])
corr_res[,1] <- gsub ("1", "",corr_res[,1])
corr_res[,2] <- gsub ("1", "",corr_res[,2])

# bind
corr_res$var3 <- paste (corr_res[,1], corr_res[,2], sep= " ~ ")
corr_res <- corr_res[which(corr_res$var3 %in% merged_cor$cor_pair),]

# bind resid corr
merged_cor$corr_res <- corr_res$value[match (merged_cor$cor_pair,corr_res$var3)] #== merged_cor$cor_pair

# plot
scatter_cor <- ggplot(data = merged_cor) +
  geom_point(mapping = aes(y = obs_cor_r, 
                           x = mean_post_r), shape = 21,
             size = 3, fill = "grey90") +
  theme_classic() +
  geom_abline(slope = 1, linetype = 3) +
  labs(y = "Observed Pearson correlation",
       x = "Posterior mean of predicted Pearson's correlation") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) +
  xlim(-0.6, 1) + ylim(-0.6, 1)

require (ggrepel)
scatter_cor + geom_text_repel(data = merged_cor, aes (y = obs_cor_r, x = mean_post_r,
                                                      label = cor_pair),
                              box.padding = 0.6, 
                              max.overlaps = Inf,
                              size=3)

ggsave(here ("output_sensitivity","scatter_cor.pdf"), 
       scatter_cor, width = 6.5, height = 6)


# resid corr
# plot
scatter_cor_res <- ggplot(data = merged_cor) +
  geom_point(mapping = aes(y = obs_cor_r, 
                           x = corr_res), 
             shape = 21,
             size = 3, fill = "grey90") +
  theme_classic() +
  geom_abline(slope = 1, linetype = 3) +
  labs(y = "Observed Pearson correlation",
       x = "Posterior mean of Pearson's residual correlation") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) +
  xlim(-0.6, 1) + ylim(-0.6, 1)

require (ggrepel)
scatter_cor_res + geom_text_repel(data = merged_cor, aes (y = obs_cor_r, x = corr_res,
                                                          label = cor_pair),
                                  box.padding = 0.6, 
                                  max.overlaps = Inf,
                                  size=3)


# predicted vs residual
scatter_cor_res <- ggplot(data = merged_cor) +
  geom_point(mapping = aes(y = mean_post_r, 
                           x = corr_res), 
             shape = 21,
             size = 3, fill = "grey90") +
  theme_classic() +
  geom_abline(slope = 1, linetype = 3) +
  labs(y = "Posterior mean of predicted Pearson's correlation",
       x = "Posterior mean of Pearson's residual correlation") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) +
  xlim(-0.6, 1) + ylim(-0.6, 1)

require (ggrepel)
scatter_cor_res + geom_text_repel(data = merged_cor, aes (y = mean_post_r, 
                                                          x = corr_res,
                                                          label = cor_pair),
                                  box.padding = 0.6, 
                                  max.overlaps = Inf,
                                  size=3)

# absolute diff

merged_cor<- cbind(merged_cor, 
                   diff = (abs(merged_cor$mean_post_r) - abs(merged_cor$corr_res)))
merged_cor<- cbind(merged_cor, 
                   diff2 = (abs(merged_cor$mean_post_r) - abs(merged_cor$obs_cor_r)))

merged_cor<-melt(merged_cor)
colnames(merged_cor)[which(colnames(merged_cor) == "variable")] <- "Correlation"
levels(merged_cor$Correlation)[which(levels(merged_cor$Correlation) == "obs_cor_r")] <- "Observed"
levels(merged_cor$Correlation)[which(levels(merged_cor$Correlation) == "corr_res")] <- "Residual"
levels(merged_cor$Correlation)[which(levels(merged_cor$Correlation) == "mean_post_r")] <- "Predicted"


# averages
merged_cor %>% group_by(Correlation) %>%
  summarise(av=mean(value),
            ra=range(value))


# organize names
# taxon
merged_cor$taxon <- gsub ("FRic", "",merged_cor$cor_pair)
merged_cor$taxon <- gsub ("Rao", "",merged_cor$taxon)
merged_cor$taxon <- firstup(merged_cor$taxon)

# order 
merged_cor$taxon<-factor (merged_cor$taxon, 
                             levels = c("Algae ~ algae",
                                        "Corals ~ corals",
                                        "Fish ~ fish",
                                        "Algae ~ corals",
                                        "Fish ~ algae",
                                        "Fish ~ corals"))

# index
merged_cor$index <- gsub ("fish", "",merged_cor$cor_pair)
merged_cor$index <- gsub ("algae", "",merged_cor$index)
merged_cor$index <- gsub ("corals", "",merged_cor$index)
merged_cor$index <- firstup(merged_cor$index)

# lollipop plot (plot of difference in correlation)

# save
pdf(file=here("output",
              "figures", 
              "diff_cor"),height=5,width=10)


# Change baseline
ggplot(merged_cor[which(merged_cor$Correlation %in% c("Observed", 
                                                      "Residual", 
                                                      "Predicted")),], 
       aes(x=value, 
           y=index,
           color = Correlation,
           fill = Correlation,
           group = cor_pair)) +
  geom_segment( aes(x=value, xend=0, 
                    y=reorder(index,  value), 
                    yend=reorder(index,  value)), 
                color="grey",
                position = position_jitter(height = 0.2, width = 0)) +
  
  geom_point(size=3,
             position = position_jitter(height = 0.2, width = 0)) + 

    scale_color_viridis_d()+
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Pearson's correlation") + 
  ylab ("Pairs of functional metrics")+
  geom_vline(aes(xintercept =0),alpha =0.5,size=1,col = "gray") +
  xlim(c(-.5,.65)) +
  theme_classic() + 
   facet_wrap(~taxon,ncol = 3,scales="free")+
  theme(legend.position = "top") 




dev.off()
