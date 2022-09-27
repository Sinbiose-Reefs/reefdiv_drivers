
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


# add 75% intervals, manually in inkscape

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




# plot
# here: https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html
pdf(here ("output","figures","fig2_temp_50p"),width=5,height=3)

color_scheme_set("blue")
mcmc_areas(posterior, pars = c( # fish
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
prob_outer = 0.95, # 95%
point_est = "median",
area_method = "equal height")

dev.off()


# plot
# here: https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html
pdf(here ("output","figures","fig2_turbidity_50p"),width=5,height=3)

color_scheme_set("blue")
mcmc_areas(posterior, pars = c( # fish
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
  area_method = "equal height") 

dev.off()


# plot
# here: https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html
pdf(here ("output","figures","fig2_SR_50p"),width=5,height=3)

color_scheme_set("blue")
mcmc_areas(posterior, pars = c( # fish
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
  area_method = "equal height")

dev.off()


# rhat
# fitting statistics
pdf(here ("output","figures", "fitting_statistics"),onefile = T)

rhat_vals <- rhat(res$best_model[[1]])
mcmc_rhat_data(rhat_vals)
mcmc_rhat(rhat_vals) + theme_bw() + yaxis_text(size = 4)


# n eff samples

neff_vals <- neff_ratio(res$best_model[[1]])
mcmc_neff_data(neff_vals)
mcmc_neff(neff_vals)  + theme_bw() + yaxis_text(size = 4)


dev.off()


# help

# explore fixed and random effects
# help here
# https://biol609.github.io/lectures/23c_brms_prediction.html


library(brms)
library(tidyverse)
load("rdata_/MCMC_selected_model.rdata")
length(res$best_model)
mod <- res$best_model[[1]]
nd <- data.frame(
  SR_fish = 0, sst_std = 0, turbidity_std = 0, SR_algae = 0, SR_corals = 0,
  region = sort(unique(mod$data$region))
)

# The rows are HMC draws, the columns correspond to the rows in `nd`
#   re_formula = NULL for incorporation of random effects.
int_preds <- brms::posterior_predict(mod, newdata = nd, re_formula = NULL)
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




# FIXED EFFECTS
require(tidybayes)

# 
fe_only <- tibble(sst_std = seq (range(res$best_model[[1]]$data$sst_std)[1],
                                 range(res$best_model[[1]]$data$sst_std)[2],
                                 0.01),
                    turbidity_std = 0,
                    SR_corals = 0,
                    SR_algae = 0,
                    SR_fish =0) %>%
  
    add_fitted_draws(res$best_model[[1]],
                     re_formula = NA,
                     scale = "response",
                     n = 1e2)
  
# summary
fe_only_mean <- fe_only %>% 
    group_by(sst_std, .category) %>%
    summarize(.value = median(.value))
  


# sst
plotSST<- ggplot(fe_only,
       aes(x = sst_std, y = .value)) +
  stat_interval(alpha = 0.1) + 
  geom_line(data = fe_only_mean, 
            color = "black", 
            lwd = 2) +
  facet_wrap (~.category, scales="free",ncol=6) +
  xlab ("SST")+theme_classic()+
  theme(legend.position = "none")


# turb
fe_only <- tibble(sst_std = 0,
                  turbidity_std = seq (range(res$best_model[[1]]$data$turbidity_std)[1],
                                       range(res$best_model[[1]]$data$turbidity_std)[2],
                                       0.01),
                  SR_corals = 0,
                  SR_algae = 0,
                  SR_fish =0) %>%
  
  add_fitted_draws(res$best_model[[1]],
                   re_formula = NA,
                   scale = "response",
                   n = 1e2)

# summary
fe_only_mean <- fe_only %>% 
  group_by(turbidity_std, .category) %>%
  summarize(.value = median(.value))


# kd490
plotKD<- ggplot(fe_only,
                 aes(x = turbidity_std, y = .value)) +
  stat_interval(alpha = 0.1) + 
  geom_line(data = fe_only_mean, 
            color = "black", 
            lwd = 2) +
  facet_wrap (~.category, scales="free",ncol=6) +
  theme(legend.position = "none")




# SR corals
fe_only <- tibble(sst_std = 0,
                  turbidity_std = 0,
                  SR_corals = seq (range(res$best_model[[1]]$data$SR_corals)[1],
                                   range(res$best_model[[1]]$data$SR_corals)[2],
                                   0.01),
                  SR_algae = 0,
                  SR_fish =0) %>%
  
  add_fitted_draws(res$best_model[[1]],
                   re_formula = NA,
                   scale = "response",
                   n = 1e2)

# summary
fe_only_mean <- fe_only %>% 
  group_by(SR_corals, .category) %>%
  summarize(.value = median(.value))




# SR corals
plotSR_corals<- ggplot(fe_only,
                aes(x = SR_corals, y = .value)) +
  stat_interval(alpha = 0.1) + 
  scale_fill_grey()+
  geom_line(data = fe_only_mean, 
            color = "black", 
            lwd = 2) +
  facet_wrap (~.category, scales="free",ncol=6) +
  theme(legend.position = "none")



# SR algae
fe_only <- tibble(sst_std = 0,
                  turbidity_std = 0,
                  SR_corals = 0,
                  SR_algae = seq (range(res$best_model[[1]]$data$SR_algae)[1],
                                  range(res$best_model[[1]]$data$SR_algae)[2],
                                  0.01),
                  SR_fish =0) %>%
  
  add_fitted_draws(res$best_model[[1]],
                   re_formula = NA,
                   scale = "response",
                   n = 1e2)

# summary
fe_only_mean <- fe_only %>% 
  group_by(SR_algae, .category) %>%
  summarize(.value = median(.value))


# SR algae
plotSR_algae<- ggplot(fe_only,
                       aes(x = SR_algae, y = .value)) +
  stat_interval(alpha = 0.1) + 
  scale_fill_grey()+
  geom_line(data = fe_only_mean, 
            color = "black", 
            lwd = 2) +
  facet_wrap (~.category, scales="free",ncol=6) +
  theme(legend.position = "none")





# SR fish
fe_only <- tibble(sst_std = 0,
                  turbidity_std = 0,
                  SR_corals = 0,
                  SR_algae = 0,
                  SR_fish =seq (range(res$best_model[[1]]$data$SR_fish)[1],
                                range(res$best_model[[1]]$data$SR_fish)[2],
                                0.01)) %>%
  
  add_fitted_draws(res$best_model[[1]],
                   re_formula = NA,
                   scale = "response",
                   n = 1e2)

# summary
fe_only_mean <- fe_only %>% 
  group_by(SR_fish, .category) %>%
  summarize(.value = median(.value))


# SR fish
plotSR_fish<- ggplot(fe_only,
                      aes(x = SR_fish, y = .value)) +
  stat_interval(alpha = 0.1) + 
  scale_fill_grey()+
  geom_line(data = fe_only_mean, 
            color = "black", 
            lwd = 2) +
  facet_wrap (~.category, scales="free",ncol=6) +
  theme(legend.position = "none")




#   
pdf (here ("output", "figures", "fig4"),
     height = 2,width=10, onefile = T)

plotSST+xlab ("SST")+theme_classic()+theme(legend.position = "none")
plotKD+xlab ("Turbidity")+theme_classic()+theme(legend.position = "none")
plotSR_algae+xlab ("Algae richness")+theme_classic()+theme(legend.position = "none")
plotSR_corals+xlab ("Coral richness")+theme_classic()+theme(legend.position = "none")
plotSR_fish+xlab ("Fish richness")+theme_classic()+theme(legend.position = "none")


dev.off()


# =================================================


# correlations & congruence inference


library(terra)
library(flexmix)
library(modeltools)
library(tidyverse)
library(brms)
library(ggdist)
require(here)
require(tidybayes)


# bayesian R2
br2<- bayes_R2(res$best_model[[1]]) # mean for each reponse
round(br2,2)
round(mean(br2 [,'Estimate']),2) # average across responses

# predicted values
all_preds <- exp(brms::posterior_predict(res$best_model[[1]]))
cor_out <- vector(mode = "list", length = nrow(all_preds))
for (i in seq_len(nrow(all_preds))) {
  cor_i <- cor(as.matrix(all_preds[i, , ]))
  ind <- which(upper.tri(cor_i, diag = TRUE), arr.ind = TRUE)
  nn <- dimnames(cor_i)
  cor_out[[i]] <- data.frame(var_a = nn[[1]][ind[, 1]],
                             var_b = nn[[2]][ind[, 2]],
                             correlation_r = cor_i[ind]) %>%
    dplyr::filter(var_a != var_b) %>%
    dplyr::mutate(.draw = i)
}

cor_out <- do.call("rbind.data.frame", cor_out) %>%
  dplyr::mutate(dplyr::across(tidyselect::starts_with("var_"),
                              ~gsub("^log", "", .x))) %>%
  dplyr::mutate(var_a = case_when(!grepl("benthos", var_a) ~ paste0(var_a, ""), TRUE ~ var_a),
                var_b = case_when(!grepl("benthos", var_b) ~ paste0(var_b, ""), TRUE ~ var_b),
                cor_pair = paste0(var_a, " ~ ", var_b))

# select correlation in the diagonal (i.e., cor of the same metric between fish and benthos)
sel_cors <- c("FRicfish1_fish ~ FRicalgae1_fish",
              "Raofish1_fish ~ Raoalgae1_fish")

# selecting interesting correlations
#cor_out_sel <- cor_out[which(cor_out$cor_pair %in% sel_cors),]
(cor_out %>%
    dplyr::group_by(cor_pair) %>%
    dplyr::summarise(ggdist::median_hdci(correlation_r),
                     prob_pos = sum(correlation_r > 0) / n(),
                     prob_neg = sum(correlation_r < 0) / n())
)

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

# summary table with median +/- 95% credible intervals
# probability of correlation being negative or positive.
(cor_out %>%
    dplyr::group_by(cor_pair) %>%
    dplyr::summarise(ggdist::median_hdci(correlation_r),
                     prob_pos = sum(correlation_r > 0) / n(),
                     prob_neg = sum(correlation_r < 0) / n())
)

# plot observed vs. mean predicted correlation
post_cor_means <- cor_out %>%
  dplyr::group_by(cor_pair) %>%
  dplyr::summarise(mean_post_r = mean(correlation_r))




# correlation from the model

cor_orig <- res$best_model[[1]]$data %>%
  dplyr::select(-sst_std,
                -turbidity_std,
                -region,
                -SR_fish,
                -SR_algae,
                -SR_corals) %>%
  dplyr::rename_with(~gsub("^log", "", .x)) %>%
  exp %>%
  cor

ind <- which(upper.tri(cor_orig, diag = TRUE), arr.ind = TRUE)
nn <- dimnames(cor_orig)
cor_orig <- data.frame(var_a = nn[[1]][ind[, 1]], var_b = nn[[2]][ind[, 2]],
                       obs_cor_r = cor_orig[ind]) %>%
  dplyr::filter(var_a != var_b) %>%
  dplyr::mutate(var_a = case_when(!grepl("benthos", var_a) ~ paste0(var_a, ""), TRUE ~ gsub("_", "", var_a)),
                var_b = case_when(!grepl("benthos", var_b) ~ paste0(var_b, ""), TRUE ~ gsub("_", "", var_b)),
                cor_pair = paste0(var_a, " ~ ", var_b))
post_cor_means$cor_pair <- gsub ("1","", post_cor_means$cor_pair)
cor_orig$cor_pair <- gsub ("_","",cor_orig$cor_pair)
mean(cor_orig$obs_cor_r)


# check original vs residual correlation
merged_cor <- merge(post_cor_means, cor_orig[, -c(1:2)], by = "cor_pair")

# bind residual correlations
require(reshape)
corr_res <- VarCorr(res$best_model[[1]] )$residual__$cor # residals coor
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

# ==============================================================
# param probability

vars_ext<-get_variables(res$best_model[[1]])[grep("sst",get_variables(res$best_model[[1]]))]


# sst SR fish
res$best_model[[1]] %>%
  spread_draws(b_logRaoalgae1_sst_std) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logRaoalgae1_sst_std),
                   prob_pos = sum(b_logRaoalgae1_sst_std > 0) / n(),
                   prob_neg = sum(b_logRaoalgae1_sst_std < 0) / n())

# sst SR benthos
model %>%
  spread_draws(b_logSRbenthos_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSRbenthos_sst),
                   prob_pos = sum(b_logSRbenthos_sst > 0) / n(),
                   prob_neg = sum(b_logSRbenthos_sst < 0) / n())

# sst Fric fish
model %>%
  spread_draws(b_logFRic_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRic_sst),
                   prob_pos = sum(b_logFRic_sst > 0) / n(),
                   prob_neg = sum(b_logFRic_sst < 0) / n())

# sst fric benthos
model %>%
  spread_draws(b_logFRicbenthos_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicbenthos_sst),
                   prob_pos = sum(b_logFRicbenthos_sst > 0) / n(),
                   prob_neg = sum(b_logFRicbenthos_sst < 0) / n())

# sst FEve fish
model %>%
  spread_draws(b_logFEve_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEve_sst),
                   prob_pos = sum(b_logFEve_sst > 0) / n(),
                   prob_neg = sum(b_logFEve_sst < 0) / n())

# sst FEve benthos
model %>%
  spread_draws(b_logFEvebenthos_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEvebenthos_sst),
                   prob_pos = sum(b_logFEvebenthos_sst > 0) / n(),
                   prob_neg = sum(b_logFEvebenthos_sst < 0) / n())

# sst Fdiv fish
model %>%
  spread_draws(b_logFDiv_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDiv_sst),
                   prob_pos = sum(b_logFDiv_sst > 0) / n(),
                   prob_neg = sum(b_logFDiv_sst < 0) / n())

# sst fdiv benthos
model %>%
  spread_draws(b_logFDivbenthos_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDivbenthos_sst),
                   prob_pos = sum(b_logFDivbenthos_sst > 0) / n(),
                   prob_neg = sum(b_logFDivbenthos_sst < 0) / n())


# reef area
vars_ext<-get_variables(model)[grep("area",get_variables(model))]

# sst SR fish
model %>%
  spread_draws(b_logSR_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSR_area),
                   prob_pos = sum(b_logSR_area > 0) / n(),
                   prob_neg = sum(b_logSR_area < 0) / n())

# sst SR benthos
model %>%
  spread_draws(b_logSRbenthos_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSRbenthos_area),
                   prob_pos = sum(b_logSRbenthos_area > 0) / n(),
                   prob_neg = sum(b_logSRbenthos_area < 0) / n())

# sst Fric fish
model %>%
  spread_draws(b_logFRic_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRic_area),
                   prob_pos = sum(b_logFRic_area > 0) / n(),
                   prob_neg = sum(b_logFRic_area < 0) / n())

# sst fric benthos
model %>%
  spread_draws(b_logFRicbenthos_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicbenthos_area),
                   prob_pos = sum(b_logFRicbenthos_area > 0) / n(),
                   prob_neg = sum(b_logFRicbenthos_area < 0) / n())

# sst FEve fish
model %>%
  spread_draws(b_logFEve_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEve_area),
                   prob_pos = sum(b_logFEve_area > 0) / n(),
                   prob_neg = sum(b_logFEve_area < 0) / n())

# sst FEve benthos
model %>%
  spread_draws(b_logFEvebenthos_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEvebenthos_area),
                   prob_pos = sum(b_logFEvebenthos_area > 0) / n(),
                   prob_neg = sum(b_logFEvebenthos_area < 0) / n())


# sst Fdiv fish
model %>%
  spread_draws(b_logFDiv_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDiv_area),
                   prob_pos = sum(b_logFDiv_area > 0) / n(),
                   prob_neg = sum(b_logFDiv_area < 0) / n())

# sst fdiv benthos
model %>%
  spread_draws(b_logFDivbenthos_area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDivbenthos_area),
                   prob_pos = sum(b_logFDivbenthos_area > 0) / n(),
                   prob_neg = sum(b_logFDivbenthos_area < 0) / n())


# sampliung area

vars_ext<-get_variables(model)[grep("samp.area",get_variables(model))]

# sst SR fish
model %>%
  spread_draws(b_logSR_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSR_samp.area),
                   prob_pos = sum(b_logSR_samp.area > 0) / n(),
                   prob_neg = sum(b_logSR_samp.area < 0) / n())

# sst SR benthos
model %>%
  spread_draws(b_logSRbenthos_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSRbenthos_samp.area),
                   prob_pos = sum(b_logSRbenthos_samp.area > 0) / n(),
                   prob_neg = sum(b_logSRbenthos_samp.area < 0) / n())

# sst Fric fish
model %>%
  spread_draws(b_logFRic_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRic_samp.area),
                   prob_pos = sum(b_logFRic_samp.area > 0) / n(),
                   prob_neg = sum(b_logFRic_samp.area < 0) / n())

# sst fric benthos
model %>%
  spread_draws(b_logFRicbenthos_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFRicbenthos_samp.area),
                   prob_pos = sum(b_logFRicbenthos_samp.area > 0) / n(),
                   prob_neg = sum(b_logFRicbenthos_samp.area < 0) / n())

# sst FEve fish
model %>%
  spread_draws(b_logFEve_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEve_samp.area),
                   prob_pos = sum(b_logFEve_samp.area > 0) / n(),
                   prob_neg = sum(b_logFEve_samp.area < 0) / n())

# sst FEve benthos
model %>%
  spread_draws(b_logFEvebenthos_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFEvebenthos_samp.area),
                   prob_pos = sum(b_logFEvebenthos_samp.area > 0) / n(),
                   prob_neg = sum(b_logFEvebenthos_samp.area < 0) / n())


# sst Fdiv fish
model %>%
  spread_draws(b_logFDiv_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDiv_samp.area),
                   prob_pos = sum(b_logFDiv_samp.area > 0) / n(),
                   prob_neg = sum(b_logFDiv_samp.area < 0) / n())

# sst fdiv benthos
model %>%
  spread_draws(b_logFDivbenthos_samp.area) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logFDivbenthos_samp.area),
                   prob_pos = sum(b_logFDivbenthos_samp.area > 0) / n(),
                   prob_neg = sum(b_logFDivbenthos_samp.area < 0) / n())




