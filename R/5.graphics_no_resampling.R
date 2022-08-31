
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



# raw correlation between taxa

cor (bind_fish_benthos[,c("SR_corals",
                     "SR_algae",
                     "SR_fish")])



# -----------------------------------------------------------

# load results of functional analyses  
load(here("output", "FD_results.RData"))




# ------------------------------------------------

# coefficients of drivers

load (here ("output_no_resampling", "MCMC_runs_multivariate_rarefied_no_aut.RData"))

# model selection

res$looic
res$param
res$best_model


# descripte statistics of fit
library("bayesplot")
library("ggplot2")
#library("rstanarm")   

# extract posterior probs
posterior <- as.array(res$best_model[[1]])
dimnames(posterior)

# plot
# here: https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html
pdf(here ("output_no_resampling","fig2_temp.pdf"),width=7,height=4)
color_scheme_set("blue")
mcmc_areas(posterior, pars = c("b_logSR_BO2_tempmean_ss_std",         
                               "b_logSRbenthos_BO2_tempmean_ss_std", 
                               "b_logFRic_BO2_tempmean_ss_std",
                               "b_logFRicbenthos_BO2_tempmean_ss_std", 
                                   "b_logFEve_BO2_tempmean_ss_std",
                               "b_logFEvebenthos_BO2_tempmean_ss_std",     
                               
                                   "b_logFDiv_BO2_tempmean_ss_std",
                                   
                               
                                   "b_logFDivbenthos_BO2_tempmean_ss_std"),
               prob = 0.5, # 85% intervals
               prob_outer = 0.95, # 95%
               point_est = "median",
           area_method = "equal height")

dev.off()

# area
pdf(here ("output_no_resampling","fig2_area.pdf"),width=7,height=4)
mcmc_areas(posterior, pars = c("b_logSR_area",         
                               "b_logSRbenthos_area", 
                               "b_logFRic_area",
                               "b_logFRicbenthos_area", 
                               "b_logFEve_area",
                               "b_logFEvebenthos_area",     
                               
                               "b_logFDiv_area",
                               
                               
                               "b_logFDivbenthos_area"),
           prob = 0.5, # 85% intervals
           prob_outer = 0.95, # 95%
           point_est = "mean",
           area_method = "equal height")
dev.off()
#sampling area
pdf(here ("output_no_resampling","fig2_samplingarea.pdf"),width=7,height=4)
mcmc_areas(posterior, pars = c("b_logSR_SamplingArea_std",         
                             "b_logSRbenthos_SamplingArea_std", 
                             "b_logFRic_SamplingArea_std",
                             "b_logFRicbenthos_SamplingArea_std", 
                             "b_logFEve_SamplingArea_std",
                             "b_logFEvebenthos_SamplingArea_std",     
                             
                               "b_logFDiv_SamplingArea_std",
                             
                               
                   "b_logFDivbenthos_SamplingArea_std"),
            prob = 0.5, # 85% intervals
            prob_outer = 0.95, # 95%
            point_est = "mean",
            area_method = "equal height")
dev.off()
# ====================================================== #

library(terra)
library(flexmix)
library(modeltools)
library(tidyverse)
library(brms)
library(ggdist)
require(here)
require(tidybayes)

# load results
#load(here ("output_sensitivity","MCMC_runs_multivariate_rarefied_no_aut.rdata"))

# recompile model with more iterations and chains
nd <- res$best_model[[1]]$data %>%
  dplyr::select(-BO2_tempmean_ss_std,
                -area,
                -SamplingArea_std) %>%
  dplyr::mutate(dplyr::across(tidyselect::everything(), log)) %>%
  dplyr::rename_with(~paste0("log", .x)) %>%
  dplyr::mutate(sst = res$best_model[[1]]$data$BO2_tempmean_ss_std,
                area = res$best_model[[1]]$data$area,
                samp.area = res$best_model[[1]]$data$SamplingArea_std)
model <- brms::brm(
  mvbind(logSR, logFRic, logFEve, logFDiv, logSR_benthos, logFRic_benthos, logFEve_benthos, logFDiv_benthos) ~ sst+area+samp.area,
  data = nd, cores = nc, chains = nc, iter = ni, warmup = nb
)

# bayesian R2
br2<- bayes_R2(model) # mean for each reponse
round(br2,2)
round(mean(br2 [,'Estimate']),2) # average across responses

# predicted values
all_preds <- exp(brms::posterior_predict(model))
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
  dplyr::mutate(var_a = case_when(!grepl("benthos", var_a) ~ paste0(var_a, "_fish"), TRUE ~ var_a),
                var_b = case_when(!grepl("benthos", var_b) ~ paste0(var_b, "_fish"), TRUE ~ var_b),
                cor_pair = paste0(var_a, " ~ ", var_b))

# select correlation in the diagonal (i.e., cor of the same metric between fish and benthos)
sel_cors <- c("SR_fish ~ SRbenthos",
              "FRic_fish ~ FRicbenthos",
              "FEve_fish ~ FEvebenthos",
              "FDiv_fish ~ FDivbenthos")

# selecting interesting correlations
cor_out_sel <- cor_out[which(cor_out$cor_pair %in% sel_cors),]
(cor_out_sel %>%
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


cor_orig <- nd %>%
  dplyr::select(-sst) %>%
  dplyr::rename_with(~gsub("^log", "", .x)) %>%
  exp %>%
  cor

ind <- which(upper.tri(cor_orig, diag = TRUE), arr.ind = TRUE)
nn <- dimnames(cor_orig)
cor_orig <- data.frame(var_a = nn[[1]][ind[, 1]], var_b = nn[[2]][ind[, 2]],
                       obs_cor_r = cor_orig[ind]) %>%
  dplyr::filter(var_a != var_b) %>%
  dplyr::mutate(var_a = case_when(!grepl("benthos", var_a) ~ paste0(var_a, "_fish"), TRUE ~ gsub("_", "", var_a)),
                var_b = case_when(!grepl("benthos", var_b) ~ paste0(var_b, "_fish"), TRUE ~ gsub("_", "", var_b)),
                cor_pair = paste0(var_a, " ~ ", var_b))
merged_cor <- merge(post_cor_means, cor_orig[, -c(1:2)])
# factors to represent indexes
#merged_cor$cor_pair <- gsub ("EstRich", "SR",merged_cor$cor_pair )
merged_cor$cor_pair <- gsub ("benthos", "Benthos",merged_cor$cor_pair )
merged_cor$cor_pair <- gsub ("fish", "Fish",merged_cor$cor_pair )
merged_cor$cor_pair <- gsub ("_", "",merged_cor$cor_pair )

# bind residual correlations
require(reshape)
corr_res <- VarCorr(model)$residual__$cor # residals coor
corr_res <- melt(corr_res[,1,]) # melt
corr_res <- corr_res[which(corr_res$value != 1),]# rm corr == 1
# edit labels
# var1
#corr_res$Var1<-gsub ("_","",corr_res$Var1)
corr_res[,1]<-gsub ("log","",corr_res[,1])
corr_res[,1] <- gsub ("EstRich", "SR",corr_res[,1])
corr_res[,1] <- gsub ("benthos", "Benthos",corr_res[,1])
corr_res[,1] <- gsub ("fish", "Fish",corr_res[,1] )
corr_res[,1] <- gsub ("_", "",corr_res[,1] )
corr_res[,1][which(corr_res[,1] == "SR")] <- "SRFish"
corr_res[,1][which(corr_res[,1] == "FRic")] <- "FRicFish"
corr_res[,1][which(corr_res[,1] == "FEve")] <- "FEveFish"
corr_res[,1][which(corr_res[,1] == "FDiv")] <- "FDivFish"
# var2
corr_res[,2]<-gsub ("_","",corr_res[,2])
corr_res[,2]<-gsub ("log","",corr_res[,2])
corr_res[,2] <- gsub ("EstRich", "SR",corr_res[,2])
corr_res[,2] <- gsub ("benthos", "Benthos",corr_res[,2])
corr_res[,2] <- gsub ("fish", "Fish",corr_res[,2])
corr_res[,2] <- gsub ("_", "",corr_res[,2])
corr_res[,2] [which(corr_res[,2] == "SR")] <- "SRFish"
corr_res[,2][which(corr_res[,2] == "FRic")] <- "FRicFish"
corr_res[,2][which(corr_res[,2] == "FEve")] <- "FEveFish"
corr_res[,2][which(corr_res[,2] == "FDiv")] <- "FDivFish"
# bind
corr_res$var3 <- paste (corr_res[,1], corr_res[,2], sep= " ~ ")
corr_res <- corr_res[which(corr_res$var3 %in% merged_cor$cor_pair),]

# bind resid corr
merged_cor$corr_res <- corr_res$value[match (merged_cor$cor_pair,corr_res$var3)] #== merged_cor$cor_pair

# plot
scatter_cor <- ggplot(data = merged_cor) +
  geom_point(mapping = aes(y = obs_cor_r, x = mean_post_r), shape = 21,
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
  geom_point(mapping = aes(y = obs_cor_r, x = corr_res), 
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
  geom_point(mapping = aes(y = mean_post_r, x = corr_res), 
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

# plot diff

a<- ggplot (merged_cor[which(merged_cor$Correlation %in% c("Observed", "Residual", "Predicted")),], 
            aes(x= value,
                y= reorder(cor_pair, - value),
                group = cor_pair,
                color = Correlation
            )) + 
  geom_point(size=2) + 
  theme_classic()  +
  xlab("Pearson's correlation") + 
  ylab ("Pairs of diversity metrics")+
  theme(legend.position = c(0.8,0.8)) + 
  geom_vline(aes(xintercept =0),alpha =0.5,size=1,col = "gray") +
  xlim(c(-1,1))


a

# ==============================================================
# param probability

vars_ext<-get_variables(model)[grep("sst",get_variables(model))]

# sst SR fish
model %>%
  spread_draws(b_logSR_sst) %>% 
  dplyr::summarise(ggdist::median_hdci(b_logSR_sst),
                   prob_pos = sum(b_logSR_sst > 0) / n(),
                   prob_neg = sum(b_logSR_sst < 0) / n())

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




