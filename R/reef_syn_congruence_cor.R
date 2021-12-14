library(terra)
library(flexmix)
library(modeltools)
library(tidyverse)
library(brms)
library(ggdist)
require(here)

# load results
load(here ("Output","MCMC_runs_multivariate_rarefied_no_aut.rdata"))
#load(here ("Output","MCMC_runs_multivariate_rarefied_no_aut_abundW.rdata"))


# recompile model with more iterations and chains
nd <- res$best_model[[1]]$data %>%
  dplyr::select(-BO2_tempmean_ss_std) %>%
  dplyr::mutate(dplyr::across(tidyselect::everything(), log)) %>%
  dplyr::rename_with(~paste0("log", .x)) %>%
  dplyr::mutate(sst = res$best_model[[1]]$data$BO2_tempmean_ss_std)
model <- brms::brm(
  mvbind(logEstRich, logFRic, logFEve, logFDiv, logEstRich_benthos, logFRic_benthos, logFEve_benthos, logFDiv_benthos) ~ sst,
  data = nd, cores = nc, chains = nc, iter = ni
)

# bayesian R2
br2<- bayes_R2(model) # mean for each reponse
mean(br2 [,'Estimate']) # average across responses

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
sel_cors <- c("EstRich_fish ~ EstRichbenthos",
              "FRic_fish ~ FRicbenthos",
              "FEve_fish ~ FEvebenthos",
              "FDiv_fish ~ FDivbenthos")

# selecting interesting correlations
cor_out_sel <- cor_out[which(cor_out$cor_pair %in% sel_cors),]

# aggregating to get averages and plot a vertical line into the histogram
#aggregate (cor_out_sel,by=list(cor_out_sel$cor_pair),
#           FUN=mean)

vline.data <- cor_out_sel %>%
  group_by(cor_pair) %>%
  summarize(z = mean(correlation_r ))

# plot posterior distribution of each pair as histogram
post_cors_plot <- ggplot(data = cor_out_sel) +
  geom_histogram(mapping = aes(x = correlation_r), fill = "dodgerblue2",
                 colour = "grey30") +
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
cor_out_sel %>%
  dplyr::group_by(cor_pair) %>%
  dplyr::summarise(ggdist::median_hdci(correlation_r),
                   prob_pos = sum(correlation_r > 0) / n(),
                   prob_neg = sum(correlation_r < 0) / n())

# plot observed vs. mean predicted correlation
post_cor_means <- cor_out %>%
  dplyr::group_by(cor_pair) %>%
  dplyr::summarise(mean_post_r = mean(correlation_r))

cor_orig <- nd %>%
  dplyr::select(-sst) %>%
  dplyr::rename_with(~gsub("^log", "", .x)) %>%
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
scatter_cor <- ggplot(data = merged_cor) +
  geom_point(mapping = aes(y = obs_cor_r, x = mean_post_r), shape = 21,
             size = 3, fill = "grey90") +
  theme_classic() +
  geom_abline(slope = 1, linetype = 3) +
  labs(y = "Observed Pearson correlation",
       x = "Mean posterior Pearson correlation") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) +
  xlim(-0.6, 1) + ylim(-0.6, 1)
scatter_cor
ggsave("scatter_cor.pdf", scatter_cor, width = 6.5, height = 6)
