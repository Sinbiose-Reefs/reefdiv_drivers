library(tidyverse)
library(brms)
library(patchwork)
require(here)
load(here ("output","MCMC_selected_model.RData"))
model <- res$best_model[[1]]
df <- model$data

# replacing standardized by natural scale vars
load(here ("data","modeling_data.RData"))
load (here ("output","FD_results.RData"))



# make plot 
make_plot_data <- function(brm_model, nd_) {
  # all posterior distributions, no random effects
  brms::posterior_epred(
    brm_model, newdata = nd_, re_formula = NA
  ) %>%
    purrr::array_tree(margin = 3) %>%
    purrr::map_dfr(function(x, nd) {
      rbind(
        mean = colMeans(x),
        apply(x, 2, quantile, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
      ) %>%
        t %>%
        cbind(nd, .)
    }, nd = nd_, .id = "rsp") %>%
    dplyr::mutate(
      rsp = dplyr::recode(rsp,
        "logFRicfish1" = "FRic Fish", "logRaofish1" = "Rao's Q Fish",
        "logFRicalgae1" = "FRic Algae", "logRaoalgae1" = "Rao's Q Algae",
        "logFRiccorals1" = "FRic Corals", "logRaocorals1" = "Rao's Q Corals"
      ),
      group = dplyr::recode(rsp,
        "FRic Fish" = "fish", "Rao's Q Fish" = "fish",
        "FRic Algae" = "algae", "Rao's Q Algae" = "algae",
        "FRic Corals" = "corals", "Rao's Q Corals" = "corals"
      ),
      rsp = forcats::fct_relevel(rsp,
        c("FRic Fish", "Rao's Q Fish", "FRic Algae", "Rao's Q Algae",
          "FRic Corals", "Rao's Q Corals")
      )
    )
}

# top plots
top_plots <- data.frame(
  sst_std = seq(min(df$sst_std), max(df$sst_std), length.out = 30),
  SR_fish = 0, SR_algae = 0, SR_corals = 0,
  turbidity_std = 0
)

# natural scale
trans_sst_std <- mean(site_covs$sst) + (sd(site_covs$sst) * top_plots$sst_std)

# mid plots
mid_plots <- data.frame(
  turbidity_std = seq(
    min(df$turbidity_std), max(df$turbidity_std), length.out = 30
  ),
  sst_std = 0,
  SR_fish = 0, SR_algae = 0, SR_corals = 0
)

# natural scale
trans_turb <- mean(site_covs$turbidity) + (sd(site_covs$turbidity) * mid_plots$turbidity_std)

# bottom plots
bot_plots <- data.frame(
  sst_std = 0,
  SR_fish = seq(min(df$SR_fish), max(df$SR_fish), length.out = 30),
  SR_algae = seq(min(df$SR_algae), max(df$SR_algae), length.out = 30),
  SR_corals = seq(min(df$SR_corals), max(df$SR_corals), length.out = 30),
  turbidity_std = 0
)



# preds

all_preds_top <- make_plot_data(model, top_plots)
all_preds_mid <- make_plot_data(model, mid_plots)
all_preds_bot <- make_plot_data(model, bot_plots) %>%
  tidyr::pivot_longer(
    cols = c("SR_fish", "SR_algae", "SR_corals"), names_to = "type",
    values_to = "richness"
  ) %>%
  dplyr::filter(
    !(type == "SR_fish" & group == "algae"),
    !(type == "SR_fish" & group == "corals"),
    !(type == "SR_algae" & group == "fish"),
    !(type == "SR_algae" & group == "corals"),
    !(type == "SR_corals" & group == "algae"),
    !(type == "SR_corals" & group == "fish"),
  )

top_pls <- ggplot(data = all_preds_top, mapping = aes(x = sst_std, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey90") +
  geom_ribbon(mapping = aes(ymin = `10%`, ymax = `90%`), fill = "grey60") +
  geom_line(mapping = aes(colour = group), lwd = 2, show.legend = FALSE) +
  scale_colour_manual(
    values = c(fish = "darkblue", algae = "darkgreen", corals = "darkorange")
  ) +
  labs(y = "", x = "Sea Surface Temperature") +
  facet_wrap(~rsp, ncol = 6, scales = "free") +
  theme_classic() + 
  scale_x_continuous(breaks= seq(min(df$sst_std), max(df$sst_std), length.out = 5),
                     labels = round(seq(min(trans_sst_std), max(trans_sst_std), length.out = 5),1)) + 
  theme (axis.text =element_text(size=7))

mid_pls <- ggplot(data = all_preds_mid, mapping = aes(x = turbidity_std, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey90") +
  geom_ribbon(mapping = aes(ymin = `10%`, ymax = `90%`), fill = "grey60") +
  geom_line(mapping = aes(colour = group), lwd = 2, show.legend = FALSE) +
  scale_colour_manual(
    values = c(fish = "darkblue", algae = "darkgreen", corals = "darkorange")
  ) +
  labs(y = "", x = "Turbidity") +
  facet_wrap(~rsp, ncol = 6, scales = "free") +
  theme_classic() + 
  scale_x_continuous(breaks= seq(min(df$turbidity_std), max(df$turbidity_std), length.out = 5),
                     labels = round(seq(min(trans_turb), max(trans_turb), length.out = 5),1)) + 
  theme (axis.text =element_text(size=7))

# transform SR in the dataset
all_preds_bot$richness_trans <- all_preds_bot$richness
all_preds_bot$richness_trans[which(all_preds_bot$group == "fish")] <- mean(FD_fish$nbsp) + (sd(FD_fish$nbsp) * all_preds_bot$richness_trans[which(all_preds_bot$group == "fish")])
all_preds_bot$richness_trans[which(all_preds_bot$group == "algae")] <- mean(FD_algae$nbsp) + (sd(FD_algae$nbsp) * all_preds_bot$richness_trans[which(all_preds_bot$group == "algae")])
all_preds_bot$richness_trans[which(all_preds_bot$group == "corals")] <- mean(df_corals$SR) + (sd(df_corals$SR) * all_preds_bot$richness_trans[which(all_preds_bot$group == "corals")])


# 
bot_pls <- ggplot(data = all_preds_bot, mapping = aes(x = richness_trans, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey90") +
  geom_ribbon(mapping = aes(ymin = `10%`, ymax = `90%`), fill = "grey60") +
  geom_line(mapping = aes(colour = group), lwd = 2, show.legend = FALSE) +
  scale_colour_manual(
    values = c(fish = "darkblue", algae = "darkgreen", corals = "darkorange")
  ) +
  labs(y = "", x = "Species richness") +
  #scale_x_continuous(breaks= seq(min(all_preds_bot$richness), max(all_preds_bot$richness), length.out = 6),
  #                   labels = round(seq(min(all_preds_bot$richness_trans), max(all_preds_bot$richness_trans), length.out = 6),2))+
  facet_wrap(~rsp, ncol = 6, scales = "free") +
  theme_classic()  + 
  theme (axis.text =element_text(size=7))
  
png(here ("output", "figures", "fig4"),width = 25, 
    height =15, res = 300,units = "cm")

#dev.new(width = 10.5, height = 6.25)
top_pls / mid_pls / bot_pls
dev.off()
