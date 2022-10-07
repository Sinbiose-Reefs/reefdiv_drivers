
## --------------------
#       FIGURES
# --------------------

## call packages
source("R/packages.R")
source("R/functions.R")



#load output

load (here ("output", 
            "MCMC_selected_model.RData"))


#---------------------

model <- res$best_model[[1]]


# correlation among predicted values

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
  dplyr::mutate(var_a = case_when(!grepl("benthos", var_a) ~ paste0(var_a, ""), TRUE ~ var_a),
                var_b = case_when(!grepl("benthos", var_b) ~ paste0(var_b, ""), TRUE ~ var_b),
                cor_pair = paste0(var_a, " ~ ", var_b))




# selecting interesting correlations
#cor_out_sel <- cor_out[which(cor_out$cor_pair %in% sel_cors),]
cor_pred <- (cor_out %>%
               dplyr::group_by(cor_pair) %>%
               dplyr::summarise(ggdist::median_hdci(correlation_r),
                                prob_pos = sum(correlation_r > 0) / n(),
                                prob_neg = sum(correlation_r < 0) / n())
)

# residual correlation
# recreate brms summary output from posterior draws
cor_res <- 
  as_draws_df(model) %>%
  dplyr::select(starts_with("rescor__")) %>%
  
  purrr::map_dfr(median_hdci, .id = "pair")





# change names
cor_res$pair <- gsub ("rescor__log", "", cor_res$pair)
cor_res$pair <- gsub ("__log", " ~ ", cor_res$pair)
cor_res$pair <- gsub ("1", "", cor_res$pair)



# correlation from the model

cor_orig <- model$data %>%
  dplyr::select(-sst_std,
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
cor_orig$cor_pair <- gsub ("_","",cor_orig$cor_pair)


# averaged corr

# obs
mean(cor_orig$obs_cor_r)
range(cor_orig$obs_cor_r)

# pred
mean(cor_pred$y)
range(cor_pred$y)

# res
mean(cor_res$y)
range(cor_res$y)


# dataframe for lollipop plot
# predicted
cor_pred_df <- data.frame (Estimate = cor_pred$y,
                           u95CI = cor_pred$ymax,
                           l95CI = cor_pred$ymin,
                           cor_pair=cor_pred$cor_pair,
                           Correlation = "Predicted")

# residual
cor_res_df <- data.frame (Estimate = cor_res$y,
                          u95CI = cor_res$ymax,
                          l95CI = cor_res$ymin,
                          cor_pair=cor_res$pair,
                          Correlation = "Residual")

# observed
cor_obs_df <- data.frame (Estimate = cor_orig$obs_cor_r,
                          u95CI = cor_orig$obs_cor_r,
                          l95CI = cor_orig$obs_cor_r,
                          cor_pair=cor_orig$cor_pair,
                          Correlation = "Observed")

# bind 
df_loll <- rbind (cor_obs_df,
                  cor_pred_df,
                  cor_res_df)

df_loll$cor_pair <- gsub ("1", "",df_loll$cor_pair)

# facets
df_loll$facetting <- gsub ("FRic","", df_loll$cor_pair)
df_loll$facetting <- gsub ("Rao","", df_loll$facetting)
df_loll$facetting <- firstup(df_loll$facetting)

# order facets
df_loll$facetting <- factor(df_loll$facetting,
                            levels  = c("Algae ~ algae",
                                        "Corals ~ corals",
                                        "Fish ~ fish",
                                        "Algae ~ corals",
                                        "Fish ~ algae",
                                        "Fish ~ corals"))

# labs Y
df_loll$cor_pair <- gsub ("fish","", df_loll$cor_pair)
df_loll$cor_pair <- gsub ("corals","", df_loll$cor_pair)
df_loll$cor_pair <- gsub ("algae","", df_loll$cor_pair)


# lollipop plot (plot of difference in correlation)
library(ggplot2)
library(ggh4x)

# save
pdf(file=here("output",
              "figures", 
              "fig5"),height=5,width=10)


# Change baseline
ggplot(df_loll, 
       aes(x=Estimate, 
           y=cor_pair,
           color = Correlation,
           fill = Correlation,
           group = cor_pair)) +
  
  
  geom_segment( aes(x=Estimate, xend=0, 
                    y=reorder(cor_pair,  Estimate), 
                    yend=reorder(cor_pair,  Estimate)), 
                color="grey",
                size=1,
                linetype='dotted',
                position = position_jitter( height=0.2,width = 0)) +
  
  
  
  
  geom_point(size=3,
             position = position_jitter(height = 0.2, width = 0)) +
  
  geom_errorbar(aes(xmin=l95CI, xmax=u95CI), 
                width=0,
                position = position_jitter(height = 0.2),
                alpha=0.5)   +
  
  scale_color_viridis_d(option ="viridis", begin = 0,end=0.75)+
  
  xlab("Pearson's correlation") + 
  ylab ("Pairs of functional metrics")+
  geom_vline(aes(xintercept =0),alpha =0.5,size=1,col = "gray") +
  xlim(c(-0.7,.9)) +
  theme_classic() + 
  facet_wrap(~facetting,ncol = 3,scales="free")+
  force_panelsizes(rows = c(0.1, 0.2)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  ) 





dev.off()


