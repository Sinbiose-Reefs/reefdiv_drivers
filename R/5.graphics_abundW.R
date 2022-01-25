
## --------------------
#       FIGURES 
# supporting information of analyzes of FD weighted by relative abundance/cover
# --------------------

## call packages
source("R/packages.R")
source("R/functions.R")


# -----------------------------------------
## correlation between variables

load(here("output", "data_to_modeling_GLM.RData"))

# ------------------------------------------ #
# Load rarefied data of fishes and benthos
# sites x transect / video x species

load (here ("output","random_composition_bentos.RData"))
load (here ("output","random_composition_fish.RData"))

# ------------------------------------------ #
#
#     Load environment data

load (here ("output","env_data.RData"))

# ------------------------------------------ #
#
#     Load results of the rarefaction considering relative abundance

# -----------------------------------------------

# MSS algorithm 

# ------------------------------------------ #

load (here ("output","FD_fish_MSS_abundW.RData"))
load (here ("output","FD_benthos_MSS_abundW.RData"))

# remove the ~15  random samples with problems (those with length==2)
FD_fish_MSS_mod<-FD_fish_MSS_abundW [-which(lapply (FD_fish_MSS_abundW,length) ==2)] # 2 is the length of problematic samples 
FD_benthos_MSS_mod<-FD_benthos_MSS_abundW [-which(lapply (FD_benthos_MSS_abundW,length) ==2)]

# also rm from composition
rdm_composition_complete_mod<-rdm_composition_complete [-which(lapply (FD_fish_MSS_abundW,length) ==2)]
rdm_composition_complete_bentos_mod<-rdm_composition_complete_bentos [-which(lapply (FD_benthos_MSS_abundW,length) ==2)]


### -------------------------
#  FIG 2
# coefficients of drivers

load (here ("output", "MCMC_runs_multivariate_rarefied_no_aut_abundW.RData"))

# model selection

res$looic
res$param
res$best_model

# extract model effects
# 95% credible interval
effects_model<- fixef(res$best_model[[1]],
      summary = TRUE,
      robust = FALSE,
      probs = c(0.025, 0.975))

#90% credible interval
effects_model90<-fixef(res$best_model[[1]],
        summary = TRUE,
        robust = FALSE,
        probs = c(0.05, 0.95))


# bind algorithm name
effects_model <- data.frame (effects_model,
                             effects_model90[,3:4])
  

# melt the list to have a df with all results
colnames(effects_model)[c(3:6)] <- c("lower95","upper95","lower90","upper90") # change colnames
# adjust these results to plot
complete_results<- data.frame (effects_model,
                               Parameter = NA,
                               Organism = NA,
                               Index = NA)
# adjust params
complete_results [grep ("Intercept",rownames (complete_results)),"Parameter"] <- "Intercept"
complete_results [grep ("BO2_tempmean_ss_std",rownames (complete_results)),"Parameter"] <- "Temperature"
# adjust organisms
complete_results [grep ("benthos",rownames (complete_results)),"Organism"] <- "Benthos"
complete_results [is.na(complete_results$Organism),"Organism"] <- "Fish"
# adjust indexes
complete_results [grep ("logEstRich",rownames (complete_results)),"Index"] <- "SR"
complete_results [grep ("logFRic",rownames (complete_results)),"Index"] <- "FRic"
complete_results [grep ("logFEve",rownames (complete_results)),"Index"] <- "FEve"
complete_results [grep ("logFDiv",rownames (complete_results)),"Index"] <- "FDiv"

# plot
dodge <- c(0.4,0.4)
pd <- position_dodge(dodge)
pdf_pt <- position_dodge(dodge)

a <- ggplot (complete_results[which (complete_results$Parameter == "Temperature"),], 
             
             aes  (y=Index, x=Estimate, 
                          colour=Organism, fill=Organism)) + 
  
  geom_errorbar(aes(xmin=lower95,xmax=upper95),width = 0.2,size=0.5,
                position=pd)  + 
  geom_errorbar(aes(xmin=lower90,xmax=upper90),width = 0.2,size=1,
                position=pd)  + 
  
  theme_classic() + 
  
  geom_point(position=(pdf_pt), 
             size=1.5)+ 
  
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray50", size=0.5)+
  
  #facet_wrap(~Algorithm+Index,scale="free",ncol=4) + 

  scale_color_manual(values=c("#f2a154", "#0e49b5")) + 
  
  xlab("Standardized effect size") + 
  
  ylab ("Metric") + 
  
  #xlim(-0.5,0.5) +
  
  theme(axis.text.x = element_text(angle = 45,size=7))

a

ggsave (file=here("output","vectorized","Fig2_abundW.pdf"),width=4,height=3)


# ---------------------------------------------
# correlation between metrics and CWM

# extract estimate of each trait level
# all levels
levels_cwm<-colnames(FD_fish_MSS_abundW[[1]]$cwm)

# remove problematic cwm estimates
prob<-(lapply (sapply(FD_fish_MSS_abundW,"[[","cwm"),nrow)) # null is prob
FD_fish_MSS_abundW <- FD_fish_MSS_abundW[which(unlist(lapply(prob,is.null)) == F)]

# extract cwm
cwm_fish <- lapply (levels_cwm, function (k) # for eeach trait level
  
  do.call(cbind, # cbind estimate of cwm for each trait level
          lapply (FD_fish_MSS_abundW, function (i) # in each rarefaction dataset
    
            i$cwm[,which(colnames(i$cwm) %in% k)] # extract estimate of CWM
      )
    )
  )

# obtain the average cwm across rarefaction datasets and trait levels
mean_cwm_fish <- do.call (cbind, 
                          lapply(cwm_fish,rowMeans)#  average across rarefied datasets
  )
colnames(mean_cwm_fish) <- levels_cwm

# the correlation between  functional metrics and cwm
correlation_metrics_cwm <- cor (cbind (cov_fish[,c("FRic","FEve", "FDiv")],
                mean_cwm_fish)
)

# the three first rows matters here
# also eliminating the first three cols (the metrics)
correlation_metrics_cwm <- correlation_metrics_cwm[1:3,-c(1:3)]
# ordering to get the highest values of correlation,if you want positive, and selecting the ntraits required

ntraits <- 2 # ntraits to plot

# names of correlated traits
corr_FRic <- correlation_metrics_cwm[1,][order(correlation_metrics_cwm[1,],decreasing=T)] [1:ntraits]
corr_FRic_neg <- correlation_metrics_cwm[1,][order(correlation_metrics_cwm[1,],decreasing=F)] [1:ntraits]
corr_FEve <- correlation_metrics_cwm[2,][order(correlation_metrics_cwm[2,],decreasing=T)] [1:ntraits]
corr_FEve_neg <- correlation_metrics_cwm[2,][order(correlation_metrics_cwm[2,],decreasing=F)] [1:ntraits]
corr_FDiv <- correlation_metrics_cwm[3,][order(correlation_metrics_cwm[3,],decreasing=T)] [1:ntraits]
corr_FDiv_neg <- correlation_metrics_cwm[3,][order(correlation_metrics_cwm[3,],decreasing=F)] [1:ntraits]

# now obtaining a df with the variation in cwm across the richness gradient (as done for each functional metric)
# first cbind richness and cwm correlated to each functional metric
# correlated to FRic
df_FRic_cwm <- data.frame (EstRich = cov_fish$BO2_tempmean_ss,
                          mean_cwm_fish[,which(colnames(mean_cwm_fish) %in% c(names(corr_FRic),names(corr_FRic_neg)))]
       )
df_FRic_cwm <- cbind(melt(df_FRic_cwm,"EstRich"),
                     Group="Fish",
                     Metric = "FRic")
                     
# correlated to FEve
df_FEve_cwm <- data.frame (EstRich = cov_fish$BO2_tempmean_ss,
                           mean_cwm_fish[,which(colnames(mean_cwm_fish) %in% c(names(corr_FEve),names(corr_FEve_neg)))]
)
df_FEve_cwm <- cbind(melt(df_FEve_cwm,"EstRich"),
                     Group="Fish",
                     Metric = "FEve")

# correlated to FDiv
df_FDiv_cwm <- data.frame (EstRich = cov_fish$BO2_tempmean_ss,
                           mean_cwm_fish[,which(colnames(mean_cwm_fish) %in% c(names(corr_FDiv),names(corr_FDiv_neg)))]
)
df_FDiv_cwm <- cbind(melt(df_FDiv_cwm,"EstRich"),
                     Group="Fish",
                     Metric = "FDiv")

# -----------------------------------------------                     
# do the same for benthos
# extract estimate of each trait level
# all levels
levels_cwm_benthos <- colnames(FD_benthos_MSS_abundW[[1]]$cwm)

# remove problematic cwm estimates
prob<-(lapply (sapply(FD_benthos_MSS_abundW,"[[","cwm"),nrow)) # null is prob
FD_benthos_MSS_abundW <- FD_benthos_MSS_abundW[which(unlist(lapply(prob,is.null)) == F)]

# extract cwm
cwm_benthos <- lapply (levels_cwm_benthos, function (k) # for eeach trait level
  
  do.call(cbind, # cbind estimate of cwm for each trait level
          lapply (FD_benthos_MSS_abundW, function (i) # in each rarefaction dataset
            
            i$cwm[,which(colnames(i$cwm) %in% k)] # extract estimate of CWM
          )
  )
)

# obtain the average cwm across rarefaction datasets and trait levels
mean_cwm_benthos <- do.call (cbind, 
                          lapply(cwm_benthos,rowMeans)#  average across rarefied datasets
)
colnames(mean_cwm_benthos) <- levels_cwm_benthos

# the correlation between  functional metrics and cwm
correlation_metrics_cwm_benthos <- cor (cbind (cov_benthos[,c("FRic","FEve", "FDiv")],
                                        mean_cwm_benthos)
)

# the three first rows matters here
# also eliminating the first three cols (the metrics)
correlation_metrics_cwm_benthos <- correlation_metrics_cwm_benthos[1:3,-c(1:3)]
# ordering to get the highest values
corr_FRic_benthos <- correlation_metrics_cwm_benthos[1,][order(correlation_metrics_cwm_benthos[1,],decreasing=T)] [1:ntraits]
corr_FRic_benthos_neg <- correlation_metrics_cwm_benthos[1,][order(correlation_metrics_cwm_benthos[1,],decreasing=F)] [1:ntraits]
corr_FEve_benthos <- correlation_metrics_cwm_benthos[2,][order(correlation_metrics_cwm_benthos[2,],decreasing=T)] [1:ntraits]
corr_FEve_benthos_neg <- correlation_metrics_cwm_benthos[2,][order(correlation_metrics_cwm_benthos[2,],decreasing=F)] [1:ntraits]
corr_FDiv_benthos <- correlation_metrics_cwm_benthos[3,][order(correlation_metrics_cwm_benthos[3,],decreasing=T)] [1:ntraits]
corr_FDiv_benthos_neg <- correlation_metrics_cwm_benthos[3,][order(correlation_metrics_cwm_benthos[3,],decreasing=F)] [1:ntraits]

# now plotting the variation in cwm across the richness gradient (as done for each functional metric)
# first cbind richness and cwm correlated to each functional metric
# correlated to FRic
df_FRic_cwm_benthos <- data.frame (EstRich = cov_benthos$BO2_tempmean_ss,
                                  mean_cwm_benthos[,which(colnames(mean_cwm_benthos) %in% c(names(corr_FRic_benthos),names(corr_FRic_benthos_neg)))]
)
df_FRic_cwm_benthos <- cbind(melt(df_FRic_cwm_benthos,"EstRich"),
                             Group="Benthos",
                             Metric = "FRic")
                             
# correlated to FEve
df_FEve_cwm_benthos <- data.frame (EstRich = cov_benthos$BO2_tempmean_ss,
                                  mean_cwm_benthos[,which(colnames(mean_cwm_benthos) %in% c(names(corr_FEve_benthos),names(corr_FEve_benthos_neg)))]
)
df_FEve_cwm_benthos <- cbind(melt(df_FEve_cwm_benthos,"EstRich"),
                             Group="Benthos",
                             Metric = "FEve")
                             
# correlated to FDiv
df_FDiv_cwm_benthos <- data.frame (EstRich = cov_benthos$BO2_tempmean_ss,
                                 mean_cwm_benthos[,which(colnames(mean_cwm_benthos) %in% c(names(corr_FDiv_benthos),names(corr_FDiv_benthos_neg)))]
)
df_FDiv_cwm_benthos <- cbind(melt(df_FDiv_cwm_benthos,"EstRich"),
                             Group="Benthos",
                             Metric = "FDiv")

# now rbind all these data
df_plot_cwm <- rbind(df_FRic_cwm,
      df_FRic_cwm_benthos,
      df_FEve_cwm,
      df_FEve_cwm_benthos,
      df_FDiv_cwm,
      df_FDiv_cwm_benthos)

# adjusting levels of the factor metric
df_plot_cwm$Metric<-factor(df_plot_cwm$Metric,
                           levels=c("FRic",
                                    "FEve",
                                    "FDiv"))
# pasting a signal to differentiate positive from negative correlations
df_plot_cwm$variable <- ifelse (df_plot_cwm$variable %in% c(names(corr_FRic_neg),
                                                            names(corr_FEve_neg),
                                                            names(corr_FDiv_neg),
                                                            names(corr_FRic_benthos_neg),
                                                            names(corr_FEve_benthos_neg),
                                                            names(corr_FDiv_benthos_neg)),
                                paste(df_plot_cwm$variable,"(-)",sep=""), 
                                paste(df_plot_cwm$variable,"(+)",sep=""))

# changing the names of trait levels
df_plot_cwm$variable<-gsub ("body_size_","Size.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Size_","Size.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Mouth_position_","Mouth.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Mobility_","Mobility.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("mobility_M","Mobile",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("mobility_S","Sessile",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Body_shape_","BShape.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("reproductive_mode_S","Sexual",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("reproductive_mode_AS","Assexual",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("growth_form_M","Massive",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("growth_form_B","Branching",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("growth_form_EF","Encr/filam",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("growth_form_E","Encrusting",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("carbonate_C","Add",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("carbonate_N","Not add",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Caudal_fin_","Fin.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("Schooling_","Schooling.",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("modularity_M","Modular",(df_plot_cwm$variable))
df_plot_cwm$variable<-gsub ("modularity_S","Sessile",(df_plot_cwm$variable))


# plot

# split data for each "facet"
leg_facet <- split(df_plot_cwm,f = df_plot_cwm$Group)


# plot
p1<-ggplot(df_plot_cwm, aes(x=EstRich,y=value,
                            group = variable,
                            fill=Group,
                            color=variable))+
  #geom_smooth(alpha=0.1,method="lm",se=T) +  #"loess"
  theme_classic()+
  theme(legend.text=element_text(size=9),
        legend.key.width= unit(2, 'cm'),
        axis.text = element_text(size=11),
        axis.title = element_text(size=14),
        strip.text = element_text(size=15),
        strip.background = element_rect(
          color="white", fill="white", size=1.5, linetype="solid"
        ))+
  facet_wrap(~Metric)+
  xlab("Species richness") + 
  ylab ("Trait frequency (CWM)") 


# do it for each "facet"
# benthos
p2 <- p1 %+% leg_facet$Benthos + 
  
  binomial_smooth(
    formula = y ~ splines::ns(x, 2),
    se=T,
    size=1,
    fill = "orange",
    alpha=0.02,
    aes(linetype=variable)) + 
  scale_color_viridis_d(option = "inferno") + 
  
  theme(axis.title.x = element_blank())


# fish
p3 <- p1 %+% leg_facet$Fish +
  
  binomial_smooth(
    formula = y ~ splines::ns(x, 2),
    se=T,
    size=1,
    fill = "blue",
    alpha=0.02,
    aes(linetype=variable)) + 
  
  scale_color_viridis_d(option = "cividis")

# here the plot
pdf(here ("output","vectorized","cwm_richness_abundW.pdf"),height=7,width=10)
grid.arrange(p2,p3)
dev.off()
