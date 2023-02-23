
load (here ("data","modeling_data_old.RData"))
load (here ("output","FD_results_old.RData"))

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

bind_fish_benthos_old<- cbind (site_covs,
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



# new data

load (here ("data","modeling_data.RData"))
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


# only fish fric did not match
bind_fish_benthos_old<-bind_fish_benthos_old[,-which(colnames(bind_fish_benthos_old) == "reef_area")]
bind_fish_benthos_old == bind_fish_benthos

bind_fish_benthos_old$FRic_fish == bind_fish_benthos$FRic_fish
round (bind_fish_benthos_old$FRic_fish,7) == round(bind_fish_benthos$FRic_fish,7)
bind_fish_benthos_old$Rao_fish == bind_fish_benthos$Rao_fish # rao is ok


# average values of metrics to present in the Results
# offshore distance was not used in the models (in the current version it was calculated using sf package)
 (round(apply(bind_fish_benthos_old[c(3:15,17:29)],2,mean,na.rm=T),3) == round(apply(bind_fish_benthos[c(3:15,17:29)],2,mean,na.rm=T),3))
(round(apply(bind_fish_benthos_old[c(3:15,17:29)],2,sd,na.rm=T),3) == round(apply(bind_fish_benthos[c(3:15,17:29)],2,sd,na.rm=T),3))

bind_fish_benthos[c(3:15,17:29)] == bind_fish_benthos_old[c(3:15,17:29)] # only fish fric differs
round (bind_fish_benthos[c(3:15,17:29)],10) == round(bind_fish_benthos_old[c(3:15,17:29)],10) # but very subtly


# data in the models  (the simplest model)
load(file=here ("output","fit_simple.RData"))
fit_simple_old <- fit_simple
load(file=here ("output","fit_sst_region.RData"))
fit_simple_new <- fit_simple

# site names are equal
rownames(fit_simple_old$data) == rownames(fit_simple_new$data)
fit_simple_new$data$FRic_fish == fit_simple_old$data$FRic_fish # not equal FRic
round(fit_simple_new$data$FRic_fish,10) == round (fit_simple_old$data$FRic_fish,10) # decimla differences
(fit_simple_old$data == fit_simple_new$data)




#  FD results
# fish fric differed
load (here ("output","FD_results.RData"))
new <- FD_fish$FRic

load (here ("output","FD_results_old.RData"))
old <- FD_fish$FRic
round(new,10) == round(old,10)
round(new,100) == round(old,100)
old == new


# checking modeling data

load (here ("data","modeling_data_old.RData"))
old_fish<- comp_fish
old_benthos<- comp_benthos
old_sitecovs<- site_covs[,-which(colnames(site_covs) == "reef_area")]
old_eff <- effort_dataframe


# new
load (here ("data","modeling_data.RData"))


table(comp_fish == old_fish) # fish composition is exactly the same
table(comp_benthos == old_benthos) #  benthic composition is exactly the same
(site_covs == old_sitecovs) # distance offshore did not match (not used) - differ because now we're using sf package
table(effort_dataframe == old_eff) # sampling data are the same





