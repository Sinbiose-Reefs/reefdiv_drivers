
# ------------------------------------------------------------------- #
#         organizacao de dados, e analises para o manuscrito 
#                 "Diversity patterns and drivers" 
#-------------------------------------------------------------------- #

## call packages
source("R/packages.R")
source("R/functions.R")
source("R/function_lomolino_richness.R")

# Load data

# ------------------------------ #
#        benthos data
# ------------------------------ #

bentos <- read.xlsx(here("data","detection","Updated_compiled_quadrats_allsites.xlsx"),
                    sheet = 1, colNames = TRUE,detectDates=F)

## converter data em data - bug do pacote openxlsx
bentos$eventDate <-convertToDate(bentos$eventDate)

# Remove benthos from dataset as below
rm_sp <- c("Areia.e.Cascalho","Desconhecido","Estrela","ourico1","ourico2",
           "Outra.ascidia","Outro.anthozoa", "Outro.echinoderma",
           "Outro.hydrozoa","Quadrado","Outro.crustaceo","Sombra")

bentos <- bentos[which(bentos$Taxon %in% rm_sp == F),]
# to lower these names 
bentos$Taxon <- tolower (bentos$Taxon)

## dados dos peixes
peixes <- read.xlsx(here("data","detection","UpdatedData_RMorais_et_al_2017.xlsx"),
                    sheet = 1, colNames = TRUE,detectDates=F)

## converter data em data - bug do pacote openxlsx
peixes$eventDate <-convertToDate(peixes$eventDate)

## filtrar os dados de peixes de acordo com o minimo de ano de bentos
# peixes <- peixes [which(peixes$eventYear >= min (bentos$eventYear) & peixes$eventYear <= max (bentos$eventYear)),]

## modificar o eventID removendo o ano, desde que escolhemos um lapso temporal que cobre 2010 a 2014, 
## de modo que possamos analisar os dados se peixes foi coletado em 2012 e bentos em 2014, por exemplo 

bentos$eventID_MOD <- substr(bentos$eventID, 1,nchar(as.character(bentos$eventID))-5) 
peixes$eventID_MOD  <- substr(peixes$eventID, 1,nchar(as.character(peixes$eventID))-5) 

## fazer um histograma pra saber o numero de eventIDS por ano

barplot_function(df1=bentos,#bentos
                 df2= peixes #peixes
                 )

## se quiser o subset de sitios que correspondem em ambos os data sets
## subset entre os datasets, desde que peixes ou bentos foram amostrados nos mesmos locais
# peixes_subset <- peixes [which(peixes$eventID_MOD %in% bentos$eventID_MOD),]
peixes_subset <- peixes

## da mesma forma, pegar o subset de ID de bentos que estao nos dados de peixes
# bentos_subset <- bentos [which(bentos$eventID_MOD %in% peixes$eventID_MOD),]
bentos_subset <- bentos
##  criar uma ID numerica para o observador
peixes_subset$ID.observer <- as.numeric(as.factor(peixes_subset$Observer))

## ID para o sitio (peixes)

peixes_subset$locality_site <- paste(peixes_subset$Locality, 
                                     peixes_subset$Site,sep=".")


## ID para o sitio (bentos)

bentos_subset$locality_site <- paste(bentos_subset$Locality, 
                                     bentos_subset$Site,sep=".")

### mapa para conferir a posicao dos sitios - com base no subset dos dados que correspondem
## tem um monte de pontos porque tem todos os eventIDs ai

initial_map_function (df1 = bentos_subset,# red
                      df2 = peixes_subset) # black

#---------------------------------------------------#
#             dados de composicao                   #
#---------------------------------------------------#

#########################################
## fish composition (site x spp)
#########################################

comp_peixes <-  cast(peixes_subset,
                     formula = Site  ~ ScientificName,
                     value= "IndCounting",
                     fun.aggregate = sum)

###########################################################################
## data for SAC (species accumulation curve, site  x transect x spp)
###########################################################################

# unuique sites
sites_fish_complete <- unique (peixes_subset$Site)
sites_fish_complete <- sites_fish_complete[order(sites_fish_complete,decreasing=F)]

# tables per site
rar_peixes <- lapply (sites_fish_complete, function (i)
    
    cast (peixes_subset [which (peixes_subset$Site == i),],
          formula = Transect_id ~ ScientificName,
          value="IndCounting",
          fun.aggregate = sum)[,-1]
    )

# SAC based on the number of samples,
# nrandom samples  (for all analysis)
niter <- 999
rarefied_richness_fish <- lapply (rar_peixes, 
                                  specaccum,
                                  method="random", 
                                  permutations=niter)

### plotting 
par (mfrow=c(1,1))
plot(NA,
     xlim=c(0,250),
     ylim=c(-5,100),
     xlab = "Number of samples",
     ylab = "Fish richness")

lapply (rarefied_richness_fish,plot,add=T)
abline(v=5,lwd=2,col="gray50",lty=2)

##---------------------------------- ## 
# richness estimate based on the minimum number of samples
## ---------------------------------- ##

min_samples <- sapply (rarefied_richness_fish,"[[","sites")
# finding the minimum number of samples across sites
min_samples<- min(unlist(sapply (min_samples,max,simplify=F)))

# finding estimated richness
est_rich <- sapply (rarefied_richness_fish,"[[","richness")
# finding the richness for min_samples
est_rich <- lapply (est_rich, function (i) i [min_samples])

# finding sd
# finding estimated richness
sd_rich <- sapply (rarefied_richness_fish,"[[","sd")
# finding the richness for min_samples
sd_rich <- lapply (sd_rich, function (i) i [min_samples])

# toatl number of samples per site
total_samples <- sapply (rarefied_richness_fish,"[[","sites")
total_samples <- unlist (lapply (total_samples, max))

# res table
res_table_samples <- data.frame (Site = sites_fish_complete,
                                 EST.Rich=unlist(est_rich),
                                 SD.Rich=unlist(sd_rich),
                                 n_samples = min_samples,
                                 n_total= total_samples)

# finding site composition based on min_samples 
# a random sample of five transects

# -------------------------------#
# 1) random sample based on min samples
# -------------------------------#

# obtain random compositions
nc <- 3
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores

# export your data and function
clusterExport(cl, c("niter", 
                    "fc_random_composition",
                    "rar_peixes",
                    "min_samples"))

rdm_composition <- parLapply (cl, seq(1,niter), function (k)
                           lapply(rar_peixes, function (i)
  
                fc_random_composition(data = i,nsamples = min_samples,
                                      replace= F)
                

                     )
                  )
                

stopCluster (cl)

# list of species to adjust cols
list_spp <- unique(unlist(lapply (rar_peixes, function (k)
   colnames(k))))

# matrix to bind (with missing names)

rdm_composition_complete <-lapply (rdm_composition, function (i)
  do.call(rbind, lapply (i, function (k) {

    bind_to_rdm_composition <- matrix(0, 
                                      ncol = length(which (list_spp %in% colnames (k)  == F)),
                                 dimnames = list (NULL,
                                                  list_spp[which (list_spp %in% colnames (k)  == F)]
                                 ))
    # cbind
    bind_rdm_composition <- cbind(k,bind_to_rdm_composition)
    # ordering by sp name to avoid problems when combining samples from diff sites
    bind_rdm_composition <- bind_rdm_composition [,order (colnames(bind_rdm_composition))];
    bind_rdm_composition

})))

rdm_composition_complete

# -------------------------------------------------------------------- #
# 2) random sample based on asymptote 
#  (Standard deviation of estimates > 1 and != 0)
#--------------------------------------------------------------------- #

## df with extracted estimates
# a site get the asymptote when richness sd was < 1 species  (no overall difference after that)
sp_accum_fish_asymptote <- lapply (seq (1,length(rarefied_richness_fish)), function (i){
  # extract statistics
  df_res <- data.frame (
    Site = sites_fish_complete[i],
    # estimated richness
    EST.Rich = rarefied_richness_fish[[i]]$richness[max(which(rarefied_richness_fish[[i]]$sd < 1 & rarefied_richness_fish[[i]]$sd > 0))],
    # standard deviation 
    SD.Rich = rarefied_richness_fish[[i]]$sd[max(which(rarefied_richness_fish[[i]]$sd < 1 & rarefied_richness_fish[[i]]$sd > 0))],
    # nsites to achieve good estimate
    n_samples = rarefied_richness_fish[[i]]$sites[max(which(rarefied_richness_fish[[i]]$sd < 1 & rarefied_richness_fish[[i]]$sd > 0))],
    # total number of samples
    tot_samples = length (rarefied_richness_fish[[i]]$sites)
  )
  # bind confidence intervals
  df_res$upper <- df_res$EST.Rich+1.96*(df_res$SD.Rich/df_res$n_samples)
  df_res$lower <- df_res$EST.Rich-1.96*(df_res$SD.Rich/df_res$n_samples )
  
  ; # return
  df_res
  
})

# list into  table
res_sp_accum_fish_asymptote <- do.call(rbind,sp_accum_fish_asymptote)

# remove those sites that did not reach asymptote
rar_peixes_asymptote <- rar_peixes[which(is.na(res_sp_accum_fish_asymptote$n_samples)!= T)]

# sample size
nsamples_asymptote_fish <- res_sp_accum_fish_asymptote$n_samples[which(is.na(res_sp_accum_fish_asymptote$n_samples) !=T)]

# obtain random compositions

cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores

# export your data and function
clusterExport(cl, c("niter", 
                    "fc_random_composition",
                    "rar_peixes_asymptote",
                    "nsamples_asymptote_fish"))

rdm_composition_asymptote <- parLapply (cl, seq(1,niter), function (k)
  
                                lapply(seq(1,length(rar_peixes_asymptote)), function (i)
    
                              fc_random_composition(data = rar_peixes_asymptote[[i]],
                                                    nsamples = nsamples_asymptote_fish[i],
                                                    replace = F)
                                  )
                              )
                          
stopCluster(cl)

## complete with missing species
rdm_composition_asymptote_complete <-lapply (rdm_composition_asymptote, function (i)
  do.call(rbind, lapply (i, function (k) {
    
    bind_to_rdm_composition <- matrix(0, 
                                      ncol = length(which (list_spp %in% colnames (k)  == F)),
                                      dimnames = list (NULL,
                                                       list_spp[which (list_spp %in% colnames (k)  == F)]
                                      ))
    # cbind
    bind_rdm_composition <- cbind(k,bind_to_rdm_composition)
    # ordering by sp name to avoid problems when combining samples from diff sites
    bind_rdm_composition <- bind_rdm_composition [,order (colnames(bind_rdm_composition))];
    bind_rdm_composition
    
  })))


############################################################################
####     rarefaction based on Lomolino function (from Diego Barneche)
#############################################################################

vegan_data <- lapply (rarefied_richness_fish, function (i)
  
  data.frame(ln_richness = log(i$richness), 
             sample = i$sites))

# MCMC settings
ni <- 30000
nb <- 20000
nt <- 20
nc <- 3

test <- run_lomolino_model(vegan_data[[14]][,-1])
exp(fixef(test)["lnxmid_Intercept", "Estimate"])

model_lomolino <- lapply(vegan_data,run_lomolino_model)
#save(model_lomolino, file = here ("output", "model_lomolino_fish.RData"))
load(here ("output", "model_lomolino_fish.RData"))

# Desse modelo voc? pode extrair a estimativa m?dia de lnasym, que representa a riqueza m?xima estimada daquele s?tio.

# asymp
asymp <- lapply (model_lomolino, function (i) 
  exp(fixef(i)["lnasym_Intercept", "Estimate"]))

# xmid = xmid is the area where half of the maximum richness is achieved
xmid <- lapply (model_lomolino, function (i) 
  exp(fixef(i)["lnxmid_Intercept", "Estimate"])
)

xmid<- lapply(xmid,round)

# df response lomolino
res_lomolino <- cbind (xmid=unlist(xmid),
                       obs_ss = unlist(lapply (rar_peixes,nrow)),
                       asymp = unlist(asymp ))

#subsetting for those sites with goood S50 estimates 

rar_peixes_lomolino <- rar_peixes [which(res_lomolino[,1] < res_lomolino[,2])]
xmid <- xmid[which(res_lomolino[,1] < res_lomolino[,2])]

# apply (res_lomolino[which(res_lomolino[,1] < res_lomolino[,2]),],2,mean)
# apply (res_lomolino[which(res_lomolino[,1] < res_lomolino[,2]),],2,sd)

# obtain composition based on lomolino function

# obtain random compositions
nc <- 3
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores

# export your data and function
clusterExport(cl, c("niter", 
                    "fc_random_composition",
                    "rar_peixes_lomolino",
                    "xmid"))


rdm_composition_lomolino <-
  
  parLapply (cl,seq(1,niter), function (k) {
    
    tryCatch (
      
      lapply(seq(1,length(rar_peixes_lomolino)), function (i)
        
        fc_random_composition(data = rar_peixes_lomolino[[i]],
                              nsamples = xmid [[i]],
                              replace =F)), 
      error = function (e)
        return (e))})


stopCluster(cl)


# list of species to adjust cols
list_spp <- unique(unlist(lapply (rar_peixes_lomolino, function (k)
  colnames(k))))

## complete with missing species
rdm_composition_lomolino_complete <-lapply (rdm_composition_lomolino, function (i)
  do.call(rbind, lapply (i, function (k) {
    
    bind_to_rdm_composition <- matrix(0, 
                                      ncol = length(which (list_spp %in% colnames (k)  == F)),
                                      dimnames = list (NULL,
                                                       list_spp[which (list_spp %in% colnames (k)  == F)]
                                      ))
    # cbind
    bind_rdm_composition <- cbind(k,bind_to_rdm_composition)
    # ordering by sp name to avoid problems when combining samples from diff sites
    bind_rdm_composition <- bind_rdm_composition [,order (colnames(bind_rdm_composition))];
    bind_rdm_composition
    
  })))

#

save (rar_peixes,# observed data
      rar_peixes_asymptote,
      rar_peixes_lomolino,
      rdm_composition_complete, # random composition based on min samples
      rdm_composition_asymptote_complete,## random composition based on SD criteria
      rdm_composition_lomolino_complete,# random composition baed on xmid
      res_table_samples, # results of estimates based on min samples
      res_sp_accum_fish_asymptote, # results of estimates based on SD
      res_lomolino, # results based on lomolino
      sites_fish_complete, # analyzed sites
      file = here("output", "random_composition_fish.RData"))


#########################################
## benthos composition (site  x spp)
#########################################

comp_bentos <-  cast(bentos_subset,
                     formula = Sites  ~ Taxon,
                     value= "Cover",
                     fun.aggregate = mean)
#
###################################################################
## data for SAC (species accumulation curve, site  x video x spp)
###################################################################
sites_bentos_complete <- unique (bentos_subset$Site)
sites_bentos_complete <- sites_bentos_complete [order(sites_bentos_complete,decreasing = F)]

rar_bentos <- lapply (sites_bentos_complete, function (i)
  
  cast (bentos_subset [which (bentos_subset$Site == i),],
        formula = Video_number ~ Taxon,
        value="Cover",
        fun.aggregate = max)[,-1])

## SAC based on the number of samples
niter<-999
rarefied_richness_bentos<- lapply (rar_bentos, 
                                  specaccum,
                                  method="random", ## method="rarefaction" - based on individuals
                                  permutations=niter)


### plotting 
par (mfrow=c(1,1))
plot(NA,
     xlim=c(0,25),
     ylim=c(-2,40),
     xlab = "Number of samples",
     ylab = "Benthos richness")

lapply (rarefied_richness_bentos,plot,add=T)
abline(v=3,lwd=2,col="gray50",lty=2)

##---------------------------------- ## 
# based on minimum number of samples
## ---------------------------------- ##
min_samples_bentos <- sapply (rarefied_richness_bentos,"[[","sites")
# finding the minimum number of samples across sites
min_samples_bentos<- min(unlist(sapply (min_samples_bentos,max,simplify=F)))

# finding estimated richness
est_rich <- sapply (rarefied_richness_bentos,"[[","richness")
# finding the richness for min_samples
est_rich <- lapply (est_rich, function (i) i [min_samples_bentos])

# finding sd
# finding estimated richness
sd_rich <- sapply (rarefied_richness_bentos,"[[","sd")
# finding the richness for min_samples
sd_rich <- lapply (sd_rich, function (i) i [min_samples_bentos])

# toatl number of samples per site
total_samples <- sapply (rarefied_richness_bentos,"[[","sites")
total_samples <- unlist (lapply (total_samples, max))

# res table
res_table_samples_bentos <- data.frame (Site = sites_bentos_complete,
                                 EST.Rich=unlist(est_rich),
                                 SD.Rich=unlist(sd_rich),
                                 n_samples = min_samples_bentos,
                                 n_total= total_samples)

# finding site composition based on min_samples 
# a random sample of five transects

# ----------------------------------------#
# 1) random sample based on min samples
# ----------------------------------------#

# obtain random compositions
nc <- 3
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores

# export your data and function
clusterExport(cl, c("niter", 
                    "fc_random_composition",
                    "rar_bentos",
                    "min_samples_bentos"))

rdm_composition_bentos <- parLapply (cl, seq(1,niter), function (k)
  lapply(rar_bentos, function (i)
    fc_random_composition(data = i[,-1],
                          nsamples = min_samples_bentos,
                          replace= F)
  ))

stopCluster (cl)

# list of species to adjust cols
list_spp <- unique(unlist(lapply (rar_bentos, function (k)
  colnames(k))))

# matrix to bind (with missing names)

rdm_composition_complete_bentos <-lapply (rdm_composition_bentos, function (i)
  do.call(rbind, lapply (i, function (k) {
    
    bind_to_rdm_composition <- matrix(0, 
                                      ncol = length(which (list_spp %in% colnames (k)  == F)),
                                      dimnames = list (NULL,
                                                       list_spp[which (list_spp %in% colnames (k)  == F)]
                                      ))
    # cbind
    bind_rdm_composition <- cbind(k,bind_to_rdm_composition)
    # ordering by sp name to avoid problems when combining samples from diff sites
    bind_rdm_composition <- bind_rdm_composition [,order (colnames(bind_rdm_composition))];
    bind_rdm_composition
    
  })))

rdm_composition_complete_bentos

# -------------------------------------------------------------------- #
# 2) random sample based on min samples
# based on asymptote (Standard deviation of estimates > 1 and != 0)
#--------------------------------------------------------------------- #

## df with extracted estimates
# a site get the asymptote when richness sd was < 1 species  (no overall difference after that)
sp_accum_bentos_asymptote <- lapply (seq (1,length(rarefied_richness_bentos)), function (i){
  # extract statistics
  df_res <- data.frame (
    Site = sites_bentos_complete[i],
    # estimated richness
    EST.Rich = rarefied_richness_bentos[[i]]$richness[max(which(rarefied_richness_bentos[[i]]$sd < 1 & rarefied_richness_bentos[[i]]$sd > 0))],
    # standard deviation 
    SD.Rich = rarefied_richness_bentos[[i]]$sd[max(which(rarefied_richness_bentos[[i]]$sd < 1 & rarefied_richness_bentos[[i]]$sd > 0))],
    # nsites to achieve good estimate
    n_samples = rarefied_richness_bentos[[i]]$sites[max(which(rarefied_richness_bentos[[i]]$sd < 1 & rarefied_richness_bentos[[i]]$sd > 0))],
    # total number of samples
    tot_samples = length (rarefied_richness_bentos[[i]]$sites)
  )
  # bind confidence intervals
  df_res$upper <- df_res$EST.Rich+1.96*(df_res$SD.Rich/df_res$n_samples)
  df_res$lower <- df_res$EST.Rich-1.96*(df_res$SD.Rich/df_res$n_samples )
  
  ; # return
  df_res
  
})

# list into  table
res_sp_accum_bentos_asymptote <- do.call(rbind,sp_accum_bentos_asymptote)

# remove those sites that did not reach asymptote
rar_bentos_asymptote <- rar_bentos[which(is.na(res_sp_accum_bentos_asymptote$n_samples)!= T)]

# sample size
nsamples_asymptote_bentos <- res_sp_accum_bentos_asymptote$n_samples[which(is.na(res_sp_accum_bentos_asymptote$n_samples) !=T)]

# obtain random compositions
nc <- 3
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores

# export your data and function
clusterExport(cl, c("niter", 
                    "fc_random_composition",
                    "rar_bentos_asymptote",
                    "nsamples_asymptote_bentos"))

rdm_composition_asymptote_bentos <- parLapply (cl, seq(1,niter), function (k)
  lapply(seq(1,length(rar_bentos_asymptote)), function (i)
    
    fc_random_composition(data = rar_bentos_asymptote[[i]][,-1],
                          nsamples = nsamples_asymptote_bentos[i],
                          replace= F)
  ))

stopCluster(cl)

## complete with missing species
rdm_composition_asymptote_complete_bentos <-lapply (rdm_composition_asymptote_bentos, function (i)
  do.call(rbind, lapply (i, function (k) {
    
    bind_to_rdm_composition <- matrix(0, 
                                      ncol = length(which (list_spp %in% colnames (k)  == F)),
                                      dimnames = list (NULL,
                                                       list_spp[which (list_spp %in% colnames (k)  == F)]
                                      ))
    # cbind
    bind_rdm_composition <- cbind(k,bind_to_rdm_composition)
    # ordering by sp name to avoid problems when combining samples from diff sites
    bind_rdm_composition <- bind_rdm_composition [,order (colnames(bind_rdm_composition))];
    bind_rdm_composition
    
  })))


############################################################################
####     rarefaction based on Lomolino function (from Diego Barneche)
#############################################################################

vegan_data <- lapply (rarefied_richness_bentos, function (i)
  
  data.frame(ln_richness = log(i$richness), 
             sample = i$sites))

# run model

model_lomolino <- lapply(vegan_data,run_lomolino_model)
#save(model_lomolino, file = here ("output", "model_lomolino_benthos.RData"))
load(here ("output", "model_lomolino_benthos.RData"))

#Desse modelo voc? pode extrair a estimativa m?dia de lnasym, que representa a riqueza m?xima estimada daquele s?tio.

# asymp
asymp <- lapply (model_lomolino, function (i) 
  exp(fixef(i)["lnasym_Intercept", "Estimate"]))

# xmid = xmid is the area where half of the maximum richness is achieved
xmid <- lapply (model_lomolino, function (i) 
  exp(fixef(i)["lnxmid_Intercept", "Estimate"])
)
xmid<- lapply(xmid,round)


# df response lomolino
res_lomolino_bentos <- cbind (xmid=unlist(xmid),
                              obs_ss = unlist(lapply (rar_bentos,nrow)),
                              asymp = unlist(asymp ))

#subsetting for those sites with goood S50 estimates 

rar_bentos_lomolino <- rar_bentos [which(res_lomolino_bentos[,1] < res_lomolino_bentos[,2])]
xmid <- xmid[which(res_lomolino_bentos[,1] < res_lomolino_bentos[,2])]

table(res_lomolino_bentos[,1] > res_lomolino_bentos[,2])

# obtain composition based on lomolino function
# basic statistics

apply (res_lomolino_bentos[which(res_lomolino_bentos[,1] < res_lomolino_bentos[,2]),],2,mean)
apply (res_lomolino_bentos[which(res_lomolino_bentos[,1] < res_lomolino_bentos[,2]),],2,sd)

# number of iterations to get niter random compositions
niter <- 999

# obtain random compositions
nc <- 3
cl <- makeCluster(nc) ## number of cores = generally ncores -1

# exportar pacote para os cores

# export your data and function
clusterExport(cl, c("niter", 
                    "fc_random_composition",
                    "rar_bentos_lomolino",
                    "xmid"))


rdm_composition_lomolino_bentos <- parLapply (cl,seq(1,niter), function (k){
  
  tryCatch (
    lapply(seq(1,length(rar_bentos_lomolino)), function (i)
      fc_random_composition(data = rar_bentos_lomolino[[i]][,-1],
                            nsamples = xmid [[i]],
                            replace=F)),
    error = function (e)
      return (e))
}

)

stopCluster(cl)


# list of species to adjust cols
list_spp <- unique(unlist(lapply (rar_bentos_lomolino, function (k)
  colnames(k))))

## complete with missing species
rdm_composition_lomolino_complete_bentos <-lapply (rdm_composition_lomolino_bentos, function (i)
  do.call(rbind, lapply (i, function (k) {
    
    bind_to_rdm_composition <- matrix(0, 
                                      ncol = length(which (list_spp %in% colnames (k)  == F)),
                                      dimnames = list (NULL,
                                                       list_spp[which (list_spp %in% colnames (k)  == F)]
                                      ))
    # cbind
    bind_rdm_composition <- cbind(k,bind_to_rdm_composition)
    # ordering by sp name to avoid problems when combining samples from diff sites
    bind_rdm_composition <- bind_rdm_composition [,order (colnames(bind_rdm_composition))];
    bind_rdm_composition
    
  })))

#

save (
      rar_bentos,
      rar_bentos_asymptote,
      rar_bentos_lomolino,
      rdm_composition_complete_bentos, # random composition based on min samples
      rdm_composition_asymptote_complete_bentos,## random composition based on SD criteria
      rdm_composition_lomolino_complete_bentos, # random composition based on lomolino 
      res_table_samples_bentos, # results of estimates based on min samples
      res_sp_accum_bentos_asymptote, # results of estimates based on SD
      res_lomolino_bentos, #results of estimates based on lomolino
      sites_bentos_complete, # analyzed sites
      file = here("output", "random_composition_bentos.RData"))


#######################################################################################
#######################################################################################
################         covariaveis de sitio                   #######################
#######################################################################################
#######################################################################################

## de sitio (aquelas que variam por sitio)

## sitios do nordeste
nord <- cast(peixes_subset,
             formula = Site  ~ Region,
             value= "IndCounting",
             fun.aggregate = sum)

## regiao 
regiao <- matrix(NA, nrow(nord),ncol = 1,
                 dimnames=list(nord$Site,
                               "Region"))

regiao [nord [,3] >0] <-  "Oceanic Islands"
regiao [nord [,2] >0] <- "Northeastern"
regiao [nord [,4] >0]<- "Southeastern"

## numero de transeccoes e videos (esforço)

#lista_sitios_peixes <- unique(comp_peixes$locality_site)

# encontrando o ID das transeccoes unicas, e vendo seu length == numero de transectos
n_transeccoes <- lapply (sites_fish_complete, function (i)
  length (unique (peixes_subset [which(peixes_subset$Site == i),"Transect_id"]))
)
# encontrando numero de videos unicos por sitio
# fish site names as comparison
n_videos <- lapply (sites_fish_complete, function (i)
  length (unique (bentos_subset [which(bentos_subset$Sites == i),"Video_number"]))
)

# montando uma tabela
tabela_esforco <- data.frame (Site = sites_fish_complete,
                              n_transectos_peixes = unlist(n_transeccoes),
                              n_videos_bentos = unlist(n_videos))
## lista de dados

dados_peixes_bentos <- list(peixes = comp_peixes,
                            bentos = comp_bentos)

## coordenadas geográficas para os mapas e analises
# peixes
coordenadas_peixes <- aggregate(peixes_subset, 
                                by= list (peixes_subset$Site), 
                                FUN=mean)[c("Group.1","Lon","Lat")]

# bentos
coordenadas_bentos <- aggregate(bentos_subset, 
                         by= list (bentos_subset$Site), 
                         FUN=mean)[c("Group.1","Lon","Lat")]

### BiO Oracle covariates
# Explore datasets in the package
# devtools::install_github("lifewatch/sdmpredictors")
layers <- list_layers()
#View (layers [grep ("Bio-ORACLE",layers$dataset_code),])

# Download specific layers to the current directory
# set prefered folder (to download data)
options(sdmpredictors_datadir=here ("data","environment"))

## chlorophil has different extent - loading and extracting in two steps         
layers_oracle <- load_layers(c("BO2_tempmean_ss","BO2_temprange_ss",
                               "BO2_ppmean_ss", "BO2_pprange_ss",
                               "BO2_salinitymean_ss", "BO2_salinityrange_ss",
                               "BO_damax","BO_damean","BO_damin"))

## these data have different extent
#layers_oracle_Chl <- load_layers (c("BO2_chlomean_ss","BO2_chlorange_ss"))

## coordinates to spatial points
# ajustando uma coordenada
coordenadas_peixes [84,2] <- as.numeric(-35.082658) ## perua preta
#
sp_points <- list(coordenadas_peixes,
                  coordenadas_bentos)

spdf <- lapply (sp_points, function (i) 
  
  SpatialPointsDataFrame(coords = i[,2:3], data = i,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")))

## extracting data

extracted_sea_data <- lapply (spdf, function (i) 
  extract (layers_oracle, i,method='simple', fun=mean))
extracted_sea_data<-lapply (seq(1,length(extracted_sea_data)), function (i) {
  rownames(extracted_sea_data[[i]]) <- sp_points[[i]]$Group.1;
  extracted_sea_data[[i]]
})

#extracted_sea_data_Chl <- extract (layers_oracle_Chl, spdf,method='simple', fun=mean)
#rownames(extracted_sea_data_Chl ) <- sp_points$Group.1

## binding these dfs
#extracted_sea_data <- cbind(extracted_sea_data,

###### distance to the coast areas

# BR coast, download from here https://mapcruzin.com/free-brazil-arcgis-maps-shapefiles.htm

BR <- readOGR(dsn=here("data", "environment","brazil-coastline"), "brazil_coastline")
crs(BR) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
BR <- spTransform(BR, CRS("+init=epsg:4326"))

# use dist2Line from geosphere - only works for WGS84 
#data bentos
sp_data <- SpatialPoints(coordenadas_bentos[,2:3])
crs(sp_data) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
sp_data <- spTransform(sp_data, CRS("+init=epsg:4326"))

# measure distance
dist_bentos <- geosphere::dist2Line(p = sp_data, 
                             line = (BR))

coordenadas_bentos <- cbind (coordenadas_bentos, dist_bentos)

ggplot(data = world) +
  geom_sf() +
  geom_point(data=coordenadas_bentos, 
             aes(x=Lon, y=Lat, col =distance)) +
  coord_sf(xlim = c(-55, -20), ylim = c(-33,0 ), expand = FALSE)


## peixes
sp_data <- SpatialPoints(coordenadas_peixes[,2:3])
crs(sp_data) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
sp_data <- spTransform(sp_data, CRS("+init=epsg:4326"))

# measure distance
dist_peixes <- geosphere::dist2Line(p = sp_data, 
                                    line = (BR))

coordenadas_peixes <- cbind (coordenadas_peixes, dist_peixes)

ggplot(data = world) +
  geom_sf() +
  geom_point(data=coordenadas_peixes, 
             aes(x=Lon, y=Lat, col =distance)) +
  coord_sf(xlim = c(-55, -20), ylim = c(-33,0 ), expand = FALSE)


## lista de covariaveis

covariates_site <- list (region = regiao,
                         site_names = sites_fish_complete,
                         coord = list(coord_bentos=coordenadas_bentos,
                                      coord_peixes = coordenadas_peixes),
                         sea_data = extracted_sea_data)

covariates_effort <- list(effort = tabela_esforco)

###############################################
# # # # # # # # #  SAVE # # # # # # # # # # # #
###############################################

### save data - para os modelos de coral

save (covariates_site , ### dados de covariaveis das bases de dados
      covariates_effort, ## variaveis de esforco
      dados_peixes_bentos, ## dados de composicao 
      file=here ("output","env_data.RData"))
