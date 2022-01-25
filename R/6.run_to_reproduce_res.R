
## --------------------
#       FIGURES
# --------------------

## call packages & functions
source("R/packages.R")
source("R/functions.R")

# define the number of iterations for rarefaction
niter <- 1000

# ncores for parallel processing
ncores <- 3

# MCMC settings (for models)
ni <- 10000 
nb <- 8000
nt <- 4
nc <- ncores

# code for organizing data
source ("R/1.codigo_organizacao_dados.R")

# code for functional analyses
source ("R/2.func_reefs.R")

# code for organize richness and functional data to modeling
source ("R/3.organize_data_to_modeling.R")

# code for modeling (multivariate linear models, GLMs)
source ("R/4.modelingGLM_rarefied.R")

# code for producing plots
source ("R/5.graphics.R")

