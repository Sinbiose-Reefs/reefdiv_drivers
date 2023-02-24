Functional diversity patterns of fish, corals and algae in the Brazilian
biogeographic province
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

Repository containing the data and scripts used in the article
“Functional diversity patterns of fish, corals and algae in the
Brazilian biogeographic province”, accepted in Journal of Biogeography
(February 2023).

Occurrence data (folder “data/detection”) are formatted according to the
Darwin Core Standards (<https://www.tdwg.org/standards/dwc/>). Trait
data follow the basic format of trait datasets: species in the rows and
traits in the columns.

The dataset of coral traits (compiled by Jessica Bleuel, PhD candidate
in the Universidade Federal do Rio Grande do Norte, supervised by
Dr. Guilherme O. Longo) has an embargo to be released (up to February
2024).

<!-- badges: start -->
<!-- badges: end -->

#### This paper was produced using the following software and associated packages:

    ## R version 4.2.2 (2022-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Portuguese_Brazil.utf8  LC_CTYPE=Portuguese_Brazil.utf8   
    ## [3] LC_MONETARY=Portuguese_Brazil.utf8 LC_NUMERIC=C                      
    ## [5] LC_TIME=Portuguese_Brazil.utf8    
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggdist_3.2.1            forcats_0.5.2           stringr_1.5.0          
    ##  [4] purrr_1.0.1             tidyr_1.3.0             tibble_3.1.8           
    ##  [7] tidyverse_1.3.2         modeltools_0.2-23       flexmix_2.3-18         
    ## [10] terra_1.7-3             RColorBrewer_1.1-3      mapdata_2.3.1          
    ## [13] plotly_4.10.1           dygraphs_1.1.1.6        xts_0.12.2             
    ## [16] zoo_1.8-11              flexdashboard_0.6.1     lubridate_1.9.1        
    ## [19] rerddap_1.0.1           readr_2.1.3             rasterVis_0.51.5       
    ## [22] viridis_0.6.2           viridisLite_0.4.1       corrplot_0.92          
    ## [25] loo_2.5.1               brms_2.18.0             Rcpp_1.0.10            
    ## [28] effects_4.2-2           carData_3.0-5           MASS_7.3-58.1          
    ## [31] emmeans_1.8.4-1         mgcv_1.8-41             nlme_3.1-160           
    ## [34] leaflet_2.1.1           sdmpredictors_0.2.14    clue_0.3-63            
    ## [37] cluster_2.1.4           FD_1.0-12.1             geometry_0.4.6.1       
    ## [40] ape_5.6-2               ade4_1.7-20             vegan_2.6-4            
    ## [43] lattice_0.20-45         permute_0.9-7           maps_3.4.1             
    ## [46] spdep_1.2-7             sf_1.0-9                spData_2.2.1           
    ## [49] raster_3.6-14           rgdal_1.6-4             abind_1.4-5            
    ## [52] dplyr_1.0.10            openxlsx_4.2.5.1        rgeos_0.6-1            
    ## [55] sp_1.6-0                scatterpie_0.1.8        ggrepel_0.9.2          
    ## [58] gridExtra_2.3           ggplot2_3.4.0           rnaturalearthdata_0.1.0
    ## [61] rnaturalearth_0.3.2     reshape2_1.4.4          reshape_0.8.9          
    ## [64] here_1.0.1             
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] estimability_1.4.1   rappdirs_0.3.3       coda_0.19-4         
    ##   [4] knitr_1.42           data.table_1.14.6    inline_0.3.19       
    ##   [7] generics_0.1.3       callr_3.7.3          proxy_0.4-27        
    ##  [10] tzdb_0.3.0           xml2_1.3.3           httpuv_1.6.8        
    ##  [13] wk_0.7.1             StanHeaders_2.26.13  assertthat_0.2.1    
    ##  [16] gargle_1.2.1         xfun_0.36            hms_1.1.2           
    ##  [19] jquerylib_0.1.4      bayesplot_1.10.0     evaluate_0.20       
    ##  [22] promises_1.2.0.1     fansi_1.0.4          readxl_1.4.1        
    ##  [25] dbplyr_2.3.0         igraph_1.3.5         DBI_1.1.3           
    ##  [28] htmlwidgets_1.6.1    tensorA_0.36.2       googledrive_2.0.0   
    ##  [31] ellipsis_0.3.2       crosstalk_1.2.0      backports_1.4.1     
    ##  [34] V8_4.2.2             insight_0.19.0       survey_4.1-1        
    ##  [37] markdown_1.4         RcppParallel_5.1.6   deldir_1.0-6        
    ##  [40] vctrs_0.5.2          cachem_1.0.6         withr_2.5.0         
    ##  [43] ggforce_0.4.1        checkmate_2.1.0      prettyunits_1.1.1   
    ##  [46] lazyeval_0.2.2       crayon_1.5.2         crul_1.3            
    ##  [49] pkgconfig_2.0.3      units_0.8-1          tweenr_2.0.2        
    ##  [52] nnet_7.3-18          rlang_1.0.6          lifecycle_1.0.3     
    ##  [55] miniUI_0.1.1.1       colourpicker_1.2.0   httpcode_0.3.0      
    ##  [58] modelr_0.1.10        cellranger_1.1.0     distributional_0.3.1
    ##  [61] rprojroot_2.0.3      polyclip_1.10-4      matrixStats_0.63.0  
    ##  [64] Matrix_1.5-1         boot_1.3-28          reprex_2.0.2        
    ##  [67] base64enc_0.1-3      processx_3.8.0       googlesheets4_1.0.1 
    ##  [70] png_0.1-8            KernSmooth_2.23-20   classInt_0.4-8      
    ##  [73] s2_1.1.2             jpeg_0.1-10          shinystan_2.6.0     
    ##  [76] scales_1.2.1         magrittr_2.0.3       plyr_1.8.8          
    ##  [79] hexbin_1.28.2        threejs_0.3.3        compiler_4.2.2      
    ##  [82] rstantools_2.2.0     lme4_1.1-31          cli_3.6.0           
    ##  [85] ps_1.7.2             Brobdingnag_1.2-9    hoardr_0.5.3        
    ##  [88] magic_1.6-1          tidyselect_1.2.0     stringi_1.7.12      
    ##  [91] mitools_2.4          yaml_2.3.7           svUnit_1.0.6        
    ##  [94] latticeExtra_0.6-30  bridgesampling_1.1-2 grid_4.2.2          
    ##  [97] sass_0.4.5           tools_4.2.2          timechange_0.2.0    
    ## [100] rstudioapi_0.14      posterior_1.3.1      farver_2.1.1        
    ## [103] digest_0.6.31        shiny_1.7.4          broom_1.0.3         
    ## [106] later_1.3.0          ncdf4_1.21           httr_1.4.4          
    ## [109] colorspace_2.1-0     rvest_1.0.3          fs_1.6.0            
    ## [112] splines_4.2.2        shinythemes_1.2.0    xtable_1.8-4        
    ## [115] jsonlite_1.8.4       nloptr_2.0.3         rstan_2.26.13       
    ## [118] ggfun_0.0.9          R6_2.5.1             pillar_1.8.1        
    ## [121] htmltools_0.5.4      mime_0.12            glue_1.6.2          
    ## [124] fastmap_1.1.0        minqa_1.2.5          DT_0.27             
    ## [127] class_7.3-20         codetools_0.2-18     pkgbuild_1.4.0      
    ## [130] mvtnorm_1.1-3        utf8_1.2.2           bslib_0.4.2         
    ## [133] arrayhelpers_1.1-0   curl_5.0.0           gtools_3.9.4        
    ## [136] zip_2.2.2            shinyjs_2.1.0        interp_1.1-3        
    ## [139] survival_3.4-0       rmarkdown_2.20       munsell_0.5.0       
    ## [142] e1071_1.7-12         haven_2.5.1          gtable_0.3.1
