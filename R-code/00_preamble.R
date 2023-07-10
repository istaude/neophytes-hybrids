
# packages ----------------------------------------------------------------

library(Matrix)
library(tidyverse)
library(data.table)
library(readxl)
library(readr)
library(sf)
library(RColorBrewer)
library(tmap)
library(effects)
library(readxl)
library(ggchicklet)
library(tidyverse)
library(cowplot)
library(hrbrthemes)
library(extrafont)
library(grid)
library(ggpubr)
library(stringr)
library(see)
library(emmeans)
library(data.table)
library(lme4)
library(stringr)
library(rnrfa)
library(EnvStats)
library(terra)
library(sf)
library(brms)
library(units)
library(WDI)
library(ggeffects)
library(patchwork)
library(tidybayes)    
library(broom)        
library(broom.mixed) 
library(V.PhyloMaker) # devtools::install_github("jinyizju/V.PhyloMaker")
library(gghalves)
library(ggridges)
library(lemon)
library(ggdist)
library(patchwork)
library(tidyr)
library(modelr)
library(paletteer)
library(ape)
library(treeio) # devtools::install_github("YuLab-SMU/treeio")
library(nodiv)
library(geiger)
library(phytools)
library(ggtree) # BiocManager::install("ggtree")
library(ggtreeExtra) # BiocManager::install("ggtreeExtra")
library(phyloseq) # BiocManager::install("phyloseq")
library(reshape2)
library(ggsci)
library(ggnewscale)
library(viridis)
library(caper)
library(rWCVP)
library(rWCVPdata)

import_public_sans()



# functions ---------------------------------------------------------------

pretty_labels <- function(x, digs = 0){
  
a <- sapply(x, 
       
  FUN = function(x, digits = digs){
  
  a <- gsub("[^0-9.-]", " ", x)
  a <- trimws(a)
  a1 <- round(as.numeric(str_split(a, " ")[[1]][1]), digits)
  a2 <- round(as.numeric(str_split(a, " ")[[1]][2]), digits)
  y <- paste(prettyNum(a1, big.mark=","), prettyNum(a2, big.mark=","), sep = " to ")
  
  y <- ifelse(is.na(x) ==  TRUE, "Missing", y)
  return(y)
}

)

factor(a,  levels = unique(a[order(x)]))
  
}



