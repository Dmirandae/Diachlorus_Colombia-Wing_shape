## Ambrosio Torres & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Research: Wing shape variation in the taxonomic recognition of species of Diachlorus Osten-Sacken (Diptera: Tabanidae) from Colombia
## Part: Mantel test for intraspecific variation Vs. geographic distances
## R version 3.1.2 & Rstudio 0.97.551
## 28th January 2015 

library(vegan)                                             ## Load 'vegan' package v. 2.0-10 ---- http://cran.r-project.org/web/packages/vegan/index.html

## D. curvipes
curvi_geogra <- read.table("curvipes_geogra.txt",header=T) ## Load the geographical distances among the areas where D. curvipes is distributed
curvi_geomet <- read.table("curvipes_geomet.txt",header=T) ## Load the geometrical distances among specimens of D. curvipes distributed in different areas
mantel(curvi_geogra, curvi_geomet)
# Mantel statistic based on Pearson's product-moment correlation 
# Call:
# mantel(xdis = curvi_geogra, ydis = curvi_geomet) 
# Mantel statistic r: 0.2771 
#       Significance: 0.138 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.277 0.345 0.359 0.359 
# Based on 999 permutations

## D. fuscistigma
fusci_geogra <- read.table("fuscistigma_geogra.txt",header=T) ## Load the geographical distances among the areas where D. fuscistigma is distributed
fusci_geomet <- read.table("fuscistigma_geomet.txt",header=T) ## Load the geometrical distances among specimens of D. fuscistigma distributed in different areas
mantel(fusci_geogra, fusci_geomet)
# Mantel statistic based on Pearson's product-moment correlation 
# Call:
# mantel(xdis = fusci_geogra, ydis = fusci_geomet) 
# Mantel statistic r: 0.05199 
#      Significance: 0.447 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.398 0.429 0.465 0.539 
# Based on 999 permutations

## D. jobbinsi
jobi_geogra <- read.table("jobbinsi_geogra.txt",header=T) ## Load the geographical distances among the areas where D. jobbinsi is distributed
jobi_geomet <- read.table("jobbinsi_geomet.txt",header=T) ## Load the geometrical distances among specimens of D. jobbinsi distributed in different areas
mantel(jobi_geogra, jobi_geomet)
# Mantel statistic based on Pearson's product-moment correlation 
# Call:
# mantel(xdis = jobi_geogra, ydis = jobi_geomet) 
# Mantel statistic r: 0.5883 
#       Significance: 0.217 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.624 0.657 0.779 0.779 
# Based on 999 permutations

## D. leucotibialis
leuco_geogra <- read.table("leucotibialis_geogra.txt",header=T) ## Load the geographical distances among the areas where D. leucotibialis is distributed
leuco_geomet <- read.table("leucotibialis_geomet.txt",header=T) ## Load the geometrical distances among specimens of D. leucotibialis distributed in different areas
mantel(leuco_geogra, leuco_geomet)
# Mantel statistic based on Pearson's product-moment correlation 
# Call:
# mantel(xdis = leuco_geogra, ydis = leuco_geomet) 
# Mantel statistic r: -0.4829 
#       Significance: 0.667 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
#     1     1     1     1 
# Based on 999 permutations
