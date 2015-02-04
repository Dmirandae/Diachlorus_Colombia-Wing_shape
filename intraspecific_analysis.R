## Ambrosio Torres & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Research: Geometric wing variation in the taxonomic recognition of species of the genus Diachlorus Osten-Sacken (Diptera: Tabanidae) from Colombia
## Part: intraspecific variation of genus Diachlorus
## R version 3.1.2 & Rstudio 0.97.551
## 28th January 2015 

######################################################################################################################################################################################
######################################################################################################################################################################################
## 1. LOAD PACKAGES AND FUNCTIONS TO THE ANALYSIS

library(colorspace)  ## Load 'colorspace' package v. 1.2-4 --- http://cran.r-project.org/web/packages/colorspace/index.html
library(cluster)     ## Load 'cluster' package v. 1.15.2 ---- http://cran.r-project.org/web/packages/cluster/index.html
library(DiscriMiner) ## Load 'Discriminer' package v. 0.1-29 ---- http://cran.r-project.org/web/packages/DiscriMiner/index.html
library(ellipse)     ## Load 'ellipse' package v. 0.3-8 ---- http://cran.r-project.org/web/packages/ellipse/index.html
library(geomorph)    ## Load 'geomorph' package v. 2.1.1 ---  http://cran.r-project.org/web/packages/geomorph/index.html
library(ggplot2)     ## Load 'ggplot2' package v. 1.0.0 ---- http://cran.r-project.org/web/packages/ggplot2/index.html
library(mclust)      ## Load 'mclust' package v. 4.3 ---- http://cran.r-project.org/web/packages/mclust/index.html
library(NbClust)     ## Load 'NbClust' package V. 2.0.3 ---- http://cran.r-project.org/web/packages/NbClust/index.html
library(shapes)      ## Load 'shapes' package v. 1.1-9 --- http://cran.r-project.org/web/packages/shapes/index.html
library(vegan)       ## Load 'vegan' package v. 2.0-10 ---- http://cran.r-project.org/web/packages/vegan/index.html
source("colLab.R")   ## Load "colLab" function, modified by us to coloring edgePar and labels of a dendrogram ---- https://stat.ethz.ch/R-manual/R-patched/library/stats/html/dendrapply.html 
                     ## and http://stackoverflow.com/

######################################################################################################################################################################################
######################################################################################################################################################################################
## 2. DEFINE THE GROUPS OF STUDY TO MAKE THE COMPARISONS IN THE ANALYSIS AND READ THE DATA (BUILD ARRAYS)

list_curvi <- list.files(pattern = "dc")     ## Make a list with all the names of landmarks-files of Diachlorus curvipes
list_fusci <- list.files(pattern = "df")     ## Make a list with all the names of landmarks-files of Diachlorus fuscistigma
list_jobi <- list.files(pattern = "dj")      ## Make a list with all the names of landmarks-files of Diachlorus jobbinsi
list_leuco <- list.files(pattern = "dl")     ## Make a list with all the names of landmarks-files of Diachlorus leucotibialis
list_nune <- list.files(pattern = "dn")      ## Make a list with all the names of landmarks-files of Diachlorus nuneztovari
list_leti <- list.files(pattern = "dt")      ## Make a list with all the names of landmarks-files of Diachlorus leticia
list_adc <- list.files(pattern = "adc")      ## Make a list with all the names of landmarks-files of Diachlorus curvipes of Amazonas
list_adf <- list.files(pattern = "adf")      ## Make a list with all the names of landmarks-files of Diachlorus fuscistigma of Amazonas
list_adn <- list.files(pattern = "adn")      ## Make a list with all the names of landmarks-files of Diachlorus nuneztovari of Amazonas
list_adl <- list.files(pattern = "adl")      ## Make a list with all the names of landmarks-files of Diachlorus leucotibialis of Amazonas
list_adj <- list.files(pattern = "adj")      ## Make a list with all the names of landmarks-files of Diachlorus jobbinsi of Amazonas
list_cdf <- list.files(pattern = "cdf")      ## Make a list with all the names of landmarks-files of Diachlorus fuscistigma of Caquetá
list_cdl <- list.files(pattern = "cdl")      ## Make a list with all the names of landmarks-files of Diachlorus leucotibialis of Caquetá
list_cdj <- list.files(pattern = "cdj")      ## Make a list with all the names of landmarks-files of Diachlorus jobbinsi of Caquetá
list_chdc <- list.files(pattern = "chdc")    ## Make a list with all the names of landmarks-files of Diachlorus curvipes of Chocó
list_chdj <- list.files(pattern = "chdj")    ## Make a list with all the names of landmarks-files of Diachlorus jobbinsi of Chocó
list_mdc <- list.files(pattern = "mdc")      ## Make a list with all the names of landmarks-files of Diachlorus curvipes of Meta
list_mdf <- list.files(pattern = "mdf")      ## Make a list with all the names of landmarks-files of Diachlorus fuscistigma of Meta
list_pdc <- list.files(pattern = "pdc")      ## Make a list with all the names of landmarks-files of Diachlorus curvipes of Putumayo
list_pdf <- list.files(pattern = "pdf")      ## Make a list with all the names of landmarks-files of Diachlorus fuscistigma of Putumayo
list_pdn <- list.files(pattern = "pdn")      ## Make a list with all the names of landmarks-files of Diachlorus nuneztovari of Putumayo
list_pdl <- list.files(pattern = "pdl")      ## Make a list with all the names of landmarks-files of Diachlorus leucotibialis of Putumayo
list_pdj <- list.files(pattern = "pdj")      ## Make a list with all the names of landmarks-files of Diachlorus jobbinsi of Putumayo
list_vdf <- list.files(pattern = "vdf")      ## Make a list with all the names of landmarks-files of Diachlorus fuscistigma of Vaupés
list_vdt <- list.files(pattern = "vdt")      ## Make a list with all the names of landmarks-files of Diachlorus leticia of Vaupés

d_curvipes <- readmulti.nts(list_curvi)      ## Build an array/object which contains all individuals of Diachlorus curvipes
d_fuscistigma <- readmulti.nts(list_fusci)   ## Build an array/object which contains all individuals of Diachlorus fuscistigma
d_jobbinsi <- readmulti.nts(list_jobi)       ## Build an array/object which contains all individuals of Diachlorus jobbinsi
d_leucotibialis <- readmulti.nts(list_leuco) ## Build an array/object which contains all individuals of Diachlorus leucotibialis
d_leticia <- readmulti.nts(list_leti)        ## Build an array/object which contains all individuals of Diachlorus leticia
d_nuneztovari <- readmulti.nts(list_nune)    ## Build an array/object which contains all individuals of Diachlorus nuneztovari
adc <- readmulti.nts(list_adc)               ## Build an array/object which contains all individuals of Diachlorus curvipes of Amazonas 
adf <- readmulti.nts(list_adf)               ## Build an array/object which contains all individuals of Diachlorus fuscistigma of Amazonas
adn <- readmulti.nts(list_adn)               ## Build an array/object which contains all individuals of Diachlorus nuneztovari of Amazonas
adl <- readmulti.nts(list_adl)               ## Build an array/object which contains all individuals of Diachlorus leucotibialis of Amazonas
adj <- readmulti.nts(list_adj)               ## Build an array/object which contains all individuals of Diachlorus jobbinsi of Amazonas
cdf <- readmulti.nts(list_cdf)               ## Build an array/object which contains all individuals of Diachlorus fuscistigma of Caquetá
cdl <- readmulti.nts(list_cdl)               ## Build an array/object which contains all individuals of Diachlorus leucotibialis of Caquetá
cdj <- readmulti.nts(list_cdj)               ## Build an array/object which contains all individuals of Diachlorus jobbinsi of Caquetá
chdc <- readmulti.nts(list_chdc)             ## Build an array/object which contains all individuals of Diachlorus curvipes of Chocó
chdj <- readmulti.nts(list_chdj)             ## Build an array/object which contains all individuals of Diachlorus jobbinsi of Chocó
mdc <- readmulti.nts(list_mdc)               ## Build an array/object which contains all individuals of Diachlorus curvipes of Meta
mdf <- readmulti.nts(list_mdf)               ## Build an array/object which contains all individuals of Diachlorus fuscistigma of Meta
pdc <- readmulti.nts(list_pdc)               ## Build an array/object which contains all individuals of Diachlorus curvipes of Putumayo
pdf <- readmulti.nts(list_pdf)               ## Build an array/object which contains all individuals of Diachlorus fuscistigma of Putumayo
pdn <- readmulti.nts(list_pdn)               ## Build an array/object which contains all individuals of Diachlorus nuneztovari of Putumayo
pdl <- readmulti.nts(list_pdl)               ## Build an array/object which contains all individuals of Diachlorus leucotibialis of Putumayo
pdj <- readmulti.nts(list_pdj)               ## Build an array/object which contains all individuals of Diachlorus jobbinsi of Putumayo
vdf <- readmulti.nts(list_vdf)               ## Build an array/object which contains all individuals of Diachlorus fuscistigma of Vaupés
vdt <- readmulti.nts(list_vdt)               ## Build an array/object which contains all individuals of Diachlorus leticia of Vaupés

######################################################################################################################################################################################
######################################################################################################################################################################################
## 3. CREATE A DIRECTORY TO PUT THE RESULTS IN A DIFERENT DIRECTORY - NOT IN THE DIRECTORY OF TE DATA (INPUT)

namedir <- paste("intraspecific_results_", date(),  sep="") ## Create a name to the new directory
dir.create(namedir, showWarnings = FALSE)                   ## Create the new directory 
setwd(paste((getwd()),  "/" ,namedir, sep=""))              ## Change the directory to the new created directory

######################################################################################################################################################################################
######################################################################################################################################################################################
## 4. GENERALIZED PROCRUSTES ANALYSIS (GPA) AND CONSESUS SHAPES

groups_intra <- list(d_curvipes=d_curvipes, d_fuscistigma=d_fuscistigma, d_jobbinsi=d_jobbinsi, d_leucotibialis=d_leucotibialis, 
                     d_leticia=d_leticia, d_nuneztovari=d_nuneztovari, adc=adc, chdc=chdc, mdc=mdc, pdc=pdc, adf=adf, cdf=cdf, 
                     mdf=mdf, pdf=pdf, vdf=cdf, adj=adj, cdj=cdj, chdj=chdj, pdj=pdj, adl=adl, cdl=cdl, pdl=pdl, adn=adn, pdn=pdn, 
                     vdt=vdt)                                                                                                         ## Make a list of the arrays/groups
alg_intra <- lapply (groups_intra, gpagen)                                                                                            ## Align all the groups in the list
alg_d_curvipes <- alg_intra$d_curvipes$coords ; writeland.tps(alg_d_curvipes, "alg_d_curvipes.tps")                                   ## Write the aligned data in .tps format to D. curvipes data
alg_d_fuscistigma <- alg_intra$d_fuscistigma$coords ; writeland.tps(alg_d_fuscistigma, "alg_d_fuscistigma.tps")                       ## Write the aligned data in .tps format to D. fuscistigma data
alg_d_jobbinsi <- alg_intra$d_jobbinsi$coords ; writeland.tps(alg_d_jobbinsi, "alg_d_jobbinsi.tps")                                   ## Write the aligned data in .tps format to D. jobbinsi data
alg_d_leucotibialis <- alg_intra$d_leucotibialis$coords ; writeland.tps(alg_d_leucotibialis, "alg_d_leucotibialis.tps")               ## Write the aligned data in .tps format to D. leucotibialis data
alg_d_leticia <- alg_intra$d_leticia$coords ; writeland.tps(alg_d_leticia, "alg_d_leticia.tps")                                       ## Write the aligned data in .tps format to D. leticia data
alg_d_nuneztovari <- alg_intra$d_nuneztovari$coords ; writeland.tps(alg_d_nuneztovari, "alg_d_nuneztovari.tps")                       ## Write the aligned data in .tps format to D. nuneztovari data
alg_adc <- alg_intra$adc$coords ; writeland.tps(alg_adc, "alg_amazonas_d_curvipes.tps")                                               ## Write the aligned data in .tps format to D. curvipes data from Amazonas
alg_adf <- alg_intra$adf$coords ; writeland.tps(alg_adf, "alg_amazonas_d_fuscistigma.tps")                                            ## Write the aligned data in .tps format to D. fuscistigma data from Amazonas
alg_adn <- alg_intra$adn$coords ; writeland.tps(alg_adn, "alg_amazonas_d_nuneztovari.tps")                                            ## Write the aligned data in .tps format to D. nuneztovari data from Amazonas
alg_adl <- alg_intra$adl$coords ; writeland.tps(alg_adl, "alg_amazonas_d_leucotibialis.tps")                                          ## Write the aligned data in .tps format to D. leucotibialis data from Amazonas
alg_adj <- alg_intra$adj$coords ; writeland.tps(alg_adj, "alg_amazonas_d_jobbinsi.tps")                                               ## Write the aligned data in .tps format to D. jobbinsi data from Amazonas
alg_cdf <- alg_intra$cdf$coords ; writeland.tps(alg_cdf, "alg_caqueta_d_fuscistigma.tps")                                             ## Write the aligned data in .tps format to D. fuscistigma data from Caqueta
alg_cdl <- alg_intra$cdl$coords ; writeland.tps(alg_cdl, "alg_caqueta_d_leucotibialis.tps")                                           ## Write the aligned data in .tps format to D. leucotibialis data from Caqueta
alg_cdj <- alg_intra$cdj$coords ; writeland.tps(alg_cdj, "alg_caqueta_d_jobbinsi.tps")                                                ## Write the aligned data in .tps format to D. jobbinsi data from Caqueta
alg_chdc <- alg_intra$chdc$coords ; writeland.tps(alg_chdc, "alg_choco_d_curvipes.tps")                                               ## Write the aligned data in .tps format to D. curvipes data from Chocó
alg_chdj <- alg_intra$chdj$coords ; writeland.tps(alg_chdj, "alg_choco_d_jobbinsi.tps")                                               ## Write the aligned data in .tps format to D. jobbinsi data from Chocó
alg_mdc <- alg_intra$mdc$coords ; writeland.tps(alg_mdc, "alg_meta_d_curvipes.tps")                                                   ## Write the aligned data in .tps format to D. curvipes data from Meta
alg_mdf <- alg_intra$mdf$coords ; writeland.tps(alg_mdf, "alg_meta_d_fuscistigma.tps")                                                ## Write the aligned data in .tps format to D. fuscistigma data from Meta
alg_pdc <- alg_intra$pdc$coords ; writeland.tps(alg_pdc, "alg_putumayo_d_curvipes.tps")                                               ## Write the aligned data in .tps format to D. curvipes data from Putumayo
alg_pdf <- alg_intra$pdf$coords ; writeland.tps(alg_pdf, "alg_putumayo_d_fuscistigma.tps")                                            ## Write the aligned data in .tps format to D. fuscistigma data from Putumayo
alg_pdn <- alg_intra$pdn$coords ; writeland.tps(alg_pdn, "alg_putumayo_d_nuneztovari.tps")                                            ## Write the aligned data in .tps format to D. nuneztovari data from Putumayo
alg_pdl <- alg_intra$pdl$coords ; writeland.tps(alg_pdl, "alg_putumayo_d_leucotibialis.tps")                                          ## Write the aligned data in .tps format to D. leucotibialis data from Putumayo
alg_pdj <- alg_intra$pdj$coords ; writeland.tps(alg_pdj, "alg_putumayo_d_jobbinsi.tps")                                               ## Write the aligned data in .tps format to D. jobbinsi data from Putumayo
alg_vdf <- alg_intra$vdf$coords ; writeland.tps(alg_vdf, "alg_vaupes_d_fuscistigma.tps")                                              ## Write the aligned data in .tps format to D. fuscistigma data from Vaupés
alg_vdt <- alg_intra$vdt$coords ; writeland.tps(alg_vdt, "alg_vaupes_d_leticia.tps")                                                  ## Write the aligned data in .tps format to D. leticia data from Vaupés

cons_names <- list(alg_d_curvipes=alg_d_curvipes, alg_d_fuscistigma=alg_d_fuscistigma, alg_d_jobbinsi=alg_d_jobbinsi, 
                   alg_d_leucotibialis=alg_d_leucotibialis, alg_d_leticia=alg_d_leticia, alg_d_nuneztovari=alg_d_nuneztovari,
                   alg_adc=alg_adc, alg_chdc=alg_chdc, alg_mdc=alg_mdc, alg_pdc=alg_pdc, alg_adf=alg_adf, alg_cdf=alg_cdf, 
                   alg_mdf=alg_mdf, alg_pdf=alg_pdf, alg_vdf=alg_cdf, alg_adj=alg_adj, alg_cdj=alg_cdj, alg_chdj=alg_chdj, 
                   alg_pdj=alg_pdj, alg_adl=alg_adl, alg_cdl=alg_cdl, alg_pdl=alg_pdl, alg_adn=alg_adn, alg_pdn=alg_pdn, 
                   alg_vdt=alg_vdt)                                                                                                   ## Make a list of the aligned arrays/groups
cons_shapes <- lapply (cons_names, mshape)                                                                                            ## Calculate the consensus shape of each array/group

######################################################################################################################################################################################
######################################################################################################################################################################################
## 5. Principal Component Analysis (PCA)

### Creating a factor to color and recognize the Diachlorus curvipes from each national natural park to do a PCA 
factors_d_curvi <-sub(".nts", "", list_curvi)
for (i in (rep(c(0,1,2,3,4,5,6,7,8,9),5))) {
     factors_d_curvi <- sub(i, "", factors_d_curvi)
}
factors_d_curvi <- factor(factors_d_curvi)                       
pca_d_curvipes <- plotTangentSpace(alg_d_curvipes, label=NULL, verbose =T, groups = factors_d_curvi, warpgrids=F) ## Perform PCA of Diachlorus curvipes
## Diachlorus curvipes PCA maked in ggplot2
label_dc <- list(expression('Amazonas'), expression('Chocó'), expression('Meta'), expression('Putumayo'))   
a <- cbind(pca_d_curvipes$pc.scores, factors_d_curvi)
a <- data.frame(a)
a$factors_d_curvi <- as.factor(a$factors_d_curvi)
df_ell <- data.frame()
for(g in levels(a$factors_d_curvi)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$factors_d_curvi==g,], ellipse(cor(PC1, PC2), 
                    scale=c(sd(PC1),sd(PC2)), 
                    centre=c(mean(PC1),mean(PC2))))),factors_d_curvi=g))
}
p_curvi <- ggplot(data=a, aes(x=PC1, y=PC2,colour=factors_d_curvi)) + geom_point(size=2, shape= 18) 
p_curvi <- p_curvi + geom_path(data=df_ell, aes(x=x, y=y,colour=factors_d_curvi), size=0.5, linetype=2) 
p_curvi <- p_curvi + scale_color_manual(name ="Areas", labels= label_dc,values=c("#1b9e77", "darkorange4", "#7570b3", "#e7298a")) 
p_curvi <- p_curvi + scale_shape_manual(name ="Areas", labels= label_dc, values=c(15, 16, 17,18)) 

pdf(file="pca_d_curvipes_ggplot2.pdf", width = 10, height = 6) 
p_curvi <- p_curvi + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) 
p_curvi + theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size=16))
dev.off()

### Creating a factor to color and recognize the Diachlorus fuscistigma from each national natural park to do a PCA 
factors_d_fusci <-sub(".nts", "", list_fusci)
for (i in (rep(c(0,1,2,3,4,5,6,7,8,9),5))) {
  factors_d_fusci <-sub(i, "", factors_d_fusci)
}
factors_d_fusci<-factor(factors_d_fusci)
pca_d_fuscistigma <- plotTangentSpace(alg_d_fuscistigma, label=NULL, verbose =T, groups = factors_d_fusci, warpgrids=F) ## Perform PCA of Diachlorus curvipes
## Diachlorus fuscistigma PCA maked in ggplot2
label_df <- list(expression('Amazonas'), expression('Caquetá'), expression('Meta'), expression('Putumayo'), expression('Vaupés') ) 
a<- cbind(pca_d_fuscistigma$pc.scores, factors_d_fusci)
a<-data.frame(a)
a$factors_d_fusci <- as.factor(a$factors_d_fusci)
df_ell <- data.frame()
for(g in levels(a$factors_d_fusci)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$factors_d_fusci==g,], ellipse(cor(PC1, PC2), 
                    scale=c(sd(PC1),sd(PC2)), 
                    centre=c(mean(PC1),mean(PC2))))),factors_d_fusci=g))
}
p_fusci <- ggplot(data=a, aes(x=PC1, y=PC2,colour=factors_d_fusci)) + geom_point(size=2, shape= 15) 
p_fusci <- p_fusci + geom_path(data=df_ell, aes(x=x, y=y,colour=factors_d_fusci), size=0.5, linetype=2) 
p_fusci <- p_fusci + scale_color_manual(name ="Areas", labels= label_df,values=c("#1b9e77", "royalblue4", "#7570b3", "#e7298a", "#e6ab02")) 
p_fusci <- p_fusci + scale_shape_manual(name ="Areas", labels= label_df, values=c(15, 7, 17, 18, 10)) 

pdf(file="pca_d_fuscistigma_ggplot2.pdf", width = 10, height = 6) 
p_fusci + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) + theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size=16))
dev.off()

### Creating a factor to color and recognize the Diachlorus jobbinsi from each national natural park to do a PCA 
factors_d_jobi <-sub(".nts", "", list_jobi)
for (i in (rep(c(0,1,2,3,4,5,6,7,8,9),5))) {
  factors_d_jobi <-sub(i, "", factors_d_jobi)
}
factors_d_jobi<-factor(factors_d_jobi)
pca_d_jobbinsi <- plotTangentSpace(alg_d_jobbinsi, label=NULL, verbose =T, groups = factors_d_jobi, warpgrids=F) ## Perform PCA of Diachlorus curvipes
## Diachlorus jobbinsi PCA maked in ggplot2
label_dj <- list(expression('Amazonas'), expression('Caquetá'), expression('Chocó'), expression('Putumayo')) 
a<- cbind(pca_d_jobbinsi$pc.scores, factors_d_jobi)
a<-data.frame(a)
a$factors_d_jobi <- as.factor(a$factors_d_jobi)
df_ell <- data.frame()
for(g in levels(a$factors_d_jobi)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$factors_d_jobi==g,], ellipse(cor(PC1, PC2), 
                    scale=c(sd(PC1),sd(PC2)), 
                    centre=c(mean(PC1),mean(PC2))))),factors_d_jobi=g))
}
p_jobi <- ggplot(data=a, aes(x=PC1, y=PC2,colour=factors_d_jobi)) + geom_point(size=2, shape= 16) 
p_jobi <- p_jobi + geom_path(data=df_ell, aes(x=x, y=y,colour=factors_d_jobi), size=0.5, linetype=2) 
p_jobi <- p_jobi + scale_color_manual(name ="Areas", labels= label_dj,values=c("#1b9e77", "royalblue4", "darkorange4", "#e7298a")) 
p_jobi <- p_jobi + scale_shape_manual(name ="Areas", labels= label_dj, values=c(15, 7, 16, 18)) 

pdf(file="pca_d_jobbinsi_ggplot2.pdf", width = 10, height = 6) 
p_jobi + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) + theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size=16))
dev.off()

### creating a factor to color and recognize the Diachlorus leucotibialis from each national natural park to do a PCA 
factors_d_leuco <-sub(".nts", "", list_leuco)
for (i in (rep(c(0,1,2,3,4,5,6,7,8,9),5))) {
  factors_d_leuco <-sub(i, "", factors_d_leuco)
}
factors_d_leuco<-factor(factors_d_leuco)
pca_d_leucotibialis <- plotTangentSpace(alg_d_leucotibialis, label=NULL, verbose =T, groups = factors_d_leuco, warpgrids=F)
## Diachlorus leucotibialis PCA maked in ggplot2
label_dl <- list(expression('Amazonas'), expression('Caquetá'), expression('Putumayo')) 
a<- cbind(pca_d_leucotibialis$pc.scores, factors_d_leuco)
a<-data.frame(a)
a$factors_d_leuco <- as.factor(a$factors_d_leuco)
df_ell <- data.frame()
for(g in levels(a$factors_d_leuco)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$factors_d_leuco==g,], ellipse(cor(PC1, PC2), 
                    scale=c(sd(PC1),sd(PC2)), 
                    centre=c(mean(PC1),mean(PC2))))),factors_d_leuco=g))
}
p_leuco <- ggplot(data=a, aes(x=PC1, y=PC2,colour=factors_d_leuco)) + geom_point(size=2, shape=17) 
p_leuco <- p_leuco + geom_path(data=df_ell, aes(x=x, y=y,colour=factors_d_leuco), size=0.5, linetype=2) 
p_leuco <- p_leuco + scale_color_manual(name ="Areas", labels= label_dl,values=c("#1b9e77", "royalblue4", "#e7298a")) 
p_leuco <- p_leuco + scale_shape_manual(name ="Areas", labels= label_dl, values=c(15, 7, 18)) 

pdf(file="pca_d_leucotibialis_ggplot2.pdf", width = 10, height = 6) 
p_leuco + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) + theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size=16))
dev.off()

### creating a factor to color and recognize the Diachlorus nuneztovari from each national natural park to do a PCA 
factors_d_nune <-sub(".nts", "", list_nune)
for (i in (rep(c(0,1,2,3,4,5,6,7,8,9),5))) {
  factors_d_nune <-sub(i, "", factors_d_nune)
}
factors_d_nune<-factor(factors_d_nune)
pca_d_nuneztovari <- plotTangentSpace(alg_d_nuneztovari, label=NULL, verbose =T, groups = factors_d_nune, warpgrids=F)
## Diachlorus nuneztovari PCA maked in ggplot2
label_dm <- list(expression('Amazonas'), expression('Putumayo'))   
a<- cbind(pca_d_nuneztovari$pc.scores, factors_d_nune)
a<-data.frame(a)
a$factors_d_nune <- as.factor(a$factors_d_nune)
df_ell <- data.frame()
for(g in levels(a$factors_d_nune)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$factors_d_nune==g,], ellipse(cor(PC1, PC2), 
                    scale=c(sd(PC1),sd(PC2)), 
                    centre=c(mean(PC1),mean(PC2))))),factors_d_nune=g))
}
p_nune <- ggplot(data=a, aes(x=PC1, y=PC2,colour=factors_d_nune, shape=factors_d_nune)) + geom_point(size=3, shape= 8) 
p_nune <- p_nune + geom_path(data=df_ell, aes(x=x, y=y,colour=factors_d_nune), size=0.5, linetype=2) 
p_nune <- p_nune + scale_color_manual(name ="Areas", labels= label_dm,values=c("#1b9e77", "#e7298a")) 
p_nune <- p_nune + scale_shape_manual(name ="Areas", labels= label_dm, values=c(15, 18)) 

pdf(file="pca_d_nuneztovari_ggplot2.pdf", width = 10, height = 6) 
p_nune + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) + theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size=16))
dev.off()

######################################################################################################################################################################################
######################################################################################################################################################################################
## 6. Hierarchical clustering
                         
datavalores_jobbinsi <- data.frame(pca_d_jobbinsi$pc.scores)            ## Make a data frame using the scores of the PCA of D. jobbinsi
distvalores_jobbinsi <- dist(datavalores_jobbinsi)                      ## Make a distance matrix of "datavalores_jobbinsi"                      
clust_valores_jobbinsi <- hclust(distvalores_jobbinsi)                  ## Perform the hierarquical cluster of D. jobbinsi
labelColors <- c("#1b9e77", "royalblue4", "darkorange4", "#e7298a")     ## Define the color of the leaves of the dendrogram
clusMember <- rep(1,length(rownames(datavalores_jobbinsi)))
clusMember[grep("a",rownames(datavalores_jobbinsi))]<-1 ##1b9e77        ## Define the color of Amazonas to #1b9e77
clusMember[grep("c",rownames(datavalores_jobbinsi))]<-2 #royalblue4     ## Define the color of Caquetá to royalblue4
clusMember[grep("ch",rownames(datavalores_jobbinsi))]<-3 #darkorange4   ## Define the color of Chocó to darkorange4
clusMember[grep("p",rownames(datavalores_jobbinsi))]<-4 ##e7298a        ## Define the color of Putumayo to #e7298a
names(clusMember) <- rownames(datavalores_jobbinsi)                     ## Put the names to the cluster member vector
clusDendro_jobbinsi <- as.dendrogram(as.hclust(clust_valores_jobbinsi)) ## Change the class of the cluster to dendrogram
clusDendro_jobbinsi <- dendrapply(clusDendro_jobbinsi, colLab)          ## Change the color of the leaves

pdf(file="hclust_D_jobbinsi.pdf", width = 10, height = 6)
par(mar=c(0.5, 4.5, 1, 1))
par(oma= c(0,0,1,0))
op = par(bg = "gray90")
par(lwd=3)
plot(clusDendro_jobbinsi,horiz=F,axes=T, ylab= "Height", cex.axis=1.3,cex.lab=1.7)
par(lwd=1)
legend("topright", pch= 21, pt.bg=c("#1b9e77", "royalblue4", "darkorange4", "#e7298a"), 
       legend=expression("Amazonas", "Caquetá", "Chocó", "Putumayo"), pt.cex = 1.7)
title(main= (expression(paste( " ", italic("Diachlorus jobbinsi")))), outer = F)
dev.off()

### Perform hierarchical cluster analysis, as above, but this time for D. leucotibialis
valores_leucotibialis <- pca_d_leucotibialis$pc.scores
datavalores_leucotibialis <- data.frame(valores_leucotibialis)
distvalores_leucotibialis <- dist(datavalores_leucotibialis)
clust_valores_leucotibialis <- hclust(distvalores_leucotibialis)
labelColors <- c("#1b9e77", "royalblue4", "#e7298a")
clusMember <- rep(1,length(rownames(datavalores_leucotibialis)))
clusMember[grep("a",rownames(datavalores_leucotibialis))]<-1 ##1b9e77
clusMember[grep("c",rownames(datavalores_leucotibialis))]<-2 #royalblue4
clusMember[grep("p",rownames(datavalores_leucotibialis))]<-3 ##e7298a
names(clusMember) <- rownames(datavalores_leucotibialis)
clusDendro_leucotibialis <- as.dendrogram(as.hclust(clust_valores_leucotibialis))
clusDendro_leucotibialis <- dendrapply(clusDendro_leucotibialis, colLab)

pdf(file="hclust_D_leucotibialis.pdf", width = 10, height = 6)
par(mar=c(0.5, 4.5, 1, 1))
par(oma= c(0,0,1,0))
op = par(bg = "gray90")
par(lwd=3)
plot(clusDendro_leucotibialis,horiz=F,axes=T, ylab= "Height", cex.axis=1.3,cex.lab=1.7)
par(lwd=1)
legend("topright", pch= 21, pt.bg=c("#1b9e77", "royalblue4", "#e7298a"), 
       legend=expression("Amazonas", "Caquetá", "Putumayo"),  pt.cex = 2, cex=1.2)
title(main= (expression(paste(" ", italic("Diachlorus leucotibialis")))), outer = F)
dev.off()

######################################################################################################################################################################################
######################################################################################################################################################################################
## 7. Define the number of clusters using model-based clustering, "Gap" statistic and Calinski-Harabasz criterion

### Model-based clustering
new_pcs_d_curvipes <- subset (pca_d_curvipes$pc.scores, select = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))         ## Make a data using the first ten components
best_model_d_curvipes <- Mclust(new_pcs_d_curvipes, G=1:5)                                                                     ## Perform model clustering of D. curvipes
summary(best_model_d_curvipes)                                                                                                 ## See the results of the model-based clustering for D. curvipes 
pdf(file="model_classification_d_curvipes.pdf", width = 10, height = 6)
plot(best_model_d_curvipes, what = c("BIC"), dimens = c(1,2))
plot(best_model_d_curvipes, what = c("classification"), dimens = c(1,2))
dev.off()

new_pcs_d_fuscistigma <- subset (pca_d_fuscistigma$pc.scores, select = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))   ## Make a data using the first ten components
best_model_d_fuscistigma <- Mclust(new_pcs_d_fuscistigma, G=1:6)                                                               ## Perform model clustering of D. fuscistigma
summary(best_model_d_fuscistigma)                                                                                              ## See the results of the model-based clustering for D. fuscistigma
pdf(file="model_classification_d_fuscistigma.pdf", width = 10, height = 6)
plot(best_model_d_fuscistigma, what = c("BIC"), dimens = c(1,2))
plot(best_model_d_fuscistigma, what = c("classification"), dimens = c(1,2))
dev.off()

new_pcs_d_jobbinsi <- subset (pca_d_jobbinsi$pc.scores, select = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))         ## Make a data using the first ten components
best_model_d_jobbinsi <- Mclust(new_pcs_d_jobbinsi, G=1:5)                                                                     ## Perform model clustering of D. jobbinsi
summary(best_model_d_jobbinsi)                                                                                                 ## See the results of the model-based clustering for D. jobbinsi 
pdf(file="model_classification_d_jobbinsi.pdf", width = 10, height = 6)
plot(best_model_d_jobbinsi, what = c("BIC"), dimens = c(1,2))
plot(best_model_d_jobbinsi, what = c("classification"), dimens = c(1,2))
dev.off()

new_pcs_d_leucotibialis <- subset (pca_d_leucotibialis$pc.scores, select = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)) ## Make a data using the first ten components
best_model_d_leucotibialis <- Mclust(new_pcs_d_leucotibialis, G=1:4)                                                             ## Perform model clustering of D. leucotibialis
summary(best_model_d_leucotibialis)                                                                                              ## See the results of the model-based clustering for D. leucotibialis 
pdf(file="model_classification_d_leucotibialis.pdf", width = 10, height = 6)
plot(best_model_d_leucotibialis, what = c("BIC"), dimens = c(1,2))
plot(best_model_d_leucotibialis, what = c("classification"), dimens = c(1,2))
dev.off()

new_pcs_d_nuneztovari <- subset (pca_d_nuneztovari$pc.scores, select = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))   ## Make a data using the first ten components
best_model_d_nuneztovari <- Mclust(new_pcs_d_nuneztovari, G=1:3)                                                               ## Perform model clustering of D. nuneztovari
summary(best_model_d_nuneztovari)                                                                                              ## See the results of the model-based clustering for D. nuneztovari 
pdf(file="model_classification_d_nuneztovari.pdf", width = 10, height = 6)
plot(best_model_d_nuneztovari, what = c("BIC"), dimens = c(1,2))
plot(best_model_d_nuneztovari, what = c("classification"), dimens = c(1,2))
dev.off()

### "Gap" statistic
gap_d_curvipes <- clusGap(new_pcs_d_curvipes, kmeans, 5, B = 1000, verbose = interactive())                                    ## Perform the gap statistic for D. curvipes
gap_d_curvipes                                                                                                                 ## See the results
pdf(file="gap_d_curvipes.pdf", width = 10, height = 6)            
plot(gap_d_curvipes)
dev.off()

gap_d_fuscistigma <- clusGap(new_pcs_d_fuscistigma, kmeans, 6, B = 1000, verbose = interactive())                              ## Perform the gap statistic for D. fuscistigma
gap_d_fuscistigma                                                                                                              ## See the results
pdf(file="gap_d_fuscistigma.pdf", width = 10, height = 6)            
plot(gap_d_fuscistigma)
dev.off()

gap_d_jobbinsi <- clusGap(new_pcs_d_jobbinsi, kmeans, 5, B = 1000, verbose = interactive())                                    ## Perform the gap statistic for D. jobbinsi
gap_d_jobbinsi                                                                                                                 ## See the results
pdf(file="gap_d_jobbinsi.pdf", width = 10, height = 6)            
plot(gap_d_jobbinsi)
dev.off()

gap_d_leucotibialis <- clusGap(new_pcs_d_leucotibialis, kmeans, 4, B = 1000, verbose = interactive())                          ## Perform the gap statistic for D. leucotibialis
gap_d_leucotibialis                                                                                                            ## See the results
pdf(file="gap_d_leucotibialis.pdf", width = 10, height = 6)            
plot(gap_d_leucotibialis)
dev.off()

gap_d_nuneztovari <- clusGap(new_pcs_d_nuneztovari, kmeans, 3, B = 1000, verbose = interactive())                              ## Perform the gap statistic for D. nuneztovari
gap_d_nuneztovari                                                                                                              ## See the results
pdf(file="gap_d_nuneztovari.pdf", width = 10, height = 6)            
plot(gap_d_nuneztovari)
dev.off()

### Calinski-Harabasz criterion
fit_curvipes <- cascadeKM(scale(new_pcs_d_curvipes, center = TRUE,  scale = TRUE), 1, 5, iter = 1000)                          ## Perform analysis using Calinski-Harabasz criterion for D. curvipes
fit_curvipes$results                                                                                                           ## See the results
pdf(file="calinski_d_curvipes.pdf", width = 10, height = 6)            
plot(fit_curvipes, sortg = TRUE, grpmts.plot = TRUE)
dev.off()

fit_fuscistigma <- cascadeKM(scale(new_pcs_d_fuscistigma, center = TRUE,  scale = TRUE), 1, 6, iter = 1000)                    ## Perform analysis using Calinski-Harabasz criterion for D. fuscistigma
fit_fuscistigma$results                                                                                                        ## See the results
pdf(file="calinski_d_fuscistigma.pdf", width = 10, height = 6)            
plot(fit_fuscistigma, sortg = TRUE, grpmts.plot = TRUE)
dev.off()

fit_jobbinsi <- cascadeKM(scale(new_pcs_d_jobbinsi, center = TRUE,  scale = TRUE), 1, 5, iter = 1000)                          ## Perform analysis using Calinski-Harabasz criterion for D. jobbinsi
fit_jobbinsi$results                                                                                                           ## See the results
pdf(file="calinski_d_jobbinsi.pdf", width = 10, height = 6)            
plot(fit_jobbinsi, sortg = TRUE, grpmts.plot = TRUE)
dev.off()

fit_leucotibialis <- cascadeKM(scale(new_pcs_d_leucotibialis, center = TRUE,  scale = TRUE), 1, 4, iter = 1000)                ## Perform analysis using Calinski-Harabasz criterion for D. leucotibialis
fit_leucotibialis$results                                                                                                      ## See the results
pdf(file="calinski_d_leucotibialis.pdf", width = 10, height = 6)            
plot(fit_leucotibialis, sortg = TRUE, grpmts.plot = TRUE)
dev.off()

fit_nuneztovari <- cascadeKM(scale(new_pcs_d_nuneztovari, center = TRUE,  scale = TRUE), 1, 3, iter = 1000)                    ## Perform analysis using Calinski-Harabasz criterion for D. nuneztovari
fit_nuneztovari$results                                                                                                        ## See the results
pdf(file="calinski_d_nuneztovari.pdf", width = 10, height = 6)            
plot(fit_nuneztovari, sortg = TRUE, grpmts.plot = TRUE)
dev.off()

### NbClust
# nb_curvi <- NbClust(new_pcs_d_curvipes, diss= NULL, distance = "euclidean", min.nc=2, max.nc=5, method = "kmeans", index = "alllong", alphaBeale = 0.1)
# nb_fusci <- NbClust(new_pcs_d_fuscistigma, diss= NULL, distance = "euclidean", min.nc=2, max.nc=5, method = "kmeans", index = "alllong", alphaBeale = 0.1)
# nb_jobi <- NbClust(new_pcs_d_jobbinsi, diss= NULL, distance = "euclidean", min.nc=2, max.nc=5, method = "kmeans", index = "alllong", alphaBeale = 0.1)
# nb_leuco <- NbClust(new_pcs_d_leucotibialis, diss= NULL, distance = "euclidean", min.nc=2, max.nc=5, method = "kmeans", index = "alllong", alphaBeale = 0.1)
# nb_nune <- NbClust(new_pcs_d_nuneztovari, diss= NULL, distance = "euclidean", min.nc=2, max.nc=5, method = "kmeans", index = "alllong", alphaBeale = 0.1)

######################################################################################################################################################################################
######################################################################################################################################################################################
## 8. Discriminant analysis (linear and quadratic)

## discriminant analysis of D. curvipes
dis_lin_d_curvipes <- linDA(new_pcs_d_curvipes, factors_d_curvi, validation="crossval")
percentage_confusion_dislineal_curvi <- NULL
for (i in 1:dim(dis_lin_d_curvipes$confusion)[1]) {
  row_new<- (dis_lin_d_curvipes$confusion[i, ]*100)/sum(dis_lin_d_curvipes$confusion[i, ])
  percentage_confusion_dislineal_curvi <- rbind(percentage_confusion_dislineal_curvi, row_new)
}
row.names(percentage_confusion_dislineal_curvi) <- dimnames(percentage_confusion_dislineal_curvi)[[2]]
write.csv(percentage_confusion_dislineal_curvi, file="dislineal_d_curvi_confusion.csv")
dis_lin_d_curvipes$error

dis_qua_d_curvipes <- quaDA(new_pcs_d_curvipes, factors_d_curvi, validation="crossval")
percentage_confusion_disqua_curvi <- NULL
for (i in 1:dim(dis_qua_d_curvipes$confusion)[1]) {
  row_new<- (dis_qua_d_curvipes$confusion[i, ]*100)/sum(dis_qua_d_curvipes$confusion[i, ])
  percentage_confusion_disqua_curvi <- rbind(percentage_confusion_disqua_curvi, row_new)
}
row.names(percentage_confusion_disqua_curvi) <- dimnames(percentage_confusion_disqua_curvi)[[2]]
write.csv(percentage_confusion_disqua_curvi, file="discua_d_curvipess_confusion.csv")
dis_qua_d_curvipes$error

## discriminant analysis of D. jobbinsi
dis_lin_d_jobbinsi <- linDA(new_pcs_d_jobbinsi, factors_d_jobi, validation="crossval")
percentage_confusion_dislineal_jobi <- NULL
for (i in 1:dim(dis_lin_d_jobbinsi$confusion)[1]) {
  row_new<- (dis_lin_d_jobbinsi$confusion[i, ]*100)/sum(dis_lin_d_jobbinsi$confusion[i, ])
  percentage_confusion_dislineal_jobi <- rbind(percentage_confusion_dislineal_jobi, row_new)
}
row.names(percentage_confusion_dislineal_jobi) <- dimnames(percentage_confusion_dislineal_jobi)[[2]]
write.csv(percentage_confusion_dislineal_jobi, file="dislineal_d_jobi_confusion.csv")
dis_lin_d_jobbinsi$error

dis_qua_d_jobbinsi <- quaDA(new_pcs_d_jobbinsi, factors_d_jobi, validation="crossval")
percentage_confusion_disqua_jobi <- NULL
for (i in 1:dim(dis_qua_d_jobbinsi$confusion)[1]) {
  row_new<- (dis_qua_d_jobbinsi$confusion[i, ]*100)/sum(dis_qua_d_jobbinsi$confusion[i, ])
  percentage_confusion_disqua_jobi <- rbind(percentage_confusion_disqua_jobi, row_new)
}
row.names(percentage_confusion_disqua_jobi) <- dimnames(percentage_confusion_disqua_jobi)[[2]]
write.csv(percentage_confusion_disqua_jobi, file="discua_d_jobbinsi_confusion.csv")
dis_qua_d_jobbinsi$error

## discriminant analysis of D. fuscistigma
dis_lin_d_fuscistigma <- linDA(new_pcs_d_fuscistigma, factors_d_fusci, validation="crossval")
percentage_confusion_dislineal_fusci <- NULL
for (i in 1:dim(dis_lin_d_fuscistigma$confusion)[1]) {
  row_new<- (dis_lin_d_fuscistigma$confusion[i, ]*100)/sum(dis_lin_d_fuscistigma$confusion[i, ])
  percentage_confusion_dislineal_fusci <- rbind(percentage_confusion_dislineal_fusci, row_new)
}
row.names(percentage_confusion_dislineal_fusci) <- dimnames(percentage_confusion_dislineal_fusci)[[2]]
write.csv(percentage_confusion_dislineal_fusci, file="dislineal_d_fusci_confusion.csv")
dis_lin_d_fuscistigma$error

dis_qua_d_fuscistigma <- quaDA(new_pcs_d_fuscistigma, factors_d_fusci, validation="crossval")
percentage_confusion_disqua_fusci <- NULL
for (i in 1:dim(dis_qua_d_fuscistigma$confusion)[1]) {
  row_new<- (dis_qua_d_fuscistigma$confusion[i, ]*100)/sum(dis_qua_d_fuscistigma$confusion[i, ])
  percentage_confusion_disqua_fusci <- rbind(percentage_confusion_disqua_fusci, row_new)
}
row.names(percentage_confusion_disqua_fusci) <- dimnames(percentage_confusion_disqua_fusci)[[2]]
write.csv(percentage_confusion_disqua_fusci, file="discua_d_fuscistigma_confusion.csv")
dis_qua_d_fuscistigma$error

## discriminant analysis of D. leucotibialis
dis_lin_d_leucotibialis <- linDA(new_pcs_d_leucotibialis, factors_d_leuco, validation="crossval")
percentage_confusion_dislineal_leuco <- NULL
for (i in 1:dim(dis_lin_d_leucotibialis$confusion)[1]) {
  row_new<- (dis_lin_d_leucotibialis$confusion[i, ]*100)/sum(dis_lin_d_leucotibialis$confusion[i, ])
  percentage_confusion_dislineal_leuco <- rbind(percentage_confusion_dislineal_leuco, row_new)
}
row.names(percentage_confusion_dislineal_leuco) <- dimnames(percentage_confusion_dislineal_leuco)[[2]]
write.csv(percentage_confusion_dislineal_leuco, file="dislineal_d_leuco_confusion.csv")
dis_lin_d_leucotibialis$error

dis_qua_d_leucotibialis <- quaDA(new_pcs_d_leucotibialis, factors_d_leuco, validation="crossval")
percentage_confusion_disqua_leuco <- NULL
for (i in 1:dim(dis_qua_d_leucotibialis$confusion)[1]) {
  row_new<- (dis_qua_d_leucotibialis$confusion[i, ]*100)/sum(dis_qua_d_leucotibialis$confusion[i, ])
  percentage_confusion_disqua_leuco <- rbind(percentage_confusion_disqua_leuco, row_new)
}
row.names(percentage_confusion_disqua_leuco) <- dimnames(percentage_confusion_disqua_leuco)[[2]]
write.csv(percentage_confusion_disqua_leuco, file="discua_d_leucotibialis_confusion.csv")
dis_qua_d_leucotibialis$error

## discriminant analysis of D. nuneztovari
dis_lin_d_nuneztovari <- linDA(new_pcs_d_nuneztovari, factors_d_nune, validation="crossval")
percentage_confusion_dislineal_nune <- NULL
for (i in 1:dim(dis_lin_d_nuneztovari$confusion)[1]) {
  row_new<- (dis_lin_d_nuneztovari$confusion[i, ]*100)/sum(dis_lin_d_nuneztovari$confusion[i, ])
  percentage_confusion_dislineal_nune <- rbind(percentage_confusion_dislineal_nune, row_new)
}
row.names(percentage_confusion_dislineal_nune) <- dimnames(percentage_confusion_dislineal_nune)[[2]]
write.csv(percentage_confusion_dislineal_nune, file="dislineal_d_nune_confusion.csv")
dis_lin_d_nuneztovari$error

dis_qua_d_nuneztovari <- quaDA(new_pcs_d_nuneztovari, factors_d_nune, validation="crossval")
percentage_confusion_disqua_nune <- NULL
for (i in 1:dim(dis_qua_d_nuneztovari$confusion)[1]) {
  row_new<- (dis_qua_d_nuneztovari$confusion[i, ]*100)/sum(dis_qua_d_nuneztovari$confusion[i, ])
  percentage_confusion_disqua_nune <- rbind(percentage_confusion_disqua_nune, row_new)
}
row.names(percentage_confusion_disqua_nune) <- dimnames(percentage_confusion_disqua_nune)[[2]]
write.csv(percentage_confusion_disqua_nune, file="discua_d_nuneztovari_confusion.csv")
dis_qua_d_nuneztovari$error

######################################################################################################################################################################################
######################################################################################################################################################################################
## 9. Riemmanian distances among species and Goodall's F test 

### D. curvipes
names_d_curvipes <- combn(names(cons_shapes[7:10]), m=2)
comb_d_curvipes <- combn(cons_shapes[7:10], m=2)
r_d_curvi <- NULL
nam <- NULL
for (i in 1:dim(comb_d_curvipes)[2]) {
     n <- paste(names_d_curvipes[,i][1], "to", names_d_curvipes[,i][2])  
     r <- riemdist(as.matrix(comb_d_curvipes[,i][1][[1]]), as.matrix(comb_d_curvipes[,i][2][[1]]))
     r_d_curvi <- rbind(r_d_curvi, r)
     nam <- c(nam, n)
}
rownames(r_d_curvi) <- nam
colnames(r_d_curvi) <- c("Riemman")
write.csv(r_d_curvi, "Riemmanian_distances_curvipes.csv")

alg_names_curvi <- combn(names(alg_intra[7:10]), m=2)
comb_alg_d_curvi <- combn(alg_intra[7:10], m=2)
name <- NULL
t_goodall_curvi <- NULL
for (i in 1:dim(comb_alg_d_curvi)[2]) {
     n <- paste(alg_names_curvi[,i][1], "and", alg_names_curvi[,i][2]) 
     t <- testmeanshapes(comb_alg_d_curvi[,i][1][1][[1]][[1]], comb_alg_d_curvi[,i][2][1][[1]][[1]], resamples = 1000, replace = FALSE, scale= TRUE)
     t_value <- c(t$G, t$G.pvalue)
     t_goodall_curvi <- rbind(t_goodall_curvi, t_value)
     name <- c(name, n)
}
rownames(t_goodall_curvi) <- name
colnames(t_goodall_curvi) <- c("G", "G_pvalue")
write.csv(t_goodall_curvi, "Goodall_curvipes.csv")

### D. fuscistigma
names_d_fuscistigma <- combn(names(cons_shapes[11:15]), m=2)
comb_d_fuscistigma <- combn(cons_shapes[11:15], m=2)
r_d_fusci <- NULL
nam <- NULL
for (i in 1:dim(comb_d_fuscistigma)[2]) {
     n <- paste(names_d_fuscistigma[,i][1], "to", names_d_fuscistigma[,i][2])  
     r <- riemdist(as.matrix(comb_d_fuscistigma[,i][1][[1]]), as.matrix(comb_d_fuscistigma[,i][2][[1]]))
     r_d_fusci <- rbind(r_d_fusci, r)
     nam <- c(nam, n)
}
rownames(r_d_fusci) <- nam
colnames(r_d_fusci) <- c("Riemman")
write.csv(r_d_fusci, "Riemmanian_distances_fuscistigma.csv")

alg_names_fusci <- combn(names(alg_intra[11:15]), m=2)
comb_alg_d_fusci <- combn(alg_intra[11:15], m=2)
name <- NULL
t_goodall_fusci <- NULL
for (i in 1:dim(comb_alg_d_fusci)[2]) {
     n <- paste(alg_names_fusci[,i][1], "and", alg_names_fusci[,i][2]) 
     t <- testmeanshapes(comb_alg_d_fusci[,i][1][1][[1]][[1]], comb_alg_d_fusci[,i][2][1][[1]][[1]], resamples = 1000, replace = FALSE, scale= TRUE)
     t_value <- c(t$G, t$G.pvalue)
     t_goodall_fusci <- rbind(t_goodall_fusci, t_value)
     name <- c(name, n)
}
rownames(t_goodall_fusci) <- name
colnames(t_goodall_fusci) <- c("G", "G_pvalue")
write.csv(t_goodall_fusci, "Goodall_fuscistigma.csv")

### D. jobbinsi
names_d_jobbinsi <- combn(names(cons_shapes[16:19]), m=2)
comb_d_jobbinsi <- combn(cons_shapes[16:19], m=2)
r_d_jobi <- NULL
nam <- NULL
for (i in 1:dim(comb_d_jobbinsi)[2]) {
     n <- paste(names_d_jobbinsi[,i][1], "to", names_d_jobbinsi[,i][2])  
     r <- riemdist(as.matrix(comb_d_jobbinsi[,i][1][[1]]), as.matrix(comb_d_jobbinsi[,i][2][[1]]))
     r_d_jobi <- rbind(r_d_jobi, r)
     nam <- c(nam, n)
}
rownames(r_d_jobi) <- nam
colnames(r_d_jobi) <- c("Riemman")
write.csv(r_d_jobi, "Riemmanian_distances_jobbinsi.csv")

alg_names_jobi <- combn(names(alg_intra[16:19]), m=2)
comb_alg_d_jobi <- combn(alg_intra[16:19], m=2)
name <- NULL
t_goodall_jobi <- NULL
for (i in 1:dim(comb_alg_d_jobi)[2]) {
     n <- paste(alg_names_jobi[,i][1], "and", alg_names_jobi[,i][2]) 
     t <- testmeanshapes(comb_alg_d_jobi[,i][1][1][[1]][[1]], comb_alg_d_jobi[,i][2][1][[1]][[1]], resamples = 1000, replace = FALSE, scale= TRUE)
     t_value <- c(t$G, t$G.pvalue)
     t_goodall_jobi <- rbind(t_goodall_jobi, t_value)
     name <- c(name, n)
}
rownames(t_goodall_jobi) <- name
colnames(t_goodall_jobi) <- c("G", "G_pvalue")
write.csv(t_goodall_jobi, "Goodall_jobbinsi.csv")

### D. leucotibialis
names_d_leucotibialis <- combn(names(cons_shapes[20:22]), m=2)
comb_d_leucotibialis <- combn(cons_shapes[20:22], m=2)
r_d_leuco <- NULL
nam <- NULL
for (i in 1:dim(comb_d_leucotibialis)[2]) {
     n <- paste(names_d_leucotibialis[,i][1], "to", names_d_leucotibialis[,i][2])  
     r <- riemdist(as.matrix(comb_d_leucotibialis[,i][1][[1]]), as.matrix(comb_d_leucotibialis[,i][2][[1]]))
     r_d_leuco <- rbind(r_d_leuco, r)
     nam <- c(nam, n)
}
rownames(r_d_leuco) <- nam
colnames(r_d_leuco) <- c("Riemman")
write.csv(r_d_leuco, "Riemmanian_distances_leucotibialis.csv")

alg_names_leuco <- combn(names(alg_intra[20:22]), m=2)
comb_alg_d_leuco <- combn(alg_intra[20:22], m=2)
name <- NULL
t_goodall_leuco <- NULL
for (i in 1:dim(comb_alg_d_leuco)[2]) {
     n <- paste(alg_names_leuco[,i][1], "and", alg_names_leuco[,i][2]) 
     t <- testmeanshapes(comb_alg_d_leuco[,i][1][1][[1]][[1]], comb_alg_d_leuco[,i][2][1][[1]][[1]], resamples = 1000, replace = FALSE, scale= TRUE)
     t_value <- c(t$G, t$G.pvalue)
     t_goodall_leuco <- rbind(t_goodall_leuco, t_value)
     name <- c(name, n)
}
rownames(t_goodall_leuco) <- name
colnames(t_goodall_leuco) <- c("G", "G_pvalue")
write.csv(t_goodall_leuco, "Goodall_leucotibialis.csv")

### D. nuneztovari
names_d_nuneztovari <- combn(names(cons_shapes[23:24]), m=2)
comb_d_nuneztovari <- combn(cons_shapes[23:24], m=2)
r_d_nune <- NULL
nam <- NULL
for (i in 1:dim(comb_d_nuneztovari)[2]) {
     n <- paste(names_d_nuneztovari[,i][1], "to", names_d_nuneztovari[,i][2])  
     r <- riemdist(as.matrix(comb_d_nuneztovari[,i][1][[1]]), as.matrix(comb_d_nuneztovari[,i][2][[1]]))
     r_d_nune <- rbind(r_d_nune, r)
     nam <- c(nam, n)
}
rownames(r_d_nune) <- nam
colnames(r_d_nune) <- c("Riemman")
write.csv(r_d_nune, "Riemmanian_distances_nuneztovari.csv")

alg_names_nune <- combn(names(alg_intra[23:24]), m=2)
comb_alg_d_nune <- combn(alg_intra[23:24], m=2)
name <- NULL
t_goodall_nune <- NULL
for (i in 1:dim(comb_alg_d_nune)[2]) {
     n <- paste(alg_names_nune[,i][1], "and", alg_names_nune[,i][2]) 
     t <- testmeanshapes(comb_alg_d_nune[,i][1][1][[1]][[1]], comb_alg_d_nune[,i][2][1][[1]][[1]], resamples = 1000, replace = FALSE, scale= TRUE)
     t_value <- c(t$G, t$G.pvalue)
     t_goodall_nune <- rbind(t_goodall_nune, t_value)
     name <- c(name, n)
}
rownames(t_goodall_nune) <- name
colnames(t_goodall_nune) <- c("G", "G_pvalue")
write.csv(t_goodall_nune, "Goodall_nuneztovari.csv")
