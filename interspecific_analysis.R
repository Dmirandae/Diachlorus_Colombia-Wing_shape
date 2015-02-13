## Ambrosio Torres & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Research: Wing shape variation in the taxonomic recognition of species of Diachlorus Osten-Sacken (Diptera: Tabanidae) from Colombia
## Part: interspecific variation of genus Diachlorus
## R version 3.1.2 & Rstudio 0.97.551
## 27th January 2015 

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
library(shapes)      ## Load 'shapes' package v. 1.1-9 --- http://cran.r-project.org/web/packages/shapes/index.html
library(vegan)       ## Load 'vegan' package v. 2.0-10 ---- http://cran.r-project.org/web/packages/vegan/index.html
source("colLab.R")   ## Load "colLab" function, modified by us to coloring edgePar and labels of a dendrogram ---- https://stat.ethz.ch/R-manual/R-patched/library/stats/html/dendrapply.html 
                     ## and http://stackoverflow.com/

######################################################################################################################################################################################
######################################################################################################################################################################################
## 2. DEFINE THE GROUPS OF STUDY TO MAKE THE COMPARISONS IN THE ANALYSIS AND READ THE DATA (BUILD ARRAYS)

list_names <- list.files(pattern = "nts")    ## Make a list with all the names of landmarks-files in .nts format, of all individuals (of all species) in study
list_curvi <- list.files(pattern = "dc")     ## Make a list with all the names of landmarks-files of Diachlorus curvipes
list_fusci <- list.files(pattern = "df")     ## Make a list with all the names of landmarks-files of Diachlorus fuscistigma
list_jobi <- list.files(pattern = "dj")      ## Make a list with all the names of landmarks-files of Diachlorus jobbinsi
list_leuco <- list.files(pattern = "dl")     ## Make a list with all the names of landmarks-files of Diachlorus leucotibialis
list_nune <- list.files(pattern = "dn")      ## Make a list with all the names of landmarks-files of Diachlorus nuneztovari
list_leti <- list.files(pattern = "dt")      ## Make a list with all the names of landmarks-files of Diachlorus leticia

diachlorus <- readmulti.nts(list_names)      ## Build an array/object which contains all individuals (of all species) in study of the genus Diachlorus, from all .nts files
d_curvipes <- readmulti.nts(list_curvi)      ## Build an array/object which contains all individuals of Diachlorus curvipes
d_fuscistigma <- readmulti.nts(list_fusci)   ## Build an array/object which contains all individuals of Diachlorus fuscistigma
d_jobbinsi <- readmulti.nts(list_jobi)       ## Build an array/object which contains all individuals of Diachlorus jobbinsi
d_leucotibialis <- readmulti.nts(list_leuco) ## Build an array/object which contains all individuals of Diachlorus leucotibialis
d_leticia <- readmulti.nts(list_leti)        ## Build an array/object which contains all individuals of Diachlorus leticia
d_nuneztovari <- readmulti.nts(list_nune)    ## Build an array/object which contains all individuals of Diachlorus nuneztovari

######################################################################################################################################################################################
######################################################################################################################################################################################
## 3. CREATE A DIRECTORY TO PUT THE RESULTS IN A DIFERENT DIRECTORY - NOT IN THE DIRECTORY OF TE DATA (INPUT)

namedir <- paste("interspecific_results_", date(),  sep="") ## Create a name to the new directory
dir.create(namedir, showWarnings = FALSE)                   ## Create the new directory 
setwd(paste((getwd()),  "/" ,namedir, sep=""))              ## Change the directory to the new created directory

######################################################################################################################################################################################
######################################################################################################################################################################################
## 4. GENERALIZED PROCRUSTES ANALYSIS (GPA) AND CONSESUS SHAPES

groups_inter <- list(diachlorus=diachlorus, d_curvipes=d_curvipes, d_fuscistigma=d_fuscistigma, d_jobbinsi=d_jobbinsi, 
               d_leucotibialis=d_leucotibialis, d_nuneztovari=d_nuneztovari, d_leticia=d_leticia)                       ## Make a list of the arrays/groups
alg_inter <- lapply(groups_inter, gpagen)                                                                               ## Align all the groups in the list                
alg_diachlorus <- alg_inter$diachlorus$coords ; writeland.tps(alg_diachlorus, "alg_diachlorus.tps")                     ## Write the aligned data in .tps format to all Diachlorus data
alg_d_curvipes <- alg_inter$d_curvipes$coords ; writeland.tps(alg_d_curvipes, "alg_d_curvipes.tps")                     ## Write the aligned data in .tps format to D. curvipes data
alg_d_fuscistigma <- alg_inter$d_fuscistigma$coords ; writeland.tps(alg_d_fuscistigma, "alg_d_fuscistigma.tps")         ## Write the aligned data in .tps format to D. fuscistigma data
alg_d_jobbinsi <- alg_inter$d_jobbinsi$coords ; writeland.tps(alg_d_jobbinsi, "alg_d_jobbinsi.tps")                     ## Write the aligned data in .tps format to D. jobbinsi data
alg_d_leucotibialis <- alg_inter$d_leucotibialis$coords ; writeland.tps(alg_d_leucotibialis, "alg_d_leucotibialis.tps") ## Write the aligned data in .tps format to D. leucotibialis data
alg_d_leticia <- alg_inter$d_leticia$coords ; writeland.tps(alg_d_leticia, "alg_d_leticia.tps")                         ## Write the aligned data in .tps format to D. leticia data
alg_d_nuneztovari <- alg_inter$d_nuneztovari$coords ; writeland.tps(alg_d_nuneztovari, "alg_d_nuneztovari.tps")         ## Write the aligned data in .tps format to D. nuneztovari data

cons_names <- list(alg_diachlorus=alg_diachlorus, alg_d_curvipes=alg_d_curvipes, alg_d_fuscistigma=alg_d_fuscistigma, 
                   alg_d_jobbinsi=alg_d_jobbinsi, alg_d_leucotibialis=alg_d_leucotibialis, alg_d_leticia=alg_d_leticia, 
                   alg_d_nuneztovari=alg_d_nuneztovari)                                                                 ## Make a list of the aligned arrays/groups
cons_shapes <- lapply(cons_names, mshape)                                                                               ## Calculate the consensus shape of each array/group

######################################################################################################################################################################################
######################################################################################################################################################################################
## 5. Principal Component Analysis (PCA)

### Creating a factor to color and recognize each species of the analysis to do a PCA of Diachlorus
factors_diachlorus <-sub(".nts", "", list_names)
for (i in (rep(c(0,1,2,3,4,5,6,7,8,9),5))) {
  factors_diachlorus <-sub(i, "", factors_diachlorus)
}
factors_diachlorus <- sub("adc", "dc", factors_diachlorus) ; factors_diachlorus <- sub("adf", "df", factors_diachlorus)
factors_diachlorus <- sub("adj", "dj", factors_diachlorus) ; factors_diachlorus <- sub("adl", "dl", factors_diachlorus)
factors_diachlorus <- sub("adn", "dn", factors_diachlorus) ; factors_diachlorus <- sub("cdf", "df", factors_diachlorus)
factors_diachlorus <- sub("cdl", "dl", factors_diachlorus) ; factors_diachlorus <- sub("cdj", "dj", factors_diachlorus)
factors_diachlorus <- sub("chdc", "dc", factors_diachlorus) ; factors_diachlorus <- sub("chdj", "dj", factors_diachlorus)
factors_diachlorus <- sub("mdc", "dc", factors_diachlorus) ; factors_diachlorus <- sub("mdf", "df", factors_diachlorus)
factors_diachlorus <- sub("pdc", "dc", factors_diachlorus) ; factors_diachlorus <- sub("pdf", "df", factors_diachlorus)
factors_diachlorus <- sub("pdj", "dj", factors_diachlorus) ; factors_diachlorus <- sub("pdl", "dl", factors_diachlorus)
factors_diachlorus <- sub("pdn", "dn", factors_diachlorus) ; factors_diachlorus <- sub("vdf", "df", factors_diachlorus)
factors_diachlorus <- sub("vdt", "dt", factors_diachlorus)
factors_diachlorus <- factor(factors_diachlorus)

pca_diachlorus <- plotTangentSpace(alg_diachlorus, label=NULL, verbose =T, groups = factors_diachlorus, warpgrids=F)          ## Perform PCA of Diachlorus

### Diachlorus PCA plotted in ggplot2
label <- list(expression(italic('D. curvipes')), expression(italic('D. fuscistigma')), expression(italic('D. jobbinsi')), 
              expression(italic('D. leucotibialis')), expression(italic('D. nuneztovari')), expression(italic('D. leticia'))) ## Label to the plot 
a <- cbind(pca_diachlorus$pc.scores, factors_diachlorus)                                                                      ## Attach the factor for the names to the data 
a <- data.frame(a) ## Convert class of a to data frame
a$factors_diachlorus <- as.factor(a$factors_diachlorus)
### Define the ellipses to the groups (95%)
df_ell <- data.frame()
for(g in levels(a$factors_diachlorus)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$factors_diachlorus==g,], ellipse(cor(PC1, PC2), 
  scale=c(sd(PC1),sd(PC2)), 
  centre=c(mean(PC1),mean(PC2))))),factors_diachlorus=g))
}
p <- ggplot(data=a, aes(x=PC1, y=PC2,colour=factors_diachlorus, shape=factors_diachlorus)) + geom_point(size=3) 
p <- p + geom_path(data=df_ell, aes(x=x, y=y,colour=factors_diachlorus), size=0.5, linetype=2) 
p <- p + scale_colour_brewer(name= "Species", labels= label, palette = "Set1") + scale_shape_discrete(name ="Species", labels= label) 
p <- p + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size=16))

pdf(file="pca_genus_Diachlorus_ggplot2.pdf", width = 10, height = 6) 
p
dev.off()

######################################################################################################################################################################################
######################################################################################################################################################################################
## 6. Hierarchical clustering

datavalores <- data.frame(pca_diachlorus$pc.scores)                                ## Make a data frame using the scores of the PCA
distvalores <- dist(datavalores)                                                   ## Make a distance matrix of "datavalores"
clust_valores <- hclust(distvalores)                                               ## Perform the hierarquical cluster
labelColors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33") ## Define the color of the leaves of the dendrogram
clusMember <- rep(1,length(rownames(datavalores)))
clusMember[grep("dc",rownames(datavalores))]<-1                                    ## Define the color of D. curvipes to #e41a1c
clusMember[grep("df",rownames(datavalores))]<-2                                    ## Define the color of D. fuscistigma to #377eb8
clusMember[grep("dj",rownames(datavalores))]<-3                                    ## Define the color of D. jobbinsi to #4daf4a
clusMember[grep("dl",rownames(datavalores))]<-4                                    ## Define the color of D. leucotibialis to #984ea3
clusMember[grep("dn",rownames(datavalores))]<-5                                    ## Define the color of D. nuneztovari to #ff7f00
clusMember[grep("dt",rownames(datavalores))]<-6                                    ## Define the color of D. leticia to #ffff33
names(clusMember) <- rownames(datavalores)                                         ## Put the names to the cluster member vector
clusDendro <- as.dendrogram(as.hclust(clust_valores))                              ## Change the class of the cluster to dendrogram
clusDendro <- dendrapply(clusDendro, colLab)                                       ## Change the color of the leaves

pdf(file="hclust_Diachlorus.pdf", width = 10, height = 6)
par(mar=c(0.5, 4.5, 1, 1))
par(oma= c(0,0,0,0))
op = par(bg = "gray90")
par(lwd=2)
plot(clusDendro,horiz=F,axes=T, ylab= "Height", cex.axis=1.3,cex.lab=1.7)
par(lwd=1)
legend("topright", pch= 21, pt.bg=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"), 
       legend=expression(italic("D. curvipes"),italic("D. fuscistigma"), italic("D. jobbinsi"), 
       italic("D. leucotibialis"),  italic("D. nuneztovari"), italic("D. leticia")), pt.cex = 2, cex=1.2)
dev.off()

######################################################################################################################################################################################
######################################################################################################################################################################################
## 7. Define the number of clusters using model-based clustering, "Gap" statistic and Calinski-Harabasz criterion

### Model-based clustering
new_pcs_diachlorus <- subset (pca_diachlorus$pc.scores, select = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8)) ## Make a data using the first eight components
best_model_diachlorus <- Mclust(new_pcs_diachlorus, G=1:7)                                                       ## Perform model clustering
summary(best_model_diachlorus)                                                                                   ## See the results of the model-based clustering
  
pdf(file="model_classification_diachlorus.pdf", width = 10, height = 6)
plot(best_model_diachlorus, what = c("BIC"), dimens = c(1,2))
plot(best_model_diachlorus, what = c("classification"), dimens = c(1,2))
dev.off()

### "Gap" statistic
gap_diachlorus <- clusGap(new_pcs_diachlorus, kmeans, 7, B = 1000, verbose = interactive())                     ## Perform the gap statistic
gap_diachlorus                                                                                                   ## See the results

pdf(file="gap_diachlorus.pdf", width = 10, height = 6)            
plot(gap_diachlorus)
dev.off()

### Calinski-Harabasz criterion
fit <- cascadeKM(scale(new_pcs_diachlorus, center = TRUE,  scale = TRUE), 1, 7, iter = 1000)                    ## Perform analysis using Calinski-Harabasz criterion
fit$results                                                                                                      ## See the results

pdf(file="calinski_diachlorus.pdf", width = 10, height = 6)            
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
dev.off()

######################################################################################################################################################################################
######################################################################################################################################################################################
## 8. Discriminant analysis (linear and quadratic)

dis_lin_diachlorus <- linDA(new_pcs_diachlorus, factors_diachlorus, validation="crossval")
percentage_confusion_dislineal <- NULL
for (i in 1:dim(dis_lin_diachlorus$confusion)[1]) {
  row_new<- (dis_lin_diachlorus$confusion[i, ]*100)/sum(dis_lin_diachlorus$confusion[i, ])
  percentage_confusion_dislineal <- rbind(percentage_confusion_dislineal, row_new)
}
row.names(percentage_confusion_dislineal) <- dimnames(percentage_confusion_dislineal)[[2]]
write.csv(percentage_confusion_dislineal, file="dislineal_diachlorus_confusion.csv")
dis_lin_diachlorus$error

dis_qua_diachlorus <- quaDA(new_pcs_diachlorus, factors_diachlorus, validation="crossval")
percentage_confusion_disqua <- NULL
for (i in 1:dim(dis_qua_diachlorus$confusion)[1]) {
  row_new<- (dis_qua_diachlorus$confusion[i, ]*100)/sum(dis_qua_diachlorus$confusion[i, ])
  percentage_confusion_disqua <- rbind(percentage_confusion_disqua, row_new)
}
row.names(percentage_confusion_disqua) <- dimnames(percentage_confusion_disqua)[[2]]
write.csv(percentage_confusion_disqua, file="discua_diachlorus_confusion.csv")
dis_qua_diachlorus$error

##Plot discriminants of Diachlorus in ggplot2
#goal_diachlorus_lin <- NULL
#for (i in 1:dim(percentage_confusion_dislineal)[1]) {
#  new_element<- percentage_confusion_dislineal[i, i]
#  goal_diachlorus_lin <- c(goal_diachlorus_lin, new_element)
#}
#goal_diachlorus_qua <- NULL
#for (i in 1:dim(percentage_confusion_disqua)[1]) {
#  new_element<- percentage_confusion_disqua[i, i]
#  goal_diachlorus_qua <- c(goal_diachlorus_qua, new_element)
#}
#SS <- c("D. curvipes", "D. fuscistigma", "D. jobbinsi", "D. leucotibialis", "D. nuneztovari", "D. leticia")
#df1 <- data.frame(Species = factor(SS), accuracy_lineal = goal_diachlorus_lin, accuracy_cuadratica = goal_diachlorus_qua) 
#pdf(file="discri_diachlorus_crossval.pdf", width = 11, height = 5.7)
#d <- ggplot(df1,aes(x=Species,y=accuracy_lineal,fill=Species)) +  stat_summary(fun.y=mean,position=position_dodge(), geom="bar") + scale_fill_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#ffff33", "#984ea3", "#ff7f00")) + scale_x_discrete("Species") +  coord_flip() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) + theme(legend.text = element_text(size = 14, face = "italic")) + theme(legend.title = element_text(size=16))
#d + geom_errorbar(aes(y=accuracy_cuadratica, ymax=accuracy_cuadratica, ymin=accuracy_cuadratica), linetype="dashed", size= 0.5, position=position_dodge()) + scale_y_continuous("Accuracy (%)", breaks = round(seq(0, 100, by = 5),1)) + theme(axis.ticks = element_blank(), axis.text.y = element_blank())
#dev.off()

######################################################################################################################################################################################
######################################################################################################################################################################################
## 9. Riemmanian distances among species and Goodall's F test 

names_diachlorus <- combn(names(cons_shapes[-1]), m=2)
comb_diachlorus <- combn(cons_shapes[-1], m=2)
r_distances <- NULL
nam <- NULL
for (i in 1:dim(comb_diachlorus)[2]) {
     n <- paste(names_diachlorus[,i][1], "to", names_diachlorus[,i][2])  
     r <- riemdist(as.matrix(comb_diachlorus[,i][1][[1]]), as.matrix(comb_diachlorus[,i][2][[1]]))
     r_distances <- rbind(r_distances, r)
     nam <- c(nam, n)
}
rownames(r_distances) <- nam
colnames(r_distances) <- c("Riemman")
write.csv(r_distances, "Riemmanian_interspecific_distances.csv")

alg_names_comb <- combn(names(alg_inter[-1]), m=2)
comb_alg_diachlorus <- combn(alg_inter[-1], m=2)
name <- NULL
t_goodall <- NULL
for (i in 1:dim(comb_alg_diachlorus)[2]) {
     n <- paste(alg_names_comb[,i][1], "and", alg_names_comb[,i][2]) 
     t <- testmeanshapes(comb_alg_diachlorus[,i][1][1][[1]][[1]], comb_alg_diachlorus[,i][2][1][[1]][[1]], resamples = 1000, replace = FALSE, scale= TRUE)
     t_value <- c(t$G, t$G.pvalue)
     t_goodall <- rbind(t_goodall, t_value)
     name <- c(name, n)
}
rownames(t_goodall) <- name
colnames(t_goodall) <- c("G", "G_pvalue")
write.csv(t_goodall, "Goodall_interspecific_significance.csv")
