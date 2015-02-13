## Ambrosio Torres & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Research: Wing shape variation in the taxonomic recognition of species of Diachlorus Osten-Sacken (Diptera: Tabanidae) from Colombia
## Part: Analysis of directional asymmetry of wing in genus Diachlorus from Colombia
## 23th January 2015
## R version 3.1.2 & Rstudio 0.97.551 

library (geomorph) ## Load 'geomorph' package v. 2.1.1 ---  http://cran.r-project.org/web/packages/geomorph/index.html
library (shapes) ## Load 'shapes' package v. 1.1-9 --- http://cran.r-project.org/web/packages/shapes/index.html

#  D. curvipes
curvi_sim<- readland.tps("iz_der_curvi.tps", specID = "ID", warnmsg = T) ## load the data of the .tps file for the species D. curvipes for both, right and left wing
curvi_sim_gpa <- gpagen(curvi_sim) ## make generalized procrustes analysis
pca_curvi_sim <- plotTangentSpace(curvi_sim_gpa$coords, label=NULL, verbose =T, warpgrids=F) ## PCA for geometric shape based on aligned configurations of landmarks
curvi_sim_der<- readland.tps("derecha_curvipes.tps", specID = "ID", warnmsg = T) ## load the data of the .tps file for the species D. curvipes for rigth wing 
curvi_sim_izq<- readland.tps("izquierda_curvipes.tps", specID = "ID", warnmsg = T) ## load the data of the .tps file for the species D. curvipes for left wing
curvi_sim__der_gpa <- gpagen(curvi_sim_der) ## make generalized procrustes analysis right wing
curvi_sim__izq_gpa <- gpagen(curvi_sim_izq) ## make generalized procrustes analysis left wing
mshape_der_curvi<- mshape(curvi_sim__der_gpa$coords) ## calculate consensus shape for right wing
mshape_izq_curvi<- mshape(curvi_sim__izq_gpa$coords) ## calculate consensus shape for left wing
riemdist(mshape_der_curvi, mshape_izq_curvi) ## Riemmanian distances between right and left wing shapes
# 0.005972782
testmeanshapes(curvi_sim__der_gpa$coords, curvi_sim__izq_gpa$coords, resamples = 1000, replace = FALSE, scale= TRUE) ## Goodall's F test to shapes 
#$G
#[1] 0.7854165
#$G.pvalue
#[1] 0.6333666


## D. fuscistigma
fusci_sim<- readland.tps("iz_der_fusci.tps", specID = "ID", warnmsg = T)
fusci_sim_gpa <- gpagen(fusci_sim)
pca_fusci_sim <- plotTangentSpace(fusci_sim_gpa$coords, label=NULL, verbose =T, warpgrids=F)
fusci_sim_der<- readland.tps("derecha_fuscistigma.tps", specID = "ID", warnmsg = T)
fusci_sim_izq<- readland.tps("izquierda_fuscistigma.tps", specID = "ID", warnmsg = T)
fusci_sim__der_gpa <- gpagen(fusci_sim_der)
fusci_sim__izq_gpa <- gpagen(fusci_sim_izq)
mshape_der_fusci<- mshape(fusci_sim__der_gpa$coords)
mshape_izq_fusci<- mshape(fusci_sim__izq_gpa$coords)
riemdist(mshape_der_fusci, mshape_izq_fusci)
# 0.004091286
testmeanshapes(fusci_sim__der_gpa$coords, fusci_sim__izq_gpa$coords, resamples = 1000, replace = FALSE, scale= TRUE)
#$G
#[1] 0.4147391
#$G.pvalue
#[1] 0.958042


## D. jobbinsi
jobi_sim<- readland.tps("iz_der_jobi.tps", specID = "ID", warnmsg = T)
jobi_sim_gpa <- gpagen(jobi_sim)
pca_jobi_sim <- plotTangentSpace(jobi_sim_gpa$coords, label=NULL, verbose =T, warpgrids=F)
jobi_sim_der<- readland.tps("derecha_jobbinsi.tps", specID = "ID", warnmsg = T)
jobi_sim_izq<- readland.tps("izquierda_jobbinsi.tps", specID = "ID", warnmsg = T)
jobi_sim__der_gpa <- gpagen(jobi_sim_der)
jobi_sim__izq_gpa <- gpagen(jobi_sim_izq)
mshape_der_jobi<- mshape(jobi_sim__der_gpa$coords)
mshape_izq_jobi<- mshape(jobi_sim__izq_gpa$coords)
riemdist(mshape_der_jobi, mshape_izq_jobi)
# 0.004099021
testmeanshapes(jobi_sim__der_gpa$coords, jobi_sim__izq_gpa$coords, resamples = 1000, replace = FALSE, scale= TRUE)
#$G
#[1] 0.3047416
#$G.pvalue
#[1] 0.98002


## D. leticia
leti_sim<- readland.tps("iz_der_leti.tps", specID = "ID", warnmsg = T)
leti_sim_gpa <- gpagen(leti_sim)
pca_leti_sim <- plotTangentSpace(leti_sim_gpa$coords, label=NULL, verbose =T, warpgrids=F)
leti_sim_der<- readland.tps("derecha_leticia.tps", specID = "ID", warnmsg = T)
leti_sim_izq<- readland.tps("izquierda_leticia.tps", specID = "ID", warnmsg = T)
leti_sim__der_gpa <- gpagen(leti_sim_der)
leti_sim__izq_gpa <- gpagen(leti_sim_izq)
mshape_der_leti <- mshape(leti_sim__der_gpa$coords)
mshape_izq_leti <- mshape(leti_sim__izq_gpa$coords)
riemdist(mshape_der_leti, mshape_izq_leti)
# 0.004920405
testmeanshapes(leti_sim__der_gpa$coords, leti_sim__izq_gpa$coords, resamples = 1000, replace = FALSE, scale= TRUE)
#$G
#[1] 0.5691968
#$G.pvalue
#[1] 0.7862138


## D. leucotibialis
leuco_sim<- readland.tps("iz_der_leuco.tps", specID = "ID", warnmsg = T)
leuco_sim_gpa <- gpagen(leuco_sim)
pca_leuco_sim <- plotTangentSpace(leuco_sim_gpa$coords, label=NULL, verbose =T, warpgrids=F)
leuco_sim_der<- readland.tps("derecha_leucotibialis.tps", specID = "ID", warnmsg = T)
leuco_sim_izq<- readland.tps("izquierda_leucotibialis.tps", specID = "ID", warnmsg = T)
leuco_sim__der_gpa <- gpagen(leuco_sim_der)
leuco_sim__izq_gpa <- gpagen(leuco_sim_izq)
mshape_der_leuco <- mshape(leuco_sim__der_gpa$coords)
mshape_izq_leuco <- mshape(leuco_sim__izq_gpa$coords)
riemdist(mshape_der_leuco, mshape_izq_leuco)
# 0.003413203
testmeanshapes(leuco_sim__der_gpa$coords, leuco_sim__izq_gpa$coords, resamples = 1000, replace = FALSE, scale= TRUE)
#$G
#[1] 0.2977132
#$G.pvalue
#[1] 0.993007


## D. nuneztovari
nune_sim<- readland.tps("iz_der_nune.tps", specID = "ID", warnmsg = T)
nune_sim_gpa <- gpagen(nune_sim)
pca_nune_sim <- plotTangentSpace(nune_sim_gpa$coords, label=NULL, verbose =T, warpgrids=F)
nune_sim_der<- readland.tps("derecha_nuneztovari.tps", specID = "ID", warnmsg = T)
nune_sim_izq<- readland.tps("izquierda_nuneztovari.tps", specID = "ID", warnmsg = T)
nune_sim__der_gpa <- gpagen(nune_sim_der)
nune_sim__izq_gpa <- gpagen(nune_sim_izq)
mshape_der_nune <- mshape(nune_sim__der_gpa$coords)
mshape_izq_nune <- mshape(nune_sim__izq_gpa$coords)
riemdist(mshape_der_nune, mshape_izq_nune)
# 0.005450411
testmeanshapes(nune_sim__der_gpa$coords, nune_sim__izq_gpa$coords, resamples = 1000, replace = FALSE, scale= TRUE)
#$G
#[1] 0.6106513
#$G.pvalue
#[1] 0.7882118
