#############################################################################################################################################################
#  Ambrosio Torres & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia) #
#  Research: Geometric wing variation in the taxonomic recognition of species of the genus Diachlorus Osten-Sacken (Diptera: Tabanidae) from Colombia       #
#  Part: README to run the analyses in R                                                                                                                    #  
#  28th January 2015                                                                                                                                        # 
#############################################################################################################################################################

### 1. Species and areas used in the analysis
*Species: D. curvipes, D. fuscistigma, D. jobbinsi, D. leucotibialis, D. leticia, D. nuneztovari
*Areas: Amazonas (PPN Amacayacu), Caquetá (PNN Chiribiquete) , Chocó (PNN Utría), Meta (PNN La Macarena), Putumayo (PNN La Paya), Vaupés (EMI Caparú).


### 2. Groups we used to do comparisons and quantify geometrical intraspecific variation
* adc -> Amazonas-Diachlorus curvipes
* adf -> Amazonas-Diachlorus fuscistigma
* adj -> Amazonas-Diachlorus jobbinsi
* adl -> Amazonas-Diachlorus leucotibialis
* adn -> Amazonas-Diachlorus nuneztovari
* cdf -> Caquetá-Diachlorus fuscistigma
* cdj -> Caquetá-Diachlorus jobbinsi
* cdl -> Caquetá-Diachlorus leucotibialis
* chdc -> Chocó-Diachlorus curvipes
* chdj -> Chocó-Diachlorus jobbinsi
* mdc -> Meta-Diachlorus curvipes
* mdf -> Meta-Diachlorus fuscistigma
* pdc -> Putumayo-Diachlorus curvipes
* pdf -> Putumayo-Diachlorus fuscistigma
* pdj -> Putumayo-Diachlorus jobbinsi
* pdl -> Putumayo-Diachlorus leucotibialis
* pdn -> Putumayo-Diachlorus nuneztovari
* vdf -> Vaupés-Diachlorus fuscistigma
* vdt -> Vaupés-Diachlorus leticia


### 3. Data and Scripts
* "interspecific_var_diachlorus.R" and "intraspecific_var_diachlorus.R" allow to run the data "datos_ala_derecha" 
* "analysis_simetry.R" runs the data "datos_alas_derecha_izquierda_asimetria"
* "mantel_diachlorus.R" runs the data "test_mantel"

### 4. Packages and functions needed to perform the analysis
* 'colorspace' package v. 1.2-4 --- http://cran.r-project.org/web/packages/colorspace/index.html
* 'cluster' package v. 1.15.2 ---- http://cran.r-project.org/web/packages/cluster/index.html
* 'Discriminer' package v. 0.1-29 ---- http://cran.r-project.org/web/packages/DiscriMiner/index.html
* 'ellipse' package v. 0.3-8 ---- http://cran.r-project.org/web/packages/ellipse/index.html
* 'geomorph' package v. 2.1.1 ---  http://cran.r-project.org/web/packages/geomorph/index.html
* 'ggplot2' package v. 1.0.0 ---- http://cran.r-project.org/web/packages/ggplot2/index.html
* 'mclust' package v. 4.3 ---- http://cran.r-project.org/web/packages/mclust/index.html
* 'shapes' package v. 1.1-9 --- http://cran.r-project.org/web/packages/shapes/index.html
* 'vegan' package v. 2.0-10 ---- http://cran.r-project.org/web/packages/vegan/index.html
* "colLab" function, modified by us to coloring edgePar and labels of a dendrogram ---- https://stat.ethz.ch/R-manual/R-patched/library/stats/html/dendrapply.html and http://stackoverflow.com/
