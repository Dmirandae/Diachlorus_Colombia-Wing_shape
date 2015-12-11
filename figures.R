## Ambrosio Torres & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Research: Geometric wing variation in the taxonomic recognition of species of the genus Diachlorus Osten-Sacken (Diptera: Tabanidae) from Colombia
## Part: intraspecific and interspecific variation in taxonomy of the species of genus Diachlorus Osten-Sacken (Diptera: Tabanidae) of Colombia
## R version 3.1.2 & Rstudio 0.97.551
## 23th January 2015 

######################################################################################################################################################################################
######################################################################################################################################################################################
## GRAPHICS FOR THE MANUSCRIPT

######################################################################################################################################################################################
## FIGURE 2
a<- cbind(pca_diachlorus$pc.scores, factors_diachlorus)
a<-data.frame(a)
a$factors_diachlorus <- as.factor(a$factors_diachlorus)
df_ell <- data.frame()
for(g in levels(a$factors_diachlorus)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$factors_diachlorus==g,], ellipse(cor(PC1, PC2), 
                    scale=c(sd(PC1),sd(PC2)), 
                    centre=c(mean(PC1),mean(PC2))))),factors_diachlorus=g))
}
fig2 <- ggplot(data=a, aes(x=PC1, y=PC2,colour=factors_diachlorus, shape=factors_diachlorus)) + geom_point(size=2) 
fig2 <- fig2 + geom_path(data=df_ell, aes(x=x, y=y,colour=factors_diachlorus), size=0.5, linetype=2) + scale_colour_brewer(palette = "Set1") 
fig2 <- fig2 + scale_shape_manual(values=c(18, 15, 16, 17, 8, 2))
fig2 <- fig2 + theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) + theme(legend.position="none") + scale_y_continuous(name="PC2 (11%)") 
fig2 <- fig2 + scale_x_continuous(name="PC1 (45%)") + theme(panel.background = element_rect(fill = "grey"))
tiff("Fig2.tiff", height = 95, width = 174, units = 'mm', type="cairo", res=1000, compression = "lzw")
fig2 + theme(text=element_text(family="Calibri"))
dev.off()
######################################################################################################################################################################################

######################################################################################################################################################################################
## FIGURE 3
tiff("Fig3.tiff", height = 95, width = 174, units = 'mm', type="cairo", res=1000, compression = "lzw")
par(mar=c(0.5, 4, 0, 0))
par(oma= c(0,0,0,0)); par(family="Calibri")
op = par(bg = "gray")
par(lwd=1.2)
plot(clusDendro,horiz=F,axes=T, ylab= "Height", cex.axis=1,cex.lab=1)
dev.off()
######################################################################################################################################################################################

######################################################################################################################################################################################
## FIGURE 4
tiff("Fig4.tiff", height = 95, width = 174, units = 'mm', type="cairo", res=1000, compression = "lzw")
multiplot((p_curvi + theme(text=element_text(family="Calibri")) + theme(panel.background = element_rect(fill = "grey")) + theme(axis.title.x = element_blank()) + theme(legend.position="none") + theme(axis.text=element_text(size=10), axis.title=element_text(size=12))),(p_fusci + theme(text=element_text(family="Calibri")) + theme(panel.background = element_rect(fill = "grey")) + theme(legend.position="none") + theme(axis.text=element_text(size=10), axis.title=element_text(size=12))), (p_jobi + theme(text=element_text(family="Calibri")) + theme(panel.background = element_rect(fill = "grey")) + theme(legend.position="none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) + theme(axis.text=element_text(size=10))), (p_leuco + theme(text=element_text(family="Calibri")) + theme(panel.background = element_rect(fill = "grey")) + theme(legend.position="none") + theme(axis.title.y = element_blank()) + theme(axis.text=element_text(size=10))), cols=2)
dev.off()
######################################################################################################################################################################################

######################################################################################################################################################################################
## FIGURE 5
tiff("Fig5.tiff", height = 75, width = 120, units = 'mm', type="cairo", res=1000, compression = "lzw")
par(mar=c(0.8, 4.5, 0, 0)); par(family="Calibri")
par(oma= c(0,0,0,0))
op = par(bg = "gray")
par(lwd=1.5)
plot(clusDendro_jobbinsi,horiz=F,axes=T, ylab= "Height", cex.axis=1.3,cex.lab=1.7)
dev.off()
######################################################################################################################################################################################

######################################################################################################################################################################################
## FIGURE 6
tiff("Fig6a.tiff", height = 90, width = 90, units = 'mm', type="cairo", res=1000, compression = "lzw")
par(oma=c(0,0,0,0))
par(mar=c(2, 2, 0, 0))
tpsgrid(cons_shapes$alg_d_leucotibialis, cons_shapes$alg_d_fuscistigma, mag= 3, ext= 0.05, ybegin= -0.25, xwidth=1, cex=0.7, ngrid=16)
dev.off()

tiff("Fig6b.tiff", height = 90, width = 90, units = 'mm', type="cairo", res=1000, compression = "lzw")
par(oma=c(0,0,0,0))
par(mar=c(2, 2, 0, 0))
tpsgrid(cons_shapes$alg_d_leucotibialis, cons_shapes$alg_d_nuneztovari, mag= 3, ext= 0.05, ybegin= -0.25, xwidth=1, cex=0.7, ngrid=16)
dev.off()

tiff("Fig6c.tiff", height = 90, width = 90, units = 'mm', type="cairo", res=1000, compression = "lzw")
par(oma=c(0,0,0,0))
par(mar=c(2, 2, 0, 0))
tpsgrid(cons_shapes$alg_vdf, cons_shapes$alg_pdf, mag= 3, ext= 0.05, ybegin= -0.25, xwidth=1, cex=0.7, ngrid=16)
dev.off()

tiff("Fig6d.tiff", height = 90, width = 90, units = 'mm', type="cairo", res=1000, compression = "lzw")
par(oma=c(0,0,0,0))
par(mar=c(2, 2, 0, 0))
tpsgrid(cons_shapes$alg_adj, cons_shapes$alg_chdj, mag= 3, ext= 0.05, ybegin= -0.25, xwidth=1, cex=0.7, ngrid=16)
dev.off()
######################################################################################################################################################################################
