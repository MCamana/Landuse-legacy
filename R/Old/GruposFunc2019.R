###########################################
#### Classificando os grupos funcionais ###
###########################################
getwd()
library(vegan)
func<-read.table("traits_funcional.txt", h=T, sep="\t", row.names = 1)
head(func)
func<-func[,!colnames(func) %in% c("Fonte.dos.dados.Morfolgicos","Fonte.de.dados.habitat" )]
peixes<-read.table('Biotic.txt', h=T)

C.peixes<-peixes[rownames(peixes)!="CCA1303",]
peixes1<-C.peixes
vecEX<-c("CGR1903","CGR1904","QUA1403","SPS1404","SAT1403","CSR1904","SPS1401","CCA1301","BAG1502","CCE1303","BAG1503","CCE1301")
peixes1<-peixes1[!rownames(peixes1) %in% vecEX,]

nomescomp<-colnames(peixes1)

vecNAMES<-make.cepnames(colnames(peixes1))
substNAMES1<-c("Astyfasc","Astyfasc.1","Asty1","Asty2","Crenscot","Crenscot.1","Gymnn","Scleminu","Scleminu.1","Scle4","Scle5","Scle3")
substNAMES2<-c("Astyfasc1","Astyfasc","Astysp1","Astysp2","Crenscot1","Crenscot","Gymnspn","Scleminu1","Scleminu","Sclesp4","Sclesp5","Sclesp3")
vecNAMES[vecNAMES %in% substNAMES1]<-substNAMES2
cbind(colnames(peixes1),as.data.frame(vecNAMES))
colnames(peixes1)<-vecNAMES


library(cluster)
library(vegan)
####Clustering
colnames(func)

##EUCLIDEAN
funcMORFO<-scale(func[,11:24])
funcMORFO<-as.data.frame(funcMORFO)
funcMORFO<-cbind(funcMORFO,func[,1:6])
rownames(funcMORFO)<-vecNAMES[!vecNAMES %in% c("Acespant","Aparaffi","Lorimela")]
#distEUCLI <- dist(funcMORFO)
distGOWER <- daisy(funcMORFO, metric = "gower")

hc2 <- hclust(distGOWER, method =  "ward.D2")
mod<-capscale(distGOWER~1, na.action=na.exclude)
summary(mod)
x11()
plot(mod, type="n")
text(mod, col="blue", cex=0.8)
## Limited output of 'summary'
head(summary(mod), tail=2)
# or correlation-based scores in PCA/RDA
scrs <- scores(capscale(distGOWER~1))

###GOWER
# rownames(func)<-make.cepnames(rownames(func))
#distGOWER<-daisy(func, metric = "gower")
# #hc2 <- hclust(distGOWER, method = "average")
# hc2 <- hclust(distGOWER, method =  "ward.D2")
# mod<-capscale(distGOWER~1)
# summary(mod)
# plot(mod, type="n")
# text(mod, col="blue", cex=0.8)
# ## Limited output of 'summary'
# head(summary(mod), tail=2)
# ## or correlation-based scores in PCA/RDA
# scrs <- scores(capscale(distGOWER~1, correlation = TRUE))

#Foram utilizados 2 critérios para separar grupos do cluster Silhueta e K
###Which number of groups
source('calinski.R')
### NUMBER of groups between 0 and 20.
ntest <- 20
res <- rep(0,ntest - 1)
for (i in 2:ntest){
  fac <- cutree(hc2, k = i)
  res[i-1] <- calinski(tab=distGOWER, fac = fac)[1]
}
x11()
par(mfrow=c(2,1))
plot(2:ntest, res, type='b', pch=20, xlab="Number of groups", ylab = "C-H index")
plot(3:ntest, diff(res), type='b', pch=20, xlab="Number of groups", ylab = "Diff in C-H index")
#Indicou 18 grupos funcionais

library(dendextend)
png("FIGS/clusterGOWER.png", width = 6, height =10, units = 'in', res = 600)
#plot(hc2, cex=0.3)
par(cex=0.3, mar=c(5, 8, 4, 1))
#plot(as.dendrogram(hc2), xlab = "", sub="", ylab = "Gower distance", horiz = T, cex.axis=1)
plot(as.dendrogram(hc2), xlab="", ylab="", main="", sub="", axes=FALSE, horiz = T)
par(cex=0.5)
title(ylab="Gower distance")
axis(2, cex.axis=1)
#groups <- cutree(hc2, k=18)
rect.dendrogram(as.dendrogram(hc2),k=18,border = "red", lwd = 0.5, lty = 3, lower_rect=-0.3, xpd=T, horiz = T)
#rect.hclust(hc2, k=18, border="red")
dev.off()
#### Quantos grupos usar? Usar método de largura da silhueta






#####################################
# sil_width <- c(NA)
# for(i in 2:12){
#   pam_fit <- pam(distGOWER,
#                  diss = TRUE,
#                  k = i)
#   sil_width[i] <- pam_fit$silinfo$avg.width
#   
# }
# 
# # # Plot sihouette width (higher is better)
# # x11(height=5,width=10)
# # plot(1:12, sil_width,
# #      xlab = "Number of clusters",
# #      ylab = "Silhouette Width")
# # lines(1:12, sil_width)
# # abline(v =18)
##########################################


#Vertical plot
# #install.packages("ape")
# library("ape")
# x11(height=10,width=6)
# par(mfrow=c(1,1), mar=c(0, 0, 0, 0))
# plot(as.phylo(hc2), cex = 0.4, label.offset = 0.002, no.margin=T)
# clus35 = cutree(hc2, 9)
# rect.hclust(hc2, k=9, border="red")
#dev.copy2pdf(device = x11, file="clusterCORTE.pdf")

spe.group2 <- as.factor(cutree(hc2, k = 18))
FGs<-rep("FG",length(levels(spe.group2)))
FGs<-paste0(FGs, levels(spe.group2))

levels(spe.group2) <- FGs
spe.group2 <- factor(spe.group2, levels=FGs)
peixesTrans<-t(peixes1)
peixesTrans<-data.frame(peixesTrans)
peixesTrans<-peixesTrans[!rownames(peixesTrans) %in% c("Acespant","Aparaffi","Lorimela"),]
length(spe.group2)
nrow(peixesTrans)
peixesTrans$sp.grupo<-spe.group2
peixesTrans[,47:48]
ncol(peixesTrans)
nomescomp<-nomescomp[!nomescomp %in% c("Acestrorhynchus.pantaneiro","Apareiodon.affinis","Loricariichthys.melanocheilus")]
tabGFS<-cbind (as.data.frame(nomescomp),peixesTrans)
write.table(tabGFS, "GFsGOWER.txt", sep=",", row.names = T)

#Agrupando soma das abundâncias por grupo funcional
library("plyr")
gruposFUN<-ddply(peixesTrans, "sp.grupo", numcolwise(sum)) 
str(gruposFUN)
rownames(gruposFUN)<-gruposFUN[,1]
gruposFUN<-gruposFUN[,!colnames(gruposFUN) %in% "sp.grupo"]
#write.table(gruposFUN, "gruposFUN1.txt")
write.table(gruposFUN, "gruposFUN2.txt")

