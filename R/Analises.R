####Libraries####
library(here) ##todas
library(iNEXT)
library(ape)
library(vegan)
require(fitdistrplus)
library(TITAN2)

##Avaliando os dados e autocorrelacao espacial
riqTOT<-specnumber(C.peixes)
mod0 <- lm (riqTOT~scale(X1993)+bacia+scale(area), data=C.dat19)
summary(mod0)
Moran.I(residuals(mod0), samples.dist.inv ,alternative="greater") #Riqueza está espacialmente autocorrelacionada

#Testando com composicao
peixeslog<-log(C.peixes+1)
PCOAfish<-capscale(peixeslog~1, distance="bray")
summary(PCOAfish)
Eixos<-scores(PCOAfish)$sites
Eixos #Tem autocorrelacao espacial para composicao
modpcoa <- lm (Eixos[,1]~X1993+bacia+log(area), data=C.dat19)
summary(modpcoa)
Moran.I(residuals(modpcoa), samples.dist.inv ,alternative="greater") #Riqueza está espacialmente autocorrelacionada

matMor<-matrix(NA,nrow (Eixos),2)
colnames(matMor)<-c("MoranI","pvalue")
rownames(matMor)<-rownames(Eixos)

#Avaliando autocoorelacao espacial
 for (i in 1:nrow(Eixos1)) {
   Rod<-Eixos1[,1]
   Rod1<-Rod[-i]
   Rod1env<-dat1[-i,]
 
   samples.dist <- as.matrix( dist(geo1[-i,]) )
   samples.dist.inv <- 1/samples.dist
   diag(samples.dist.inv) <- 0
   
   modpcoa <- lm (Rod1~X1993+bacia+log(area), data= Rod1env)
   MorR<-Moran.I(residuals(modpcoa), samples.dist.inv ,alternative="greater") #Riqueza está espacialmente autocorrelacionada
   matMor[i,]<-c(MorR$observed,MorR$p.value)
     }
 matMor<-as.data.frame(matMor)
 teste<-matMor[order(matMor$MoranI),]
 
 plot(teste$MoranI, 1:nrow(teste))
 text(teste$MoranI, 1:nrow(teste), label=rownames(teste), cex = 0.8)

#Criando objetos para novos dados sem autocorrelacao espacial
Eixos1<-Eixos
dat1<-C.dat19
geo1<-geo
peixes1<-C.peixes

vecEX<-c("CGR1903","CGR1904","QUA1403","SPS1404","SAT1403","CSR1904","SPS1401","CCA1301","BAG1502","CCE1303","BAG1503","CCE1301")
peixes1<-peixes1[!rownames(peixes1) %in% vecEX,]
dat1<-dat1[!dat1$Pontos %in% vecEX,]
geo1<-geo[!rownames(geo) %in% vecEX,]
riqTOT1<-specnumber(peixes1)

peixeslog1<-log(peixes1+1)
PCOAfish<-capscale(peixeslog1~1, distance="bray")
summary(PCOAfish)
Eixos1<-scores(PCOAfish)$sites

samples.dist1 <- as.matrix( dist(geo1) )
samples.dist.inv1 <- 1/samples.dist1
diag(samples.dist.inv1) <- 0

modpcoa1 <- lm (Eixos1[,1]~X1993+bacia+log(area), data=dat1)
summary(modpcoa1)
Moran.I(residuals(modpcoa1), samples.dist.inv1 ,alternative="greater") #Riqueza está espacialmente autocorrelacionada

modriq <- lm (riqTOT1~X1993+bacia+log(area), data=dat1)
Moran.I(residuals(modriq), samples.dist.inv1 ,alternative="greater") #Riqueza está espacialmente autocorrelacionada

plot(geo1[,1],geo1[,2], pch=19, cex=0.5)
text(geo1[,1],geo1[,2], label=rownames(geo1), cex = 0.8)

matMor<-matrix(NA,nrow (Eixos1),2)
colnames(matMor)<-c("MoranI","pvalue")
rownames(matMor)<-rownames(Eixos1)

#####FINAL FICAMOS COM peixes1 e dat1
#Riqueza total
riqTOT<-specnumber(peixes1)
plot (riqTOT~dat1$area)
abline(lm(riqTOT~(dat1$area)))

#####Rodando os modelos lineares do efeito da vegetação anual sobre a riqueza de espécies

#Avaliar normalidade dos dados de Riqueza
 normal=fitdist(riqTOT,"norm")
 lnormal=fitdist(riqTOT,"lnorm")
 par(mfrow=c(1,2), mar=c(4,4,2,2))
 cdfcomp(list(normal,lnormal),horizontals=F, lwd=2,addlegend=T,legendtext=c("Normal","LNormal"))
 qqcomp(list(normal,lnormal),addlegend=T,legendtext=c("Normal","Lnormal"))###a Lognormal foi a melhor
 par(mfrow=c(1,1), cex.lab=1, cex.axis=1)

## Conferindo se os nomes das linhas (sitios amostrais) sao iguais
rownames(dat1)<-dat1$Pontos
dat1<-dat1[,!colnames(dat1) %in% "Pontos"]
rownames(dat1)==rownames(peixes1)
str(dat1)

## Agora, criamos modelos de cada gravar o coef. de inclinação padronizado (tamanho do efeito) so da cobertura vegetal isolada
#RIQUEZA TOTAL
RIQmod<-matrix(NA,29,2)
colnames(RIQmod)
colnames(RIQmod)<-c("Effect.size","pvalue")
rownames(RIQmod)<- colnames(dat1)[1:29]

for (i in 1:29){
  mod0 <- lm (scale(riqTOT)~scale(dat1[,i])+bacia+scale(area), data=dat1)
  RIQmod[i,1] <- coef(mod0)[2]
  pval<-anova(mod0)[5][1,]
  RIQmod[i,2] <- pval
  }

####RDA do efeito da vegetacao em diferentes anos sobre a composicao####
#Matriz de peixes SEM os pontos de 2019 
peixeslog<-log(peixes1+1)

##LOOP para obter R2 da RDA
COMPmod<-matrix(NA, 29,2)
colnames(COMPmod)
colnames(COMPmod)<-c("Effect.size","pvalue")
rownames(COMPmod)<- colnames(dat1)[1:29]

for (i in 1:29){
  mod0 <- rda(peixeslog~dat1[,i]+area+Condition(bacia), data=dat1)
  COMPmod[i,1] <- RsquareAdj(mod0)$r.squared
  pval<-anova(mod0, by="terms")$`Pr(>F)`[1]
  COMPmod[i,2] <- pval
  names(COMPmod)[i] <- colnames(dat19)[i]
}

COMPmod<-as.data.frame(COMPmod)
COMPmod
mod2001 <- rda(peixeslog~dat1$X2001+area+Condition(bacia), data=dat1)
anova(mod2001, by="terms")

mod2013 <- rda(peixeslog~dat1$X2013+area+Condition(bacia), data=dat1)
anova(mod2013, by="terms")

##Avaliando o efeito das trajetorias de land use nas comunidades
traject<-read.table(here("data", "processed", "predit_titan.txt"), h=T) # Contem dados de 2019
head(traject)
nrow(traject)
vecEX<-c("CGR1903","CGR1904","QUA1403","SPS1404","SAT1403","CSR1904","SPS1401","CCA1301","BAG1502","CCE1303","BAG1503","CCE1301")
C.traject<-traject[rownames(traject)!="CCA1303",]
traject.1<-C.traject[!rownames(C.traject) %in% vecEX,]
nrow(traject.1)

##Carregando grupos funcionais. Estes dados foram gerados no script D
peixesFUN<-read.table(here("data", "processed", "gruposFUN.txt"), h=T)
peixesFUN<- peixesFUN[,colSums(ifelse(peixesFUN>0,1,0))>=3]
peixesFUN<-t(peixesFUN)

#Excluindo spp e grupos funcionais com menos de 3 ocorrencias
colSums(peixes)
peixesTAX19 <- peixes1[,colSums(ifelse(peixes1>0,1,0))>=3]
peixesFUN19 <- peixesFUN[,colSums(ifelse(peixesFUN>0,1,0))>=3]

ncol(peixes1)
ncol(peixesTAX19)

ncol(peixesFUN19)
ncol(peixesFUN)

#Aqui começamos a conferencia dos dados para a TITAN
LOGpeixesTAX19<-log(peixesTAX19+1)
LOGpeixesFUN19<-log(peixesFUN19+1)

##############
##INSERIR DADOS DO SCRIPT GRUPOS FUNCIONAIS###
#############



#checando os nomes das matrizes
rownames(traject)
rownames(peixes1)
colnames(peixes1)
rownames(peixesFUN)
rownames(dat1)

#Exportando as abreviaturas dos peixes
vecNAMES<-make.cepnames(colnames(peixes1))
substNAMES1<-c("Astyfasc","Astyfasc.1","Asty1","Asty2","Crenscot","Crenscot.1","Gymnn","Scleminu","Scleminu.1","Scle4","Scle5","Scle3")
substNAMES2<-c("Astyfasc1","Astyfasc","Astysp1","Astysp2","Crenscot1","Crenscot","Gymnspn","Scleminu1","Scleminu","Sclesp4","Sclesp5","Sclesp3")
vecNAMES[vecNAMES %in% substNAMES1]<-substNAMES2
cbind(colnames(peixes1),as.data.frame(vecNAMES))
colnames(peixes1)<-vecNAMES

#####Taxonomico####
#%Native cover 1993
plot (degrass1993, LOGpeixesTAX19$Scleminu)
plot (degrass2013, LOGpeixesTAX19$Scleminu)
nrow(LOGpeixesTAX19)
titanNAT93 <- titan(degrass1993, LOGpeixesTAX19, ivTot=F, nBoot = 1000)
rownames(titanNAT93$sppmax)<-colnames(LOGpeixesTAX19)
titanNAT93$sumz.cp
sppNAT93<-titanNAT93$sppmax
plotTaxa(titanNAT93)
sppNAT93[sppNAT93[,"filter"]>0,]


#%Native cover 2013
titanNAT13 <- titan(degrass2013, LOGpeixesTAX19, ivTot=F, nBoot = 1000)
rownames(titanNAT13$sppmax)<-colnames(LOGpeixesTAX19)
titanNAT13$sumz.cp
sppNAT13<-titanNAT13$sppmax
plotTaxa(titanNAT13)
sppNAT13[sppNAT13[,"filter"]>0,]


#MAGNITUDE
magnitude<-(-1*traject.1$mag) #Transformando magnitude em positivo
titanMAG <- titan(magnitude, LOGpeixesTAX19, ivTot=F, nBoot = 1000)
rownames(titanMAG$sppmax)<-colnames(LOGpeixesTAX19)
titanMAG$sumz.cp
sppMAG<-titanMAG$sppmax
plotTaxa(titanMAG, z2=F)
sppMAG[sppMAG[,"filter"]>0,]

#FREQUENCIA
titanFREQ <- titan(traject.1$freq, LOGpeixesTAX19, ivTot=F, nBoot = 1000)
titanFREQ$sumz.cp
sppFREQ<-titanFREQ$sppmax
plotTaxa(titanFREQ)
sppFREQ[sppFREQ[,"filter"]>0,]

#DURACAO T40
titanT40 <- titan(traject.1$T60, LOGpeixesTAX19, ivTot=F, nBoot = 1000)
rownames(titanT40$sppmax)<-colnames(LOGpeixesTAX19)
titanT40$sumz.cp
sppT40<-titanT40$sppmax
plotTaxa(titanT40, z1=F)
sppT40[sppT40[,"filter"]>0,]

#DURACAO T20
titanT20 <- titan(traject.1$T80, LOGpeixesTAX19, ivTot=F, nBoot = 1000)
rownames(titanT20$sppmax)<-colnames(LOGpeixesTAX19)
titanT20$sumz.cp
plotTaxa(titanT20, z1=F)
sppT20<-titanT20$sppmax
plotTaxa(titanT20, z1=F)
sppT20[sppT20[,"filter"]>0,]


##### Funcional######
FUN.titanNAT93 <- titan(degrass1993, LOGpeixesFUN19, ivTot=F, nBoot = 1000)
FUN.titanNAT93$sumz.cp
plotTaxa(FUN.titanNAT93, legend = T, z2=F)
plotSumz(FUN.titanNAT93, filter = T)
FUNnat93<-FUN.titanNAT93$sppmax
plotTaxa(FUN.titanNAT93,z2=F)
FUNnat93[FUNnat93[,"filter"]>0,]

#%Native cover 2013
FUN.titanNAT13 <- titan(degrass2013, LOGpeixesFUN19, ivTot=F, nBoot = 1000)
FUN.titanNAT13$sumz.cp
plotTaxa(FUN.titanNAT13, z2=F)
FUNnat13<-FUN.titanNAT13$sppmax
plotTaxa(FUN.titanNAT13)
FUNnat13[FUNnat13[,"filter"]>0,]

#MAGNITUDE
magnitude<-(-1*traject.1$mag) #Transformando magnitude em positivo
FUN.titanMAG <- titan(magnitude, LOGpeixesFUN19, ivTot=F, nBoot = 1000)
FUN.titanMAG$sumz.cp
plotTaxa(FUN.titanMAG, z2=F)
plotSumz(FUN.titanMAG, filter = F)
FUNMAG<-FUN.titanMAG$sppmax
plotTaxa(FUN.titanMAG,z2=F)
FUNMAG[FUNMAG[,"filter"]>0,]

#FREQUENCIA
FUN.titanFREQ <- titan(traject.1$freq, LOGpeixesFUN19, ivTot=F, nBoot = 1000)
names(FUN.titanFREQ)
plotTaxa(FUN.titanFREQ, z2=F)
plotSumz(FUN.titanFREQ, filter = T)

#DURACAO T40
FUN.titanT40 <- titan(traject.1$T60, LOGpeixesFUN19, ivTot=F, nBoot = 1000)
FUN.titanT40$sumz.cp
plotTaxa(FUN.titanT40)
plotSumz(FUN.titanT40, filter = T)
FUNT40<-FUN.titanT40$sppmax
plotTaxa(FUN.titanT40,z2=F)
FUNT40[FUNT40[,"filter"]>0,]

#DURACAO T20
FUN.titanT20 <- titan(traject.1$T80, LOGpeixesFUN19, ivTot=F, nBoot = 1000)
FUN.titanT20$sumz.cp
plotTaxa(FUN.titanT20, z2=F)
plotSumz(FUN.titanT20, filter = T)
FUNT20<-FUN.titanT20$sppmax
plotTaxa(FUN.titanT20,z2=F)
FUNT20[FUNT20[,"filter"]>0,]