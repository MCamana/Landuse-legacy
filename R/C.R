####Libraries####
library(here) 
library(iNEXT)
library(ape)

##Carregando os dados
here()
dat19 <- read.table (here::here("data","processed","dados_cobertura19.txt"), header=T)
str(dat19)
t(dat19)
# #ponto QUA1403 tem uma área muito discrepante (134km2), por isso foi retirado das analises
dat19<-dat19[rownames(dat19)!="QUA1403",]
hist(dat19$area)

#Vou carregar a matriz de abund?ncia para a an?lise
peixes<-read.table (here::here("data","processed","biotic.txt"), header=T)
#peixes<-peixes[rownames(peixes)!="QUA1403",]
head(peixes)

#Avaliando o nível de completness das comunidades
out2 <- iNEXT(t(peixes), datatype="abundance")
ggiNEXT(out2, type=2, se=TRUE)
out2$DataInfo #Usamos ponto de corte de nível de sample coverage (completness SC) = 0.95 (95%)
NomeRemover<-out2$DataInfo[out2$DataInfo$SC<0.94,]$site
NomeRemover<-as.character(NomeRemover)

###Filtrando todos pontos menos 1 que nao tem amostragem completa
C.dat19<-dat19[dat19$Pontos!=NomeRemover,]
C.peixes<-peixes[rownames(peixes)!=NomeRemover,]

### Avaliando a correlaçao espacial de Moran ####
coords<-read.csv (here::here("data","processed","ptos_coordenadas.csv"), header=T,  sep=",")
C.coords<-coords[coords$Codigo_Ca!=NomeRemover,]

##Avaliando distribuicao dos pontos
plot(coords$X, coords$Y, cex=0.5, pch=19)
text(coords$X, coords$Y, labels=coords$Codigo_Ca)

#Conferindo o ponto excluido
all(C.coords$Codigo_Ca==C.dat19$Pontos)
all(C.coords$Codigo_Ca==rownames(C.peixes))

##Novo objeto para coordenadas
geo<-cbind(C.coords$X, C.coords$Y)
rownames(geo)<-C.coords$Codigo_Ca
samples.dist <- as.matrix( dist(geo) )
samples.dist.inv <- 1/samples.dist
diag(samples.dist.inv) <- 0

##Dados de morfologia
func<-read.table(here::here("data","processed", "traits_funcional.txt"), header=T, sep="\t", row.names = 1)
func<-func[,!colnames(func) %in% c("Fonte.dos.dados.Morfolgicos","Fonte.de.dados.habitat" )]

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

