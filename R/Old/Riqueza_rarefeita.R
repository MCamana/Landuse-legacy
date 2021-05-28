##Neste script vou gerar a riqueza rarefeita dos pontos amostrados nos projetos PPBIO e NEXUS

#selecionamos para onde ir? o os arquivos do R
setwd ('C:/Users/Ecologia/Documents/Usu?rios/Camana/Mestrado/Resultados/An?lises')

#Caminho padr?o dos arquivos
getwd()

#Vou carregar a matriz de abund?ncia para a an?lise
dados.riqueza<-read.table('rarefeita.txt', h=T)

names(dados.riqueza)

dados.riqueza<-t(dados.riqueza)

#carregamos pacote vegan para rerefa??o
library (vegan)

#Riqueza total
riqTOT<-specnumber(dados.riqueza)

#Riqueza rarefeita
rarecurve (dados.riqueza, label=F)
min(rowSums(dados.riqueza))
riqRAR<-rarefy (dados.riqueza, 70)

#Agora vou exportar o documento como um txt
#write.table(riqueza.final, ('C:/Users/Ecologia/Documents/Usu?rios/Camana/Mestrado/Resultados/An?lises/riqueza_rarefeita.txt'))
