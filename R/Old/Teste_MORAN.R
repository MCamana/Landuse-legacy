coords<-read.csv("ptos_coordenadas.csv", h=T, sep=",")
head(coords)

#Carregar antes os arquivos em Analises.R
dat19<-read.table ('dados_cobertura19.txt', h=T) ##Buscamos os dados brutos na pasta
str(dat19)
#Vou carregar a matriz de abund?ncia para a an?lise
peixes<-read.table('Biotic.txt', h=T)
#carregamos pacote vegan para rerefa??o


#Removendo outlier em área
dat19<-dat19[rownames(dat19)!="QUA1403",]
peixes<-peixes[rownames(peixes)!="QUA1403",]
coords<-coords[coords$Codigo_Ca!="QUA1403",]

library (vegan)
#Riqueza total
riqTOT<-specnumber(peixes)

nrow(coords)
nrow(dat19)

#conferindo de os nomes são iguais
coords$Codigo_Ca==rownames(dat19)
coords$Codigo_Ca==rownames(peixes)

#mandando dados para a mesma planilha
dat19$RIQ<-riqTOT
dat19$x<-coords$X
dat19$y<-coords$Y

# datTESTE<-dat19[rownames(dat19)!="QUA1403",]
# coordsTESTE<-coords[coords$Codigo_Ca!="QUA1403",]

#Gerando um modelo exemplo pra testar autocorrelaçao espacial
mod0 <- lm (RIQ~scale(X1993)+bacia+scale(area), data=dat19)
summary(mod0)

#Testando se tem autocorrelacao com o indice de Moran
library(ape)
geo<-cbind(coords$X, coords$Y)
samples.dist <- as.matrix( dist(geo) )
samples.dist.inv <- 1/samples.dist
diag(samples.dist.inv) <- 0
Moran.I(residuals(mod0), samples.dist.inv ,alternative="two.sided")
#Sim, temos autocorrelaco positiva espacial na riqueza

################ Selecionando eixos de PCNM espaciais que explicam riqueza
library(vegan) # will be used for PCNM
colnames(geo)<-c("x","y")
# all paiwise euclidean distances between the cells
xy.dist <- dist(geo)

# PCNM axes of the dist. matrix (from 'vegan' package)
pcnm.axes <- pcnm(xy.dist)$vectors
pcnm.axes<-as.data.frame(pcnm.axes)

#Selecionando eixos significativos
modSPA<-lm(dat19$RIQ~., data=pcnm.axes)

########Criando funcoes para a selecao backward de p valor ##
has.interaction <- function(x,terms){
  out <- sapply(terms,function(i){
    sum(1-(strsplit(x,":")[[1]] %in% strsplit(i,":")[[1]]))==0
  })
  return(sum(out)>0)
}
model.select <- function(model,keep,sig=0.05,verbose=F){
  counter=1
  # check input
  if(!is(model,"lm")) stop(paste(deparse(substitute(model)),"is not an lm object\n"))
  # calculate scope for drop1 function
  terms <- attr(model$terms,"term.labels")
  if(missing(keep)){ # set scopevars to all terms
    scopevars <- terms
  } else{            # select the scopevars if keep is used
    index <- match(keep,terms)
    # check if all is specified correctly
    if(sum(is.na(index))>0){
      novar <- keep[is.na(index)]
      warning(paste(
        c(novar,"cannot be found in the model",
          "\nThese terms are ignored in the model selection."),
        collapse=" "))
      index <- as.vector(na.omit(index))
    }
    scopevars <- terms[-index]
  }
  
  # Backward model selection : 
  
  while(T){
    # extract the test statistics from drop.
    test <- drop1(model, scope=scopevars,test="F")
    
    if(verbose){
      cat("-------------STEP ",counter,"-------------\n",
          "The drop statistics : \n")
      print(test)
    }
    
    pval <- test[,dim(test)[2]]
    
    names(pval) <- rownames(test)
    pval <- sort(pval,decreasing=T)
    
    if(sum(is.na(pval))>0) stop(paste("Model",
                                      deparse(substitute(model)),"is invalid. Check if all coefficients are estimated."))
    
    # check if all significant
    if(pval[1]<sig) break # stops the loop if all remaining vars are sign.
    
    # select var to drop
    i=1
    while(T){
      dropvar <- names(pval)[i]
      check.terms <- terms[-match(dropvar,terms)]
      x <- has.interaction(dropvar,check.terms)
      if(x){i=i+1;next} else {break}              
    } # end while(T) drop var
    
    if(pval[i]<sig) break # stops the loop if var to remove is significant
    
    if(verbose){
      cat("\n--------\nTerm dropped in step",counter,":",dropvar,"\n--------\n\n")              
    }
    
    #update terms, scopevars and model
    scopevars <- scopevars[-match(dropvar,scopevars)]
    terms <- terms[-match(dropvar,terms)]
    
    formul <- as.formula(paste(".~.-",dropvar))
    model <- update(model,formul)
    
    if(length(scopevars)==0) {
      warning("All variables are thrown out of the model.\n",
              "No model could be specified.")
      return()
    }
    counter=counter+1
  } # end while(T) main loop
  return(model)
}
##############################################################

#Selecionando os eixos
eixosSIG<-model.select (modSPA)
summary(eixosSIG)
nomesEIXOS<-names(eixosSIG$coefficients)[-1]

#eixosSPAT<-pcnm.axes[,c("PCNM1","PCNM3","PCNM11")]
eixosSPAT<-pcnm.axes[,nomesEIXOS]

#Adicionando a uma matrix nova
novaDATA<-cbind(dat19,eixosSPAT)

#Testando se os residuos do modelo com os eixos espaciais não sao mais correlacionados espacialmente 
mod0 <- lm (scale(RIQ)~scale(X1993)+scale(area)+bacia+scale(PCNM1)+scale(PCNM3)+scale(PCNM11), data=novaDATA)
anova(mod0)
summary(mod0)
Moran.I(residuals(mod0), samples.dist.inv ,alternative="two.sided")

#RIQUEZA TOTAL
#library(nlme)
RIQmod<-list()
for (i in 1:29){
  a<-colnames(novaDATA)[i]
  dataROD<-novaDATA[,c("RIQ",a,"area","bacia","x","y","PCNM1","PCNM3","PCNM11")]
  colnames(dataROD)[2]<-"LU"
  
  mod0 <- lm(scale(RIQ)~scale(LU)+scale(area)+bacia+scale(PCNM1)+scale(PCNM3)+scale(PCNM11),
              data=dataROD)

  
  #Teste usando modelos GLS
  # skip_to_next <- FALSE
  # tryCatch(mod0 <- gls(RIQ~LU+area+bacia,
  #                      correlation = corExp(form = ~x + y, nugget = TRUE),
  #                      data=dataROD), error = function(e) { skip_to_next <<- TRUE})
  # if(skip_to_next) { next }   
  
  RIQmod[i] <- coef(mod0)[2]
  names(RIQmod)[i] <- colnames(novaDATA)[i]
}


str(RIQmod)
RIQmod<-data.frame(RIQmod)
RIQmod<-t(RIQmod)
colnames(RIQmod)<-"Effect.size"

x11(height=5,width=6)
par(mfrow=c(1,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
plot(RIQmod, type="l", xlab = "", ylab="Effect size (richness)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.2)
axis(side = 1, at = seq(1,29,1), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
mtext("Years", las=0, cex=1.2, 1, 2.8, outer=F)

mod0 <- lm(scale(RIQ)~scale(X1993)+scale(area)+bacia+scale(PCNM1)+scale(PCNM3)+scale(PCNM11),
           data=novaDATA)
anova(mod0)


#Matriz de peixes SEM os pontos de 2019 
peixeslog<-log(peixes+1)
##LOOP para obter R2 da RDA
geo<-cbind(coords$X, coords$Y)
samples.dist <- as.matrix( dist(geo) )
samples.dist.inv <- 1/samples.dist
diag(samples.dist.inv) <- 0
nrow(coords)
nrow(dat19)

# Matrix of PCNM variables
# Run PCNM
library(vegan) # will be used for PCNM
colnames(geo)<-c("x","y")
# all paiwise euclidean distances between the cells
xy.dist <- dist(geo)

# PCNM axes of the dist. matrix (from 'vegan' package)
pcnm.axes <- pcnm(xy.dist)$vectors
pcnm.axes<-as.data.frame(pcnm.axes)

PCNM <- rda(peixeslog,pcnm.axes)
anova.cca(PCNM)
R2a <- RsquareAdj(PCNM)$adj.r.squared
library(adespatial)
PCNM.fwd <- forward.sel(peixeslog, as.matrix(pcnm.axes))
PCNMsel<-pcnm.axes[,PCNM.fwd$variables]


COMPmod<-list()
for (i in 1:29){
  nomeROD<-colnames(dat19)[i]
  envROD<-dat19[,c(nomeROD,"area","bacia")]
  colnames(envROD)<-c("LU","area","bacia")
  #mod0 <- rda(peixeslog,envROD, PCNMsel)
  #mm1 <- model.matrix(~ LU + area, envROD)
  #mm2 <- model.matrix(~ bacia, dat19)
  #as<-varpart(peixeslog, envROD, PCNMsel)
  #as2<-as$part$fract[3]
  #as3<-as2$Adj.R.squared[3]
  #COMPmod[i] <-as3
  rda.result <- rda(peixeslog ~ LU + Condition (area) + Condition (bacia) + Condition (as.matrix(PCNMsel)), data = envROD)
  COMPmod[i] <- RsquareAdj(rda.result)$r.squared
  names(COMPmod)[i] <- colnames(dat19)[i]
}

COMPmod<-t(data.frame(COMPmod))

# mod1993 <- rda(peixeslog~dat19$X1993+area+Condition(bacia), data=dat19)
# anova(mod1993, by="terms")
# 
# mod2013 <- rda(peixeslog~dat19$X2013+area+Condition(bacia), data=dat19)
# anova(mod2013, by="terms")
x11(height=5,width=6)
par(mfrow=c(1,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
plot(COMPmod, type="l", xlab = "", ylab="Effect size (composition)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.2)
axis(side = 1, at = seq(1,29,1), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
mtext("Years", las=0, cex=1.2, 1, 2.8, outer=F)


