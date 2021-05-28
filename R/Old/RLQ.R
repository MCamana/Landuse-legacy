#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#                                                            #
#               RLQ AND FOURTH-CORNER ANALYSES               #
#                                                            #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

#Packages
library (vegan)
library(ade4)

#Traits matrix
traits.sp <-read.table ('traits_sp.txt', h=T)
rownames(traits.sp)

### Species list matrix with abundances
spp_ab <- read.table ('species_list.txt', h=T)
spp_ab <- t(spp_ab)  #Row = sites

colnames(spp_ab) == rownames(traits.sp)

#### Making abbreviated species names  
names1 <- make.cepnames(rownames(traits.sp))
colnames(spp_ab)    <- names1
rownames(traits.sp) <- names1
all(rownames (traits.sp)==colnames(spp_ab))

#Presence/absence
spp_ab[spp_ab>0]<-1

#Environmental variables
env <- read.table ('ENV.csv', h=T)

#@@@@@@@@@@@@@@@@@@@@#
#### RLQ Analysis ####
#@@@@@@@@@@@@@@@@@@@@#

#### As dataframe
spp_ab<-as.data.frame(spp_ab)
traits.sp<-as.data.frame(scale(traits.sp))
env<-as.data.frame(scale(env))

#RLQ matrices
coa1       <- dudi.coa(spp_ab, scannf = FALSE)
pca.traits <- dudi.pca(traits.sp, row.w = coa1$cw, scannf = FALSE)
pca.env    <- dudi.pca(env, row.w = coa1$lw, scannf = FALSE)

rlq1 <- rlq(pca.env, coa1, pca.traits, scannf = FALSE)
summary(rlq1)
plot(rlq1)

scatter (rlq1, posieig = "bottomright")
scatter (coa1)
scatter (pca.traits)
scatter (pca.env)

##Individual Plots
#s.arrow(rlq1$l1)
#s.arrow(rlq1$c1)
#s.label(rlq1$lQ, boxes = FALSE)


# Percentage of co-Inertia for each axis
100*rlq1$eig/sum(rlq1$eig) #71.9% fist and 22.5% second axes 


## weighted correlations axes / env.
t(pca.env$tab)%*%(diag(pca.env$lw))%*%as.matrix(rlq1$mR)


## weighted correlations axes / traits.
t(pca.traits$tab)%*%(diag(pca.traits$lw))%*%as.matrix(rlq1$mQ)

#Permutation test
rnda<-randtest(rlq1, modeltype = 6, nrepet = 999)
rnda
plot (rnda)

####Fourthcorner Permutarion test type=6, reducing poterian type error I, plus p-value correction
Srlq <- fourthcorner2(env, spp_ab, traits.sp, modeltype = 6, p.adjust.method.G = "bonferroni", nrepet = 999)

#results
Srlq$trRLQ


####Clustering
hc2 <- hclust(dist(rlq1$lQ, method = "euclidean"), method = "ward.D2")
plot(hc2, cex=.8)


###Which number of groups?
source('calinski.R')

### NUMBER of groups between 0 and 12.
ntest <- 12
res <- rep(0,ntest - 1)
for (i in 2:ntest){
  fac <- cutree(hc2, k = i)
  res[i-1] <- calinski(tab=rlq1$lQ, fac = fac)[1]
}
#par(mfrow=c(1,2))
plot(2:ntest, res, type='b', pch=20, xlab="Number of groups", ylab = "C-H index")
plot(3:ntest, diff(res), type='b', pch=20, xlab="Number of groups", ylab = "Diff in C-H index")


#### 4 GROUPS SEEMS A BETTER OPTIONS ####
hc2 <- hclust(dist(rlq1$lQ, method = "euclidean"), method =  "ward.D2")
plot(hc2, cex=1)
groups <- cutree(hc2, k=4) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(hc2, k=4, border="red")

spe.group2 <- as.factor(cutree(hc2, k = 4))
levels(spe.group2) <- c("C","A","B","D")
spe.group2 <- factor(spe.group2, levels=c("C","A","B","D"))

s.class(rlq1$lQ, spe.group2, col= 1:nlevels(spe.group2))
#s.arrow(rlq1$c1, add.plot = T,clab=1, boxes=T)
s.arrow(rlq1$l1*3, add.plot = T,clab=.6, boxes=F)


### PLOTING DIFFERENCES BETWEEN GROUPS FOR EACH TRAIT
source("corratio.R")
eta2 <- cor.ratio(traits.sp, data.frame(spe.group2), weights = rep(1, length(spe.group2)))
#par(mfrow=n2mfrow(ncol(traits.sp)))
#plot(table(spe.group2,traits.sp[,1]), main =names(traits.sp)[1])
for(i in 2:ncol(traits.sp)){
  label <- paste(names(traits.sp)[i], "(cor.ratio =", round(eta2[i-1],3), ")")
  plot(traits.sp[,i]~spe.group2, main = label, border = 1:nlevels(spe.group2))
}


####PLOTING CLUSTERS AGAINST PREDICTORS
colnames(spp_ab)==rownames(rlq1$lQ)
colnames(spp_ab)

A<-t(spp_ab[,spe.group2=='A'])
B<-t(spp_ab[,spe.group2=='B'])
C<-t(spp_ab[,spe.group2=='C'])
D<-t(spp_ab[,spe.group2=='D'])

as.matrix(rownames(A))
as.matrix(rownames(B))
as.matrix(rownames(C))
as.matrix(rownames(D))

Afreq<-rowSums(spp_ab[,spe.group2=='A'])/rowSums(spp_ab)
Bfreq<-rowSums(spp_ab[,spe.group2=='B'])/rowSums(spp_ab)
Cfreq<-rowSums(spp_ab[,spe.group2=='C'])/rowSums(spp_ab)
Dfreq<-rowSums(spp_ab[,spe.group2=='D'])/rowSums(spp_ab)
s.class(rlq1$lQ, spe.group2, col= 1:nlevels(spe.group2))
