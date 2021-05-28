####Libraries####
library(here) ##todas
library(iNEXT)
library(ape)
library(plotTaxa2.R) ###output

##Gráfico com o resultado das regressoes
png(here("output", "figures", "LM.png"))
par(mfrow=c(1,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
vecpoints<-ifelse(RIQmod[,2]<0.05,19,1)
plot(RIQmod[,1]~c(1985:2013), xlab = "", ylab="Effect size (richness)", xaxt="n", pch=vecpoints, mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.2)
summary(lm(RIQmod[,1]~c(1985:2013)))
axis(side = 1, at = c(1985:2013), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
mtext("Years", las=0, cex=1.2, 1, 2.8, outer=F)
dev.off()


##Efeito do land use em 1993 e 2013
degrass1993<-100-dat1$X1993
degrass2013<-100-dat1$X2013
png(here("output", "supp", "LM_year.png"))
par(mfrow=c(1,1), mar=c(3.2, 3.2, .1, 0.1), cex=1, las=0, tcl=-0.3)
plot(riqTOT~degrass1993, xlab='Native vegetation conversion (%)',ylab = 'Species richness', mgp=c(2, 0.5, 0),tcl=-0.3, lwd=1, bg="orange", pch=21, cex=1.5, cex.lab=1.2)
points(riqTOT~degrass2013, xlab='Native vegetation conversion (%)',ylab = 'Species richness',  mgp=c(2, 0.5, 0),tcl=-0.3, lwd=1, bg="blue", pch=21, cex=1.5, cex.lab=1.2)
abline (lm(riqTOT~degrass1993),lwd=3, col="orange")
abline (lm(riqTOT~degrass2013),lwd=3, col="blue", lty = 2)
legend("topright", legend=c("1993", "2013"),col=c("black", "black"), cex=1, pt.bg=c("orange","blue"), pch=c(21,21),pt.cex=1.2, bty="o")
dev.off()


##Grafico RDA
png(here("output", "figures", "RDA.png"))
par(mfrow=c(1,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
vecpoints1<-ifelse(COMPmod[,2]<0.05,19,1)
plot(COMPmod[,1]~c(1985:2013), xlab = "", ylab="Effect size (composition)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.2, pch=vecpoints1)
axis(side = 1, at = c(1985:2013), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
mtext("Years", las=0, cex=1.2, 1, 2.8, outer=F)
summary(lm(COMPmod[,1]~c(1985:2013)))
dev.off()


####### PLOT COMPOSTO RIQ E COMP ########
png(here("output", "figures", "Plots.png"))
par(mfrow=c(2,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
vecpoints<-ifelse(RIQmod[,2]<0.05,19,1)
plot(RIQmod[,1], xlab = "", ylab="Effect size (richness)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.4, ylim=c(0.15,0.56), pch=vecpoints)
axis(side = 1, at = seq(1,29,1), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
legend("topright", legend=c("A"),cex=1.5, bty="n", text.font=2)
vecpoints1<-ifelse(COMPmod[,2]<0.05,19,1)
plot(COMPmod[,1], xlab = "", ylab="Effect size (composition)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.4, ylim=c(0.05,0.072), pch=vecpoints1)
axis(side = 1, at = seq(1,29,1), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
mtext("Years", las=0, cex=1.4, 1, 3, outer=F)
legend("topright", legend=c("B"),cex=1.5, bty="n",text.font=2)
dev.off()

##GRAFICOS TAX TITAN

png (png(here("output", "figures", "titan_tax.png")))
par(mfrow=c(3,2), mar=c(2.6, 4, 0.5, 3.5), cex=1, las=0, tcl=-0.3)
plotTaxa2(titanNAT93, xlabel = "NVL 1993 (%)",cex.axis = 0.9, cex.taxa = 0.9, xmin=0, xmax=80)
legend("topright", legend=c("A"),cex=1.5, bty="n",text.font=2)
plotTaxa2(titanNAT13, xlabel = "NVL 2013 (%)",cex.axis = 0.9, cex.taxa = 0.9, xmin=0, xmax=80)
titanNAT13$sumz.cp
legend("topright", legend=c("B"),cex=1.5, bty="n",text.font=2)
plotTaxa2(titanMAG, z2=F, xlabel = "Magnitude (accumulated%)",cex.axis = 0.9, cex.taxa = 0.9, xmin=0, xmax=80)
legend("topright", legend=c("C"),cex=1.5, bty="n",text.font=2)
plotTaxa2(titanT20, z1=F, xlabel = "Duration T20 (years)",cex.axis = 0.9, cex.taxa = 0.9)
legend("topright", legend=c("D"),cex=1.5, bty="n",text.font=2)
plotTaxa2(titanT40, z1=F, xlabel = "Duration T40 (years)",cex.axis = 0.9, cex.taxa = 0.9)
titanT40$sumz.cp
legend("topright", legend=c("E"),cex=1.5, bty="n",text.font=2)
dev.off()

##GRAFICOS FUN TITAN

png (png(here("output", "figures", "titan_fun.png")))
par(mfrow=c(2,2), mar=c(2.6, 4, 0.5, 3.5), cex=1, las=0, tcl=-0.3)
plotTaxa2(FUN.titanNAT93, z2=F, xlabel = "NVL 1993 (%)",cex.axis = 0.9, cex.taxa = 0.9, xmin=0, xmax=80)
legend("topright", legend=c("A"),cex=1.5, bty="n",text.font=2)
plotTaxa2(FUN.titanNAT13, xlabel = "NVL 2013 (%)",cex.axis = 0.9, cex.taxa = 0.9, z2=F, xmin=0, xmax=80)
legend("topright", legend=c("B"),cex=1.5, bty="n",text.font=2)
plotTaxa2(FUN.titanMAG, z2=F, xlabel = "Magnitude (accumulated%)",cex.axis = 0.9, cex.taxa = 0.9, xmin=0, xmax=80)
legend("topright", legend=c("C"),cex=1.5, bty="n",text.font=2)
plotTaxa2(FUN.titanT20, z2=F,xlabel = "Duration T20 (years)",cex.axis = 0.9, cex.taxa = 0.9)
legend("topright", legend=c("D"),cex=1.5, bty="n",text.font=2)
dev.off()

