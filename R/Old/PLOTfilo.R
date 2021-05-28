func<-read.table("funcional.txt", h=T)
head(func)

library(vegan)
func2<-func
rownames(func2)<-make.cepnames(rownames(func))

library(cluster)
####Clustering
distGOWER<-daisy(func2, metric = "gower")
hc2 <- hclust(distGOWER, method = "average")

#install.packages(c("factoextra", "dendextend"))
library(factoextra)
par(mfrow=c(1,1), mar=c(3.2, 3.2, .1, 0.1), cex=1, las=0, tcl=-0.3)
x11(height = 10, width = 5)
asd<-fviz_dend(hc2, k = 36, cex = 0.3, horiz = TRUE, main ="", rect = T, lower_rect=-0.035, lwd=0.2, phylo_layout=T)

asd + labs(subtitle = "hjust = 1 and vjust = 0 place tottom\nright cornor of the text to (x, y)")


dev.copy2pdf(device = x11, file="PHYLO.pdf")

text(6, 2, "the text is CENTERED around (x,y) = (6,2) by default",
     cex = .8)

?fviz_dend
