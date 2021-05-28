library(plyr)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(ClustOfVar)
library(ggdendro)

# Cluster analysis
fit = hc2

# Labels data frames
df2 <- data.frame(cluster = cutree(fit, k =18), 
                  states = factor(fit$labels, levels = fit$labels[fit$order]))

spe.group2 <- as.factor(cutree(hc2, k = 18))
FGs<-rep("FG",length(spe.group2))
FGs<-paste0(FGs, spe.group2)

df2$cluster<-FGs

df3 <- ddply(df2, .(cluster), summarise, pos = mean(as.numeric(states)))


library(cluster); library(ggdendro); library(ggplot2)

hc2
distGOWER
funcMORFO


hcaf   <- hc2
k     <- 18
clustf <- cutree(hcaf,k=k)  # k clusters

dendrf    <- dendro_data(hcaf, type="rectangle") # convert for ggplot
clust.dff <- data.frame(label=rownames(funcMORFO), 
                        cluster=factor(clustf)) 
dendrf[["labels"]]   <- merge(dendrf[["labels"]],clust.dff, by="label")
rectf <- aggregate(x~cluster,label(dendrf),range)
rectf <- data.frame(rectf$cluster,rectf$x)
ymax <- mean(hcaf$height[length(hcaf$height)-((k-2):(k-1))])

png("FIGS/TESTE.png", width = 6, height =10, units = 'in', res = 600)
ggplot() + 
  geom_segment(data=segment(dendrf), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=label(dendrf), aes(x, y, label= label, hjust=1.1,
                                    color=cluster),
            size=1.5) +
  geom_rect(data=rectf, aes(xmin=X1-0.5, xmax=X2+.5, ymin=0, ymax=-0.25),
            color="red", fill=NA, lwd=0.1)+
  geom_text(data = df3, aes(x = pos, y=-0.2,label = cluster),size=2.5)+
  coord_flip() +
  scale_y_continuous(name="Gower distance",expand = c(0, 0)) + 
  scale_x_discrete(expand = c(0.01, 0.01)) +
  scale_color_discrete(name="Cluster") +
theme(legend.position="none",
      axis.line.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.y=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x = element_text(),  # show x-axis labels
      axis.ticks.x = element_line(), # show x-axis tick marks
      axis.line.x = element_line(),
      panel.background=element_rect(fill="white"),
      panel.grid=element_blank())
dev.off()

