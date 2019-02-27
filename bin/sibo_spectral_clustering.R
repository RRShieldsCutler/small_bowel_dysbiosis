library(kernlab)
library(ggplot2)

#############################################################
##### Comparing the CLOUD method to spectral clustering #####

# Starts with the aitchison distance matrix with all samples from the SIBO study
# Uses metadata file 

spec.clus.sibo.all <- sibo.aitch

sc <- specc(spec.clus.sibo.all, centers=2)
plot(spec.clus.sibo.all, col=sc, pch=20)            # estimated classes (x)
points(spec.clus.sibo.all, col=sibo.map$healthy_binary, pch="O") # true classes (<>)
legend(x=5, legend =c("H", "S"),
       col=sibo.map$healthy_binary)

#sc@.Data
spec.map <- sibo.map.cloud
spec.map$spec.clust <- as.character(sc@.Data)
table(spec.map$spec.clust, spec.map$dys_class)


# With just symptomatics
spec.clus.sibo.sx <- sibo.aitch[sibo.alpha.cloud$Row.names, sibo.alpha.cloud$Row.names]
spec.map.sx <- sibo.map.cloud[sibo.map.cloud$Row.names %in% sibo.alpha.cloud$Row.names, ]

sc.sx <- specc(spec.clus.sibo.sx, centers=2)
# plot(spec.clus.sibo.all, col=sc, pch=20)            # estimated classes (x)
# points(spec.clus.sibo.all, col=sibo.map$healthy_binary, pch="O") # true classes (<>)
# legend(x=5, legend =c("H", "S"),
#        col=sibo.map$healthy_binary)

#sc@.Data
spec.map.sx$spec.clust <- as.character(sc.sx@.Data)
# Compare the two methods classifications
table(spec.map.sx$spec.clust, spec.map.sx$dys_class)


# Show the SC distribution on the PCA plot
sc.pca <- as.data.frame(cmdscale(spec.clus.sibo.all, k=5))

pcnames = c()
for(i in 1:ncol(sc.pca)){
  pcnames[i] <- paste0("PC",i)
}
colnames(sc.pca) = pcnames

sc.plot.map <- merge(spec.map, sc.pca, by.x = 1, by.y = 0)

sc.plot <- ggplot(sc.plot.map, aes(x=PC1, y=PC2)) +
  geom_point(aes(shape=spec.clust, color=dys_class), size=4, alpha=0.75) +
  theme_classic() + theme(axis.text = element_text(color='black'),
                          aspect.ratio = 1)
sc.plot
