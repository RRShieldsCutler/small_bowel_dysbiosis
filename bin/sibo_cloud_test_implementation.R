library(tibble)
library(ggplot2)
library(vegan)
library(gridExtra)
library(randomForest)
library(gplots)
library(heatmap3)
library(stringr)
library(Boruta)

#################################################################################
#### Apply the CLOUD test to classify dysbiosis relative to healthy controls ####
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6080375/
# Code adapted from:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6080375/bin/40168_2018_514_MOESM1_ESM.docx

# Uses metadata, and Aitchison distances based on OTU tables as calculated in accompanying scripts

#### CLOUD test for symptomatic patients ####

source(file = '~/bin/cloud_test.R')
# Start with the healthy distribution to figure out the k
k.neighbors <- 37  # We want the entire healthy group to represent a healthy neighborhood

healthy.ids <- rownames(sibo.map)[sibo.map$healthy_pos_neg == 'healthy']
healthy.cloud <- NULL
test.df <- sibo.aitch[rownames(sibo.aitch) %in% healthy.ids,
                      colnames(sibo.aitch) %in% healthy.ids]
healthy.cloud <- piecewise_kn_V1(d=test.df, test.ix=c(1:nrow(test.df)), k=k.neighbors)
sum(healthy.cloud$pvals < 0.05)

sx.ids <- rownames(sibo.map)[sibo.map$healthy_binary != 'healthy']
sibo.cloud <- NULL
for (p in sx.ids) {
  test.df <- sibo.aitch[rownames(sibo.aitch) %in% p | rownames(sibo.aitch) %in% healthy.ids,
                        colnames(sibo.aitch) %in% p | colnames(sibo.aitch) %in% healthy.ids]
  cloud.result <- piecewise_kn_V1(d=test.df, test.ix=c(1:nrow(test.df)), k=k.neighbors)
  cloud.result <- cloud.result[p,]
  sibo.cloud <- rbind(sibo.cloud, cloud.result)
}
table(sibo.cloud$pvals < 0.05)

# Repeat for western diet
# healthy.ids <- rownames(all.map)[all.map$healthy_binary == 'healthy' & all.map$study == 'sibo']
# healthy.cloud <- NULL
# test.df <- all.aitch[rownames(all.aitch) %in% healthy.ids,
#                       colnames(all.aitch) %in% healthy.ids]
# healthy.cloud <- piecewise_kn_V1(d=test.df, test.ix=c(1:nrow(test.df)), k=k.neighbors)
# sum(healthy.cloud$pvals < 0.05)

west.ids <- rownames(all.map)[all.map$study == 'westd_study']
west.cloud <- NULL
for (p in west.ids) {
  test.df <- all.aitch[rownames(all.aitch) %in% p | rownames(all.aitch) %in% healthy.ids,
                        colnames(all.aitch) %in% p | colnames(all.aitch) %in% healthy.ids]
  cloud.result <- piecewise_kn_V1(d=test.df, test.ix=c(1:nrow(test.df)), k=k.neighbors)
  cloud.result <- cloud.result[p,]
  west.cloud <- rbind(west.cloud, cloud.result)
}
# table(west.cloud$pvals < 0.05)

# Cutoff for classification, at 2xSD from the log2stat
shapiro.test(log2(healthy.cloud$stats))
healthy.cloud$log2stats <- log2(healthy.cloud$stats)
sibo.cloud$log2stats <- log2(sibo.cloud$stats)
west.cloud$log2stats <- log2(west.cloud$stats)
cutoff.val <- (mean(healthy.cloud$log2stats) + 2*sd(healthy.cloud$log2stats))
table(healthy.cloud$log2stats > cutoff.val)
table(sibo.cloud$log2stats > cutoff.val)
table(west.cloud$log2stats > cutoff.val)


# Doing some analysis with these
west.map.cloud <- merge(west.cloud, west.map, by=0)
plot(log2(west.map.cloud$stats) ~ west.map.cloud$Vist_number)

# Write the westd map with CLOUD scores
west.map.cloud$dys_class <- west.map.cloud$log2stats > cutoff.val
west.map.cloud$dys_class[west.map.cloud$dys_class == T] <- 'dysbiotic'
west.map.cloud$dys_class[west.map.cloud$dys_class == F] <- 'healthy-like'

west.map.cloud.w <- west.map.cloud
west.map.cloud.w$stats <- round(west.map.cloud.w$stats, digits = 5)
west.map.cloud.w$pvals <- round(west.map.cloud.w$pvals, digits = 5)
west.map.cloud.w$log2stats <- round(west.map.cloud.w$log2stats, digits = 5)
colnames(west.map.cloud.w)[1] <- 'sampleID'
write.table(west.map.cloud.w, file = 'metadata_cloud_westd_180406.txt', sep = '\t', quote = F, row.names = F)


sibo.map.cloud <- merge(sibo.map, rbind(sibo.cloud, healthy.cloud), by=0)
plot(log2(sibo.map.cloud$stats) ~ sibo.map.cloud$DI)
boxplot(log2(sibo.map.cloud$stats) ~ sibo.map.cloud$Abx_4wks_before_EGD)


### These groups were also classified using a similar approach to the symptomatic index,
### but it involved using the previous classification as part of the input data, so the 
### CLOUD method avoids the potential bias by using a fresh approach. However, we wanted
### to compare whether these two approaches gave similar classifications. It was saved
### as the DI in the imported metadata table. But the log2stats and dys_class are the 
### scores and classifications used in the final analysis.

# RF index and log2stat agree well
cor.test(sibo.map.cloud$DI, sibo.map.cloud$log2stats)
plot(sibo.map.cloud$DI, sibo.map.cloud$log2stats)

sibo.map.cloud$dys_class <- sibo.map.cloud$log2stats > cutoff.val
sibo.map.cloud$dys_class[sibo.map.cloud$dys_class == T] <- 'dysbiotic'
sibo.map.cloud$dys_class[sibo.map.cloud$dys_class == F] <- 'healthy-like'

# Check out the associations with Alpha diversity (calculated from same OTU table using QIIME)
table(sibo.map.cloud$healthy_or_sb_positive_or_sb_negative, sibo.map.cloud$dys_class)
sibo.alphadiv <- read.delim('rarefied5000_alphadivs.txt', sep = '\t', check.names = F, row.names = 1)
sibo.alphadiv$Row.names <- rownames(sibo.alphadiv)
sibo.alpha.cloud <- merge(sibo.map.cloud, sibo.alphadiv, by='Row.names', all = F)
wilcox.test(sibo.alpha.cloud$PD_whole_tree ~ sibo.alpha.cloud$dys_class)

# Write the CLOUD stat table
sibo.map.cloud.w <- sibo.map.cloud
sibo.map.cloud.w$stats <- round(sibo.map.cloud.w$stats, digits = 5)
sibo.map.cloud.w$pvals <- round(sibo.map.cloud.w$pvals, digits = 5)
sibo.map.cloud.w$log2stats <- round(sibo.map.cloud.w$log2stats, digits = 5)
colnames(sibo.map.cloud.w)[1] <- 'sampleID'
write.table(sibo.map.cloud.w, file = '../../data/metadata_cloud_SIBO_study.txt', quote = F, row.names = F, sep = '\t')

# Just sx
sibo.alpha.cloud.w <- rbind(sibo.map.cloud, west.map.cloud)
sibo.alpha.cloud.w$stats <- round(sibo.alpha.cloud.w$stats, digits = 5)
sibo.alpha.cloud.w$pvals <- round(sibo.alpha.cloud.w$pvals, digits = 5)
sibo.alpha.cloud.w$log2stats <- round(sibo.alpha.cloud.w$log2stats, digits = 5)
colnames(sibo.alpha.cloud.w)[1] <- 'sampleID'
write.table(sibo.alpha.cloud.w, file = '../../data/metadata_cloud_alpha_symptomatic.txt', quote = F, row.names = F, sep = '\t')

# Add these values to the plot
sibo.pca.2 <- tibble::column_to_rownames(sibo.pca, var = 'Row.names')
sibo.pca.cloud <- merge(sibo.pca.2, sibo.cloud, by=0)
healthy.pca.cloud <- merge(sibo.pca.2, healthy.cloud, by=0)

ggplot() +
  geom_point(aes(x=sibo.pca.cloud$PC1, y=sibo.pca.cloud$PC2,
                 color=sibo.pca.cloud$stats), size=2.8, alpha=0.75) +
  geom_point(aes(x=healthy.pca.cloud$PC1, y=healthy.pca.cloud$PC2,
                 color=healthy.pca.cloud$stats), size = 1.5, alpha = 0.5, shape=17) +
  scale_color_gradient(low="#FFCC99", high="#660000", guide = 'colourbar') +
  theme_classic() + theme(
    aspect.ratio = 1, axis.title.y = element_text(size=9),
    axis.text = element_text(color='black', size=8))
  
# Tukey boxplot for healthy vs sx CLOUD scores
score.vs.healthy <- ggplot(sibo.map.cloud, aes(x=healthy_binary, y=log2stats, color=healthy_binary)) +
  geom_boxplot(outlier.colour = 'white') + geom_jitter(width = .27, aes(color=healthy_binary), alpha=0.75) +
  theme_classic() + theme(axis.text = element_text(color='black')) +
  # scale_color_manual(values = c("healthy" = "blue", "sibo_patient"='red')) +
  scale_color_manual(values = c("healthy" = "red", "sibo_patient"='blue'),
                    labels=c("healthy" = "healthy", "sibo_patient"='symptomatic')) +
  scale_x_discrete(labels=c("healthy" = "healthy", "sibo_patient"='symptomatic')) +
  labs(x=NULL,y='Dysbiosis CLOUD score') +
  guides(color=F)
# ggsave(score.vs.healthy, file='../../results_dysbiosis_index/cloud_symp-vs-healthy_boxes_v2.png', height = 3, width = 3, dpi=300)



####### Run random forest model on the new classifications ######
####### and Boruta feature selection to yield significantly differentiating OTUs

# Read in the genus-level table, generated from the OTU table using QIIME2 
genus <- read.delim('sibo_prok_otu_L6.txt', header=1, check.names = F, row.names = 1, sep = '\t')

# CLR transform
genus.c <- t(genus); eps <- 0.5
genus.c <- genus.c * (1 - rowSums(genus.c==0) * eps / rowSums(genus.c))
genus.c[genus.c == 0] <- eps
genus.c <- sweep(genus.c, 1, rowSums(genus.c), '/')
ls <- log(genus.c)
genus.c <- t(ls - rowMeans(ls))
genus.c <- genus.c[, !is.nan(colSums(genus.c))]
genus.clr <- as.data.frame(genus.c)

sibo.rf.map <- tibble::column_to_rownames(sibo.map.cloud, var = 'Row.names')
rf.ids <- intersect(colnames(genus.clr), rownames(sibo.rf.map))
rd.ids <- sort(rf.ids)
genus.clr <- t(genus.clr[, rf.ids])
sibo.rf.map <- sibo.rf.map[rf.ids, ]
dim(genus)
dim(sibo.rf.map)

# Random forest model
rfc.model <-randomForest(x = genus.clr, y=factor(sibo.rf.map$dys_class), ntree = 1000, importance = T, keep.forest = T)
rfc.model$importance
rfc.prob <- predict(rfc.model, x=genus.clr, type = 'prob')
sibo.rf.map$cloud_prob <- rfc.prob[,1]
pairs(cbind(sibo.rf.map$cloud_prob, sibo.rf.map$log2stats, sibo.rf.map$DI))
plot(sibo.rf.map$cloud_prob, sibo.rf.map$log2stats)
plot(sibo.rf.map$cloud_prob, sibo.rf.map$DI)

# Boruta feature selection
boruta.rfimp <- Boruta(x=genus.clr, y=factor(sibo.rf.map$dys_class))
length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed'])
sig.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed']),
       FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
sig.taxa
length(sig.taxa)

# Plot the features as heatmap
sig.genera.clr <- genus.clr[, colnames(genus.clr) %in% sig.taxa]
colnames(sig.genera.clr) <- lapply(X=colnames(sig.genera.clr),
                                   FUN = function(xx) strsplit(as.character(xx),'__',fixed=T)[[1]][7])

fullmap <- read.delim('../../data/metadata_cloud_SIBO_study_180406.txt', sep = '\t', check.names = F, row.names = 1, header=1)
dys.ids <- as.character(rownames(fullmap[fullmap$dys_class=='dysbiotic',]))
healthylike.ids <- as.character(rownames(fullmap[fullmap$dys_class=='healthy-like', ]))

dys_healthy_order <- c(dys.ids, healthylike.ids)

sig.genera.clr <- sig.genera.clr[dys_healthy_order, ]
fullmap2 <- fullmap[dys_healthy_order, ]

hm.meta <- data.frame(fullmap2[,'dys_class'], row.names = rownames(fullmap2), col=fullmap2$dys_class=='healthy-like')
hm.meta$col[hm.meta$col==TRUE] <- '#E69F00'
hm.meta$col[hm.meta$col==FALSE] <- '#009E73'

png(paste0("../../results_dysbiosis_index/Boruta_heatmap_sibo_v4.png"),  # create PNG for the heat map
          width = 8.5*300,                        # 5 x 300 pixels
          height = 6*300,
          res = 300,                              # 300 pixels per inch
          pointsize = 6)                          # smaller font size
heatmap3(file='../../results_dysbiosis_index/test.pdf', x=t(sig.genera.clr), showColDendro = F, showRowDendro = F,
         ColSideColors = cbind(`Dysbiotic = purple`=hm.meta$col), cexRow = 2.8,
         margins=c(2,7), Colv=NA)
         # ColSideAnn = hm.meta2,
         # ColSideFun = function(x) showAnn(x),ColSideWidth=1.2,ColSideColors=, cexRow = 1)
dev.off()


#### The SIBO vs DI test ####
wilcox.test(sx.map.cloud$log2stats ~ sx.map.cloud$micro_binary)$p.value
fisher.test(table(sx.map.cloud$dys_class, sx.map.cloud$micro_binary))
