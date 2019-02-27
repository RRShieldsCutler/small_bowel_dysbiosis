library(tibble)
library(vegan)

# First, standard filtering of the OTU table for singletons/rare OTUs/low depth:

sibo.otus = read.delim(file = 'sibo_ampotu.txt', header = 1, row.names = 1, sep = '\t', check.names = F)
sum(sibo.otus)
# Total number of reads is ~6.8 million, so our threshold for low counts/singleton will be 7 counts.
dim(sibo.otus)
colnames(sibo.otus) <- lapply(colnames(sibo.otus), FUN = function(xx) paste0('S', xx))
rownames(sibo.otus) <- lapply(rownames(sibo.otus), FUN = function(xx) gsub(x = xx, pattern = ' ', replacement = '_'))
sibo.otus <- sibo.otus[(rowSums(sibo.otus) >= 7), ]  # Remove singletons at 7 based on counts threshold
sibo.otus = sibo.otus[rowMeans(sibo.otus > 0) >= 0.01, ]  # Remove rare sibo.otus in less than 1% of samples
# Check for and remove low depth samples
depths = colSums(sibo.otus)
sort(depths)[1:50]
max(depths)
dim(sibo.otus)
sibo.otus = sibo.otus[, depths >= 1000]
dim(sibo.otus)

# Now import the OTU table from diet intervention study
west.otus = read.delim(file = 'westd_ampotu.txt',
                       header = 1, row.names = 1, sep = '\t', check.names = F)

# Total number of reads is 3.1 million, so our threshold for low counts will be 3 counts.
dim(west.otus)
rownames(west.otus) <- lapply(rownames(west.otus), FUN = function(xx) gsub(x = xx, pattern = ' ', replacement = '_'))
west.otus <- west.otus[(rowSums(west.otus) >= 3), ]  # Remove singletons based on the total read depth

# Get the aspirate samples only
west.otus.a <- west.otus[, grep('aspirate', colnames(west.otus))]
west.otus.a = west.otus.a[rowMeans(west.otus.a > 0) >= 0.01, ]
west.otus.a <- west.otus.a[(rowSums(west.otus.a) >= 3), ]
# Check for and remove low depth samples
depths = colSums(west.otus.a)
sort(depths)[1:20]
max(depths)
dim(west.otus.a)
west.otus.a = west.otus.a[, depths >= 1000]
dim(west.otus.a)

## Get the metadata, yo.
sibo.map <- read.csv('../../data/meta.data.with.DI.score20180214.csv', header=1, check.names = F, row.names = 1)
# sibo.map <- tibble::column_to_rownames(sibo.map, var = '')

# Organize tables
sample.ids.sibo <- intersect(colnames(sibo.otus), rownames(sibo.map))
sample.ids.sibo <- sort(sample.ids.sibo)
sibo.map <- sibo.map[sample.ids.sibo, ]
sibo.otus <- sibo.otus[, sample.ids.sibo]

# Get the diet study map
map <- read.delim('/Users/Robin/Box Sync/knights_box/mayo_small_bowel/amplicon_project/westd/data/westd_map2.txt',
                  sep = '\t', header=1, row.names = 1, check.names = F)
# map <- tibble::column_to_rownames(map, var = 'SampleID')
sample.ids.west <- intersect(rownames(map), colnames(west.otus.a))
sample.ids.west <- sort(sample.ids.west)
map <- map[sample.ids.west, ]
west.otus.a <- west.otus.a[, sample.ids.west]

# Generate some useful new variables
map$pre_post_fiber[map$Vist_number==1] <- 'high_fiber'
map$pre_post_fiber[map$Vist_number==2] <- 'low_fiber'
map$pre_post_fiber <- factor(map$pre_post_fiber,
                             levels=c('high_fiber','low_fiber'),
                             ordered=T)

map$micro_binary <- map$Aspirate_micro_result
map$micro_binary[map$Aspirate_Micro_result==0] <- 'positive'
map$micro_binary[map$Aspirate_Micro_result==2] <- 'positive'
map$micro_binary[map$Aspirate_Micro_result==3] <- 'negative'
map$micro_binary <- as.factor(map$micro_binary)

map$healthy_binary <- 'healthy'

map$healthy_pos_neg <- map$Aspirate_micro_result
map$healthy_pos_neg[map$Aspirate_Micro_result==0] <- 'healthy_positive'
map$healthy_pos_neg[map$Aspirate_Micro_result==2] <- 'healthy_positive'
map$healthy_pos_neg[map$Aspirate_Micro_result==3] <- 'healthy_negative'
map$healthy_pos_neg <- as.factor(map$healthy_pos_neg)

west.map <- map

west.otus.a[is.na(west.otus.a)] <- 0

# Multiplicative replacement and CLR
west.otus.a.c <- t(west.otus.a); eps <- 0.5
west.otus.a.c <- west.otus.a.c * (1 - rowSums(west.otus.a.c==0) * eps / rowSums(west.otus.a.c))
west.otus.a.c[west.otus.a.c == 0] <- eps
west.otus.a.c <- sweep(west.otus.a.c, 1, rowSums(west.otus.a.c), '/')
ls <- log(west.otus.a.c)
west.otus.a.c <- t(ls - rowMeans(ls))
west.otus.a.c <- west.otus.a.c[, !is.nan(colSums(west.otus.a.c))]
west.otus.a.clr <- as.data.frame(west.otus.a.c)

west.aitch <- as.data.frame(as.matrix(vegdist(t(west.otus.a.clr), method = 'euclidean')))


#### Generating the aitchison distance matrix ####

sibo.otus[is.na(sibo.otus)] <- 0

# Multiplicative replacement and CLR
sibo.otus.c <- t(sibo.otus); eps <- 0.5
sibo.otus.c <- sibo.otus.c * (1 - rowSums(sibo.otus.c==0) * eps / rowSums(sibo.otus.c))
sibo.otus.c[sibo.otus.c == 0] <- eps
sibo.otus.c <- sweep(sibo.otus.c, 1, rowSums(sibo.otus.c), '/')
ls <- log(sibo.otus.c)
sibo.otus.c <- t(ls - rowMeans(ls))
sibo.otus.c <- sibo.otus.c[, !is.nan(colSums(sibo.otus.c))]
sibo.otus.clr <- as.data.frame(sibo.otus.c)

sibo.aitch <- as.matrix(vegdist(t(sibo.otus.clr), method = 'euclidean'))

# Save the table as a file
# sibo.aitch.w <- tibble::rownames_to_column(data.frame(sibo.aitch), var = 'sampleID')
# write.table(sibo.aitch.w, file = '../sibo-westd/sibo_aitchison_table.txt', sep = '\t', quote = F, col.names = T, row.names = T)

# Combine the two studies' metadata and OTU tables
west.map$study <- 'westd_study'
sibo.map$healthy_pos_neg <- sibo.map$healthy_or_sb_positive_or_sb_negative
sibo.map$study <- 'sibo_study'
all.map <- rbind(west.map[, c('healthy_pos_neg','healthy_binary','micro_binary','study')], sibo.map[, c('healthy_pos_neg','healthy_binary','micro_binary','study')])
all.otus <- merge(sibo.otus, west.otus.a, by = 0, all = T)
all.otus <- tibble::column_to_rownames(all.otus, var = 'Row.names')
all.otus[is.na(all.otus)] <- 0

# Organize names
all.ids <- intersect(colnames(all.otus), rownames(all.map))
all.ids <- sort(all.ids)
all.otus <- all.otus[,all.ids]
all.map <- all.map[all.ids,]

# Multiplicative replacement and CLR
all.otus.c <- t(all.otus); eps <- 0.5
all.otus.c <- all.otus.c * (1 - rowSums(all.otus.c==0) * eps / rowSums(all.otus.c))
all.otus.c[all.otus.c == 0] <- eps
all.otus.c <- sweep(all.otus.c, 1, rowSums(all.otus.c), '/')
ls <- log(all.otus.c)
all.otus.c <- t(ls - rowMeans(ls))
all.otus.c <- all.otus.c[, !is.nan(colSums(all.otus.c))]
all.otus.clr <- as.data.frame(all.otus.c)

all.aitch <- as.data.frame(as.matrix(vegdist(t(all.otus.clr), method = 'euclidean')))