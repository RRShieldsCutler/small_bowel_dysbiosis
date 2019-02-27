require(pROC)
require(randomForest)
require(ggplot2)

load(file = 'Data.wk.RData')
#############################################################################
# Random Forest-based prediction to form sIndex
set.seed(123)
prop <- data.obj$abund.list[['Species']]
prop <- t(prop) / colSums(prop)
response <- data.obj$meta.dat$healthy_binary
rf.obj <- randomForest(x=prop, y=response, sampsize=table(response))


# OOB probability of being symptomatic (sIndex)
sIndex <- predict(rf.obj,  type='prob')[, 'sibo_patient']
data <- data.frame(sIndex = sIndex, Group = response)

pdf('Boxplot_OOB_prediction_sIndex.pdf', height = 4, width = 4)
ggplot(data, aes(x = Group, y = sIndex , col = Group)) +
		geom_boxplot() +
		geom_jitter(width = 0.2, alpha = 0.75, size = 2) +
		theme_bw() + 
		theme(legend.position="none") +
		scale_colour_manual(values=c("#CC6666", "#9999CC"))
dev.off()

pdf('ROC_plot_OOB_prediction_sIndex .pdf', height = 5, width = 5)
roc2 <- roc(response, sIndex , 
		plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
		print.auc=TRUE, show.thres=TRUE)

ci(roc2)  # 95% CI: 0.8435-0.9489 (DeLong)

sens.ci <- ci.se(roc2, specificities=seq(0, 1, len=20))
plot(sens.ci, type="shape", col="lightblue")
plot(sens.ci, type="bars")
dev.off()

#############################################################################
# Association of SIndex with host factors (linear regression)
sIndex2 <- sIndex 
sIndex2[sIndex2 == 1] <- 0.995
y <- binomial()$linkfun(sIndex2)   # Logit scale
data <- data.obj$meta.dat
ind <- data$healthy_binary == 'sibo_patient'

data <- data[ind, ]
y <- y[ind]
sIndex <- sIndex[ind]

# Marignal correlation with host factors (logit scale)
data$Sex <- data$sex_0M_1F
data$SB_status <- factor(data$sb_status)
data$GI_Surgery <- as.character(data$GI_Surgery)
data$GI_Surgery[data$GI_Surgery != '0'] <- '1'
data$GI_Surgery <- factor(data$GI_Surgery)
levels(data$GI_Surgery) <- c('No', 'Yes')
data$Antibiotics <- as.character(data$Antibiotic_rx)
data$Antibiotics[data$Antibiotics %in% c('', '5')] <- 'No'
data$Antibiotics[data$Antibiotics != 'No'] <- 'Yes'
data$Antibiotics <- factor(data$Antibiotics)

data$PPI <- as.character(data$PPI)
data$PPI[data$PPI == '5'] <- 'No'
data$PPI[data$PPI != 'No' & !is.na(data$PPI)] <- 'Yes'
data$PPI <- factor(data$PPI)

summary(lm(y ~ Age, data = data))
summary(lm(y ~ Sex, data = data))
summary(lm(y ~ BMI, data = data))
summary(lm(y ~ SB_status, data = data))
summary(lm(y ~ GI_Surgery, data = data))
summary(lm(y ~ Antibiotics, data = data))
summary(lm(y ~ PPI, data = data))

# Joint association with host factors (logit scale)
summary(lm(y ~ Age + GI_Surgery + Antibiotics + PPI, data = data))

#############################################################################
# Boxplots - both logit scale and original probability scale
data$sIndex <- y
data$AgeGroup <- as.character(data$Age)
data$AgeGroup[data$Age >= 50] <- 'Age>=50'
data$AgeGroup[data$Age < 50] <- 'Age<50'
obj1 <- ggplot(data, aes(x = AgeGroup, y = sIndex, col = AgeGroup)) +
		geom_boxplot() +
		geom_jitter(width = 0.2, alpha = 0.75, size = 2) +
		theme_bw() + 
		theme(legend.position="none") +
		scale_colour_manual(values=c("#CC6666", "#9999CC"))

obj2 <- ggplot(data, aes(x = Antibiotics, y = sIndex, col = Antibiotics)) +
		geom_boxplot() +
		geom_jitter(width = 0.2, alpha = 0.75, size = 2) +
		theme_bw() + theme(legend.position="none") +
		scale_colour_manual(values=c("#CC6666", "#9999CC"))

obj3 <- ggplot(data, aes(x = GI_Surgery, y = sIndex, col = GI_Surgery)) +
		geom_boxplot() +
		geom_jitter(width = 0.2, alpha = 0.75, size = 2) +
		theme_bw() + theme(legend.position="none") +
		scale_colour_manual(values=c("#CC6666", "#9999CC"))

obj4 <- ggplot(subset(data, !is.na(PPI)), aes(x = PPI, y = sIndex, col = PPI)) +
		geom_boxplot() +
		geom_jitter(width = 0.2, alpha = 0.75, size = 2) +
		theme_bw() + theme(legend.position="none") +
		scale_colour_manual(values=c("#CC6666", "#9999CC"))

pdf('HostFactorSIndex_logit.pdf', height = 6, width = 6)
multiplot(obj1, obj2, obj3, obj4, cols=2)
dev.off()

data$sIndex <- sIndex
data$AgeGroup <- as.character(data$Age)
data$AgeGroup[data$Age >= 50] <- 'Age>=50'
data$AgeGroup[data$Age < 50] <- 'Age<50'
obj1 <- ggplot(data, aes(x = AgeGroup, y = sIndex, col = AgeGroup)) +
		geom_boxplot() +
		geom_jitter(width = 0.2, alpha = 0.75, size = 2) +
		theme_bw() + 
		theme(legend.position="none") +
		scale_colour_manual(values=c("#CC6666", "#9999CC"))


obj2 <- ggplot(data, aes(x = Antibiotics, y = sIndex, col = Antibiotics)) +
		geom_boxplot() +
		geom_jitter(width = 0.2, alpha = 0.75, size = 2) +
		theme_bw() + theme(legend.position="none") +
		scale_colour_manual(values=c("#CC6666", "#9999CC"))

obj3 <- ggplot(data, aes(x = GI_Surgery, y = sIndex, col = GI_Surgery)) +
		geom_boxplot() +
		geom_jitter(width = 0.2, alpha = 0.75, size = 2) +
		theme_bw() + theme(legend.position="none") +
		scale_colour_manual(values=c("#CC6666", "#9999CC"))

obj4 <- ggplot(subset(data, !is.na(PPI)), aes(x = PPI, y = sIndex, col = PPI)) +
		geom_boxplot() +
		geom_jitter(width = 0.2, alpha = 0.75, size = 2) +
		theme_bw() + theme(legend.position="none") +
		scale_colour_manual(values=c("#CC6666", "#9999CC"))

pdf('HostFactorSIndex_probability.pdf', height = 6, width = 6)
multiplot(obj1, obj2, obj3, obj4, cols=2)
dev.off()
#############################################################################

