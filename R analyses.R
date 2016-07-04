# For any questions, I am available by e-mail: ilkkasipila@yahoo.com

# Required libraries
library(BayesLCA)
library(e1071)
library(coda)
library(pheatmap)
library(ggplot2)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(grid)
library(clinfun)
library(logisticPCA)
library(reshape2)
library(devtools)

# Setup: choose the working directory by replacing the text within the quotes. 
# This should be the same folder as where the Presence_absence.txt file is located. 
# By default, it is set to be in the same folder as the R script and workspace
setwd("C:/Users/Ilkka/OneDrive/RMA thesis/Excel/CSV files")

# Import the tab limited text file. Requires that the faunal data in the original database 
# is first transformed into ones and zeros.
df <- read.csv("Presence_absence.txt", header = T, sep="\t", na.strings = "")

# Set seed to keep results consistent
set.seed(123)

# Color palette created with colorspace
pal <- c("#008599", "#0097A8", "#00A8B7", "#28BAC6", "#7CCCD5", "#ADDDE3", "#D4ECF0", "#F4FAFB", "#FCF8F5", "#F2E6D8", "#E7D1B7", "#D9BA91", "#CAA366", "#BA8B29", "#A77300", "#935B00")

##################################################################################################################################
##################################################################################################################################
# Testing the normality assumption

# QQ-plots for the mobility variables
qqnorm(df$Overall.Group.Mobility)
qqnorm(df$Mean.Effort.Per.Raw.Material)
qqnorm(df$Qnt.of.lithics)
qqnorm(df$MaxDist)
qqnorm(df$Sources)

# Shapiro-Wilk test for normality
shapiro.test(df$Overall.Group.Mobility)
shapiro.test(df$Mean.Effort.Per.Raw.Material)
shapiro.test(df$Qnt.of.lithics)
shapiro.test(df$MaxDist)
shapiro.test(df$Sources)

##################################################################################################################################
##################################################################################################################################
# Analysis 1: Graphing of mobility variables 

# Remove cases without climatic stages
dfadj <- df[complete.cases(df[,17]),]

# Re-configure the climatic stages into an ordered factor
dfadj$reorder <- factor(dfadj$Stage,
                        levels = c('Saalian','Eemian', 'Early Weichselian', 'Late Weichselian'),ordered = TRUE)

# Get midpoints for sources, quantity of lithics, and MaxDist
sourcemidpoint <- (max(dfadj$Sources)+min(dfadj$Sources))/2
quantitymidpoint <- (max(log(dfadj$Qnt.of.lithics))+min(log(dfadj$Qnt.of.lithics)))/2
maxdistmidpoint <- (max(dfadj$MaxDist)+min(dfadj$MaxDist))/2

# Linear models
# Subsetting data to time period
saalianlm <- subset(dfadj, reorder == "Saalian")
eemianlm <- subset(dfadj, reorder == "Eemian")
eweichlm <- subset(dfadj, reorder == "Early Weichselian")
lweichlm <- subset(dfadj, reorder == "Late Weichselian")

# Linear models
saalianmod <- lm(Mean.Effort.Per.Raw.Material~Overall.Group.Mobility, saalianlm)
eemianmod <- lm(Mean.Effort.Per.Raw.Material~Overall.Group.Mobility, eemianlm)
eweichmod <- lm(Mean.Effort.Per.Raw.Material~Overall.Group.Mobility, eweichlm)
lweichmod <- lm(Mean.Effort.Per.Raw.Material~Overall.Group.Mobility, lweichlm)

# Summarise models
summary(saalianmod)
summary(eemianmod)
summary(eweichmod)
summary(lweichmod)

# Plot residuals for each time interval
par(mfrow=c(2,2))
plot(saalianmod)

par(mfrow=c(2,2))
plot(eemianmod)

par(mfrow=c(2,2))
plot(eweichmod)

par(mfrow=c(2,2))
plot(lweichmod)

# Plots
# source_gist call provides the function stat_smooth_func which pastes the regression function and R2-value on the graph
source_gist("524eade46135f6348140")

# Linear regression models
linear <- ggplot(dfadj, aes(Overall.Group.Mobility, Mean.Effort.Per.Raw.Material)) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, xpos=30, ypos=1, col=2, size=2.5) +
  geom_smooth(method="lm",se=FALSE, color="red") +
  geom_point(size=1.2) +
  scale_colour_gradient2(low="#80cdc1", midpoint = sourcemidpoint, mid="#0868ac", high="#dfc27d",
                         name ="Number \nof sources")+
  facet_grid(~reorder)+
  xlab("Overall mobility")+
  ylab("Mean effort \nper raw material")+
  theme(legend.margin=unit(0,"cm"),
        plot.margin = unit(x = c(0, 3.7, 0.25, 0), units = "cm"),
        plot.title = element_text(size=24),
        legend.title=element_text(size=14),
        legend.text = element_text(size=12),
        axis.title=element_text(size=12))

# Sources against mobility
mobsource <- ggplot(dfadj)+
  aes(x = Overall.Group.Mobility, y = Mean.Effort.Per.Raw.Material, colour = Sources)+
  geom_smooth()+
  geom_point(size=1.2)+
  scale_colour_gradient2(low="#80cdc1", midpoint = sourcemidpoint, mid="#0868ac", high="#dfc27d",
                         name ="Number \nof sources")+
  facet_grid(~reorder)+
  xlab("Overall mobility")+
  ylab("Mean effort \nper raw material")+
  theme(legend.margin=unit(0,"cm"),
        plot.margin = unit(x = c(0, .95, 0.25, 0), units = "cm"),
        plot.title = element_text(size=24),
        legend.title=element_text(size=14),
        legend.text = element_text(size=12),
        axis.title=element_text(size=12))

# Quantity of lithics against mobility
mobqnt <- ggplot(dfadj)+
  aes(x = Overall.Group.Mobility, y = Mean.Effort.Per.Raw.Material, colour = log(Qnt.of.lithics))+
  geom_smooth()+
  geom_point(size=1.2)+
  scale_colour_gradient2(low="#80cdc1", midpoint = quantitymidpoint, mid="#0868ac", high="#dfc27d",
                         name ="(Log) Quantity \nof lithics")+
  facet_grid(~reorder)+
  xlab("Overall mobility")+
  ylab("Mean effort \nper raw material")+
  theme(legend.margin=unit(0,"cm"),
        plot.margin = unit(x = c(0, 0.04, 0.25, 0), units = "cm"),
        legend.title=element_text(size=14),
        legend.text = element_text(size=12),
        axis.title=element_text(size=12))

# MaxDist against mobility
mobdist <- ggplot(dfadj)+
  aes(x = Overall.Group.Mobility, y = Mean.Effort.Per.Raw.Material, colour = MaxDist)+
  geom_smooth()+
  geom_point(size=1.2)+
  scale_colour_gradient2(low="#80cdc1", midpoint = maxdistmidpoint, mid="#0868ac", high="#dfc27d",
                         name ="Maximum \ntransport \ndistance")+
  facet_grid(~reorder)+
  xlab("Overall mobility")+
  ylab("Mean effort \nper raw material")+
  theme(legend.margin=unit(0,"cm"),
        plot.margin = unit(x = c(0, 1.05, 0.25, 0), units = "cm"),
        legend.title=element_text(size=14),
        legend.text = element_text(size=12),
        axis.title=element_text(size=12))

# Take all of the above plots and place them in one. Code for multiplot function taken from: 
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot(linear, mobsource, mobqnt, mobdist, cols=1)

##################################################################################################################################
##################################################################################################################################
# Analysis 2: East-west cline tests

# Analysis 2.1.: Kendall's tau
# Change the format of lat & long to numeric
df$Longitude <- as.numeric(df$Longitude)
df$Latitude <- as.numeric(df$Latitude)

# First extract the required columns, then execute correlation
tau <- merge(df[,5:6], df[,27:31], by=,0)
tau <- tau[,2:ncol(tau)]
cor(tau, method="kendall")

# Extract p-value for correlation, code copied from Tyler Rinker: 
# http://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

cor.test.p(tau)


# K-means clustering for latitude and longitude
# Create a dataframe with lat & long columns
cluster <- data.frame(df$Longitude)
cluster$latitude <- df$Latitude

# Determine number of clusters - code copied from: 
# http://www.statmethods.net/advstats/cluster.html
wss <- (nrow(cluster)-1)*sum(apply(cluster,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(cluster, 
                                     centers=i)$withinss)
wssdf <- data.frame(wss)
wssdf$cluster <- c(1:15)

# Create a plot with within groups of sum of squares and number of clusters for elbow/joint analysis
clusterselect <- ggplot(wssdf)+
  aes(x = cluster, y=wss)+
  geom_point(size=3)+
  geom_line()+
  theme(legend.title=element_text(size=20),
        legend.text = element_text(size=14),
        plot.title = element_text(size=24),
        axis.line= element_line(colour="black"),
        axis.title.x= element_text(size= 20, colour="black"),
        axis.title.y= element_text(size= 20, colour="black"),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=10))+
  scale_x_continuous(breaks = c(1:16))+ 
  ylab("Within groups sum of squares")+
  xlab("Number of clusters")+
  ggtitle("Cluster selection for K-means cluster analysis")

# K-Means Cluster Analysis
fit <- kmeans(cluster, 2) # 2 cluster solution

# Get cluster means 
aggregate(cluster,by=list(fit$cluster),FUN=mean)

# Append cluster assignment
cluster_result <- data.frame(cluster, fit$cluster)
cluster_result <- merge(df[,1:4], cluster_result, by=,0)
cluster_result <- cluster_result[,2:ncol(cluster_result)]
cluster_result <- merge(cluster_result, df[,27:31], by=,0)
cluster_result <- cluster_result[,2:ncol(cluster_result)]

# Convex hull, function copied from user Andy W who in turn based it on Gota Morota's work which I could not locate: 
# http://stats.stackexchange.com/questions/22805/how-to-draw-neat-polygons-around-scatterplot-regions-in-ggplot2/22855
find_hull <- function(cluster_result) cluster_result[chull(cluster_result$df.Longitude, cluster_result$latitude), ]
hulls <- ddply(cluster_result, "fit.cluster", find_hull)

# Cluster plot
clusterplot <- ggplot(cluster_result)+
  aes(x = df.Longitude, y = latitude, colour = factor(fit.cluster), fill = factor(fit.cluster))+
  geom_polygon(data=hulls, alpha=0.5)+
  geom_point(size = 4, show.legend = F)+
  scale_fill_discrete(name="Cluster",
                      breaks=c(1,2),
                      labels=c("West", "East"))+
  theme(legend.title=element_text(size=20),
        legend.text = element_text(size=14),
        plot.title = element_text(size=24),
        axis.line= element_line(colour="black"),
        axis.title.x= element_text(size= 20, colour="black"),
        axis.title.y= element_text(size= 20, colour="black"),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=10))+
  guides(colour=F)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("K-means cluster membership")

# Show plot
clusterplot

# Analysis 2.2.: Mann-Whitney U test
wilcox.test(cluster_result$MaxDist ~ cluster_result$fit.cluster, mu = 0, alt="two.sided", conf.int=T, conf.level=0.95, paired=F, exact =T, correct=T)
wilcox.test(cluster_result$Overall.Group.Mobility ~ cluster_result$fit.cluster, mu = 0, alt="two.sided", conf.int=T, conf.level=0.95, paired=F, exact =T, correct=T)
wilcox.test(cluster_result$Mean.Effort.Per.Raw.Material ~ cluster_result$fit.cluster, mu = 0, alt="two.sided", conf.int=T, conf.level=0.95, paired=F, exact =T, correct=T)
wilcox.test(cluster_result$Sources ~ cluster_result$fit.cluster, mu = 0, alt="two.sided", conf.int=T, conf.level=0.95, paired=F, exact =T, correct=T)
wilcox.test(cluster_result$Qnt.of.lithics ~ cluster_result$fit.cluster, mu = 0, alt="two.sided", conf.int=T, conf.level=0.95, paired=F, exact =T, correct=T)

# Analysis 2.3.: T-test for mean effort
t_cluster1 <- cluster_result[cluster_result$fit.cluster == 1,]
t_cluster2 <- cluster_result[cluster_result$fit.cluster == 2,]
t.test(x = t_cluster1$Mean.Effort.Per.Raw.Material, y = t_cluster2$Mean.Effort.Per.Raw.Material, alternative= "two.sided", mu=0, conf.level = 0.95)

##################################################################################################################################
##################################################################################################################################
# Analysis 3: Temporal analysis

# Jonckheere-Terpstra test for trends

# This is the same code as in the first analysis, and to avoid doing this twice it is commented. However, should you wish to only execute J-T test,
# remove the first hashtages from the following two lines.
# dfadj <- df[complete.cases(df[,17]),] # Remove cases without climatic stages
# dfadj$reorder <- factor(dfadj$Stage, levels = c('Saalian','Eemian', 'Early Weichselian', 'Late Weichselian'),ordered = TRUE) # Re-configure the climatic stages into an ordered factor

# Tests:
attach(dfadj)
jonckheere.test(Qnt.of.lithics, reorder, alternative ="increasing")
jonckheere.test(Overall.Group.Mobility, reorder, alternative ="increasing")
jonckheere.test(Mean.Effort.Per.Raw.Material, reorder, alternative ="increasing")
jonckheere.test(MaxDist, reorder, alternative ="increasing")
jonckheere.test(Sources, reorder, alternative ="increasing")
detach(dfadj)

# Graphs of boxplots for quantity of lithics, overall mobility, maxdist, and N of sources
JT_lithics <- ggplot(dfadj) + 
  aes(reorder, log(Qnt.of.lithics), fill=reorder)+
  geom_boxplot()+
  scale_fill_discrete(name="Climatic stage")+
  ylab("(Ln) Quantity of lithics")+
  xlab("Climatic stage")+
  ggtitle("Temporal change in lithic quantities")+
  theme(legend.title=element_text(size=20),
        legend.text = element_text(size=14),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        panel.border= element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=24),
        axis.line= element_line(colour="black"),
        axis.title.x= element_text(size= 20, colour="black"),
        axis.title.y= element_text(size= 20, colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=10))

JT_mobility <- ggplot(dfadj) + 
  aes(reorder, Overall.Group.Mobility, fill=reorder)+
  geom_boxplot()+
  scale_fill_discrete(name="Climatic stage")+
  ylab("Overall mobility")+
  xlab("Climatic stage")+
  ggtitle("Temporal change in overall mobility")+
  theme(legend.title=element_text(size=20),
        legend.text = element_text(size=14),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        panel.border= element_blank(),
        plot.title = element_text(size=24),
        axis.line= element_line(colour="black"),
        axis.title.x= element_text(size= 20, colour="black"),
        axis.title.y= element_text(size= 20, colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=10))

JT_effort <- ggplot(dfadj) + 
  aes(reorder, Mean.Effort.Per.Raw.Material, fill=reorder)+
  geom_boxplot()+
  scale_fill_discrete(name="Climatic stage")+
  ylab("Mean effort per raw material")+
  xlab("Climatic stage")+
  ggtitle("Temporal change in mean effort per raw material")+
  theme(legend.title=element_text(size=20),
        legend.text = element_text(size=14),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        panel.border= element_blank(),
        plot.title = element_text(size=24),
        axis.line= element_line(colour="black"),
        axis.title.x= element_text(size= 20, colour="black"),
        axis.title.y= element_text(size= 20, colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=10))

JT_dist <- ggplot(dfadj) + 
  aes(reorder, MaxDist, fill=reorder)+
  geom_boxplot()+
  scale_fill_discrete(name="Climatic stage")+
  ylab("Maximum transport distance")+
  xlab("Climatic stage")+
  ggtitle("Temporal change in maximum transport distances")+
  theme(legend.title=element_text(size=20),
        legend.text = element_text(size=14),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        panel.border= element_blank(),
        plot.title = element_text(size=24),
        axis.line= element_line(colour="black"),
        axis.title.x= element_text(size= 20, colour="black"),
        axis.title.y= element_text(size= 20, colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=10))

JT_sources <- ggplot(dfadj) + 
  aes(reorder, Sources, fill=reorder)+
  geom_boxplot()+
  scale_fill_discrete(name="Climatic stage")+
  ylab("Sources")+
  xlab("Climatic stage")+
  ggtitle("Temporal change in the number of utilised sources")+
  theme(legend.title=element_text(size=20),
        legend.text = element_text(size=14),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        panel.border= element_blank(),
        plot.title = element_text(size=24),
        axis.line= element_line(colour="black"),
        axis.title.x= element_text(size= 20, colour="black"),
        axis.title.y= element_text(size= 20, colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=10))

# Show plots
JT_lithics
JT_mobility
JT_effort
JT_dist
JT_sources

##################################################################################################################################
##################################################################################################################################
# Analysis 4: Palaeoenvironmental grouping and analysis of mobility variables

# Setup:
# Reduce the number of assemblages to only include assemblages with two or more ungulate species
subset1 <- df[df$N_Ungulates >= 2, ]

# Choose only ungulates as variables
subset2 <- subset1[c(1:35,62:90)]

# Choose only ungulates that occur more than twice
subset4 <- subset2[,36:64]
subset4 <- subset4[colSums(subset4) > 2]
subset4 <- merge(subset2[,1:35], subset4, by=,0)
subset4 <- subset4[,2:ncol(subset4)]

# Analysis 4.1.: Logistic PCA 
# A worked example of the logistic PCA available here:
# https://cran.r-project.org/web/packages/logisticPCA/vignettes/logisticPCA.html

# Select only ungulates for analysis
binPCA <- subset4[c(36:ncol(subset4))]

# Cross-validation to get the tuning parameter m
logpca_cv = cv.lpca(binPCA, ks = 2:8, ms = 1:10)

# Create a similar (but modifiable) cross-validation plot as with command: plot(logpca_cv)
lpca_cv_plot <- data.frame(logpca_cv)
lpca_cv_plot <- t(lpca_cv_plot)
lpca_cv_plot <- data.frame(lpca_cv_plot)
m <- c(1:10)
lpca_cv_plot <- cbind(lpca_cv_plot, m = m)
lpca_cv_plot <- melt(lpca_cv_plot, id="m")

ggplot(lpca_cv_plot)+
  aes(x=m, y=value, colour=variable)+
  geom_line(size=2)+
  ylab("Negative Log-likelihood")+
  scale_colour_discrete(name="No. of PCs",
                        labels=c(2:8))+  
  scale_x_continuous(breaks = c(0:12),
                     labels = c(0:12))+
  ggtitle("Cross-validation of Logistic PCA models")+
  theme(legend.title=element_text(size=20),
        legend.text = element_text(size=14),
        panel.grid.major= element_line(size=1),
        panel.grid.minor= element_blank(),
        panel.border= element_blank(),
        plot.title = element_text(size=24),
        axis.line= element_line(colour="black"),
        axis.title.x= element_text(size= 20, colour="black"),
        axis.title.y= element_text(size= 20, colour="black"),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=10))

# Logistic PCA grouping
logpca_model = logisticPCA(binPCA, k = 5, m=4, main_effects = T)
logpca_model2 <- data.frame(logpca_model$PCs) 

# Component loadings and transposing them
loadings <- data.frame(logpca_model$U)
names(loadings) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
loadings_t <- t(loadings) 
loadings_t <- data.frame(loadings_t)

# Rebuilding the dataframe
logPCA_merge <- merge(subset4[,1:35], logpca_model2, by=,0)
logPCA_merge <- logPCA_merge[,2:ncol(logPCA_merge)]
logPCA_merge <- logPCA_merge[complete.cases(logPCA_merge[,17]),]
logPCA_merge$reorder <- factor(logPCA_merge$Stage, c("Saalian", "Eemian", "Early Weichselian", "Late Weichselian"))

# Explain PCs by their component loadings with heatmap
# Label columns
names(loadings_t)[1:21] <- c("S. scrofa", "C. elaphus", "C. simplicidens", "D. dama", "M. giganteus", "C. capreolus", "R. tarandus", "S. tatarica",
                             "R. rupicapra", "R. pyrenaica", "C. ibex", "C. pyrenaica", "H. bonali", "B. priscus", "B. primigenius", "M. primigenius",
                             "P. antiquus", "S. hemitoechus", "C. antiquitatis", "E. caballus", "E. hydruntinus")
# Heat map
pheatmap(as.matrix(loadings_t), color = pal, border_color = "black", cellwidth = 30, cellheight = 30, display_numbers = T, scale="none", cluster_rows = F, legend = T, show_rownames = T, show_colnames = T, main = "PCA grouping based on component loadings", labels_row = c("PC 1", "PC 2","PC 3","PC 4","PC 5"))

##################################################################################################################################
##################################################################################################################################
# Analysis 4.2.: LCA analyses

# Bayesian LCA with bootstrapped Expectation Maximization algorithm for ungulates in subset4 dataframe
EM2 <- blca.boot(subset4[,36:ncol(subset4)], 2, start.vals = "across", iter = 200, B = 1000)
EM3 <- blca.boot(subset4[,36:ncol(subset4)], 3, start.vals = "across", iter = 200, B = 1000)
EM4 <- blca.boot(subset4[,36:ncol(subset4)], 4, start.vals = "across", iter = 200, B = 1000)
EM5 <- blca.boot(subset4[,36:ncol(subset4)], 5, start.vals = "across", iter = 200, B = 1000)
EM6 <- blca.boot(subset4[,36:ncol(subset4)], 6, start.vals = "across", iter = 200, B = 1000)
EM7 <- blca.boot(subset4[,36:ncol(subset4)], 7, start.vals = "across", iter = 200, B = 1000)
EM8 <- blca.boot(subset4[,36:ncol(subset4)], 8, start.vals = "across", iter = 200, B = 1000)

# Bayesian LCA with Gibbs sampling
gibbs2 <- blca.gibbs(subset4[,36:ncol(subset4)], 2, relabel = T, burn.in = 150, thin = 1/10, iter = 100000)
gibbs3 <- blca.gibbs(subset4[,36:ncol(subset4)], 3, relabel = T, burn.in = 150, thin = 1/10, iter = 100000)
gibbs4 <- blca.gibbs(subset4[,36:ncol(subset4)], 4, relabel = T, burn.in = 150, thin = 1/10, iter = 100000)
gibbs5 <- blca.gibbs(subset4[,36:ncol(subset4)], 5, relabel = T, burn.in = 150, thin = 1/10, iter = 100000)
gibbs6 <- blca.gibbs(subset4[,36:ncol(subset4)], 6, relabel = T, burn.in = 150, thin = 1/10, iter = 100000)
gibbs7 <- blca.gibbs(subset4[,36:ncol(subset4)], 7, relabel = T, burn.in = 150, thin = 1/10, iter = 100000)
gibbs8 <- blca.gibbs(subset4[,36:ncol(subset4)], 8, relabel = T, burn.in = 150, thin = 1/10, iter = 100000)

# AIC model fits for EM models
c(EM2$AIC, EM3$AIC, EM4$AIC, EM5$AIC, EM6$AIC, EM7$AIC, EM8$AIC) 

# Gibbs sampling model fit
c(gibbs2$DIC, gibbs3$DIC, gibbs4$DIC, gibbs5$DIC, gibbs6$DIC, gibbs7$DIC, gibbs8$DIC)

# Heat map plot for EM5 model
df1 <- data.frame(EM5$itemprob)
names(df1)[1:21] <- c("S. scrofa", "C. elaphus", "C. simplicidens", "D. dama", "M. giganteus", "C. capreolus", "R. tarandus", "S. tatarica",
                "R. rupicapra", "R. pyrenaica", "C. ibex", "C. pyrenaica", "H. bonali", "B. priscus", "B. primigenius", "M. primigenius",
                "P. antiquus", "S. hemitoechus", "C. antiquitatis", "E. caballus", "E. hydruntinus")
pheatmap(as.matrix(df1), color = pal, border_color = "black", cellwidth = 30, cellheight = 30, scale="none", display_numbers = T, show_rownames = T, show_colnames = T, main = "Group membership probabilities \n(EM algorithm)", labels_row = c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"))

# Heat map plot for gibbs4 model
df2 <- data.frame(gibbs4$itemprob)
names(df2)[1:21] <- c("S. scrofa", "C. elaphus", "C. simplicidens", "D. dama", "M. giganteus", "C. capreolus", "R. tarandus", "S. tatarica",
                      "R. rupicapra", "R. pyrenaica", "C. ibex", "C. pyrenaica", "H. bonali", "B. priscus", "B. primigenius", "M. primigenius",
                      "P. antiquus", "S. hemitoechus", "C. antiquitatis", "E. caballus", "E. hydruntinus")
pheatmap(as.matrix(df2), color = pal, border_color = "black", cellwidth = 30, cellheight = 30, scale="none", display_numbers = T, legend = T, legend_breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), show_rownames = T, show_colnames = T, main = "Group membership probabilities \n(Gibbs sampling)", labels_row = c("Group 1", "Group 2", "Group 3", "Group 4"))


# Create a dataframe with four columns showing the probability of group membership
gibbsdf <- data.frame(gibbs4$Z)
gibbsdf <- merge(subset4[,1:35], gibbs4$Z, by=,0)

# Drop sites without climatic stages
gibbsdf1 <- gibbsdf[complete.cases(gibbsdf[,18]),]
gibbsdf1 <- gibbsdf1[,2:ncol(gibbsdf1)]

# Reorder by stage
gibbsdf1$reorder <- factor(gibbsdf1$Stage, c("Saalian", "Eemian", "Early Weichselian", "Late Weichselian"))


##################################################################################################################################
##################################################################################################################################
# Kendall's tau for LCA and mobility 
gibbstau <- merge(gibbsdf[,28:32], gibbsdf[,37:40], by=,0)
gibbstau <- gibbstau[,2:ncol(gibbstau)]
gibbstau_result <- cor(gibbstau, method="kendall")

# cor.test.p same function as in the first correlation analysis
gibbstau_p <- cor.test.p(gibbstau_result)

##################################################################################################################################
##################################################################################################################################
# Kendall's tau for PCA and mobility 
kent <- merge(logPCA_merge[,27:31], logPCA_merge[,36:40], by=,0)
kent <- kent[,2:11]
kent_result <- cor(kent, method="kendall")

# cor.test.p same function as in the first correlation analysis
kent_p <- cor.test.p(kent_result)

##################################################################################################################################
##################################################################################################################################
# Scatterplots of PCs against mobility variables
# PCs & overall mobility
pcplot_mob_1 <- ggplot(logPCA_merge)+
  aes(x=Overall.Group.Mobility, y=X1)+
  geom_point()+
  xlab("Overall \nmobility")+
  ylab("PC 1")

pcplot_mob_2 <- ggplot(logPCA_merge)+
  aes(x=Overall.Group.Mobility, y=X2)+
  geom_point()+
  xlab("Overall \nmobility")+
  ylab("PC 2")

pcplot_mob_3 <- ggplot(logPCA_merge)+
  aes(x=Overall.Group.Mobility, y=X3)+
  geom_point()+
  xlab("Overall \nmobility")+
  ylab("PC 3")

pcplot_mob_4 <- ggplot(logPCA_merge)+
  aes(x=Overall.Group.Mobility, y=X4)+
  geom_point()+
  xlab("Overall \nmobility")+
  ylab("PC 4")

pcplot_mob_5 <- ggplot(logPCA_merge)+
  aes(x=Overall.Group.Mobility, y=X5)+
  geom_point()+
  xlab("Overall \nmobility")+
  ylab("PC 5")

# PCs and mean effort
pcplot_eff_1 <- ggplot(logPCA_merge)+
  aes(x=Mean.Effort.Per.Raw.Material, y=X1)+
  geom_point()+
  xlab("Mean effort \nper raw material")+
  ylab("PC 1")

pcplot_eff_2 <- ggplot(logPCA_merge)+
  aes(x=Mean.Effort.Per.Raw.Material, y=X2)+
  geom_point()+
  xlab("Mean effort \nper raw material")+
  ylab("PC 2")

pcplot_eff_3 <- ggplot(logPCA_merge)+
  aes(x=Mean.Effort.Per.Raw.Material, y=X3)+
  geom_point()+
  xlab("Mean effort \nper raw material")+
  ylab("PC 3")

pcplot_eff_4 <- ggplot(logPCA_merge)+
  aes(x=Mean.Effort.Per.Raw.Material, y=X4)+
  geom_point()+
  xlab("Mean effort \nper raw material")+
  ylab("PC 4")

pcplot_eff_5 <- ggplot(logPCA_merge)+
  aes(x=Mean.Effort.Per.Raw.Material, y=X5)+
  geom_point()+
  xlab("Mean effort \nper raw material")+
  ylab("PC 5")

# PCs & transport distance
pcplot_dist_1 <- ggplot(logPCA_merge)+
  aes(x=MaxDist, y=X1)+
  geom_point()+
  xlab("Maximum \ntransport distance")+
  ylab("PC 1")

pcplot_dist_2 <- ggplot(logPCA_merge)+
  aes(x=MaxDist, y=X2)+
  geom_point()+
  xlab("Maximum \ntransport distance")+
  ylab("PC 2")

pcplot_dist_3 <- ggplot(logPCA_merge)+
  aes(x=MaxDist, y=X3)+
  geom_point()+
  xlab("Maximum \ntransport distance")+
  ylab("PC 3")

pcplot_dist_4 <- ggplot(logPCA_merge)+
  aes(x=MaxDist, y=X4)+
  geom_point()+
  xlab("Maximum \ntransport distance")+
  ylab("PC 4")

pcplot_dist_5 <- ggplot(logPCA_merge)+
  aes(x=MaxDist, y=X5)+
  geom_point()+
  xlab("Maximum \ntransport distance")+
  ylab("PC 5")

# PCs & number of sources
pcplot_src_1 <- ggplot(logPCA_merge)+
  aes(x=Sources, y=X1)+
  geom_point()+
  xlab("No. of \nsources")+
  ylab("PC 1")

pcplot_src_2 <- ggplot(logPCA_merge)+
  aes(x=Sources, y=X2)+
  geom_point()+
  xlab("No. of \nsources")+
  ylab("PC 2")

pcplot_src_3 <- ggplot(logPCA_merge)+
  aes(x=Sources, y=X3)+
  geom_point()+
  xlab("No. of \nsources")+
  ylab("PC 3")

pcplot_src_4 <- ggplot(logPCA_merge)+
  aes(x=Sources, y=X4)+
  geom_point()+
  xlab("No. of \nsources")+
  ylab("PC 4")

pcplot_src_5 <- ggplot(logPCA_merge)+
  aes(x=Sources, y=X5)+
  geom_point()+
  xlab("No. of \nsources")+
  ylab("PC 5")

# PCs & log(quantity of lithics)
pcplot_qnt_1 <- ggplot(logPCA_merge)+
  aes(x=log(Qnt.of.lithics), y=X1)+
  geom_point()+
  xlab("(Ln) quantity \nof lithics")+
  ylab("PC 1")

pcplot_qnt_2 <- ggplot(logPCA_merge)+
  aes(x=log(Qnt.of.lithics), y=X2)+
  geom_point()+
  xlab("(Ln) quantity \nof lithics")+
  ylab("PC 2")

pcplot_qnt_3 <- ggplot(logPCA_merge)+
  aes(x=log(Qnt.of.lithics), y=X3)+
  geom_point()+
  xlab("(Ln) quantity \nof lithics")+
  ylab("PC 3")

pcplot_qnt_4 <- ggplot(logPCA_merge)+
  aes(x=log(Qnt.of.lithics), y=X4)+
  geom_point()+
  xlab("(Ln) quantity \nof lithics")+
  ylab("PC 4")

pcplot_qnt_5 <- ggplot(logPCA_merge)+
  aes(x=log(Qnt.of.lithics), y=X5)+
  geom_point()+
  xlab("(Ln) quantity \nof lithics")+
  ylab("PC 5")

# Multiplot
multiplot(pcplot_mob_1, pcplot_mob_2, pcplot_mob_3, pcplot_mob_4, pcplot_mob_5, 
          pcplot_eff_1, pcplot_eff_2, pcplot_eff_3, pcplot_eff_4, pcplot_eff_5,
          pcplot_dist_1, pcplot_dist_2, pcplot_dist_3, pcplot_dist_4, pcplot_dist_5,
          pcplot_src_1, pcplot_src_2, pcplot_src_3, pcplot_src_4, pcplot_src_5,
          pcplot_qnt_1, pcplot_qnt_2, pcplot_qnt_3, pcplot_qnt_4, pcplot_qnt_5,
          cols=5)


##################################################################################################################################
##################################################################################################################################
# Scatterplots of Latent Classes against mobility variables
# LCs & overall mobility
lcplot_mob_1 <- ggplot(gibbsdf)+
  aes(x=Overall.Group.Mobility, y=`Group 1`)+
  geom_point()+
  xlab("Overall \nmobility")

lcplot_mob_2 <- ggplot(gibbsdf)+
  aes(x=Overall.Group.Mobility, y=`Group 2`)+
  geom_point()+
  xlab("Overall \nmobility")
  
lcplot_mob_3 <- ggplot(gibbsdf)+
  aes(x=Overall.Group.Mobility, y=`Group 3`)+
  geom_point()+
  xlab("Overall \nmobility")
  
lcplot_mob_4 <- ggplot(gibbsdf)+
  aes(x=Overall.Group.Mobility, y=`Group 4`)+
  geom_point()+
  xlab("Overall \nmobility")
  
# LCs and mean effort
lcplot_eff_1 <- ggplot(gibbsdf)+
  aes(x=Mean.Effort.Per.Raw.Material, y=`Group 1`)+
  geom_point()+
  xlab("Mean effort \nper raw material")

lcplot_eff_2 <- ggplot(gibbsdf)+
  aes(x=Mean.Effort.Per.Raw.Material, y=`Group 2`)+
  geom_point()+
  xlab("Mean effort \nper raw material")

lcplot_eff_3 <- ggplot(gibbsdf)+
  aes(x=Mean.Effort.Per.Raw.Material, y=`Group 3`)+
  geom_point()+
  xlab("Mean effort \nper raw material")

lcplot_eff_4 <- ggplot(gibbsdf)+
  aes(x=Mean.Effort.Per.Raw.Material, y=`Group 4`)+
  geom_point()+
  xlab("Mean effort \nper raw material")

# LCs & transport distance
lcplot_dist_1 <- ggplot(gibbsdf)+
  aes(x=MaxDist, y=`Group 1`)+
  geom_point()+
  xlab("Maximum \ntransport distance")

lcplot_dist_2 <- ggplot(gibbsdf)+
  aes(x=MaxDist, y=`Group 2`)+
  geom_point()+
  xlab("Maximum \ntransport distance")

lcplot_dist_3 <- ggplot(gibbsdf)+
  aes(x=MaxDist, y=`Group 3`)+
  geom_point()+
  xlab("Maximum \ntransport distance")

lcplot_dist_4 <- ggplot(gibbsdf)+
  aes(x=MaxDist, y=`Group 4`)+
  geom_point()+
  xlab("Maximum \ntransport distance")

# LCs & number of sources
lcplot_src_1 <- ggplot(gibbsdf)+
  aes(x=Sources, y=`Group 1`)+
  geom_point()+
  xlab("No. of \nsources")

lcplot_src_2 <- ggplot(gibbsdf)+
  aes(x=Sources, y=`Group 2`)+
  geom_point()+
  xlab("No. of \nsources")

lcplot_src_3 <- ggplot(gibbsdf)+
  aes(x=Sources, y=`Group 3`)+
  geom_point()+
  xlab("No. of \nsources")

lcplot_src_4 <- ggplot(gibbsdf)+
  aes(x=Sources, y=`Group 4`)+
  geom_point()+
  xlab("No. of \nsources")

# LCs & log(quantity of lithics)
lcplot_qnt_1 <- ggplot(gibbsdf)+
  aes(x=log(Qnt.of.lithics), y=`Group 1`)+
  geom_point()+
  xlab("(Ln) quantity \nof lithics")

lcplot_qnt_2 <- ggplot(gibbsdf)+
  aes(x=log(Qnt.of.lithics), y=`Group 2`)+
  geom_point()+
  xlab("(Ln) quantity \nof lithics")

lcplot_qnt_3 <- ggplot(gibbsdf)+
  aes(x=log(Qnt.of.lithics), y=`Group 3`)+
  geom_point()+
  xlab("(Ln) quantity \nof lithics")

lcplot_qnt_4 <- ggplot(gibbsdf)+
  aes(x=log(Qnt.of.lithics), y=`Group 4`)+
  geom_point()+
  xlab("(Ln) quantity \nof lithics")


# Multiplot
multiplot(lcplot_mob_1, lcplot_mob_2, lcplot_mob_3, lcplot_mob_4, 
          lcplot_eff_1, lcplot_eff_2, lcplot_eff_3, lcplot_eff_4,
          lcplot_dist_1, lcplot_dist_2, lcplot_dist_3, lcplot_dist_4,
          lcplot_src_1, lcplot_src_2, lcplot_src_3, lcplot_src_4,
          lcplot_qnt_1, lcplot_qnt_2, lcplot_qnt_3, lcplot_qnt_4,
          cols=5)



##################################################################################################################################
##################################################################################################################################
# Generalized Linear Model for LCA
### Overall mobility and environment
quasigm1 <- glm(`Group 1`~Overall.Group.Mobility,data=gibbsdf, family=quasibinomial(link="logit"))
quasigm2 <- glm(`Group 2`~Overall.Group.Mobility,data=gibbsdf, family=quasibinomial(link="logit"))
quasigm3 <- glm(`Group 3`~Overall.Group.Mobility,data=gibbsdf, family=quasibinomial(link="logit"))
quasigm4 <- glm(`Group 4`~Overall.Group.Mobility,data=gibbsdf, family=quasibinomial(link="logit"))

# Summaries 
summary(quasigm1)
summary(quasigm2)
summary(quasigm3)
summary(quasigm4)

# Chi-square baseline
chitest1 <- glm(`Group 1`~1,data=gibbsdf, family=quasibinomial(link="logit"))
chitest2 <- glm(`Group 2`~1,data=gibbsdf, family=quasibinomial(link="logit"))
chitest3 <- glm(`Group 3`~1,data=gibbsdf, family=quasibinomial(link="logit"))
chitest4 <- glm(`Group 4`~1,data=gibbsdf, family=quasibinomial(link="logit"))

# Chi-square tests 
anova(chitest1,quasigm1, test="Chisq")
anova(chitest2,quasigm2, test="Chisq")
anova(chitest3,quasigm3, test="Chisq")
anova(chitest4,quasigm4, test="Chisq")

# Return the increase in the odds (& CI 95%) of the site belonging to the particular environment
# with one unit increase in overall mobility
exp(cbind(OR = coef(quasigm1), confint(quasigm1)))
exp(cbind(OR = coef(quasigm2), confint(quasigm2)))
exp(cbind(OR = coef(quasigm3), confint(quasigm3)))
exp(cbind(OR = coef(quasigm4), confint(quasigm4)))


### Mean effort and environment
quasigmean1 <- glm(`Group 1`~Mean.Effort.Per.Raw.Material,data=gibbsdf, family=quasibinomial(link="logit"))
quasigmean2 <- glm(`Group 2`~Mean.Effort.Per.Raw.Material,data=gibbsdf, family=quasibinomial(link="logit"))
quasigmean3 <- glm(`Group 3`~Mean.Effort.Per.Raw.Material,data=gibbsdf, family=quasibinomial(link="logit"))
quasigmean4 <- glm(`Group 4`~Mean.Effort.Per.Raw.Material,data=gibbsdf, family=quasibinomial(link="logit"))

# Summaries
summary(quasigmean1)
summary(quasigmean2)
summary(quasigmean3)
summary(quasigmean4)

# Chi-square test
anova(chitest1,quasigmean1, test="Chisq")
anova(chitest2,quasigmean2, test="Chisq")
anova(chitest3,quasigmean3, test="Chisq")
anova(chitest4,quasigmean4, test="Chisq")

# Return the increase in the odds (& CI 95%) of the site belonging to the particular environment
# with one unit increase in mean effort per raw material
exp(cbind(OR = coef(quasigmean1), confint(quasigmean1)))
exp(cbind(OR = coef(quasigmean2), confint(quasigmean2)))
exp(cbind(OR = coef(quasigmean3), confint(quasigmean3)))
exp(cbind(OR = coef(quasigmean4), confint(quasigmean4)))


### Maximum transport distance and environment
quasigdist1 <- glm(`Group 1`~MaxDist,data=gibbsdf, family=quasibinomial(link="logit"))
quasigdist2 <- glm(`Group 2`~MaxDist,data=gibbsdf, family=quasibinomial(link="logit"))
quasigdist3 <- glm(`Group 3`~MaxDist,data=gibbsdf, family=quasibinomial(link="logit"))
quasigdist4 <- glm(`Group 4`~MaxDist,data=gibbsdf, family=quasibinomial(link="logit"))

# Summaries
summary(quasigdist1)
summary(quasigdist2)
summary(quasigdist3)
summary(quasigdist4)

# Chi-square test
anova(chitest1,quasigdist1, test="Chisq")
anova(chitest2,quasigdist2, test="Chisq")
anova(chitest3,quasigdist3, test="Chisq")
anova(chitest4,quasigdist4, test="Chisq")

# Return the increase in the odds (& CI 95%) of the site belonging to the particular environment
# with one unit increase in maximum transport distance
exp(cbind(OR = coef(quasigdist1), confint(quasigdist1)))
exp(cbind(OR = coef(quasigdist2), confint(quasigdist2)))
exp(cbind(OR = coef(quasigdist3), confint(quasigdist3)))
exp(cbind(OR = coef(quasigdist4), confint(quasigdist4)))

### Sources and environment
quasisource1 <- glm(`Group 1`~+Sources, data=gibbsdf, family=quasibinomial(link="logit"))
quasisource2 <- glm(`Group 2`~+Sources,data=gibbsdf, family=quasibinomial(link="logit"))
quasisource3 <- glm(`Group 3`~+Sources,data=gibbsdf, family=quasibinomial(link="logit"))
quasisource4 <- glm(`Group 4`~+Sources,data=gibbsdf, family=quasibinomial(link="logit"))

# Summaries
summary(quasisource1)
summary(quasisource2)
summary(quasisource3)
summary(quasisource4)

# Chi-square test
anova(chitest1,quasisource1, test="Chisq")
anova(chitest2,quasisource2, test="Chisq")
anova(chitest3,quasisource3, test="Chisq")
anova(chitest4,quasisource4, test="Chisq")

# Return the increase in the odds (& CI 95%) of the site belonging to the particular environment
# with one unit increase in number of sources
exp(cbind(OR = coef(quasisource1), confint(quasisource1)))
exp(cbind(OR = coef(quasisource2), confint(quasisource2)))
exp(cbind(OR = coef(quasisource3), confint(quasisource3)))
exp(cbind(OR = coef(quasisource4), confint(quasisource4)))

### Quantity of lithics and environment 
quasiqnt1 <- glm(`Group 1`~ log(Qnt.of.lithics), data=gibbsdf, family=quasibinomial(link="logit"))
quasiqnt2 <- glm(`Group 2`~ log(Qnt.of.lithics), data=gibbsdf, family=quasibinomial(link="logit"))
quasiqnt3 <- glm(`Group 3`~ log(Qnt.of.lithics), data=gibbsdf, family=quasibinomial(link="logit"))
quasiqnt4 <- glm(`Group 4`~ log(Qnt.of.lithics), data=gibbsdf, family=quasibinomial(link="logit"))

# Summaries
summary(quasiqnt1)
summary(quasiqnt2)
summary(quasiqnt3)
summary(quasiqnt4)

# Chi-square test
anova(chitest1,quasiqnt1, test="Chisq")
anova(chitest2,quasiqnt2, test="Chisq")
anova(chitest3,quasiqnt3, test="Chisq")
anova(chitest4,quasiqnt4, test="Chisq")

# Return the increase in the odds (& CI 95%) of the site belonging to the particular environment
# with one unit increase in quantity of lithics
exp(cbind(OR = coef(quasiqnt1), confint(quasiqnt1)))
exp(cbind(OR = coef(quasiqnt2), confint(quasiqnt2)))
exp(cbind(OR = coef(quasiqnt3), confint(quasiqnt3)))
exp(cbind(OR = coef(quasiqnt4), confint(quasiqnt4)))
