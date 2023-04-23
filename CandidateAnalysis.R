rm(list=ls())
setwd("M:/w2k/BigDataBio/Brassica Report")
load("brassicaData.rda")
library(ggplot2)
#####

# Creating list of lm- of correlation between genes and trait

lmlist <- lapply(geneNames, function(gene){
  lm(merge.trait$Trait~merge.trait[,gene])
})
names(lmlist) <- geneNames

#Creating a list of ANOVA results of lms

anovalist <- lapply(lmlist, function(gene){
  anova(gene)
})
names(anovalist) <- geneNames

#Creating a df of ANOVA results

results.trait <- as.data.frame(matrix(nrow = 0, ncol = 8))
colnames(results.trait) <- c("Df", "Sum.Sq", "Mean.Sq", "F.value", "P.value", "R2", "Intercept", "Gradient")

for (gene in 1:length(geneNames)){
  anova <- as.data.frame((anovalist[[gene]])[1,])
  intercept <- as.data.frame(coefficients(lmlist[[gene]]))[1,1]
  gradient <- as.data.frame(coefficients(lmlist[[gene]]))[2,1]
  R2 <- as.data.frame(summary(lmlist[[gene]])$r.squared)
  nextrow <- as.data.frame(c(anova, R2, intercept, gradient))
  colnames(nextrow) <- colnames(results.trait)
  results.trait <- rbind(results.trait, nextrow)
}
rownames(results.trait) <- geneNames

#removing unnecessary variables from image
rm(anova, R2, intercept, gradient, nextrow, gene)


#removing unexpressed genes.

results.trait <- results.trait[-which(is.na(results.trait[,"Gradient"])==T),]

#One degrees of freedom therefore remove first column

results.trait <- results.trait[,-1]

#order data frame in ascending p value

results.trait <- results.trait[order(results.trait$P.value, decreasing = F),]

#####
#SUBSETTING SIGNIFICANT GENES#
pval <- 0.0001
merge.trait.candidates <- subset(results.trait, P.value<pval)
geneNames.candidates <- rownames(merge.trait.candidates)

#####
#Function to plot lm

plot.trait <- function(gene){
  plot(merge.trait[,gene], merge.trait$Trait,
       xlab = paste(gene, "RPKM"),
       ylab = "Trait"
       )
  abline(lm(merge.trait$Trait~merge.trait[,gene]), col = "red")
}

#####  
##GRAPHING CANDIDATE GENES TO TRAIT.
dev.new(width = 4.13, height = 5.83, unit = "in")
par(mfrow=c(4,4))
par(mar = c(1, 0.8, 0.8, 0.4) + 0.02)
for (i in 1:16){
  plot(merge.trait[,geneNames.candidates[i]], merge.trait$Trait,
       xlab = geneNames.candidates[i], ylab = NA,
       xaxt = "n", yaxt = "n",
       mgp = c(0, 1, 0))
  abline(lm(merge.trait$Trait~merge.trait[,geneNames.candidates[i]]), col ="red")
}

multiplot.candidate <- recordPlot()

pdf(file = "candidates_trait_multiplot.pdf")
multiplot.candidate
dev.off()

#####

qqnorm(results.trait$P.value, pch = 1, frame = FALSE)
qqline(results.trait$P.value, col = "steelblue", lwd = 2)
library("car")
qqPlot(results.trait$P.value)


boxplot(x = merge.trait.candidates$Gradient)
#####

pairs(merge.trait[,geneNames.candidates])

pairs

#####

save(results.trait, geneNames.candidates, file = "candidateAnalysis.rda")

load("candidateAnalysis.rda")
