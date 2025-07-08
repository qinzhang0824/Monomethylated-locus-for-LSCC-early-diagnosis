library("ggpubr")
library(ROCR)
library(randomForest)
library(caret)
library(pROC)
library(dplyr)
library(varSelRF)
library(glmnet)
library(RColorBrewer)
args=commandArgs(T)

setwd(args[3])

Data<-read.table(args[1],sep='\t',header=T)
Data$Type <- factor(Data$Type,levels=c("Normal","Tumor"))

my_comparisons <- list(c("Normal","Tumor"))

pdf(args[2],width=3,height=5)
p <- ggboxplot(Data, x = "Type", y = "cg10364040",color = "Type",add = "jitter",palette = c("#ED0000","blue"),ylab = "TCGA cg10364040 beta value")+stat_compare_means(comparisons=my_comparisons,label = "p.format")+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
p
dev.off()
