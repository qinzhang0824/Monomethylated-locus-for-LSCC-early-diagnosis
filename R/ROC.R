library(ROCR)
library(randomForest)
library(caret)
library(pROC)
library(dplyr)
library(varSelRF)
library(glmnet)

da <- read.table("TCGA.HNSC.GSE40279.352s.437s.picked.beta_boxplot_input.xls",header=TRUE,sep="\t")

data <-as.data.frame(da)
data$Type <-factor(data$Type,levels=c("Normal","Tumor"))

glm.model <- glm(Group ~ cg10364040, family = "binomial", data = data)

prediction <- predict(glm.model, data, type="response")

proc<-roc(data$Type,data$cg10364040,ci=T)

pdf('TCGA_GEO_ROC.pdf',width=8,height=8)
plot(proc,
     print.auc=TRUE, 
     print.thres=TRUE,
     main="ROC",
     col="Red")

dev.off()
