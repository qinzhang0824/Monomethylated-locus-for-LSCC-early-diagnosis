library(ChAMP)
library(ROCR)
library(randomForest)
library(caret)
library(pROC)
library(dplyr)
library(varSelRF)
library(glmnet)


setwd('./')
hnsc_info=read.table("TCGA.HNSC.HumanMethylation450.all.picked.sample.xls",header=T,sep="\t")
rownames(hnsc_info)<-hnsc_info$Sample_Name
hnsc_info$Primary_Tumor_Site <- as.vector(hnsc_info$Primary_Tumor_Site)
hnsc_beta <- read.table("TCGA.HNSC.HumanMethylation450.all.picked.beta.xls",header=T,row.names=1,sep="\t")
hnsc_beta <- as.matrix(hnsc_beta)

myLoad<-champ.filter(beta=hnsc_beta,pd=hnsc_info,fixOutlier=FALSE)
myimpute <- champ.impute(beta=myLoad$beta,pd=myLoad$pd)
myimpute <- champ.impute(beta=myLoad$beta,pd=myLoad$pd,ProbeCutoff=0.7,SampleCutoff=0.3)
champ.QC(beta=hnsc_beta,pheno=hnsc_info$Primary_Tumor_Site,Rplot=FALSE)

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
QC.GUI(beta=myNorm,pheno=myLoad$pd$Primary_Tumor_Site,arraytype="450k")

myCombat <- champ.runCombat(beta = myNorm,
                            pd = myLoad$pd,
                            batchname = c("Slide"),
                            variablename = c("Sample_Type"))

myDMPPBC <- champ.DMP(beta=myNormPBC,pheno=myimpute$pd$Primary_Tumor_Site)
myDMPBMIQ <- champ.DMP(beta=myNormBMIQ,pheno=myimpute$pd$Primary_Tumor_Site)

logfc_cutoff <- 0.3
df_dmp_pbc <- (myDMPPBC$Tumor_to_Normal[,1:9])
df_dmp_bmiq <- (myDMPBMIQ$Tumor_to_Normal[,1:9])

df_dmp_pbc$change <- ifelse(df_dmp_pbc$adj.P.Val < 0.05 & abs(df_dmp_pbc$logFC) > logfc_cutoff,ifelse(df_dmp_pbc$logFC > logfc_cutoff ,'UP','DOWN'),'NOT')
df_dmp_bmiq$change <- ifelse(df_dmp_bmiq$adj.P.Val < 0.05 & abs(df_dmp_bmiq$logFC) > logfc_cutoff,ifelse(df_dmp_bmiq$logFC > logfc_cutoff ,'UP','DOWN'),'NOT')

id_pbc <- rownames(df_dmp_pbc[1])
write.table(data.frame(id_pbc,df_dmp_pbc), file="Total_DMP_PBCNormalized_tumor_vs_Para_addAnnote.xls",col.names=T,row.names=F,quote=F,sep="\t")

id_bmiq <- rownames(df_dmp_bmiq[1])
write.table(data.frame(id_bmiq,df_dmp_bmiq), file="Total_DMP_BMIQNormalized_tumor_vs_Para_addAnnote.xls",col.names=T,row.names=F,quote=F,sep="\t")

delete <- read.table('Delete_sample_list',header=T)
af_delete <- myNormBMIQ[,(!(myimpute$pd$Sample_Name %in% delete$DeleteID))]
mygroup <- myimpute$pd[!(myimpute$pd$Sample_Name %in% delete$DeleteID),]

tumor_sample <- mygroup[mygroup$Primary_Tumor_Site=='Tumor',]
normal_sample <- mygroup[mygroup$Primary_Tumor_Site=='Normal',]

mytumor <-af_delete[,(colnames(af_delete) %in% tumor_sample$Sample_Name)]
mynormal <-af_delete[,(colnames(af_delete) %in% normal_sample$Sample_Name)]

set.seed(111)
trainindex <- sample(1:528,353,replace = FALSE)
train_tumor <- tumor[trainindex,]
validation_tumor <- tumor[-trainindex,]
set.seed(222)
trainindex <- sample(1:656,437,replace = FALSE)
train_normal <- normal[trainindex,]
validation_normal <- normal[-trainindex,]

train_clinical <- rbind(train_tumor,train_normal)
validation_clinical <- rbind(validation_tumor,validation_normal)


mytrain_tumor <- af_delete[,(colnames(af_delete) %in% train_tumor$Sample_Name)]
mytrain_normal <- af_delete[,(colnames(af_delete) %in% train_normal$Sample_Name)]
train_all <- cbind(mytrain_tumor,mytrain_normal)


myval_tumor <- af_delete[,(colnames(af_delete) %in% validation_tumor$Sample_Name)]
myval_normal <- af_delete[,(colnames(af_delete) %in% validation_normal$Sample_Name)]
val_all <- cbind(myval_tumor,myval_normal)


mytrainDMP <- champ.DMP(beta=train_all,pheno=train_clinical$Primary_Tumor_Site)

logfc_cutoff <- 0.4
df_dmp <- (mytrainDMP$Tumor_to_Normal[,1:9])
df_dmp$change <- ifelse(df_dmp$adj.P.Val < 0.05 & abs(df_dmp$logFC) > logfc_cutoff,ifelse(df_dmp$logFC > logfc_cutoff ,'UP','DOWN'),'NOT')

df_dmp_para <- (mytrainDMP$Tumor_to_Normal[,1:9])
df_dmp_para$change <- ifelse(df_dmp_para$adj.P.Val < 0.05 & abs(df_dmp_para$logFC) > logfc_cutoff,ifelse(df_dmp_para$logFC > logfc_cutoff ,'UP','DOWN'),'NOT')
id<- rownames(df_dmp[1])

#######################
pdf('DMProbe_volcano_plot.pdf',height=10,width=10)
this_tile <- paste0('Cutoff for logFC is 0.4&FDR<0.05','\nThe number of up gene is ',nrow(df_dmp[df_dmp$change =='UP',]),'\nThe number of down gene is ',nrow(df_dmp[df_dmp$change =='DOWN',]))
g <- ggplot(data=df_dmp,aes(x=logFC, y=-log10(adj.P.Val),color=change)) +ylim(0,60)+geom_point(alpha=0.4, size=1) +theme_set(theme_set(theme_bw(base_size=10)))+xlab("log2 fold change") + ylab("-log10 adj.p-value") +ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+scale_colour_manual(values = c('blue','black','red'))
g
dev.off()
######################DMP Heatmap 
library(pheatmap)
choose_cpg <- rownames(df_dmp[df_dmp$change != "NOT",])
plot_matrix <- train_all[choose_cpg,]
annotation_col <- data.frame(Sample=train_clinical$Primary_Tumor_Site)
rownames(annotation_col) <- colnames(plot_matrix)
ann_colors = list(Sample = c(Normal="green", Tumor="red"))
pdf('DMProbe_heatmap_plot.pdf',height=20,width=20)
p <- pheatmap(plot_matrix,show_colnames = F,annotation_col = annotation_col,border_color=NA,color = colorRampPalette(colors = c("blue","white","red"))(100),annotation_colors = ann_colors)
p
dev.off()
######################PCA 
library("FactoMineR")
library("factoextra")
dat <- t(plot_matrix)
group_list <-train_clinical$Primary_Tumor_Site
dat.pca <- PCA(dat, graph = FALSE)
pdf('DMProbe_PCA_plot.pdf',height=10,width=10)
c <- fviz_pca_ind(dat.pca,geom.ind = "point",col.ind = group_list,addEllipses = TRUE,legend.title = "Groups")
c
dev.off()
#############################挑选前adj.p.value 1000的cpg
DEG <- df_dmp[df_dmp$change %in% c("UP", "DOWN"),]
DEG <- DEG %>%
  arrange(adj.P.Val)

DEG <- DEG[1:1000,]
data <- train_all[rownames(train_all) %in% rownames(DEG) ,]
data <- data.frame(t(data))
data$sampleID <- rownames(data)
train_clinical$sampleID <-train_clinical$Sample_Name
dat <- merge(train_clinical,data,by = "sampleID")

valdata <- val_all[rownames(val_all) %in% rownames(DEG) ,]
valdata <- data.frame(t(valdata))
valdata$sampleID <- rownames(valdata)
validation_clinical$sampleID <-validation_clinical$Sample_Name
valdat <- merge(validation_clinical,valdata,by = "sampleID")
tbvaldat <- merge(hnsc_info,tb_val_data,by = "sampleID")

######################### LASSO
dat_dis <- dat[,c(5:length(colnames(dat)))]
dat$Primary_Tumor_Site <- factor(dat$Primary_Tumor_Site)
for (i in 1:500){
  model.lasso <- cv.glmnet(
    as.matrix(dat_dis), dat$Primary_Tumor_Site,
    family="binomial",
    alpha=1,
    nfolds=10
  )
  # plot(model.lasso)
  idx <- with(model.lasso, which.min(abs(glmnet.fit$lambda - lambda.1se)))
  coefs.lasso <- with(model.lasso$glmnet.fit, c(
    `(Intercept)`=unname(a0[idx]),
    beta[,idx]
  ))
  
  coefs.lasso <- coefs.lasso[abs(coefs.lasso) > .Machine$double.eps]
  if(i == 1){
    candidates <- names(coefs.lasso)
  }else{
    candidates <- c(candidates, names(coefs.lasso))
  }
  if (i %% 50 == 0){
    cat(paste(i, "predictions accomplished!\n", sep = " "))
  } 
}

# pick CpGs that was picked over 450 times in 500 iterations
predictors <- unique(candidates)
predictors <- predictors[which(predictors != "(Intercept)")]
cat(paste("there are", length(predictors), "predictors\n"))
hnsc_freq <- table(candidates)
predictors <- names(hnsc_freq)[which(hnsc_freq >= 450)[-1]]
tplassocpg <- names(hnsc_freq)[which(hnsc_freq >= 450)[-1]]
cat(paste(length(cpg), "predictors are picked\n"))


glm_mat <- dat_dis[,lassocpg]
glm_data <- as.data.frame(glm_mat)
colnames(glm_data) <- lassocpg
glm_data$group <- factor(dat$Primary_Tumor_Site)
f <- paste("group ~ ", lassocpg[1], sep = "")
if(length(lassocpg) != 1){
  for (j in 2:length(lassocpg)){
    f <- paste(f, " + ", lassocpg[j], sep = "")
  }
}
glm.model <- glm(f, family = "binomial", data = glm_data)
glm_mat <- valdata[,lassocpg]
glm_data <- as.data.frame(glm_mat)
colnames(glm_data) <- lassocpg
prediction <- predict(glm.model, glm_data, type="response")
perf <- performance(prediction(prediction, validation_clinical$Primary_Tumor_Site), "tpr", "fpr")
# ROC plot
plot(perf)
# AUC
print(performance(prediction(prediction, validation_clinical$Primary_Tumor_Site), "auc")@y.values)

dis_overlap <- dat[,overlap_cpg]
val_overlap <- valdat[,overlap_cpg]

dat$Primary_Tumor_Site <- factor(dat$Primary_Tumor_Site)
valdat$Primary_Tumor_Site <- factor(valdat$Primary_Tumor_Site)

dis_overlap$group <- factor(dat$Primary_Tumor_Site)
val_overlap$group <- factor(valdat$Primary_Tumor_Site)

glm.model <- glm(group ~ cg10364040, family = "binomial", data = dis_overlap)
prediction <- predict(glm.model, dis_overlap, type="response")
######Youden index 找cutoff
proc<-roc(val_overlap$group,prediction,ci=T)
plot(proc,
     print.auc=TRUE,  #显示AUC
     print.thres=TRUE,#显示最佳cutoff
     main="ROC",
     col="Red")
	 
pre <- c()
for (i in 1:378){
   pre <- append(pre,ifelse(prediction[[i]] >0.5,'Tumor','Normal'))
  }
caret::confusionMatrix(as.factor(pre),as.factor(dis_overlap$group),positive="Tumor")

#################################################################################################### RF
set.seed(1234)
core <- detectCores()
forkCL <- makeForkCluster(round(0.85*core, 0))
hnsc_RFBoot <- varSelRFBoot(as.matrix(dat_dis), factor(dat$Primary_Tumor_Site, ordered = F), c.sd = 1,
                            mtryFactor = 1, ntree = 5000, ntreeIterat = 2000,
                            vars.drop.frac = 0.3, bootnumber = 1000,
                            whole.range = FALSE,
                            recompute.var.imp = FALSE, usingCluster = TRUE,
                            srf = hnsc_RF, TheCluster = forkCL)
# extract results the variables chosen by RFBoot
randomForestcpg <- hnsc_RFBoot$all.data.vars
# the prediction compared to observation
hnsc_rf <- hnsc_RFBoot$all.data.randomForest
table(hnsc_rf$predicted, hnsc_rf$y)

load("tcga_tumor_paracancerous_Lasso_14marker_trainning.Rdata")
load("tcga_tumor_paracancerous_RandomForest_15marker_trainning.Rdata")

load("hnsc_352s_437s_filtered_imputed_Normed_overlap1139_clinical.Rdata")
load("TCGA_GEO_176s_219s_Validation_overlap1139_Validation_beta_clinical.Rdata")

lasso_13marker <- tplassocpg[-6]
rf_14marker <- tprandomforestcpg[-8]

lasso_train_matrix <- dat[,lasso_13marker]
lasso_train_matrix <- t(lasso_train_matrix)
colnames(lasso_train_matrix) <- dat$sampleID

annotation_col <- data.frame(Sample=dat$Primary_Tumor_Site)
rownames(annotation_col) <- colnames(lasso_train_matrix)

ann_colors = list(Sample = c(Blood="green", Tumor="red"))
pdf('lasso_13marker_Train_heatmap_plot.pdf',height=6,width=12)
p <- pheatmap(lasso_train_matrix,show_colnames = F,annotation_col = annotation_col,border_color=NA,color = colorRampPalette(colors = c("blue","white","red"))(100),annotation_colors = ann_colors,cluster_cols = FALSE,main="Training(n=789)", fontsize_row = 12)
p
dev.off()

lasso_val_matrix <- tbvaldat[,lasso_13marker]
lasso_val_matrix <- t(lasso_val_matrix)
colnames(lasso_val_matrix) <- tbvaldat$sampleID

annotation_col <- data.frame(Sample=tbvaldat$Primary_Tumor_Site)
rownames(annotation_col) <- colnames(lasso_val_matrix)

ann_colors = list(Sample = c(Blood="green", Tumor="red"))
pdf('lasso_13marker_Validation_heatmap_plot.pdf',height=6,width=12)
p <- pheatmap(lasso_val_matrix,show_colnames = F,annotation_col = annotation_col,border_color=NA,color = colorRampPalette(colors = c("blue","white","red"))(100),annotation_colors = ann_colors,cluster_cols = FALSE,main="Validation(n=395)", fontsize_row = 12)
p
dev.off()

#################################################
rf_train_matrix <- dat[,rf_14marker]
rf_train_matrix <- t(rf_train_matrix)
colnames(rf_train_matrix) <- dat$sampleID

annotation_col <- data.frame(Sample=dat$Primary_Tumor_Site)
rownames(annotation_col) <- colnames(rf_train_matrix)

ann_colors = list(Sample = c(Blood="green", Tumor="red"))
pdf('rf_14marker_Train_heatmap_plot.pdf',height=6,width=12)
p <- pheatmap(rf_train_matrix,show_colnames = F,annotation_col = annotation_col,border_color=NA,color = colorRampPalette(colors = c("blue","white","red"))(100),annotation_colors = ann_colors,cluster_cols = FALSE,main="Training(n=789)", fontsize_row = 12)
p
dev.off()

rf_val_matrix <- tbvaldat[,rf_14marker]
rf_val_matrix <- t(rf_val_matrix)
colnames(rf_val_matrix) <- tbvaldat$sampleID

annotation_col <- data.frame(Sample=tbvaldat$Primary_Tumor_Site)
rownames(annotation_col) <- colnames(rf_val_matrix)

ann_colors = list(Sample = c(Blood="green", Tumor="red"))
pdf('rf_14marker_Validation_heatmap_plot.pdf',height=6,width=12)
p <- pheatmap(rf_val_matrix,show_colnames = F,annotation_col = annotation_col,border_color=NA,color = colorRampPalette(colors = c("blue","white","red"))(100),annotation_colors = ann_colors,cluster_cols = FALSE,main="Validation(n=395)", fontsize_row = 12)
p
dev.off()

#################################################




