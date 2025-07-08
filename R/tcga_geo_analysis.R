library(ChAMP)
library(ROCR)
library(randomForest)
library(caret)
library(pROC)
library(dplyr)
library(varSelRF)
library(glmnet)


setwd('./')
load('hnsc_346s_32s_filter_impute_NormBMIQ_Trainning.Rdata')
load('hnsc_175_16s_filter_impute_NormBMIQ_Validation.Rdata')
