setwd("H:/feature_selection")
library(readxl)
library(glmnet)
library(plotmo)
library(MASS) 
library(caret)
##############DSTB & LTB
dstb_ltb<-read_excel("fold_change_DS-TB_vs_LTB.xlsx",sheet=3)
class(dstb_ltb$padj)= "numeric"
dstb_ltb_genes<- dstb_ltb[which(dstb_ltb$padj<0.05),]

####extracting normalized values of significant genes
norm <- read.delim("vst_normalized_count.txt" ,header=TRUE,sep = "\t",dec =".",fill = TRUE)
rownames(norm)<-norm$X
norm<- norm[,-1]
norm_dstb_ltb <- norm[,c(8:21)]
norm_dstb_ltb <- norm_dstb_ltb[which(rownames(norm_dstb_ltb)%in%dstb_ltb_genes$`Ensembl ID`),]
norm_dstb_ltb<-as.data.frame(t(norm_dstb_ltb))

###lasso feature selection using glmnet R package (0 - ltb, 1 - dstb)
set.seed(123)
norm_dstb_ltb$cohort <- c(rep(0, 8), rep(1, 6))
norm_dstb_ltb$cohort <- factor(norm_dstb_ltb$cohort)

##change data in matrix form
x <- as.matrix(norm_dstb_ltb[, -1082])
class=as.matrix(norm_dstb_ltb[,1082])
class(class)="numeric"

lasso.mod <- glmnet(x, y = class, alpha = 1, family = 'binomial',nlambda = 200)
plot(lasso.mod, xvar = 'lambda',label = TRUE)
plot_glmnet(lasso.mod, label=TRUE,main="Lasso")                          

#k-fold cross-validation to find the best lambda
cvfit <- cv.glmnet(x, class, alpha = 1,nfolds =3,type.measure="deviance",family = 'binomial')

#Best lambda value that minimizes mean squared error
best_lambda <- cvfit$lambda.min
#Display mean squared error by lambda value
plot(cvfit)
text(-4, 1.4, "N=8", col = "red", cex = 1.2) 
text(-4, 1.35, "Optimal λ = 0.0048", col = "red", cex = 1.2) 


lasso_best <- glmnet(x, class, alpha = 1, family = 'binomial', lambda = best_lambda)

library(coefplot)

coefplot(lasso_best, intercept = F, title="Coefficient Plot of lasso selected miRNAs",
         xlab="Coefficient Value", ylab="miRNA", cex=5)
#Display shrunk and eliminated coefficients
imp_feature<-coef(lasso_best)
# c(imp_feature@i+1)
lasso_imp_features <- c()
coeffecient <- c()
for (i in 1:length(imp_feature@i+1)) {
  lasso_imp_features[i] <- imp_feature@Dimnames[[1]][[c(imp_feature@i+1)[i]]]
  coeffecient[i] <- imp_feature[c(imp_feature@i+1)[i]]
}

lasso_sel_features<- data.frame(lasso_imp_features,coeffecient)[-1,]


#### randomForest feature selection using  randomForest R package  ####ran in manju's system
####https://www.geeksforgeeks.org/random-forest-approach-in-r-programming/
###https://www.geeksforgeeks.org/random-forest-approach-for-classification-in-r-programming/
##ref: https://www.youtube.com/watch?v=6EXPYzbfLCE
library("randomForest")
set.seed(123)
repeat_cv <- trainControl(method='repeatedcv', number=10)
# 
mtry <- tuneRF(x,factor(class), ntreeTry=500,
               plot=TRUE)
best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
print(mtry)
print(best.m) 
grid <- expand.grid(.mtry=c(32)) 
RF<-train(cohort ~ .,data = norm_dstb_ltb, method='rf', trControl=repeat_cv, tuneLength=15,tuneGrid=grid, metric='Accuracy')
classifier_RF<- RF$finalModel
var_imp <- varImp(classifier_RF, scale=FALSE)$importance
classifier_RF 
# # # Plotting model 
plot(classifier_RF)

# # Importance plot 
rf_imp_feature<-importance(classifier_RF) 
# # Variable importance plot  
varImpPlot(classifier_RF,main="Random Forest (Top 30 features)" )
rf_sel_features <- data.frame(rownames(rf_imp_feature),rf_imp_feature)


#####Xgboost using Xgboost R pacakge
##ref:https://www.projectpro.io/recipes/visualise-xgboost-feature-importance-r
##https://www.analyticsvidhya.com/blog/2016/01/xgboost-algorithm-easy-steps/
###https://www.hackerearth.com/practice/machine-learning/machine-learning-algorithms/beginners-tutorial-on-xgboost-parameter-tuning-r/tutorial/

library(caret) 
library(tidyverse)
library(xgboost)
set.seed(123)
xgb_train = xgb.DMatrix(data = data.matrix(x), label = class)
#parameters
params <- list(subsample=0.5, colsample_bytree=0.5)
xgb_cv<- xgb.cv( params = params, data = xgb_train, nrounds = 1000, nfold = 10,early_stopping_rounds = 1000)
xgb_cv$best_iteration
xgb <- xgboost(data = xgb_train, 
               nround=xgb_cv$best_iteration,nfolds=10,
               subsample = 0.5,
               colsample_bytree = 0.5) ##nround= no. of trees, subsample suggested by help page
summary(xgb)
# Compute feature importance matrix
importance_matrix = xgb.importance(colnames(train), model = xgb)
xgboost_sel_features<-as.data.frame(importance_matrix)
xgb.ggplot.importance(importance_matrix,main="Xgboost")


#########svm feature selection 
###ref:https://github.com/johncolby/SVM-RFE
####https://dataaspirant.com/support-vector-machine-classifier-implementation-r-caret-package/
###https://topepo.github.io/caret/using-your-own-model-in-train.html#illustrative-example-1-svms-with-laplacian-kernels
###https://copyprogramming.com/howto/using-an-svm-for-feature-selection#using-an-svm-for-feature-selection
###https://topepo.github.io/caret/recursive-feature-elimination.html#the-selectsize-function
##ref: https://www.appsloveworld.com/r/100/102/seed-object-for-reproducible-results-in-parallel-operation-in-caret?expand_article=1

library(doParallel) 
set.seed(123)
cl <- makeCluster(detectCores(), type='PSOCK') 
registerDoParallel(cl) 
set.seed(123)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]]<-sample.int(n=1081, 17)
#for the last model
seeds[[11]]<-sample.int(1081, 1)
rfecontrol <- rfeControl(functions = caretFuncs,
                         method="cv",number=10,seeds=seeds)
svm<-rfe(x=norm_dstb_ltb[,1:1081],y=as.numeric(class) ,
    sizes = c(1:10,50,80,100,300,400,500), rfeControl = rfecontrol ,
    method = "svmRadial")
svm_sel_features <- predictors(svm)
plot(svm, type = c("g", "o"),main="SVM")
data.frame(svm$optVariables)



########LDA      
###ref:https://topepo.github.io/caret/recursive-feature-elimination.html#rfe
###ref:https://www.kaggle.com/code/mubashirsultan/linear-discriminant-analysis
##https://rdrr.io/cran/caret/src/R/rfe.R
#lda <- lda(cohort ~ .,data = train) 

ibrary(doParallel) 
set.seed(123)
cl <- makeCluster(detectCores(), type='PSOCK') 
registerDoParallel(cl) 
set.seed(123)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]]<-sample.int(n=1081, 17)
seeds[[11]]<-sample.int(1081, 1)
rfecontrol <- rfeControl(functions = ldaFuncs,
                         method="cv",number=10,seeds=seeds)
lda<-rfe(x=norm_dstb_ltb[,1:1081],y=factor(class) ,
         sizes = c(1:10,50,80,100,300,400,500), rfeControl = rfecontrol ,
         method = NULL)
lda_sel_features <- predictors(lda)
plot(lda, type = c("g", "o"),main="LDA")
data.frame(lda$optVariables)


###plsda
###https://rdrr.io/github/xia-lab/MetaboAnalystR/man/PLSDA.CV.html
###https://mixomicsteam.github.io/mixOmics-Vignette/id_05.html
cl <- makeCluster(detectCores(), type='PSOCK') 
registerDoParallel(cl) 
library(mixOmics)
set.seed(123)
tune.splsda <- tune.splsda(x, factor(class), validation = 'Mfold', ncomp =10,
                           folds = 10, nrepeat = 10,dist = 'max.dist',test.keepX = c(1:10,50,100,200,300,400,500),progressBar = T,cpus=11)
plot(tune.splsda, sd = TRUE)
# # The optimal number of components according to our one-sided t-tests
ncomp <- tune.splsda$choice.ncomp$ncomp
# # Optimal number of variables to select
select.keepX <- tune.splsda$choice.keepX[1:ncomp]
select.keepX
splsda <- splsda(x, factor(class), ncomp = ncomp, keepX = select.keepX)
select.name <- selectVar(splsda, comp = 1)$name