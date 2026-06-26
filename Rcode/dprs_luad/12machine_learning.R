setwd('E:/端粒课题/Rproject/12_机器学习建模/')
#这里主要使用上面多组学的数据进行机器学习建模-----------------------------------
#基因选择每个亚型的前20的基因
rm(list = ls())
options(stringsAsFactors = F)
# 1. 根据TERT基因 无监督一致性聚类
#### 
library(GSVA)
library(ggplot2)
library(limma)
library(pheatmap)
library(reshape2)
library(mclust)
# 使用15个RNA进行建模-----------------------------------------------------------
load('../11_movics包分型/exp_cli_clusterGroup.rdata')
#去除生存时间为0的样本
cli3 <- cli3_unique[cli3_unique$time != 0, ]
exp3 <- exp[,colnames(exp)%in%cli3$SampleID]

mygene <- c('TERT','TNKS','KLF4','TERF2','TINF2',
            'MIR3682','MIR140','MIR3605','MIR324','MIR25',
            'MALAT1','MEG3','SNHG5','FOXD2-AS1','SNHG1')

# 整理机器学习输入矩阵----------------------------------------------------------
# 列为生存时间，生存状态和基因名，行为样本名
exp1 <- exp3[rownames(exp3) %in% mygene,]
exp1 <- t(exp1)
table(rownames(exp1)==cli3$SampleID)
# TRUE 
# 524
sig_gene_exp <- cbind(cli3[,c(10,9)],exp1)
rownames(sig_gene_exp) <- cli3$PatientD
colnames(sig_gene_exp)[1] <- 'OS.time'
colnames(sig_gene_exp)[2] <- 'OS'
sig_gene_exp$OS.time <- sig_gene_exp$OS.time / 365
write.csv(sig_gene_exp,file = 'sig_gene_exp.csv')

# 加载机器学习组合算法使用的R包-------------------------------------------------
library(ggplot2)
library(ggsci)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(mixOmics)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library("superpc")
library(devtools)
library("devtools")

##10种方法，127种组合----------------------------------------------------------- 
##################################
#### 准备工作 ####
##################################
# mm是包含训练集和验证集的list

#将总数据集分为训练集（占 70%）和测试集（占 30%）
iris <- sig_gene_exp
seed <- 1234
set.seed(seed)
colnames(iris)

ind <- sample(2,nrow(iris),replace = TRUE,prob = c(0.7,0.3))
train_data <- iris[ind==1,] #训练数据集
test_data <- iris[ind==2,] #测试数据集
#①这里修改了，加上了总体队列
mm <- list(Training_Dataset = train_data,
           Testing_Dataset = test_data,
           Total_Dataset = sig_gene_exp)

val_data_list <- mm

est_data <- mm$Training_Dataset
val_data_list <- mm

pre_var <- colnames(est_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',pre_var)]})
a <- val_dd_list[[1]]

rm(mm)
result <- c()

# 挑选出最佳rf_nodesize---------------------------------------------------------
for (rf_nodesize  in 1:30) {
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
               ntree = 3000,nodesize = rf_nodesize,##该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- 'RSF'
  cc$rfNodeSize = rf_nodesize
  result <- rbind(result,cc)
}
result <- result[order(result$Cindex,decreasing = T),]
result1 <- result

rf_nodesize <- 1

# 挑选出最佳的树----------------------------------------------------------------
result <- c()
for (bestntree  in seq(500,3000,200)) {
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
               ntree = bestntree,nodesize = rf_nodesize,##该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- 'RSF'
  cc$rfbestntree = bestntree
  result <- rbind(result,cc)
}

result2 <- result
result2 <- result2[order(result2$Cindex,decreasing = T),]
result2

rf_nodesize <- 1
rfbestntree <- 1100 #这里原来是2100 我调整为1100
seed <- 1234


# 机器学习组合算法开始----------------------------------------------------------
##################################
#### 1-1.RSF ####
##################################
result <- c()
set.seed(seed)
# fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
#              ntree = 1000,nodesize = rf_nodesize,##该值建议多调整
#              splitrule = 'logrank',
#              importance = T,
#              proximity = T,
#              forest = T,
#              seed = seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
             ntree = rfbestntree,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
# best <- which.min(fit$err.rate)
# set.seed(seed)
# fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
#              ntree = best,nodesize = rf_nodesize,##该值建议多调整  
#              splitrule = 'logrank',
#              importance = T,
#              proximity = T,
#              forest = T,
#              seed = seed)

rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result,cc)


##################################
#### 1-2.rsf+enet ####
##################################

vi <- data.frame(imp=vimp.rfsrc(fit)$importance)
vi$imp <- (vi$imp-min(vi$imp))/(max(vi$imp)-min(vi$imp))
vi$ID <- rownames(vi)

ggplot(vi,aes(imp,reorder(ID,imp)))+
  geom_bar(stat = 'identity',fill='#FF9933',color='black',width=0.7)+
  geom_vline(xintercept = 0.4,color='grey50',linetype=2)+
  labs(x='Relative importance by Random Forest',y=NULL)+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.x = element_text(size = 11,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.title = element_text(size=13,color='black'),
        legend.text = element_text(size=12,color='black'),
        legend.title = element_text(size=13,color='black'))+
  scale_y_discrete(expand = c(0.03,0.03))+
  scale_x_continuous(expand = c(0.01,0.01))

rid <- rownames(vi)[vi$imp>0.25] ## 0.4 can be adjust
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}


##################################
#### 1-3.RSF + Ridge ####
##################################
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 0)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Ridge')
result <- rbind(result, cc)


##################################
#### 1-4.RSF + Lasso ####
##################################
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Lasso')
result <- rbind(result, cc)


##################################
#### 1-5.rsf+CoxBoost ####
##################################
library(snowfall)
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
plot(fit)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result,cc)


##################################
#### 1-6.rsf+plsRcox ####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=3)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + plsRcox')
result <- rbind(result,cc)

##################################
#### 1-7.rsf+superpc ####
##################################

data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + SuperPC')
result <- rbind(result,cc)


##################################
#### 1-8.rsf+gbm ####
##################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
plot(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + GBM')
result <- rbind(result,cc)

##################################
#### 1-9.rsf+survivalsvm ####
##################################

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + survival-SVM')
result <- rbind(result,cc)


##################################
#### 1-10.rsf+stepcox ####
##################################

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + StepCox','[',direction,']')
  result <- rbind(result,cc)
}

##################################
#### 2-1.Enet ####
##################################

x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time,est_dd$OS))

for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}

##################################
#### 3-1.Lasso####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = 'cox', alpha = 1)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso')
result <- rbind(result, cc)

##################################
#### 3-2.Lasso+RSF####
##################################

set.seed(seed)
fit = cv.glmnet(x1, x2,family = "cox",alpha=1,nfolds = 10)
coef.min = coef(fit, s = "lambda.min") 
active.min = which(as.numeric(coef.min)!=0)
rid <- colnames(x1)[active.min]

est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 3000,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = best,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'Lasso + RSF'
result <- rbind(result,cc)


##################################
#### 3-3.Lasso+enet ####
##################################

est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}


##################################
#### 3-4.Lasso + Ridge ####
##################################
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 0)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'Ridge')
result <- rbind(result, cc)



##################################
#### 3-5.Lasso+CoxBoost ####
##################################

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result,cc)

##################################
#### 3-6.Lasso+plsRcox ####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + plsRcox')
result <- rbind(result,cc)

##################################
#### 3-7.Lasso+superpc ####
##################################

data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + SuperPC')
result <- rbind(result,cc)

##################################
#### 3-8.Lasso+gbm ####
##################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + GBM')
result <- rbind(result,cc)

##################################
#### 3-9.Lasso+survivalsvm ####
##################################

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + survival-SVM')
result <- rbind(result,cc)


##################################
#### 3-10.Lasso+StepCox ####
##################################

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + StepCox','[',direction,']')
  result <- rbind(result,cc)
}


# 此部分是将StepCox算法和其他算法组合放在一个for循环执行的，直接运行，如果报错则执行下面注释掉的分开运行的StepCox
##################################
#### 3.StepCox ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd), direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

for (direction in c("both", "backward", "forward")) {
  ##当基因名有-时，需要将基因名中的-改成.(因为My.stepwise.coxph函数识别不了-)
  colnames(est_dd) = gsub("-","\\.",colnames(est_dd))
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  rid = gsub("\\.","-",rid)
  colnames(est_dd) = gsub("\\.","-",colnames(est_dd))
  
  est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
  val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                              trace=TRUE, start.penalty = 500, parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                        maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
  fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                  stepno = cv.res$optimal.step, penalty = pen$penalty)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
  result <- rbind(result, cc)
  
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
    rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
    cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
    result <- rbind(result, cc)
  }
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
  cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold=10, #例文描述：10-fold cross-validation
                  family = "cox", alpha = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
  result <- rbind(result, cc)
  set.seed(seed)
  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time,
                 event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "cox", alpha = 0)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
  result <- rbind(result, cc)
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
               ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
  result <- rbind(result, cc)
  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
               censoring.status = est_dd2$OS,
               featurenames = colnames(est_dd2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                       n.fold = 10,
                       n.components = 3,
                       min.features = 5,
                       max.features = nrow(data$x),
                       compute.fullcv = TRUE,
                       compute.preval = TRUE)
  rs <- lapply(val_dd_list2, function(w){
    test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2], RS = rr)
    return(rr2)
  })
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
  result <- rbind(result, cc)
  fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
  result <- rbind(result, cc)
}



##################################
#### 4-1.CoxBoost ####
##################################

set.seed(seed)
colnames(est_dd)[1] <- "OS.time"
pen <- optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)

##################################
#### 4-2.CoxBoost+Enet ####
##################################

rid <- names(coef(fit)[which(coef(fit)!=0)])
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}

##################################
#### 4-3.CoxBoost + Ridge####
##################################

set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 0)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Ridge')
result <- rbind(result, cc)


##################################
#### 4-4.CoxBoost + lasso####
##################################

set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Lasso')
result <- rbind(result, cc)


##################################
#### 4-4.CoxBoost+RSF ####
##################################

set.seed(seed)
# fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
#              ntree = 1000,nodesize = rf_nodesize,##该值建议多调整
#              splitrule = 'logrank',
#              importance = T,
#              proximity = T,
#              forest = T,
#              seed = seed)
# best <- which.min(fit$err.rate)
# set.seed(seed)
# fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
#              ntree = best,nodesize = rf_nodesize,##该值建议多调整
#              splitrule = 'logrank',
#              importance = T,
#              proximity = T,
#              forest = T,
#              seed = seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 330,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'CoxBoost + RSF'
result <- rbind(result,cc)


##################################
#### 4-5.CoxBoost+plsRcox ####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + plsRcox')
result <- rbind(result,cc)



##################################
#### 4-6.CoxBoost+superpc ####
##################################

data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + SuperPC')
result <- rbind(result,cc)

##################################
#### 4-7.CoxBoost+gbm ####
##################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + GBM')
result <- rbind(result,cc)

##################################
#### 4-8.CoxBoost+survivalsvm ####
##################################

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + survival-SVM')
result <- rbind(result,cc)

##################################
#### 4-9.CoxBoost+stepcox ####
##################################

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + StepCox','[',direction,']')
  result <- rbind(result,cc)
}


##################################
#### 5.plsRcox####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd[,pre_var],time=est_dd$OS.time,status=est_dd$OS),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(est_dd[,pre_var],time=est_dd$OS.time,event=est_dd$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)

##################################
#### 6.superpc####
##################################

data <- list(x=t(est_dd[,-c(1,2)]),y=est_dd$OS.time,censoring.status=est_dd$OS,featurenames=colnames(est_dd)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=3, 
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result,cc)

##################################
#### 7.GBM ####
##################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = 1000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)

##################################
#### 8.survivalsvm ####
##################################

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd, gamma.mu = 1)

rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('survival-SVM')
result <- rbind(result,cc)

####结束####--------------------------------------------------------------------
result2 <- result 

result2%>%
  ggplot(aes(Cindex,reorder(Model,Cindex)))+
  geom_bar(width = 0.7,stat = 'summary',fun='mean',fill='orange2')+
  theme_classic()+
  labs(y=NULL)


dd <- result2%>%  group_by(Model)
dd <- as.data.frame(dd)
dd1 <- aggregate(dd$Cindex, by=list(type=dd$Model),mean)

library(ggbreak)
dd1$x <- as.numeric(format(dd1$x, digits = 3, nsmall = 3))
dd1 <- dd1[order(dd1$x,decreasing = T),]

ggplot(dd1,aes(dd1$x,reorder(dd1$type,dd1$x)))+
  geom_bar(width=0.7,stat = 'identity',fill='orange')+
  scale_x_break(c(0.05,0.53),scales = 20)

library(tidyr)
dd2 <- pivot_wider(result2,names_from = 'ID',values_from = 'Cindex')%>%as.data.frame()
# dd2$mean_Cindex = rowMeans(select(dd2,c(TCGA,GSE31210,GSE72094,GSE26939)))
# dd2$mean_TCGA_GSE31210 = rowMeans(select(dd2,c(TCGA,GSE31210)))

# 仅绘制GEO验证集的C指数热图
#②这里修改了，把2：3改成了2：4
dt <- dd2[, 2:4]
rownames(dt) <- dd2$Model
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

##热图绘制----------------------------------------------------------------------
Cindex_mat=dt
avg_Cindex <- apply(Cindex_mat, 1, mean)           # 计算每种算法在所有队列中平均C-index
avg_Cindex <- sort(avg_Cindex, decreasing = T)     # 对各算法C-index由高到低排序
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]      # 对C-index矩阵排序

avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数

col_ha = columnAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                             gp = gpar(fill = "steelblue", col = NA),
                                             add_numbers = T, numbers_offset = unit(-10, "mm"),
                                             axis_param = list("labels_rot" = 0),
                                             numbers_gp = gpar(fontsize = 9, col = "white"),
                                             width = unit(3, "cm")),
                          show_annotation_name = F)

#CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # 设置队列颜色
CohortCol <- rainbow( ncol(Cindex_mat))
names(CohortCol) <- colnames(Cindex_mat)
row_ha = rowAnnotation("Cohort" = colnames(Cindex_mat),
                       col = list("Cohort" = CohortCol),
                       show_annotation_name = F)



Cindex_mat = t(Cindex_mat)
cellwidth = 1
cellheight=0.3
pdf("cindex111.pdf",14,20)
Heatmap(as.matrix(Cindex_mat), name = "C-index",
        #left_annotation = row_ha, 
        #top_annotation = col_ha,
        # col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 黄绿配色
        col = c("#4195C1", "#FFFFFF", "#CB5746"), # 红蓝配色
        rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
        cluster_columns = FALSE, cluster_rows = FALSE, # 不进行聚类，无意义
        show_column_names = TRUE, 
        show_row_names = TRUE,
        row_names_side = "left",
        width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
        #height = unit(cellheight * nrow(Cindex_mat), "cm"),
        column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
        column_title = NULL,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                    x, y, gp = gpar(fontsize = 5))
        }
)

dev.off()

#### 筛选出RSF------------------------------------------------------------------
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
             ntree = rfbestntree,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
plot(fit)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
data <- list(x = t(est_dd[, -c(1, 2)]), y = est_dd$OS.time,
             censoring.status = est_dd$OS,
             featurenames = colnames(est_dd)[-c(1,2)])

save(fit, val_dd_list,file = 'my_mlresult_RSF15.rdata')
#输出重要基因列表---------------------------------------------------------------
importantGenes <- data$featurenames

aa.gene.1 <- c()
for (i in 1:length(importantGenes)) {
  aa.gene <- importantGenes[i] 
  aa.gene.1 <- paste(aa.gene,aa.gene.1,sep = "+")
}
aa.gene.1
est_dd3 <- est_dd
colnames(est_dd3)[3] <- 'FOXD2.AS1'
res.cox <- coxph(Surv(OS.time,OS)~TNKS+TINF2+TERT+TERF2+SNHG5+SNHG1+MIR3682+MIR3605+MIR324+MIR25+MIR140+MEG3+MALAT1+KLF4+FOXD2.AS1,
                 est_dd3)
summary(res.cox)
res.cox <- as.data.frame(res.cox$coefficients)

res.cox$gene <- rownames(res.cox)
colnames(res.cox)[1] <- "coefficients"

res.cox

##绘制基因系数图----------------------------------------------------------------

library(ggplot2)

# 排序数据框
df <- as.data.frame(res.cox)
# 排序数据框
df <- df[order(df$coefficients, decreasing = TRUE), ]

# 绘制图形
ggplot(df, aes(x = reorder(gene, coefficients), y = coefficients)) +
  geom_segment(aes(x = reorder(gene, coefficients), xend = reorder(gene, coefficients), 
                   y = 0, yend = coefficients), color = "skyblue", size = 1) +  # 先绘制垂直线段
  geom_point(color = "orange", size = 3) +  # 后绘制红色的点
  geom_hline(yintercept = 0, linetype = "dashed") +  # 绘制x轴的0位置的虚线
  coord_flip() +  # 横向条形图
  labs(title = "",
       x = "",
       y = "Coefficient") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),  # 添加黑色边框
    axis.ticks = element_line(color = "black", size = 0.5)  # 添加坐标轴刻度
  )
write.table(res.cox,"res.cox.coefficients.0630.txt",quote = F,sep = "\t",row.names = F)
#手动保存8*6--------------------------------------------------------------------

# 我属实有些蒙圈了 但我觉得这是将系数和基因名相乘作为模型的公式-----------------
res.cox$risk <- paste0("(",res.cox$coef,")","*","est_dd3$",res.cox$gene)

aa.gene.2 <- c()
for (i in 1:length(res.cox$gene)) {
  aa.gene <- res.cox$risk[i] 
  aa.gene.2 <- paste0(aa.gene,"+",aa.gene.2)
}
aa.gene.2 

est_dd3$riskscore <- (0.369460901946597)*est_dd3$FOXD2.AS1+(-0.328247928957762)*est_dd3$KLF4+(0.0631004437118076)*est_dd3$MALAT1+(0.24138664181875)*est_dd3$MEG3+(-0.168615135604993)*est_dd3$MIR140+(-0.0523561558725742)*est_dd3$MIR25+(-0.194900224040818)*est_dd3$MIR324+(0.202295290454816)*est_dd3$MIR3605+(0.14456861112338)*est_dd3$MIR3682+(0.123634166937939)*est_dd3$SNHG1+(-0.174220644981062)*est_dd3$SNHG5+(0.674072978757864)*est_dd3$TERF2+(0.121916046405159)*est_dd3$TERT+(-0.440793835176886)*est_dd3$TINF2+(-0.204409953502085)*est_dd3$TNKS
est_dd3$group <- ifelse(est_dd3$riskscore > median(est_dd3$riskscore),"high_risk","low_risk" )

median(est_dd3$riskscore)

head(est_dd3)
sel.exp <- est_dd3
head(sel.exp)
write.table(sel.exp,"riskscore-survival.groups.0630.txt",quote = F,sep = "\t",row.names = F)
write.csv(sel.exp,"riskscore-survival.groups.0630.csv")

#绘制训练集KM和ROC曲线-----------------------------------------------------------------
require(rms) #加载rms包
require(survival) #加载survival包
require(survminer) #加载survminer包
#用surv_fit函数对生存数据对象拟合生存函数，创建K-M曲线
graphGRADE<-surv_fit(Surv(OS.time,OS)~group,data=sel.exp) 
pdf("riskscore-survival.genes-plot.0630.pdf",5,5.5)

ggsurvplot(graphGRADE,
           data = sel.exp,
           risk.table = TRUE,
           pval = T)
dev.off()
#绘制累积风险函数。变为累积风险函数后，不能添加中位生存线，y轴的标题一定要正确，若是更改务必注意

Survival_ROC_input<-sel.exp[,c('OS','OS.time','riskscore','group')]
#Survival_ROC_input$OS.time<-Survival_ROC_input$OS.time/365
colnames(Survival_ROC_input) <- c('status','time','riskscore','group')
library(timeROC)

time_ROC_input<-Survival_ROC_input
time_ROC<-timeROC(T=time_ROC_input$time, 
                  delta=time_ROC_input$status, 
                  marker=time_ROC_input$riskscore,
                  cause=1, 
                  weighting = "marginal", 
                  times = c(1,3,5), 
                  ROC=TRUE,
                  iid=TRUE
)
time_ROC 

summary(time_ROC) 
time_ROC$AUC
time_ROC.res<-data.frame(TP_3year=time_ROC$TP[,1],
                         FP_3year=time_ROC$FP[,1],  
                         TP_5year=time_ROC$TP[,2],  
                         FP_5year=time_ROC$FP[,2], 
                         TP_10year=time_ROC$TP[,3],
                         FP_10year=time_ROC$FP[,3]) 
time_ROC$AUC 
p <- ggplot() +
  geom_line(data = time_ROC.res, aes(x = FP_3year, y = TP_3year), size = 1, color = "#EA6433") +
  geom_line(data = time_ROC.res, aes(x = FP_5year, y = TP_5year), size = 1, color = "#16A5D9") +
  geom_line(data = time_ROC.res, aes(x = FP_10year, y = TP_10year), size = 1, color = "#FFC21A") +
  geom_line(aes(x = c(0, 1), y = c(0, 1)), color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text", x = 0.75, y = 0.25, size = 5,
           label = paste0("AUC at 1 year = ", round(time_ROC$AUC[[1]], 3)), color = "#EA6433") +
  annotate("text", x = 0.75, y = 0.15, size = 5,
           label = paste0("AUC at 3 years = ", round(time_ROC$AUC[[2]], 3)), color = "#16A5D9") +
  annotate("text", x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", round(time_ROC$AUC[[3]], 3)), color = "#FFC21A") +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次网格线
  )

print(p)

ggsave(p, filename="riskscore.ROC.pdf", width=5,
       height=5.5)  

write.csv(sel.exp,"riskscore.ROC.csv")

##绘制验证集KM生存曲线和ROC曲线-------------------------------------------------

Testing_Dataset <- val_data_list[["Testing_Dataset"]]
os_data <- Testing_Dataset[,c(1:2)]
Testing_Dataset <- Testing_Dataset[,colnames(Testing_Dataset)%in%importantGenes]
Testing_Dataset <- cbind(os_data,Testing_Dataset)

Testing_Dataset$riskscore <- (0.369460901946597)*Testing_Dataset$`FOXD2-AS1`+(-0.328247928957762)*Testing_Dataset$KLF4+(0.0631004437118076)*Testing_Dataset$MALAT1+(0.24138664181875)*Testing_Dataset$MEG3+(-0.168615135604993)*Testing_Dataset$MIR140+(-0.0523561558725742)*Testing_Dataset$MIR25+(-0.194900224040818)*Testing_Dataset$MIR324+(0.202295290454816)*Testing_Dataset$MIR3605+(0.14456861112338)*Testing_Dataset$MIR3682+(0.123634166937939)*Testing_Dataset$SNHG1+(-0.174220644981062)*Testing_Dataset$SNHG5+(0.674072978757864)*Testing_Dataset$TERF2+(0.121916046405159)*Testing_Dataset$TERT+(-0.440793835176886)*Testing_Dataset$TINF2+(-0.204409953502085)*Testing_Dataset$TNKS
Testing_Dataset$group <- ifelse(Testing_Dataset$riskscore > median(Testing_Dataset$riskscore),"high_risk","low_risk" )

median(Testing_Dataset$riskscore)

head(Testing_Dataset)

write.table(Testing_Dataset,"riskscore-survival.groups.test.0630.txt",quote = F,sep = "\t",row.names = F)
write.csv(Testing_Dataset,"riskscore-survival.groups.test.0630.csv")

####KM####
#用surv_fit函数对生存数据对象拟合生存函数，创建K-M曲线
graphGRADE<-surv_fit(Surv(OS.time,OS)~group,data=Testing_Dataset)

pdf("riskscore-survival.genes-plot.test.0630.pdf",5,5.5)
ggsurvplot(graphGRADE,
           data = Testing_Dataset,
           risk.table = TRUE,
           pval = T)
dev.off()

####ROC####
Survival_ROC_input<-Testing_Dataset[,c('OS','OS.time','riskscore','group')]
#Survival_ROC_input$OS.time<-Survival_ROC_input$OS.time/365
colnames(Survival_ROC_input) <- c('status','time','riskscore','group')
library(timeROC)

time_ROC_input<-Survival_ROC_input
time_ROC<-timeROC(T=time_ROC_input$time, 
                  delta=time_ROC_input$status, 
                  marker=time_ROC_input$riskscore,
                  cause=1, 
                  weighting = "marginal", 
                  times = c(1,3,5), 
                  ROC=TRUE,
                  iid=TRUE
)
time_ROC 

summary(time_ROC) 
time_ROC$AUC
time_ROC.res<-data.frame(TP_3year=time_ROC$TP[,1],
                         FP_3year=time_ROC$FP[,1],  
                         TP_5year=time_ROC$TP[,2],  
                         FP_5year=time_ROC$FP[,2], 
                         TP_10year=time_ROC$TP[,3],
                         FP_10year=time_ROC$FP[,3]) 
time_ROC$AUC 
p <- ggplot() +
  geom_line(data = time_ROC.res, aes(x = FP_3year, y = TP_3year), size = 1, color = "#EA6433") +
  geom_line(data = time_ROC.res, aes(x = FP_5year, y = TP_5year), size = 1, color = "#16A5D9") +
  geom_line(data = time_ROC.res, aes(x = FP_10year, y = TP_10year), size = 1, color = "#FFC21A") +
  geom_line(aes(x = c(0, 1), y = c(0, 1)), color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text", x = 0.75, y = 0.25, size = 5,
           label = paste0("AUC at 1 year = ", round(time_ROC$AUC[[1]], 3)), color = "#EA6433") +
  annotate("text", x = 0.75, y = 0.15, size = 5,
           label = paste0("AUC at 3 years = ", round(time_ROC$AUC[[2]], 3)), color = "#16A5D9") +
  annotate("text", x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", round(time_ROC$AUC[[3]], 3)), color = "#FFC21A") +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次网格线
  )

print(p)

ggsave(p, filename="riskscore.ROC.test.pdf", width=5,
       height=5.5)  

##绘制总体队列KM生存曲线和ROC曲线-----------------------------------------------
Total_Dataset <- sig_gene_exp
os_data <- Total_Dataset[,c(1:2)]
Total_Dataset <- Total_Dataset[,colnames(Total_Dataset)%in%importantGenes]
Total_Dataset <- cbind(os_data,Total_Dataset)

Total_Dataset$riskscore <- (0.369460901946597)*Total_Dataset$`FOXD2-AS1`+(-0.328247928957762)*Total_Dataset$KLF4+(0.0631004437118076)*Total_Dataset$MALAT1+(0.24138664181875)*Total_Dataset$MEG3+(-0.168615135604993)*Total_Dataset$MIR140+(-0.0523561558725742)*Total_Dataset$MIR25+(-0.194900224040818)*Total_Dataset$MIR324+(0.202295290454816)*Total_Dataset$MIR3605+(0.14456861112338)*Total_Dataset$MIR3682+(0.123634166937939)*Total_Dataset$SNHG1+(-0.174220644981062)*Total_Dataset$SNHG5+(0.674072978757864)*Total_Dataset$TERF2+(0.121916046405159)*Total_Dataset$TERT+(-0.440793835176886)*Total_Dataset$TINF2+(-0.204409953502085)*Total_Dataset$TNKS
Total_Dataset$group <- ifelse(Total_Dataset$riskscore > median(Total_Dataset$riskscore),"high_risk","low_risk" )

median(Total_Dataset$riskscore)

head(Total_Dataset)

write.table(Total_Dataset,"riskscore-survival.groups.total.0630.txt",quote = F,sep = "\t",row.names = F)
write.csv(Total_Dataset,"riskscore-survival.groups.total.0630.csv")
####KM####
#用surv_fit函数对生存数据对象拟合生存函数，创建K-M曲线
graphGRADE<-surv_fit(Surv(OS.time,OS)~group,data=Total_Dataset) 

pdf("riskscore-survival.genes-plot.total.0630.pdf",5,5.5)
ggsurvplot(graphGRADE,
           data = Total_Dataset,
           risk.table = TRUE,
           pval = T)
dev.off()

####ROC####
Survival_ROC_input<-Total_Dataset[,c('OS','OS.time','riskscore','group')]
#Survival_ROC_input$OS.time<-Survival_ROC_input$OS.time/365
colnames(Survival_ROC_input) <- c('status','time','riskscore','group')
library(timeROC)

time_ROC_input<-Survival_ROC_input
time_ROC<-timeROC(T=time_ROC_input$time, 
                  delta=time_ROC_input$status, 
                  marker=time_ROC_input$riskscore,
                  cause=1, 
                  weighting = "marginal", 
                  times = c(1,3,5), 
                  ROC=TRUE,
                  iid=TRUE
)
time_ROC 

summary(time_ROC) 
time_ROC$AUC
time_ROC.res<-data.frame(TP_3year=time_ROC$TP[,1],
                         FP_3year=time_ROC$FP[,1],  
                         TP_5year=time_ROC$TP[,2],  
                         FP_5year=time_ROC$FP[,2], 
                         TP_10year=time_ROC$TP[,3],
                         FP_10year=time_ROC$FP[,3]) 
time_ROC$AUC 
p <- ggplot() +
  geom_line(data = time_ROC.res, aes(x = FP_3year, y = TP_3year), size = 1, color = "#EA6433") +
  geom_line(data = time_ROC.res, aes(x = FP_5year, y = TP_5year), size = 1, color = "#16A5D9") +
  geom_line(data = time_ROC.res, aes(x = FP_10year, y = TP_10year), size = 1, color = "#FFC21A") +
  geom_line(aes(x = c(0, 1), y = c(0, 1)), color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text", x = 0.75, y = 0.25, size = 5,
           label = paste0("AUC at 1 year = ", round(time_ROC$AUC[[1]], 3)), color = "#EA6433") +
  annotate("text", x = 0.75, y = 0.15, size = 5,
           label = paste0("AUC at 3 years = ", round(time_ROC$AUC[[2]], 3)), color = "#16A5D9") +
  annotate("text", x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", round(time_ROC$AUC[[3]], 3)), color = "#FFC21A") +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次网格线
  )

print(p)

ggsave(p, filename="riskscore.ROC.total.pdf", width=5,
       height=5.5)  
#保存机器学习得到的分组信息
save(est_dd3,cli3,file = 'ml_result.rdata')
