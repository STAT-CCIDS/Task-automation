

##############################################################################
#### Feature Selection
#### -------------------------------------------------------------------------
# This is the function for feature selection.
# It returns the selected features and coefficients using Lasso, Stability, Elastic net, and Random forest
# which is used for exploring significant features especially when n < p.
#
# To get proper result, there are few things that you should follow and it is written below.
#
# 1. The format for input data must be a data frame.
#
# 2. The type of variable must be clarified.
#
# 3. It returns file via "xlsx" format. 
#
# 4. Parameter description
#    -data : data (data.frame format) which contains variables
#    -out.type : outcome type, "binary" or "survival"
#    -out1 : column name of outcome variable (time variable for survival outcome)
#    -out2 : column name of event variable (only for survival outcome)
#    -clinical : clinical feature name
#    -feature1: first feature name
#    -feature2: last feature name
#    -method: name of selection method.
#             Only allowable choice at present is "LASSO", "stability", "elastic net", or "random forest".
#             If you don't specify it, the results will be returned using all methods.
#    -nfold: number of folds (dafult 10)
#    -seed : seed number.
#            If you don't specify it, it will be automatically assigned to your current date (YYYYMMDD).
#    -B : number of replicates (default 100)
#    -filename : name of the xlsx file
##############################################################################

featureSel <- function(data, out.type, out1, out2=NA, clinical=NA, feature1, feature2, method=NA, nfold=10, seed=NA, B=100, filename="Output") {
  tStart <- Sys.time()
  
  ## NA's are not allowed
  pMiss <- function(x) {sum(is.na(x))/length(x)*100}
  if(length(which(apply(data,2,pMiss)!=0))!=0) {stop("NA's are not allowed")}
  
  
  ## Install packages
  list.of.packages <- c("survival", 
                        "caret",
                        "glmnet",
                        "glmnetUtils",
                        "stabs",
                        "rms",
                        "randomForestSRC",
                        "ggRandomForests",
                        "randomForest",
                        "progress",
                        "xlsx",
                        "doParallel")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {install.packages(new.packages)}
  
  library(survival, quietly=TRUE)
  library(caret, quietly=TRUE)
  library(glmnet, quietly=TRUE)
  library(glmnetUtils, quietly=TRUE)
  library(stabs, quietly=TRUE)
  library(rms, quietly=TRUE)
  library(randomForestSRC, quietly=TRUE)
  library(ggRandomForests, quietly=TRUE)
  library(randomForest, quietly=TRUE)
  library(progress, quietly=TRUE)
  library(xlsx, quietly=TRUE)
  library(doParallel, quietly=TRUE)
  registerDoParallel(detectCores())
  detectCores() #the number of logical CPU cores on Window
  options(rf.cores=detectCores(), mc.cores=detectCores()) #Cores for parallel processing
  
  
  ## Functions for feature selection
  ##############################################################################
  #### featureLSE
  #### -------------------------------------------------------------------------
  # rs.id: resampled id
  # data.b: repeated data
  # features.b: repeated features 
  # type.measure: "deviance" for Cox, "auc" for two-class logistic
  ##############################################################################
  featureLSE <- function(data, out.type, out1, out2=NA, seed, B, method, nfold) {
    
    data.b <- list() 
    features.b <- list()
    rs.id <- matrix(0, nrow=nrow(data), ncol=B)
    fit.cv.b <- list()
    min.coef.b <- list()
    
    set.seed(seed)
    pb <- txtProgressBar(min=0, max=(B-1), style=3)
    features <- model.matrix(~.,data[,-which(names(data) %in% c(out1,out2))])[,-1]
    for(i in 1:B) {
      rs.id[,i] <- sample(1:nrow(features), nrow(features), replace=T)
      data.b[[i]] <- data[rs.id[,i],]
      features.b[[i]] <- features[rs.id[,i],]
      
      ## Lasso
      if(method=="Lasso") {
        if(out.type=="survival") {
          fit.cv.b[[i]] <- cv.glmnet(features.b[[i]], Surv(data.b[[i]][,out1], data.b[[i]][,out2]), family="cox", type.measure="deviance", nfolds=nfold)
        }
        if(out.type=="binary") {
          fit.cv.b[[i]] <- cv.glmnet(features.b[[i]], data.b[[i]][,out1], family="binomial", type.measure="auc", nfolds=nfold)
        }
        min.coef.b[[i]] <- coef(fit.cv.b[[i]], s="lambda.min")
        min.coef <- matrix(NA, ncol=B, nrow=nrow(min.coef.b[[1]]), dimnames=list(rownames(min.coef.b[[1]]), 1:B))
      }
      
      ## Stability
      if(method=="Stability") {
        if(out.type=="survival") {
          return("Stability selection is not available for the survival outcome")
        }
        if(out.type=="binary") {
          fit.cv.b[[i]] <- stabsel(features.b[[i]], data.b[[i]][,out1], fitfun=glmnet.lasso, PFER=1, cutoff=0.75)
          min.coef.b[[i]] <- fit.cv.b[[i]]$max
          min.coef <- matrix(NA, ncol=B, nrow=length(min.coef.b[[1]]), dimnames=list(names(min.coef.b[[1]]), 1:B))
        }
      }
      
      ## Elastic net
      if(method=="Elastic net") {
        if(out.type=="survival") {
          elnet.cv <- list()
          alpha <- seq(0, 1, len = 11)^3
          elnet.dev <- rep(0, length(alpha))
          elnet.cv[[i]] <- cva.glmnet(features.b[[i]], Surv(data.b[[i]][,out1], data.b[[i]][,out2]), family="cox", type.measure="deviance", nfolds=nfold)
          for(k in 1:11) { elnet.dev[k] <- min(elnet.cv[[i]]$modlist[[k]]$cvm) }
          fit.cv.b[[i]] <- glmnetUtils::glmnet(features.b[[i]], Surv(data.b[[i]][,out1], data.b[[i]][,out2]), family="cox",
                                               alpha=alpha[which.min(elnet.dev)], lambda=elnet.cv[[i]]$modlist[[which.min(elnet.dev)]]$lambda.min)
        }
        if(out.type=="binary") {
          grid <- expand.grid(.alpha=seq(0.1,0.9,by=0.1), .lambda=seq(0,1,by=0.01))
          elnet.cv <- train(features.b[[i]], data.b[[i]][,out1], method="glmnet", trControl=trainControl(method="repeatedcv",number=nfold), tuneGrid=grid)
          alpha <- elnet.cv$bestTune[[1]]
          lambda <- elnet.cv$bestTune[[2]]
          fit.cv.b[[i]] <- glmnet(features.b[[i]], data.b[[i]][,out1], family="binomial", alpha=alpha, lambda=lambda)
        }
        min.coef.b[[i]] <- fit.cv.b[[i]]$beta
        min.coef <- matrix(NA, ncol=B, nrow=nrow(min.coef.b[[1]]), dimnames=list(rownames(min.coef.b[[1]]), 1:B))
      }
      
      Sys.sleep(0.01)
      setTxtProgressBar(pb,i)
    }
    close(pb)
    
    min.coef <- as.data.frame(min.coef)
    for(i in 1:B) {min.coef[,i] <- min.coef.b[[i]]}
    
    min.coef <- as.matrix(min.coef)
    m <- c(0); abs_m <- c(0); se <- c(0); rsd <- c(0); abs_rsd <- c(0)
    for(i in 1:nrow(min.coef)) {
      m[i] <- mean(as.numeric(min.coef[i,]))
      abs_m[i] <- abs(m[i])
      se[i] <- sd(as.numeric(min.coef[i,]))
      rsd[i] <- (se[i]/m[i])
      abs_rsd[i] <- (se[i]/abs_m[i])
    }
    
    return(cbind(min.coef, mean=m, abs_m=abs_m, se=se, rsd=rsd, abs_rsd=abs_rsd))
  }
  
  
  ##############################################################################
  #### Random Survival Forests
  #### -------------------------------------------------------------------------
  # ntree: total number of trees
  # mtry: number of variables entering in each division as candidates
  #       (sqrt(p) in classification, survival / p/3 in regression)
  # nodesize: minimum number of observations in a terminal node (3 in survival)
  # nsplit: number of random points to explore in continous predictors
  # importance=T: prints out variable importance ranking
  # proximity=T: compute this metric
  ##############################################################################
  featureRF <- function(data, out.type, out1, out2=NA, seed) {
    if(out.type=="survival") {
      form <- as.formula(paste0("Surv(",out1,",",out2,")~."))
      o <- tune(form, data, seed)
      selected.var <- var.select(rfsrc(form, data, seed, ntree=o$rf$ntree, mtry=o$rf$mtry, nodesize=o$rf$nodesize, importance=TRUE), verbose=FALSE)$varselect
      return(selected.var[order(match(rownames(selected.var),names(data)[-which(names(data) %in% c(out1,out2))])),])
    }
    if(out.type=="binary") {
      form <- as.formula(paste0(out1,"~."))
      o <- tune(form, data, seed)
      selected.var <- var.select(rfsrc(form, data, seed, ntree=o$rf$ntree, mtry=o$rf$mtry, nodesize=o$rf$nodesize, importance=TRUE), verbose=FALSE)$varselect
      return(selected.var[order(match(rownames(selected.var),names(data)[-which(names(data)==out1)])),])
    }
  }
  
  
  ## Dataset for feature selection
  dd <- data[,c(which(names(data) %in% c(out1,out2,clinical)), which(names(data)==feature1):which(names(data)==feature2))]
  
  ## Set seed
  if(is.na(seed)) {seed=format(Sys.Date(), "%Y%b%d")}; set.seed(seed) #20190218 수정

  ## Output
  if(is.na(method)) {
    print("Lasso")
    f.LASSO <- featureLSE(dd, out.type, out1, out2, seed, B, method="Lasso", nfold)
    write.xlsx(f.LASSO, file=paste(filename,"_",out1,"_all.xlsx",sep=""), sheetName="Lasso")
    
    print("Stability")
    f.STAB <- featureLSE(dd, out.type, out1, out2, seed, B, method="Stability", nfold=NA)
    #if out.type="survival", print "Stability selection is not available for the survival outcome"
    write.xlsx(f.STAB, file=paste(filename,"_",out1,"_all.xlsx",sep=""), sheetName="Stability", append=TRUE)
    
    print("Elastic net")
    f.ELNET <- featureLSE(dd, out.type, out1, out2, seed, B, method="Elastic net", nfold)
    write.xlsx(f.ELNET, file=paste(filename,"_",out1,"_all.xlsx",sep=""), sheetName="Elastic net", append=TRUE)
    
    print("Random forest")
    f.RF <- featureRF(dd, out.type, out1, out2, seed)
    write.xlsx(f.RF, file=paste(filename,"_",out1,"_all.xlsx",sep=""), sheetName="Random forest", append=TRUE)
    
  } else {
    print(method)
    if(method=="Random forest") {
      write.xlsx(featureRF(dd, out.type, out1, out2, seed), file=paste(filename,"_",out1,"_",method,".xlsx",sep=""), sheetName=method)
    } else {
      write.xlsx(featureLSE(dd, out.type, out1, out2, seed, B, method, nfold), file=paste(filename,"_",out1,"_",method,".xlsx",sep=""), sheetName=method)
    }
  }
  
  return(Sys.time()-tStart)
}






##############################################################################
#### Examples 
#### -------------------------------------------------------------------------
# check and set up the variable type first
# NA's are not allowed in the feature variables
##############################################################################
setwd("C:/Users/gradbios/Google 드라이브/Automation/")

## survival outcome
dat <- read.csv("Example_survival.csv")
dat <- dat[complete.cases(dat),]
featureSel(dat, out.type="survival", "OS2", "Death", feature1="major_axis", feature2="Calvarial_remodeling", method="Random forest", seed=20181120, B=100) 
# -> Result using defined method only


## binary outcome
library(mlbench)
data("breast")
table(breast$status)
breast <- breast[complete.cases(breast),]
sum(is.na(breast))
featureSel(breast, out.type="binary", "status", clinical="tsize", feature1="mean_radius", feature2="worst_fractaldim", seed=20181120, B=100) 
# -> Time difference of XXXX secs

