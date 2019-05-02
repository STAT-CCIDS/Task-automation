###############################
##### internal validation #####
###############################

# This is the function for doing internal validation
# It returns a bootstrap-validated estimate of the AUC, original AUC and multivariable regression table 
# which is used for reporting summary statstics in conventional papers.

# To get proper result from it, there are few things that 
# you should follow and it is written below.

# 1. The format for input data must be a data frame
#    - It can be done by simple command, "as.data.frame"

# 2. The type of variable must be clarified.
#    - If you want to treat the variable as "categorical" variable
#      you should transform it as factor and it accomplished by R command, "as.factor(x)".
#    - In case of continuous, use "as.numeric(x)" for trasnformation.

# 3. It returns file via "csv" format. 

# 4. Parameter description : 
#    -out.type : outcome type, "binary" or "survival"
#    -val.type : type of internal validation, "bootstrap" or "cross"
#    -data : data (data.frame format) which contains independent variables
#    -var.list : There are two opions are available.
#                1) var.list = c("var1 : var2") : summary up var1 to var2
#                2) var.list = c("var1","var2","var3"...) - more than 2var : summary up to only desginated var.
#    -event.colname : "event" variable name (dependent variable)
#    -time.colname : "follow-up time" variable name (only available for survival analysis)
#    -save.dir : specify working directory. 
#                If you don't, it will be automatically assigned to your current path, "setwd()". 
#    -filename : "filename", in case of specfied working directory  or full path for the file.
#    -set.seed : specify seed number. 
#                If you don't, it will be automatically assigned to your current date
#    -B : bootstrap number(only available for val.type "bootstrap"). default 1000
#    -k : cross validation number(only available for val.type "cross"). default 10

#===========================================================================================================================================    
#===========================================================================================================================================    

internal.validation=function( out.type=c("binary","survival"),val.type=c("bootstrap","cross"),
                              data,
                              var.list=NULL,
                              event.colname=NULL,time.colname=NULL,
                              save.dir=NA, filename="Internal Validation result",
                              set.seed=NA,
                              B=1000,k=10){
  
  list.of.packages = c("brglm","progress","glmnet","survival","pROC","cvAUC","coxphf","stringr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {install.packages(new.packages)}
  
  library(brglm)
  library(progress)
  library(glmnet)
  library(survival)
  library(pROC)
  library(ROCR)
  library(cvAUC)
  library(coxphf)
  library(stringr)
  
  #===========================================================================================================================================    
  #Data handling============================================================================================================================== 
  #data reordering (dependent variable, independent variables)================================================================================    
  
  
  if(regexpr(":",var.list)[1] != -1){
    splitstr=strsplit(var.list,split = ":")
    start.colname=str_trim(splitstr[[1]][1],side="both")
    end.colname=str_trim(splitstr[[1]][2],side="both")
    start.colnum=which(colnames(data)==start.colname)
    end.colnum=which(colnames(data)==end.colname)
    colnum.list=c(start.colnum:end.colnum)
  }else{
    colnum.list=c()
    for(i in 1:length(var.list)){
      colnum.list[i]= which(colnames(data)==var.list[i])
    }
  } 
  
  if(out.type=="binary"){
    
    event.colnum=which(colnames(data)==event.colname)
    data=data[,c(event.colnum,colnum.list)]
    event=1
    
  }else if(out.type=="survival"){
    
    #event column number 
    event.colnum=which(colnames(data)==event.colname)
    
    time.colnum=which(colnames(data)==time.colname)
    data=data[,c(time.colnum,event.colnum,colnum.list)]
    event=2
  }
  
  #delete NA  
  data=na.omit(data)
  
  #test complete separation for Firth penalized regression
  com.sep=rep(0,ncol(data)-1)
  i=4
  for(i in (event+1):ncol(data)){
    if(is.factor(data[,i])){
      com.sep[i]=ifelse(sum(table(data[,event],data[,i])==0)>0,1,0)
    }
  }
  #Firth peneralized regression option 
  firth.op=ifelse(sum(com.sep)>0,TRUE,FALSE)
  
  #===========================================================================================================================================    
  #Start bootstrap validation=============================================================================================================== 
  
  set.seed = ifelse(is.na(set.seed),format(Sys.Date(),"%Y%m%d"),set.seed)
  
  if(val.type=="bootstrap"){
    
    b.cindex<-rep(0,B) #c index vector for bootstrap sample
    o.cindex<-rep(0,B)  #c index vector fot original data
    diff.cindex=rep(0,B)
    pb = txtProgressBar(min=0,max=B-1,style=3)
    
    #outcome type: binary (logistic regression)===================================================================================================   
    if(out.type=="binary"){  
      i=1
      while(i < B){
        
        tryCatch({ 
          
          #bootstrap sample
          data.b=data[sample(dim(data)[1],dim(data)[1],replace=TRUE),] 
          
          #fit trained model using bootstrap sample
          if(firth.op){
            fit=suppressWarnings(brglm(data.b[,1]~., family="binomial", data=data.b[,-1]))#Firth penalized logistic regression
          }else{
            fit=glm(data.b[,1]~.,data =data.b[,-1],family="binomial" )
          }
          
          pred=fit$fitted.values
          roc=roc(data.b[,1],pred)
          b.cindex[i]=as.vector(roc$auc) # c index of bootstrap sample using trained model 
          
          o=predict(fit,newdata=data, type="response") 
          roc.o=roc(data[,1],o) 
          o.cindex[i]=as.vector(roc.o$auc)  # c index of original data using trained model 
          
          diff.cindex[i]=b.cindex[i]-o.cindex[i] # calculate the optimism
          
          i=i+1
        }
        , warning = function(w) {}
        , error = function(e) {})
        Sys.sleep(0.01)
        setTxtProgressBar(pb,i)
      }; close(pb)
      
      # fit model using total data 
      if(firth.op){
        fit=brglm(data[,1]~.,data =data[,-1],family="binomial",pl=F)
      }else{
        fit=glm(data[,1]~.,data =data[,-1],family="binomial" )
      }
      
      pred=fit$fitted.values
      original.cindex=roc(data[,1] , pred,ci=T)
      original.AUC=original.cindex$auc # c index of original data using model
      
      #95% confidence interval of original auc 
      original.AUC_25=original.cindex$ci[1]
      original.AUC_975=original.cindex$ci[2]
      
      
      b.validated.est=as.vector(original.AUC)-mean(diff.cindex) # bootstrap-validated estimate of the AUC 
      
      #95% confidence interval of b.validated.est 
      quan_25=quantile(diff.cindex,probs=c(0.025))
      quan_975=quantile(diff.cindex,probs=c(0.975))
      b.validated.est_975=original.cindex$auc-quan_25   
      b.validated.est_25=original.cindex$auc-quan_975 
      
      
      #outcome type: survival (Cox proportional hazard regression)=================================================================================   
    }else if(out.type=="survival"){
      i=1
      while(i < B){
        
        tryCatch({ 
          
          #bootstrap sample
          dat.b=data[sample(dim(data)[1],dim(data)[1],replace=TRUE),]
          
          t.data=data.frame(data[,-c(1,2)]) #total data 
          names(t.data)=names(data)[-c(1,2)]
          t.surv=Surv(data[,1],data[,2]) #survival time of total data 
          surv=Surv(dat.b[,1],dat.b[,2]) #bootstrap sample 
          data.b=data.frame(dat.b[,-c(1,2)]) #survival time of bootstrap sample 
          names(data.b)=names(dat.b)[-c(1,2)]
          
          #fit trained model using bootstrap sample
          if(firth.op){
            fit=coxphf(surv~.,data =data.b)#Firth penalized survival regression
          }else{
            fit=coxph(surv~.,data =data.b)
          }
          
          concor=survConcordance(surv ~ predict(fit), data=data.b) 
          b.cindex[i]=concor$concordance # c index of bootstrap sample using trained model 
          
          o=survConcordance(t.surv ~predict(fit,newdata =t.data), data=t.data)
          o.cindex[i]=o$concordance  # c index of original data using trained model 
          
          diff.cindex[i]=b.cindex[i]-o.cindex[i] # calculate the optimism 
          
          i=i+1
        }
        , warning = function(w) {}
        , error = function(e) {})
        Sys.sleep(0.01)
        setTxtProgressBar(pb,i)
      }; close(pb)
      
      # fit model using total data 
      if(firth.op){
        fit=coxphf(t.surv ~.,data =t.data,pl=F)
      }else{
        fit=coxph(t.surv ~.,data =t.data)
      }
      
      original.cindex=survConcordance(t.surv ~predict(fit), data=t.data) 
      original.AUC=original.cindex$concordance # c index of original data using model
      
      #95% confidence interval of original auc 
      original.AUC_25=original.AUC-1.96*original.cindex$std.err
      original.AUC_975=original.AUC+1.96*original.cindex$std.err
      
      b.validated.est=original.AUC-mean(diff.cindex) # bootstrap-validated estimate of the AUC 
      
      #95% confidence interval of original auc 
      quan_25=quantile(diff.cindex,probs=c(0.025))
      quan_975=quantile(diff.cindex,probs=c(0.975))
      b.validated.est_975=original.cindex$concordance-quan_25 
      b.validated.est_25=original.cindex$concordance-quan_975 
      
    } 
    originalAUC=paste(round(original.AUC,3),"(",round(original.AUC_25,3),",",round(original.AUC_975,3),")",sep = "")
    result=paste(round(b.validated.est,3),"(",round(b.validated.est_25,3),",",round(b.validated.est_975,3),")",sep = "")
    
    AUC.result=data.frame(conf=c(originalAUC,result),p_value=NA)
    rownames(AUC.result)=c("Original AUC","Bootstrap validation AUC")
    
    #===========================================================================================================================================    
    #Start cross validation===================================================================================================================== 
    
  } else if(val.type=="cross"){
    
    cindex = rep(0,k)
    data$id = sample(1:k, nrow(data), replace = TRUE) #k-fold cross validation
    list = 1:k
    pred.list=list();y.list=list()
    
    #outcome type: binary (logistic regression)==================================================================================================   
    
    #fit trained model using k fold sample  
    if(out.type=="binary"){
      
      for(i in 1:k){
        
        # remove rows with id i from dataframe to create training set
        # select rows with id i to create test set
        trainingset = subset(data, id %in% list[-i])
        testset = subset(data, id %in% c(i))
        
        trainingset = trainingset[,-which(colnames(trainingset)=="id")]
        testset = testset[,-which(colnames(testset)=="id")]
        
        
        if(firth.op){
          fit=suppressWarnings(brglm(trainingset[,1]~., family="binomial", data=trainingset[,-1]))#Firth penalized logistic regression
        }else{
          fit=glm(trainingset[,1]~.,data =trainingset[,-1],family="binomial" )
        }
        pred=predict(fit,newdata=testset, type="response")
        pred.list[[i]]=as.vector(pred)
        y.list[[i]]=testset[,1]
        
        cv.auc=ci.cvAUC(pred.list,y.list) 
        # c index of cross validation
        #=============================================================
        #Erin LeDell,Maya L. Petersen(2015)
        #Computationally Efficient Confidence Intervals for Cross-validated Area Under the ROC Curve Estimates
        #Electron J Stat. 2015 ; 9(1): 1583–1607
        #=============================================================
        
        if(firth.op){
          fit=brglm(data[,1]~.,data =data[,-1],family="binomial",pl=F)
        }else{
          fit=glm(data[,1]~.,data =data[,-1],family="binomial" )
        }
        pred=fit$fitted.values
        original.cindex=roc(data[,1] , pred,ci=T) # c index of original data
        
        original.AUC=original.cindex$auc
        
        #95% confidence interval of original auc 
        original.AUC_25=original.cindex$ci[1]
        original.AUC_975=original.cindex$ci[2]
        
        
        result=paste(round(cv.auc$cvAUC,3),"(",round(cv.auc$ci[2],3),",",round(cv.auc$ci[1],3),")",sep = "")
        originalAUC=paste(round(original.AUC,3),"(",round(original.AUC_25,3),",",round(original.AUC_975,3),")",sep = "")
        
      }
      
      #outcome type: survival (Cox proportional hazard regression)=================================================================================   
      
    }else if(out.type=="survival"){
      
      
      for(i in 1:k){
        
        # remove rows with id i from dataframe to create training set
        # select rows with id i to create test set
        trainingset = subset(data, id %in% list[-i])
        testset = subset(data, id %in% c(i))
        
        trainingset = trainingset[,-which(colnames(trainingset)=="id")]
        testset = testset[,-which(colnames(testset)=="id")]
        
        t.data=data[,-c(1,2)]
        t.surv=Surv(data[,1],data[,2])
        testset=testset[,-c(1,2)]
        test.surv=Surv(testset[,1],testset[,2])
        training.surv=Surv(trainingset[,1],trainingset[,2])
        trainingset=trainingset[,-c(1,2)]
        
        if(firth.op){
          fit=coxphf(training.surv~.,data =trainingset)
        }else{
          fit=coxph(training.surv~.,data =trainingset)
          
        }
        
        pred=survConcordance(test.surv ~predict(fit,newdata =testset), data=testset) # original data에 대한 c index
        cindex[i]=pred$concordance
        result=mean(cindex)
        
        # c index about original data
        if(firth.op){
          fit=coxphf(t.surv ~.,data =t.data)
        }else{
          fit=coxph(t.surv ~.,data =t.data)
        }
        
        original.cindex=survConcordance(t.surv ~predict(fit), data=t.data) # c index about original data
        original.AUC=original.cindex$auc
      }
    }; 
    
    AUC.result=data.frame(conf=c(originalAUC,result),p_value=NA)
    rownames(AUC.result)=c("Original AUC","Cross validation AUC")
    
  }
  
  #===========================================================================================================================================    
  #multivariate regression =================================================================================================================== 
  
  if(out.type=="binary"){
    
    if(firth.op){
      cat("Logistic regression with Firth’s penalized maximum likelihood estimation\n")
    }else{
      cat("Logistic regression\n")
    }
    ret=summary(fit)
    p.value=ret$coefficients[,4]
    odds=as.vector(exp(coef(fit)))
    conf=as.data.frame(exp(confint(fit)))
    conf=paste(round(odds[-1],3)," (",round(conf$`2.5 %`[-1],3),", ",round(conf$`97.5 %`[-1],3),")",sep = "")
    regression.result=data.frame(conf=conf,p_value=round(p.value[-1],3))
    
  }else if(out.type=="survival"){
    
    if(firth.op){
      ret=summary(fit)
      p.value=ret$prob
      odds=exp(coef(fit))
      conf=as.data.frame(exp(confint(fit)))
      conf=paste(round(odds,3)," (",round(conf$`2.5 %`,3),", ",round(conf$`97.5 %`,3),")",sep = "")
      cat("Cox proportional hazard regression with Firth’s penalized maximum likelihood estimation\n")
    }
    else{
      ret=summary(fit)
      p.value=ret$coefficients[,5]
      odds=as.vector(exp(coef(fit)))
      conf=as.data.frame(exp(confint(fit)))
      conf=paste(round(odds,3)," (",round(conf$`2.5 %`,3),", ",round(conf$`97.5 %`,3),")",sep = "")
      cat("Cox proportional hazard regression\n")
    }
    regression.result=data.frame(conf=conf,p_value=round(p.value,3))
  }
  
  result=rbind(regression.result,AUC.result)
  
  dir = ifelse(is.na(save.dir),getwd(),save.dir)
  filename=paste(dir,"/",out.type,filename,".xlsx",sep = "")
  tryCatch({ 
    xlsx::write.xlsx(file=filename,result)}
    , error = function(e) {cat("Remove ",filename," file and try again")})
  
  return(result)
}
