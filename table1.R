##############################################################################
#### Data description: Table 1
#### -------------------------------------------------------------------------
# This is the function for making Table 1.
# It returns summarized table which is used for reporting summary statstics in conventional papers.
#
# To get proper result, there are few things that you should follow and it is written below.
#
# 1. The format for input data must be a data frame.
#
# 2. The type of variable must be clarified.
#
# 3. It returns file via "csv" format so that the working directory must be designated
#    or state full path to "filename" parameter.
#
# 4. Parameter description
#    -data : data (data.frame format)
#    -out1 : column name of outcome variable (Only if there are two outcomes)
#    -clinical : clinical feature name
#    -filename: name of the csv file
#
# Example codes will be offered below.
##############################################################################



#    -var.list : There are three opions are available.
#                1) var.list = "var1" : summary up var1 to end column of dataset
#                2) var.list = c("var1","var2) : summary up var1 to var2
#                3) var.list = c("var1","var2","var3"...) - more than 2var : summary up to only desginated var.
#    -group : group variable 


table1 <- function (data, out1, clinical, filename="Table 1") {
  
  ## Install packages
  list.of.packages <- c("plyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {install.packages(new.packages)}
  library(plyr)
  
  group <- as.factor(data[,which(names(data)==out1)])
  glev <- levels(group)
  dd <- data[,which(names(data) %in% clinical)]
  
  rownum <- 0
  index <- c()
  #g.name <- group
  #g.index<-grep(group,colnames(dat))
  
  if(length(clinical)==1){
    index <- seq(grep(var.list,colnames(dat)),ncol(dat))
  }else if(length(clinical)==2){
    index<- seq(grep(var.list[1],colnames(dat)),grep(var.list[2],colnames(dat)))
  }else{
    for(i in 1:length(var.list)){
      index<-c(index,grep(var.list[i],colnames(dat)))
    }
  }
  
  for(i in index){
    if(colnames(dat)[i]==g.name){next}
    if(is.factor(dat[,i])){
      rownum<-rownum+length(levels(dat[,i]))
    }else{
      rownum<-rownum+1
    }
  }
  table <- matrix(NA,ncol = 4,nrow= rownum)
  colnames(table)<-c("variable",levels(as.factor(group)),"p-value")
  j<-1 #row indicator.
  
  for(i in index){
    if(colnames(dat)[i]==g.name){next}
    
    lev<-levels(dat[,i])
    
    if(is.factor(dat[,i])){ # when the variable is a categorical var.
      for(l in j : (j+length(lev)-1)){
        table[l,1]<-paste(colnames(dat)[i]," ",lev[l-j+1])
      }
      contingency<-table(dat[,i],dat[,g.index])
      if(sum(suppressWarnings(chisq.test(contingency)$expected<5))>=1){ #expected value <5 , exact test.
        table[j,"p-value"]<-round(fisher.test(contingency)$p.value,3)
      }else{ #expected value >5, chisq test.
        table[j,"p-value"]<-round(chisq.test(contingency)$p.value,3)
      }
      
      prop<-round(prop.table(contingency,margin=2)*100,1)
      
      for(k in j:(j+length(lev)-1)){
        table[k,2]<-paste(contingency[(k-j+1),1]," (",prop[(k-j+1),1],"%)")
        table[k,3]<-paste(contingency[(k-j+1),2]," (",prop[(k-j+1),2],"%)")
      }          
      
      j<-j+length(lev)
      
    }
    
    if(!is.factor(dat[,i])){
      
      #chcek the normality
      dat[,i]<-as.numeric(dat[,i])
      nom.p <- ddply(data, .(out1), function(data) {data.frame(shapiro.test = shapiro.test(data[,i])$p.value)})
      
      
      tgr1<-as.numeric(subset(dat[,i],group==as.numeric(glev[1])))
      tgr2<-as.numeric(subset(dat[,i],group==as.numeric(glev[2])))
      if(sum(nom.p[,2][-3] > 0.05)==2){#normality assumption is valid
        table[j,1]<-colnames(dat)[i]
        table[j,"p-value"]<-round(t.test(tgr1,tgr2)$p.value,3)
        table[j,2]<-paste(round(mean(tgr1,na.rm=T))," ± ",round(sd(tgr1,na.rm=T),3))
        table[j,3]<-paste(round(mean(tgr2,na.rm=T))," ± ",round(sd(tgr2,na.rm=T),3))
      }else{#not valid
        table[j,1]<-colnames(dat)[i]
        table[j,"p-value"]<-round(wilcox.test(tgr1,tgr2)$p.value,3)
        table[j,2]<-paste(round(median(tgr1,na.rm=T),3)," (",round(quantile(tgr1,na.rm=T)[2],3),", ",round(quantile(tgr1,na.rm=T)[4],3),")")
        table[j,3]<-paste(round(median(tgr2,na.rm=T),3)," (",round(quantile(tgr2,na.rm=T)[2],3),", ",round(quantile(tgr2,na.rm=T)[4],3),")")
      }
      j<-j+1
    }
    
  }
  write.csv(table,paste(filename,".csv",sep=""),row.names = F)
  return(table)
}
