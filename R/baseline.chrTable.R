#' @title Create Baseline Characteristics Table by 4 Groups
#'
#' @description The \code{chrTable.by.4group} function takes data and creates a table for baseline characteristics stratified by 4 groups. Both categorical and numerical variables are included and analyzed together. P values will be bolded in the Markdown document if it is less than 0.05.
#' @param variable a vector of variables
#' @param type a vector of types either categorical or numerical (\code{"cat"} or \code{"num"})
#' @param label a vector of labels shown in the table
#' @param data data frame dataset
#' @param digits number of decimal points for numerical variable
#' @param summary three options are available for the numerical variable to be displayed; If \code{median}, median with min and max is calculated. If \code{median.iq}, median with 25 and 75 percentiles is calculated. If \code{mean}, mean with its standard deviation is calculated. Default is \code{"median"}.
#' @param vargroup a variable by which data is staratified
#' @param vargroup.lab a vector of grouping labels for \code{vargroup} in the table
#' @param displaytot logical; If \code{TRUE}, additional column for the whole cohort is included and default is \code{TRUE}.
#' @param chisqtest logical for categorical variable; If \code{TRUE}, Chisquare test is performed. If \code{FALSE}, Fisher's exact test is performed and default is \code{FALSE}.
#' @param anova logical for numerical variable; If \code{TRUE}, anova test is performed. If \code{FALSE}, Kruskal-Wallis test is performed and default is \code{FALSE}.
#' @param useNA option whether or not \code{NULL} value is displayed for the categorical variable in the table. If \code{TRUE}, character in argument \code{na} is displayed for \code{NULL} value with its percentage. Default is \code{FALSE}.
#' @param na character to be displayed for \code{NULL} value if \code{useNA} is \code{TRUE}. Default is \code{"N/A"}.
#' @param col1st character label displayed in the 1st cell in the table
#' @importFrom coin wilcox_test
#' @export
#' @return summary matrix with calculated baseline characteristics by group with P values from selected statistical tests
#' @examples
#' group <- c(rep(1,25), rep(2,25), rep(3,25), rep(4,25))
#' var1 <- sample(1:100, 100, replace=TRUE)
#' var2 <- c(sample(1:50, 25, replace=TRUE), sample(25:75, 25, replace=TRUE), sample(30:80, 25, replace=TRUE), sample(50:100, 25, replace=TRUE))
#' df <- cbind(group, var1, var2)
#'
#' chrTable.by.4group(variable=c("var1","var2"),type=c("num","num"),
#'                    label=c("var1", "var2"), data=df, summary="median",
#'                    vargroup="group", chisqtest=FALSE, anova=FALSE,
#'                    useNA=FALSE, col1st="")
#'



chrTable.by.4group <- function(variable,type,label,data,digits,summary="median",vargroup,vargroup.lab=NULL,
                               displaytot=T,transpose=F,chisqtest=F, anova=F,
                               useNA=F, na="N/A", col1st=""){

  data <- na.handle.2(data, useNA, na)

  numvar<-length(variable)
  grand.table<-NULL
  col<-NULL

  for (k in 1:length(vargroup)){

    var.group<-factor(data[,which(colnames(data)==vargroup[k])])
    group.level<-levels(var.group)
    dat1=subset(data,data[,which(colnames(data)==vargroup[k])]==group.level[1])  # subset of group level[1]
    dat2=subset(data,data[,which(colnames(data)==vargroup[k])]==group.level[2])  # subset of group level[2]
    dat3=subset(data,data[,which(colnames(data)==vargroup[k])]==group.level[3])  # subset of group level[3]
    dat4=subset(data,data[,which(colnames(data)==vargroup[k])]==group.level[4])  # subset of group level[3]

    #Create Table ===================================================================
    #Group 1

    myTable.g1<-NULL
    for(i in 1:numvar){

      if (type[i]=="num"){          # IF numeric variable

        if (summary=="median"){
          row.tab<-c(paste(label[i],', median (min,max)',sep=""),
                     getMedian(dat1[,which(colnames(dat1)==variable[i])],digits=digits))

          row.tab[2]<-getMedian(data[,which(colnames(data)==variable[i])],digits=digits)[1]       #Assign non-missing N total
          myTable.g1<-rbind(myTable.g1,c(row.tab,rep(i,1)))

        }else if (summary=="median.mid"){
          row.tab<-c(paste(label[i],', median (25%,75%)',sep=""),
                     getMedian.mid(dat1[,which(colnames(dat1)==variable[i])],digits=digits))

          row.tab[2]<-getMedian.mid(data[,which(colnames(data)==variable[i])],digits=digits)[1] #Assign non-missing N total
          myTable.g1<-rbind(myTable.g1,c(row.tab,rep(i,1)))

        }else if (summary=="mean"){
          row.tab<-c(paste(label[i],', mean (SD)',sep=""),
                     getMean(dat1[,which(colnames(dat1)==variable[i])],digits=digits))
          row.tab[2]<-getMedian(data[,which(colnames(data)==variable[i])],digits=digits)[1]     #Assign non-missing N total
          myTable.g1<-rbind(myTable.g1,c(row.tab,rep(i,1)))
        }

      }else if (type[i]=="cat"){  # IF categorical variable

        ncat=length(levels(as.factor(data[,which(colnames(data)==variable[i])])))

        row.tab<-c(paste(label[i],', No. (%)',sep=""),NA,NA)
        catvar <- dat1[,which(colnames(dat1)==variable[i])]
        catvar<-factor(catvar)
        numCat<-length(levels(catvar))
        cat.tab=NULL

        if (numCat==0){
          cat.tab<-rbind(cat.tab,c(NA,levels(catvar)[1],getNP(catvar,levels(catvar)[1],digits=digits)))
        } else {
          for(j in 1:numCat){
            cat.tab<-rbind(cat.tab,c(NA,levels(catvar)[j],getNP(catvar,levels(catvar)[j],digits=digits)))
          }
        }

        if (numCat!=ncat){cat.tab<-rbind(cat.tab,matrix(c(NA,NA,NA),nrow=ncat-numCat,ncol=3))}

        cat.var=rbind(row.tab,cat.tab)
        myTable.g1<-rbind(myTable.g1,cbind(cat.var,rep(i,nrow(cat.var))))
      }
    }
    rownames(myTable.g1) <- NULL
    myTable.g1<-as.data.frame(myTable.g1)


    #Group 2
    myTable.g2<-NULL
    for(i in 1:numvar){

      if (type[i]=="num"){        # IF numeric variable

        if (summary=="median"){
          row.tab<-c(paste(label[i],', median (min,max)',sep=""),
                     getMedian(dat2[,which(colnames(dat2)==variable[i])],digits=digits))

          row.tab[2]<-getMedian(data[,which(colnames(data)==variable[i])],digits=digits)[1]     #Assign non-missing N total
          myTable.g2<-rbind(myTable.g2,c(row.tab,rep(i,1)))

        }else if (summary=="median.mid"){
          row.tab<-c(paste(label[i],', median (25%,75%)',sep=""),
                     getMedian.mid(dat2[,which(colnames(dat2)==variable[i])],digits=digits))

          row.tab[2]<-getMedian.mid(data[,which(colnames(data)==variable[i])],digits=digits)[1] #Assign non-missing N total
          myTable.g2<-rbind(myTable.g2,c(row.tab,rep(i,1)))

        }else if (summary=="mean"){
          row.tab<-c(paste(label[i],', mean (SD)',sep=""),
                     getMean(dat2[,which(colnames(dat2)==variable[i])],digits=digits))
          row.tab[2]<-getMedian(data[,which(colnames(data)==variable[i])],digits=digits)[1]     #Assign non-missing N total
          myTable.g2<-rbind(myTable.g2,c(row.tab,rep(i,1)))
        }

      }else if (type[i]=="cat"){ # IF categorical variable

        ncat=length(levels(as.factor(data[,which(colnames(data)==variable[i])])))

        row.tab<-c(paste(label[i],', No. (%)',sep=""),NA,NA)
        catvar <- dat2[,which(colnames(dat2)==variable[i])]
        catvar<-factor(catvar)
        numCat<-length(levels(catvar))
        cat.tab=NULL

        if (numCat==0){
          cat.tab<-rbind(cat.tab,c(NA,levels(catvar)[1],getNP(catvar,levels(catvar)[1],digits=digits)))
        } else {
          for(j in 1:numCat){
            cat.tab<-rbind(cat.tab,c(NA,levels(catvar)[j],getNP(catvar,levels(catvar)[j],digits=digits)))
          }
        }

        if (numCat!=ncat){cat.tab<-rbind(cat.tab,matrix(c(NA,NA,NA),nrow=ncat-numCat,ncol=3))}

        cat.var=rbind(row.tab,cat.tab)
        myTable.g2<-rbind(myTable.g2,cbind(cat.var,rep(i,nrow(cat.var))))
      }
    }
    rownames(myTable.g2) <- NULL
    myTable.g2<-as.data.frame(myTable.g2)


    #Group 3

    myTable.g3<-NULL
    for(i in 1:numvar){

      if (type[i]=="num"){        # IF numeric variable

        if (summary=="median"){
          row.tab<-c(paste(label[i],', median (min,max)',sep=""),
                     getMedian(dat3[,which(colnames(dat3)==variable[i])],digits=digits))

          row.tab[2]<-getMedian(data[,which(colnames(data)==variable[i])],digits=digits)[1]     #Assign non-missing N total
          myTable.g3<-rbind(myTable.g3,c(row.tab,rep(i,1)))

        }else if (summary=="median.mid"){
          row.tab<-c(paste(label[i],', median (25%,75%)',sep=""),
                     getMedian.mid(dat3[,which(colnames(dat3)==variable[i])],digits=digits))

          row.tab[2]<-getMedian.mid(data[,which(colnames(data)==variable[i])],digits=digits)[1] #Assign non-missing N total
          myTable.g3<-rbind(myTable.g3,c(row.tab,rep(i,1)))

        }else if (summary=="mean"){
          row.tab<-c(paste(label[i],', mean (SD)',sep=""),
                     getMean(dat3[,which(colnames(dat3)==variable[i])],digits=digits))
          row.tab[2]<-getMedian(data[,which(colnames(data)==variable[i])],digits=digits)[1]     #Assign non-missing N total
          myTable.g3<-rbind(myTable.g3,c(row.tab,rep(i,1)))
        }

      }else if (type[i]=="cat"){ # IF categorical variable

        ncat=length(levels(as.factor(data[,which(colnames(data)==variable[i])])))

        row.tab<-c(paste(label[i],', No. (%)',sep=""),NA,NA)
        catvar <- dat3[,which(colnames(dat3)==variable[i])]
        catvar<-factor(catvar)
        numCat<-length(levels(catvar))
        cat.tab=NULL

        if (numCat==0){
          cat.tab<-rbind(cat.tab,c(NA,levels(catvar)[1],getNP(catvar,levels(catvar)[1],digits=digits)))
        } else {
          for(j in 1:numCat){
            cat.tab<-rbind(cat.tab,c(NA,levels(catvar)[j],getNP(catvar,levels(catvar)[j],digits=digits)))
          }
        }

        if (numCat!=ncat){cat.tab<-rbind(cat.tab,matrix(c(NA,NA,NA),nrow=ncat-numCat,ncol=3))}

        cat.var=rbind(row.tab,cat.tab)
        myTable.g3<-rbind(myTable.g3,cbind(cat.var,rep(i,nrow(cat.var))))
      }
    }
    rownames(myTable.g3) <- NULL
    myTable.g3<-as.data.frame(myTable.g3)



    #Group 4
    myTable.g4<-NULL
    for(i in 1:numvar){

      if (type[i]=="num"){        # IF numeric variable

        if (summary=="median"){
          row.tab<-c(paste(label[i],', median (min,max)',sep=""),
                     getMedian(dat4[,which(colnames(dat4)==variable[i])],digits=digits))

          row.tab[2]<-getMedian(data[,which(colnames(data)==variable[i])],digits=digits)[1]     #Assign non-missing N total
          myTable.g4<-rbind(myTable.g4,c(row.tab,rep(i,1)))

        }else if (summary=="median.mid"){
          row.tab<-c(paste(label[i],', median (25%,75%)',sep=""),
                     getMedian.mid(dat4[,which(colnames(dat4)==variable[i])],digits=digits))

          row.tab[2]<-getMedian.mid(data[,which(colnames(data)==variable[i])],digits=digits)[1] #Assign non-missing N total
          myTable.g4<-rbind(myTable.g4,c(row.tab,rep(i,1)))

        }else if (summary=="mean"){
          row.tab<-c(paste(label[i],', mean (SD)',sep=""),
                     getMean(dat4[,which(colnames(dat4)==variable[i])],digits=digits))
          row.tab[2]<-getMedian(data[,which(colnames(data)==variable[i])],digits=digits)[1]     #Assign non-missing N total
          myTable.g4<-rbind(myTable.g4,c(row.tab,rep(i,1)))
        }

      }else if (type[i]=="cat"){ # IF categorical variable

        ncat=length(levels(as.factor(data[,which(colnames(data)==variable[i])])))

        row.tab<-c(paste(label[i],', No. (%)',sep=""),NA,NA)
        catvar <- dat4[,which(colnames(dat4)==variable[i])]
        catvar<-factor(catvar)
        numCat<-length(levels(catvar))
        cat.tab=NULL

        if (numCat==0){
          cat.tab<-rbind(cat.tab,c(NA,levels(catvar)[1],getNP(catvar,levels(catvar)[1],digits=digits)))
        } else {
          for(j in 1:numCat){
            cat.tab<-rbind(cat.tab,c(NA,levels(catvar)[j],getNP(catvar,levels(catvar)[j],digits=digits)))
          }
        }

        if (numCat!=ncat){cat.tab<-rbind(cat.tab,matrix(c(NA,NA,NA),nrow=ncat-numCat,ncol=3))}

        cat.var=rbind(row.tab,cat.tab)
        myTable.g4<-rbind(myTable.g4,cbind(cat.var,rep(i,nrow(cat.var))))
      }
    }
    rownames(myTable.g4) <- NULL
    myTable.g4<-as.data.frame(myTable.g4)



    #Total
    myTable.tot<-NULL

    for(i in 1:numvar){
      if (type[i]=="num"){

        if (summary=="median"){
          row.tab<-c(paste(label[i],', median (min,max)',sep=""),
                     getMedian(data[,which(colnames(data)==variable[i])],digits=digits))
          myTable.tot<-rbind(myTable.tot,c(row.tab,rep(i,1)))

        }else if (summary=="median.mid"){
          row.tab<-c(paste(label[i],', median (25%,75%)',sep=""),
                     getMedian.mid(data[,which(colnames(data)==variable[i])],digits=digits))
          myTable.tot<-rbind(myTable.tot,c(row.tab,rep(i,1)))

        }else if (summary=="mean"){
          row.tab<-c(paste(label[i],', mean (SD)',sep=""),
                     getMean(data[,which(colnames(data)==variable[i])],digits=digits))
          myTable.tot<-rbind(myTable.tot,c(row.tab,rep(i,1)))
        }

      }else if (type[i]=="cat"){
        row.tab<-c(paste(label[i],', No. (%)',sep=""),NA,NA)
        catvar <- data[,which(colnames(data)==variable[i])]
        catvar<-factor(catvar)
        numCat<-length(levels(catvar))
        cat.tab=NULL
        for(j in 1:numCat){
          cat.tab<-rbind(cat.tab,c(NA,levels(catvar)[j],getNP(catvar,levels(catvar)[j],digits=digits)))
        }
        cat.var=rbind(row.tab,cat.tab)
        myTable.tot<-rbind(myTable.tot,cbind(cat.var,rep(i,nrow(cat.var))))
      }
    }
    rownames(myTable.tot) <- NULL

    #Add a column with numbers 1 to nrow(myTable.tot) to do sorting afterwards
    index=as.numeric(seq(1,nrow(myTable.tot)))
    myTable.tot<-as.data.frame(myTable.tot)
    myTable.tot<-cbind(myTable.tot,index)

    # P value ===============================================================================
    #Compute column of p-values comparing groups
    mypval=NULL
    for(i in 1:numvar){
      if (type[i]=="num"){

        group <- data[,which(colnames(data)==vargroup[k])]
        var=data[,which(colnames(data)==variable[i])]
        nonmis.var=sum(!is.na(var))
        xg <- split(var, group)

        if (anova==T){
          lm=lm(data[,which(colnames(data)==variable[i])] ~ data[,which(colnames(data)==vargroup[k])], data = data)
          anova=anova(lm)
          pval= round(anova[1,5],4)

        } else {
          kt=kruskal.test(data[,which(colnames(data)==variable[i])] ~ data[,which(colnames(data)==vargroup[k])], data = data)
          pval=round(kt$p.value,4)
        }

        mypval<-c(mypval,pval)

      }else if (type[i]=="cat"){

        if (chisqtest==T){

          #If there is missing or unknown we omit this category to do the test
          tempvar=as.character(data[,which(colnames(data)==variable[i])])
          if (sum((tempvar==na|tempvar==""),na.rm=T)!=0){                # HH: Added additional condition for NA
            tempvar[(tempvar==na|tempvar=="")]<-NA
            chisq=chisq.test(data[,which(colnames(data)==vargroup[k])] , tempvar)
          } else if (sum((tempvar==na|tempvar==""),na.rm=T)==0){
            chisq=chisq.test(data[,which(colnames(data)==vargroup[k])] , data[,which(colnames(data)==variable[i])])
          }
          catvar <- data[,which(colnames(data)==variable[i])]
          catvar<-factor(catvar)
          pval=round(chisq$p.value,4)
          pval.col=c("",pval,rep("",length(levels(catvar))-1))
          mypval<-c(mypval,pval.col)


        } else if (chisqtest==F){

          #If there is missing or unknown we omit this category to do the test
          tempvar=as.character(data[,which(colnames(data)==variable[i])])
          if (sum((tempvar==na|tempvar==""),na.rm=T)!=0){                    # HH: Added addtional condition for NA such as ""

            tempvar[(tempvar==na|tempvar=="")]<-NA     # calculate after transforming na to real NULL value

            if (nlevels(as.factor(tempvar))<2){      # If number of levels is just 1, than no statistical test is conducted
              ft<-list(p.value=NA)
              ft$p.value <- NA
            }else{
              ft=fisher.test(data[,which(colnames(data)==vargroup[k])] , tempvar, workspace=2e8)
            }

          } else if (sum(tempvar==na|tempvar=="",na.rm=T)==0){

            if (nlevels(as.factor(tempvar))<2){     # If number of levels is just 1, than no statistical test is conducted
              ft<-list(p.value=NA)
              ft$p.value <- NA
            }else{
              ft=fisher.test(data[,which(colnames(data)==vargroup[k])] , data[,which(colnames(data)==variable[i])], workspace=2e8)
            }
          }
          catvar <- data[,which(colnames(data)==variable[i])]
          catvar<-factor(catvar)
          pval=round(ft$p.value,4)
          pval.col=c("",pval,rep("",length(levels(catvar))-1)) # p value is located from the second row
          mypval<-c(mypval,pval.col)

        }
      }

    }                                             # End of p value calculation

    # HH: P values NEJM format -----------------------------------------------
    tmp <- tmpc <- as.numeric(mypval)
    tmp[tmpc >= 0.005 & !is.na(tmpc)] <- formatC(tmpc[tmpc >= 0.005 & !is.na(tmpc)], digit=2, format="f")
    tmp[tmpc < 0.005 & tmpc >= 0.001 & !is.na(tmpc)] <- formatC(tmpc[tmpc < 0.005 & tmpc >= 0.001 & !is.na(tmpc)], digit=3, format="f")
    tmp[tmpc >=0.995] <- ">0.99"
    tmp[tmpc < 0.001] <- "<.001"
    tmp[tmpc < 0.05 & !is.na(tmpc)]  <- as.character(paste("**",tmp[tmpc < 0.05 & !is.na(tmpc)],"**",sep=""))  ## added 01/04/2019
    mypval <- tmp

    #Combine table groups -----------------------------------------------------
    colnames(myTable.tot)<-c("variable","level","measure.tot","number","index")
    colnames(myTable.g1) <-c("variable","level","measure.g1","number")
    colnames(myTable.g2) <-c("variable","level","measure.g2","number")
    colnames(myTable.g3) <-c("variable","level","measure.g3","number")
    colnames(myTable.g4) <-c("variable","level","measure.g4","number")

    # Merge tables ------------------------------------------------------------
    #table <- plyr::join_all(list(myTable.tot, myTable.g1, myTable.g2), type="left")
    table <- merge(myTable.tot,myTable.g1,by=c("variable","level","number"),all.x=T)
    table <- merge(table,myTable.g2,by=c("variable","level","number"),all.x=T)
    table <- merge(table,myTable.g3,by=c("variable","level","number"),all.x=T)
    table <- merge(table,myTable.g4,by=c("variable","level","number"),all.x=T)

    #myTable1=merge(myTable.tot,myTable.g1,by=c("variable","level","number"),all.x=T)
    #myTable2=merge(myTable1,myTable.g2,by=c("variable","level","number"),all.x=T)

    #myTable2=myTable2[order(myTable2$index),c(1,2,4,6,7)]
    # table <- table[order(table$index), c(1,2,4,6,7,8)]
    table <- table[order(table$index), c(1,2,4,6,7,8,9)]
    myTable=cbind(table,mypval)
    myTable=as.matrix(myTable)
    rownames(myTable) <- NULL

    # Fill out empty frequency in either group ---------------------------------
    myTable[,'measure.g1'][!is.na(myTable[,'level']) & is.na(myTable[,'measure.g1'])] <- "0 (0.00%)"
    myTable[,'measure.g2'][!is.na(myTable[,'level']) & is.na(myTable[,'measure.g2'])] <- "0 (0.00%)"
    myTable[,'measure.g3'][!is.na(myTable[,'level']) & is.na(myTable[,'measure.g3'])] <- "0 (0.00%)"
    myTable[,'measure.g4'][!is.na(myTable[,'level']) & is.na(myTable[,'measure.g4'])] <- "0 (0.00%)"

    N =  nrow(data)
    N1 = nrow(dat1)
    N2 = nrow(dat2)
    N3 = nrow(dat3)
    N4 = nrow(dat4)

    colnames(myTable) <- c(col1st,
                           '',
                           paste('Total'," ",'(N',"=",N,')', sep=""),
                           paste(ifelse(length(vargroup.lab)!=0, vargroup.lab[1], group.level[1])," ",'(N',"=",N1,')', sep=""),
                           paste(ifelse(length(vargroup.lab)!=0, vargroup.lab[2], group.level[2])," ",'(N',"=",N2,')', sep=""),
                           paste(ifelse(length(vargroup.lab)!=0, vargroup.lab[3], group.level[3])," ",'(N',"=",N3,')', sep=""),
                           paste(ifelse(length(vargroup.lab)!=0, vargroup.lab[4], group.level[4])," ",'(N',"=",N4,')', sep=""),
                           "P-value")

    #Combine all Tables
    grand.table=cbind(grand.table,myTable)

  }                                             ## End of first FOR LOOP

  if (displaytot==F){grand.table<-grand.table[,-c(3)]}

  grand.table[is.na(grand.table)]<-rep('',sum(is.na(grand.table)))

  return(grand.table)
}

