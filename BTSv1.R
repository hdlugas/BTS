list.of.packages <- c("shiny","shinythemes","Hmisc","ggplot2","psych","car","qqplotr","ggpubr","survival","survminer","plyr","Rmisc","nlme","ISwR","DT")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(shinythemes)
library(Hmisc)
library(ggplot2)
library(psych)
library(car)
library(qqplotr)
library(ggpubr)
library(survival)
library(survminer)
library(plyr)
library(Rmisc)
library(nlme)
library(DT)
library(ISwR)

paircomp <- function(y,x1,g,data,phm="holm",nonparm=FALSE){
  if(F){
    y="size"
    x1="trt"
    g="time"
    data=data
    phm="holm"
    nonparm=FALSE
  }
  xlev = levels(factor(data[,x1]))
  tlen = length(xlev)
  aout = c()
  for(i in 1:tlen){
    tpos = which(data[,x1]==xlev[i])
    if(length(tpos)>0){
      tdata = data[tpos,]
      if(nonparm){
        tmp = pairwise.wilcox.test(tdata[,y],tdata[,g],p.adjust.method="none")
      }else{
        tmp = pairwise.t.test(tdata[,y],tdata[,g],p.adjust.method="none")
      }
      tout = tmp$p.value
      trow = rownames(tout)
      tcol = colnames(tout)
      trn = length(trow)
      tcn = length(tcol)
      tmpout = c()
      for(si in 1:trn){
        for(sj in 1:tcn){
          toutnum = tout[si,sj]
          if(!is.na(toutnum)){
            tmpout = rbind(tmpout,c(xlev[i],trow[si],tcol[sj],toutnum))
          }
        }
      }
      aout = rbind(aout,tmpout)
    }
  }
  dout = as.data.frame(aout)
  dimnames(dout)[[2]] = c(x1,paste0(g,"1"),paste0(g,"2"),"unadjusted.p")
  dout$adjusted.p = p.adjust(as.numeric(as.character(dout$unadjusted.p)),method=phm)
  dout[,4] = sprintf("%.4f",as.numeric(dout[,4]))
  dout[,5] = sprintf("%.4f",as.numeric(dout[,5]))
  dout
}

medianCI <- function(X,conf.level=0.95,na.rm=T)
{
  if(F){
    X = 1:10 #data[which(data$time==0 & data$trt=="Control"),"size"]
    conf.level = 0.95
    na.rm = T
  }
  alpha=1-conf.level
  if(na.rm)
    X = na.omit(X)
  n <- length(X)
  Xsort=sort(X)
  CIl=Xsort[max(1,ceiling(-qnorm(1-alpha/2)*sqrt(n)/2+n/2))]
  CIr=Xsort[min(n,floor(qnorm(1-alpha/2)*sqrt(n)/2+(n+1)/2))]
  c(CIl,CIr)
}

summarySEM <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                       conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     median = median (xx[[col]], na.rm=na.rm),
                     mlow = medianCI (xx[[col]], na.rm=na.rm)[1],
                     mup = medianCI (xx[[col]], na.rm=na.rm)[2]
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

sigp <- function(p){
  tp = sprintf("%.3f",p)
  if(tp=="0.000"){
    p = "p < 0.001"
  }else if(tp=="1.000"){
    p = "p > 0.99"
  }else{
    p = paste("p = ",tp,sep="")
  }
  p
}

sigt <- function(p){
  tp = sprintf("%.3f",p)
  if(tp=="0.000"){
    p = "<0.001"
  }else if(tp=="1.000"){
    p = ">0.99"
  }else{
    p = tp
  }
  p
}

panel.cor <- function(x, y, digits=2, cex.cor,method="spearman")
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1),xaxt="n",yaxt="n")
  #r <- (cor(x, y, method=method, use="pairwise.complete.obs"))
  #txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y, method=method)
  tcor <- c(as.numeric(test$estimate))
  tp <- c(as.numeric(test$p.value))
  txt <- format(c(tcor, 0.123456789), digits=digits)[1]
  Signif <- ifelse(round(tp,3)<0.001,"p<0.001",paste("p=",round(tp,3)))  
  text(0.5, 0.75, paste("r=",txt), font=2)
  if(tp<0.05)
    text(.5, .25, Signif, font=4, col="blue")
  else
    text(.5, .25, Signif, font=2, col="black")
}

panel.lm <- function (x, y,  pch = 18, cex=0.8, col="red",
                      col.lm = "black", method="pearson") 
{   
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch,ylim = ylim, xlim= xlim, col=col,xaxt="n",yaxt="n") #,...)
  ok <- is.finite(x) & is.finite(y)
  test <- cor.test(x,y, method=method)
  if(test$p.value<0.05){
    col.lm = "blue"
  }
  if (any(ok)){
    if(method=="pearson"){
      abline(lm(y[ok]~ x[ok]),col = col.lm, lwd=2) #, ...)
    }
  }
}

cor.mest <- function(mat, method, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method=method, ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$estimate
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

cor.mp <- function(mat, method, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method=method, ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

extract_help <- function(pkgfn, csvfile, to = c("txt", "html", "latex", "ex"))
{
  
  desc.data <- NULL
  if(pkgfn == ""){
    desc.data = csvfile$datapath
  }else if(pkgfn != ""){
    tname = strsplit(pkgfn,"::")[[1]]
    tpos.pkg = which(rdlist$results[,1]==tname[1])
    tpos.dat = which(rdlist$results[,3]==tname[2])
    if(length(tpos.pkg)>0 & length(tpos.dat)>0){
      
      pkg = tname[1]
      fn = tname[2]
      
      to <- match.arg(to)
      rdbfile <- file.path(find.package(pkg), "help", pkg)
      rdb <- tools:::fetchRdDB(rdbfile, key = fn)
      convertor <- switch(to, 
                          txt   = tools::Rd2txt, 
                          html  = tools::Rd2HTML, 
                          latex = tools::Rd2latex, 
                          ex    = tools::Rd2ex
      )
      f <- function(x) capture.output(convertor(x))
      if(is.null(fn)){ 
        desc.data = lapply(rdb, f) 
      }else{
        desc.data = f(rdb)
      }
    }
  }
  desc.data
}

# list of all R data
rdlist = data(package = .packages(all.available = TRUE))

#shinyUI(
ui =
  fluidPage(
    
    tags$head(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      )
    ),
    
    theme = shinytheme("united"),
    #titlePanel("Basic Testing System (BTS), Version CB7600"),
    #titlePanel(title=div(img(src="BTS_image.jpg",height="5%",width="5%"), "Basic Testing System (BTS), Version 1.0")),
    titlePanel("BTS",title=div(img(src="BTS_image.jpg",height="3.7%",width="3.7%",align="top"), "Basic Testing System (BTS), Version 1.0",style="text-align: center; color: #00aed1;")),
    #br(),
    #br(),
    br(),

    sidebarPanel(
      
      #wellPanel(
      
      # SliderbarPanel for BTS tab
      conditionalPanel(#condition = "$('li.active a').first().html()==='About T-test'",
        'input.method === "About BTS"',
        helpText("Basic Testing System (BTS)",align="center",style="color: #2d6a4f;"),
        div(img(src='BTS_image.jpg',height="70%",width="70%"),style="text-align: center;") #,align="center")
      ),
      
      # SlidebarPanel for file upload tab
      conditionalPanel(#condition = "$('li.active a').first().html()==='Data View'",
        'input.method === "Upload Data"',
        
        fileInput('file1', 'Choose CSV file',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        
        checkboxInput('header', 'Header', TRUE),
        radioButtons('sep', 'Separator',
                     c(Comma=',',
                       Semicolon=';',
                       Tab='\t'),
                     ','),
        radioButtons('quote', 'Quote',
                     c(None='',
                       'Double Quote'='"',
                       'Single Quote'="'"),
                     '"'),
        
        tags$hr(),
        textInput("rname.data", "Choose R data (Type package_name::data_name; e.g., survival::aml)", ""),
        tags$hr()
        #,
        #radioButtons("warn",
        #             "Do you want to turn R warnings on?",
        #             choices = c("No" = "-1",
        #                         "Yes" = "0"))
      ),
      
      # SlidebarPanel for variable view tab
      conditionalPanel(#condition = "$('li.active a').first().html()==='Data View'",
        'input.method === "Variable View (Continuous)"',
        selectInput("var1_var", 
                    label = "Please Select a Numerical Variable",
                    ""
        ),
        sliderInput("bins.var",
                    "Numer of bins:",
                    min = 1,
                    max = 50,
                    value = 2
        ),
        radioButtons("transform.var",
                     "Please choose a type of transformation:",
                     choices = c("No transformation" = "none",
                                 "log" = "logtrans", 
                                 "sqaure root" = "srtrans")),
        radioButtons("isgroup_var",
                     "Do you want to look at by a group?",
                     choices = c("No" = "no",
                                 "Yes" = "yes")),
        
        conditionalPanel(condition = "input.isgroup_var === 'yes'",
                         selectInput("group.var", 
                                     label = "Please Select a group Variable",
                                     ""
                         ),
        ),
        
        helpText("Note: This is only for numerical variables")
        
      ),
      
      # SliderbarPanel for t-test tab
      conditionalPanel(#condition = "$('li.active a').first().html()==='T-test'",
        'input.method === "T-test"',
        radioButtons("nonparm_ttest",
                     "Nonparametric test?",
                     choices = c("No" = "no",
                                 "Yes" = "yes")),
        radioButtons("transform.ttest",
                     "Please choose a type of transformation:",
                     choices = c("No transformation" = "none",
                                 "log" = "logtrans", 
                                 "sqaure root" = "srtrans")),
        radioButtons("sample_ttest",
                     "Please choose one sample t test or two sample t test:",
                     choices = c("One sample" = "oneSamp",
                                 "Paired sample" = "pairSamp",
                                 "Two sample" = "twoSamp")),
        selectInput("var1_ttest", 
                    label = "Please Select a Numerical Variable (X1 or X)",
                    ""
        ),
        conditionalPanel(condition = "input.sample_ttest == 'pairSamp'",
                         selectInput("pair2.ttest", 
                                     label = "Please Select a Numerical Variable (X2)",
                                     ""
                         )
        ),
        
        conditionalPanel(condition = "input.sample_ttest == 'twoSamp'",
                         radioButtons("isgroup_ttest",
                                      "Do you want to look at by a group?",
                                      choices = c("No" = "no",
                                                  "Yes" = "yes")),
                         
                         conditionalPanel(condition = "input.isgroup_ttest == 'no'",
                                          selectInput("var2.ttest", 
                                                      label = "Please Select a Numerical Variable (X2)",
                                                      ""
                                          ),
                         ),
                         
                         conditionalPanel(condition = "input.isgroup_ttest == 'yes'",
                                          selectInput("group2.ttest", 
                                                      label = "Please Select a Group Variable (G)",
                                                      ""
                                          ),
                         ),
                         
                         radioButtons("varequal.ttest",
                                      "Are the two samples have equal variance:",
                                      choices = c("Yes" = "yes",
                                                  "No" = "no"))
        ),
        
        selectInput("tail.ttest",
                    label = "Please Select an alternative hypothesis you want to test:",
                    choices = c("Not Equal (X1 != X2)" = "two.sided", 
                                "Less (X1 < X2)" = "less",
                                "Greater (X1 > X2)" = "greater")),
        conditionalPanel(condition = "input.sample_ttest == 'oneSamp'",
                         numericInput("null.ttest",
                                      "Mean value You Want to Test (X2)",
                                      value = 0
                                      
                         )
        ),
        numericInput("conf.ttest",
                     label = "Please Select a confidence level:",
                     value = 0.95,
                     min = 0.8,
                     max = 0.99),
        helpText("Note: Please assign a number between 0 and 1 in the numeric Input")
        
      ),
      
      # SliderbarPanel for one-way ANOVA
      conditionalPanel(#condition = "$('li.active a').first().html()==='About T-test'",
        'input.method === "One-Way ANOVA"',
        radioButtons("nonparm_anova",
                     "Nonparametric test?",
                     choices = c("No" = "no",
                                 "Yes" = "yes")),
        selectInput("var1_anova", 
                    label = "Please Select a Numerical Variable (X1)",
                    ""
        ),
        radioButtons("transform.anova",
                     "Please choose a type of transformation for X1:",
                     choices = c("No transformation" = "none",
                                 "log" = "logtrans", 
                                 "sqaure root" = "srtrans")
        ),
        selectInput("var2_anova", 
                    label = "Please Select a group Variable (G)",
                    ""
        ),
        radioButtons("varequal.anova",
                     "Are the samples have equal variance:",
                     choices = c("Yes" = "yes",
                                 "No" = "no")
        ),
        
        radioButtons("posthoc.anova",
                     "Please choose a correction method for pairwise post-hoc p-values",
                     choices = c("Holm's" = "holm",
                                 "Bonferroni" = "bonferroni", 
                                 "Benjamini-Hochberg" = "BH",
                                 "None" = "none")
        )
        #,
        #helpText("Note: Please assign a number between 0 and 1 in the numeric Input")
        
      ),
      
      # SliderbarPanel for Two-way ANOVA
      conditionalPanel(#condition = "$('li.active a').first().html()==='About T-test'",
        'input.method === "Two-Way ANOVA"',
        selectInput("var1_tanova", 
                    label = "Please select an outcome variable (Y)",
                    ""
        ),
        radioButtons("transform.tanova",
                     "Please choose a type of transformation for Y only for the hypothesis testing:",
                     choices = c("No transformation" = "none",
                                 "log" = "logtrans", 
                                 "sqaure root" = "srtrans",
                                 "rank" = "ranktrans")
        ),
        selectInput("var2_tanova", 
                    label = "Please select the first factor (X)",
                    ""
        ),
        selectInput("var3_tanova", 
                    label = "Please select the second factor (G)",
                    ""
        ),
        
        textInput("xlab.tanova", "X-axis label (X):", "Time (in days)"),
        textInput("ylab.tanova", "Y-axis label (Y):", "Tumor volume"),
        textInput("tlab.tanova", "Legend label (G):", "Treatment"),
        numericInput("pdodge.tanova", label = "Gap of boxes between groups", value = 1),
        
        radioButtons("posthoc.tanova",
                     "Please choose a correction method for pairwise post-hoc p-values",
                     choices = c("Holm's" = "holm",
                                 "Bonferroni" = "bonferroni", 
                                 "Benjamini-Hochberg" = "BH",
                                 "None" = "none")
        )
        #,
        #helpText("Note: Please assign a number between 0 and 1 in the numeric Input")
        
      ),
      
      # SliderbarPanel for Tumor growth curve
      conditionalPanel(#condition = "$('li.active a').first().html()==='About T-test'",
        'input.method === "Tumor Growth Analysis"',
        selectInput("var1_tmg", 
                    label = "Please select a variable for individual ids (I)",
                    ""
        ),
        selectInput("var2_tmg", 
                    label = "Please select an outcome variable (Y)",
                    ""
        ),
        selectInput("var3_tmg", 
                    label = "Please select a time variable (T)",
                    ""
        ),
        selectInput("var4_tmg", 
                    label = "Please select a group variable (G)",
                    ""
        ),
        
        textInput("xlab.tmg", "X-axis label (T):", "Time (in days)"),
        textInput("ylab.tmg", "Y-axis label (Y):", "Tumor volume"),
        radioButtons("bar.tmg",
                     "Please choose a type of error bar in the mean plot:",
                     choices = c("mean+-SD" = "msd",
                                 "mean+-SE" = "mse", 
                                 "mean+-CI" = "mci",
                                 "median+-CI" = "mdci")
        ),
        radioButtons("transform.tmg",
                     "Please choose a type of transformation for Y only for the hypothesis testing:",
                     choices = c("No transformation" = "none",
                                 "log" = "logtrans", 
                                 "sqaure root" = "srtrans",
                                 "rank" = "ranktrans")
        ),
        numericInput("pdodge.tmg", label = "Gap of error bars between groups", value = 0)
        #,
        #helpText("Note: Please assign a number between 0 and 1 in the numeric Input")
        
      ),
      
      # SliderbarPanel for Categorical Test
      conditionalPanel(#condition = "$('li.active a').first().html()==='About T-test'",
        'input.method === "Categorical Test"',
        radioButtons("test_cat",
                     "Please choose a type of test:",
                     choices = c("Chi-square" = "chisq",
                                 "Fisher's exact" = "fisher", 
                                 "McNemar's" = "mcnemar",
                                 "Trend" = "trend"
                     )
        ),
        selectInput("var1_cat", 
                    label = "Please select a categorical variable (X1)",
                    ""
        ),
        conditionalPanel(condition = "input.test_cat !== 'trend'",
                         selectInput("var2_cat", 
                                     label = "Please select a categorical variable (X2)",
                                     ""
                         )
        ),
        conditionalPanel(condition = "input.test_cat === 'trend'",
                         selectInput("time_cat", 
                                     label = "Please select a time or group variable (G)",
                                     ""
                         )
        )
        #,
        
        #helpText("Note: Please assign a number between 0 and 1 in the numeric Input")
        
      ),
      
      # SliderbarPanel for Regression
      conditionalPanel(#condition = "$('li.active a').first().html()==='About T-test'",
        'input.method === "Regression"',
        radioButtons("model_reg",
                     "Please choose a type of regression model:",
                     choices = c("Linear" = "linear",
                                 "Logistic" = "logistic", 
                                 "Cox" = "cox")
        ),
        conditionalPanel(condition = "input.model_reg === 'linear'",
                         selectInput("res.reg", 
                                     label = "Please select a response variable",
                                     "",
                                     multiple=FALSE
                         ),
                         radioButtons("restrans_reg",
                                      "Please choose a type of transformation for a response variable:",
                                      choices = c("No transformation" = "none",
                                                  "log" = "logtrans", 
                                                  "sqaure root" = "srtrans")
                         )
        ),
        conditionalPanel(condition = "input.model_reg == 'logistic'",
                         selectInput("logit.reg", 
                                     label = "Please select a response variable",
                                     "",
                                     multiple=FALSE
                         )
        ),
        conditionalPanel(condition = "input.model_reg === 'cox'",
                         selectInput("ct.reg", 
                                     label = "Please select a time duration variable",
                                     "",
                                     multiple=FALSE
                         ),
                         selectInput("cs.reg", 
                                     label = "Please select a status variable",
                                     "",
                                     multiple=FALSE
                         )
        ),
        selectInput("indep_reg", 
                    label = "Please select covariates",
                    "",
                    multiple=TRUE
        ),
        selectInput("indepcat_reg", 
                    label = "Please indicate categorical covariates among the selected covariates",
                    "",
                    multiple=TRUE
        ),
        selectInput("scale.reg", 
                    label = "Please select continuous covariates for scaling",
                    "",
                    multiple=TRUE
        )
        #,
        #textInput("int.reg", "Please type the interaction terms  with separating with a space", ""),
        
        #helpText("Note: Please assign a number between 0 and 1 in the numeric Input")
        
      ),
      
      # SliderbarPanel for Correlation
      conditionalPanel(#condition = "$('li.active a').first().html()==='About T-test'",
        'input.method === "Correlation"',
        selectInput("var_cor", 
                    label = "Please Select Numerical Variable(s)",
                    "",
                    multiple=TRUE
        ),
        selectInput("log.cor", 
                    label = "Please Select Numerical Variable(s) for log-transformation",
                    "",
                    multiple=TRUE
        ),
        selectInput("sqr.cor", 
                    label = "Please Select Numerical Variable(s) for squre root-transformation",
                    "",
                    multiple=TRUE
        ),
        radioButtons("cmet_cor",
                     "Please choose a type of correlation:",
                     choices = c("Pearson's" = "pearson",
                                 "Spearman's" = "spearman", 
                                 "Kendall's" = "kendall")
        )
        #,
        
        #helpText("Note: Please assign a number between 0 and 1 in the numeric Input")
        
      ),
      
      # SliderbarPanel for Survival Curve
      conditionalPanel(#condition = "$('li.active a').first().html()==='About T-test'",
        'input.method === "Survival Curve"',
        selectInput("t_surv", 
                    label = "Please Select a survival time variable",
                    ""
        ),
        selectInput("s_surv", 
                    label = "Please Select a survival status variable",
                    ""
        ),
        
        radioButtons("unit1.surv",
                     "Input Unit:",
                     choices = c("Days" = "days",
                                 "Months" = "months",
                                 "Years" = "years"
                     )
        ),
        radioButtons("unit2.surv",
                     "Ouput Unit:",
                     choices = c("Days" = "days",
                                 "Months" = "months",
                                 "Years" = "years"
                     )
        ),
        
        textInput("xlab.surv", "X-axis label:", "Days After Diagnosis"),
        textInput("ylab.surv", "Y-axis label:", "Overall Survival"),
        
        numericInput("rate_surv", label = "Time point for a survival rate", value = 1),
        
        radioButtons("group_surv",
                     "Group comparison?",
                     choices = c("No" = "no",
                                 "Yes" = "yes")
        ),
        
        conditionalPanel(condition = "input.group_surv == 'yes'",
                         
                         selectInput("g.surv", 
                                     label = "Please Select a group variable",
                                     ""
                         ),
                         
                         radioButtons("pval.surv",
                                      "Add the p-value in the plot?",
                                      choices = c("Yes" = "yes",
                                                  "No" = "no")
                         ),
                         
                         radioButtons("grouporder_surv",
                                      "Change the group order",
                                      choices = c("No" = "no",
                                                  "Yes" = "yes")
                         ),
                         conditionalPanel(condition = "input.grouporder_surv == 'yes'",
                                          textInput("g.ref.surv", "Group order (starting from a reference group with separating with a space):", "")
                         )
        )
        #,        
        
        #helpText("Note: Please assign a number between 0 and 1 in the numeric Input")
        
      ),
      
      # SliderbarPanel for contingency table
      conditionalPanel(#condition = "$('li.active a').first().html()==='T-test'",
        'input.method === "Contingency Table"',
        #wellPanel(
        radioButtons("type_ctab",
                     "Please choose a type of test:",
                     choices = c("Chi-square test" = "chisq",
                                 "Fisher's exact test" = "fisher", 
                                 "McNemar's test" = "mcnemar",
                                 #"Binominal test" = "binom",
                                 "Trend test" = "trend"
                     )),
        checkboxInput('header.ctab', 'Header', FALSE),
        textAreaInput("data.ctab", "Data:", "39 27\n66 153"),
        
        conditionalPanel(condition = "input.type_ctab == 'fisher'",
                         
                         selectInput("tail.ctab",
                                     label = "Please Select an alternative hypothesis you want to test:",
                                     choices = c("Not Equal" = "two.sided", 
                                                 "Less" = "less",
                                                 "Greater" = "greater")),
                         numericInput("conf.ctab",
                                      label = "Please Select a confidence level:",
                                      value = 0.95,
                                      min = 0.8,
                                      max = 0.99)
        )
      ),
      
      # SliderbarPanel for multiple comparison correction
      conditionalPanel(#condition = "$('li.active a').first().html()==='T-test'",
        'input.method === "Multiple Comparison Correction"',
        #wellPanel(
        radioButtons("holm_mcorrect",
                     "Holm's (1979)",
                     choices = c("No" = "no",
                                 "Yes" = "yes")),
        radioButtons("bonferroni_mcorrect",
                     "Bonferroni correction",
                     choices = c("No" = "no",
                                 "Yes" = "yes")),
        radioButtons("hommel_mcorrect",
                     "Hommel (1988)",
                     choices = c("No" = "no",
                                 "Yes" = "yes")),
        radioButtons("hochberg_mcorrect",
                     "Hochberg (1988)",
                     choices = c("No" = "no",
                                 "Yes" = "yes")),
        radioButtons("BH_mcorrect",
                     "Benjamini-Hochberg (1995; BH; aka FDR)",
                     choices = c("No" = "no",
                                 "Yes" = "yes")),
        radioButtons("BY_mcorrect",
                     "Benjamini-Yekutieli (2001; BY)",
                     choices = c("No" = "no",
                                 "Yes" = "yes")),
        textAreaInput("data.mcorrect", "Insert p-values:", "0.01 0.5 0.7 0.03\n0.023 0.0001")
      )
      #,submitButton("Submit") # View", icon("refresh"))
    ),
    mainPanel(
      tabsetPanel(
        id = 'method',
        tabPanel('About BTS',
                 #h2("The Basic Testing System (BTS) is an R/Shiny based statistical software package for performing basic statistical analyses."),
                 #h2("It must be used only for CB7600-related activities, and distribution outside of class is prohibited."),
                 #h2("Developer: Seongho Kim"),
                 #h2("Questions?"),
                 #h2("Email: kimse at karmanos dot org"),
                 #br(),
                 #br()
                 #,div(img(src='BTS_image.jpg',height="30%",width="30%"),style="text-align: center;") #,align="center")
                 
                 h4("Welcome to Basic Testing System (BTS), an R/Shiny-based statistical software app for teaching, learning, and performing basic statistical/biostatistical analyses. It can be utilized as user-friendly, supportive statistical software for easily teaching and learning basic statistical or biostatistical concepts within biological or non-statistical disciplines. Besides these supportive learning and teaching capacity, for non-statistical scientists, it can be easily utilized to conduct basic statistical or biostatistical analyses for their research. The following statistical analyses can be conducted using the app:",style="color: #194b78;"),
                
                 tags$ul(   
                   tags$li("The contingency table, including Chi-square test, Fisher’s exact test, McNemar’s test, and trend test.",style="color: #00aed1;"),
                   tags$li("Multiple comparison correction adjustments to p-values, including Holm’s, Bonferroni, Hommel, Hochberg, Benjamini-Hochberg, and Benjamini-Yekutieli adjustments.",style="color: #006aa2;"),
                   tags$li("Evaluation of numerical and graphical summary statistics of datasets.",style="color: #00aed1;"),
                   tags$li("Non-parametric and parametric tests for one, two, and paired sample comparisons, including t-test and Wilcoxon rank sum test.",style="color: #006aa2;"),
                   tags$li("One-Way and Two-Way ANOVA.",style="color: #00aed1;"), 
                   tags$li("Tumor growth Analysis.",style="color: #006aa2;"),
                   tags$li("Tests for categorical data, including Chi-square test, Fisher’s exact test, McNemar’s test, and trend test.",style="color: #00aed1;"), 
                   tags$li("Correlation calculation between variables, including Pearson, Spearman, and Kendall correlations.",style="color: #006aa2;"), 
                   tags$li("Regression analysis, including linear and logistic regression.",style="color: #00aed1;"),
                   tags$li("Survival analysis, including Cox regression.",style="color: #006aa2;")
                 ),
                 
                 h4("For many of these analyses, an example has been provided by default. By studying and manipulating these examples, users can easily learn how to conduct valid statistical analyses. Some approaches will activate upon uploading data to the app.",style="color: #194b78;"),
                 
                 h4("The user tutorial for BTS Version 1.0 is available below:",style="color: #ff8500;"),
                 tags$a(href="BTSv1.0_User Tutorial_v07262024.pdf","User Tutorial Verision 1.0",target="_blank"),
                 h5("Please direct your questions to Hunter Dlugas (Email: dlugash at karmanos.org), Janaka Liyanage (Email: liyanagej at karmanos.org) and Seongho Kim (Email: kimse at karmanos.org; founder and creator), the developers of the software.",style="color: #00aed1;"),
                 br()
        ),
        
        tabPanel('Contingency Table',
                 h2("Data"),
                 verbatimTextOutput('tab.ctab'),
                 conditionalPanel(condition = "input.type_ctab == 'chisq'",
                                  h2("Expected values"),
                                  verbatimTextOutput('etab.ctab')
                 ),
                 #h2("Summary statistics"),
                 #verbatimTextOutput('summ.ctab'),
                 h2("Hypothesis test"),
                 verbatimTextOutput('out.ctab')
        ),
        tabPanel('Multiple Comparison Correction',
                 h2("p-values"),
                 verbatimTextOutput('tab.mcorrect'),
                 #h2("Summary statistics"),
                 #verbatimTextOutput('summ.ctab'),
                 h2("Corrected p-values"),
                 verbatimTextOutput('out.mcorrect'),
                 h2("Details from R help"),
                 h6("The Holm (1979), Bonferroni, Hommel (1988), and Hochberg (1988) methods are designed to give strong control of the family-wise error rate. There seems no reason to use the unmodified Bonferroni correction because it is dominated by Holm's method, which is also valid under arbitrary assumptions"),
                 h6("Hochberg's and Hommel's methods are valid when the hypothesis tests are independent or when they are non-negatively associated (Sarkar, 1998; Sarkar and Chang, 1997). Hommel's method is more powerful than Hochberg's, but the difference is usually small and the Hochberg p-values are faster to compute."),
                 h6("The BH and BY methods of Benjamini, Hochberg, and Yekutieli control the false discovery rate (FDR), the expected proportion of false discoveries amongst the rejected hypotheses. The false discovery rate is a less stringent condition than the family-wise error rate, so these methods are more powerful than the others.")
        ),
        
        tabPanel("Upload Data", 
                 #"This panel is intentionally left blank",
                 navbarPage(
                   title = 'DataTable Options',
                   tabPanel('Data',     DT::dataTableOutput('tab.data')),
                   tabPanel('Summary',verbatimTextOutput('summ.data')),
                   tabPanel('Description',htmlOutput('desc.data')),
                   tabPanel('R Data list', DT::dataTableOutput('rdata.data'))
                 )
        ),
        
        tabPanel('Variable View (Continuous)',
                 conditionalPanel(condition = "(input.var1_var !== 'undefined' && input.var1_var.length > 0)",
                 
                   fluidRow(column(10, offset = 1,
                                   h2("Histogram"),
                                   plotOutput('hisgraph.var'))),
                   fluidRow(column(10, offset = 1,
                                   h2("Density plot"),
                                   plotOutput('dengraph.var'))),
                   fluidRow(column(10, offset = 1,
                                   h2("Boxplot"),
                                   plotOutput('boxgraph.var'))),
                   fluidRow(column(10, offset = 1,
                                   h2("Q-Q plot"),
                                   plotOutput('qqgraph.var'))),
                   fluidRow(column(10, offset = 1,
                                   h2("Sharpiro-Wilk's normality test"),
                                   verbatimTextOutput('swtest.var')
                   )
                   ),
                   fluidRow(column(10, offset = 1,
                                   h2("Summary"),
                                   verbatimTextOutput('summ.var')
                   )
                   )
                 
                 )
        ),           
        tabPanel('T-test',
                 conditionalPanel(condition = "(input.var1_ttest !== 'undefined' && input.var1_ttest.length > 0)",

                   fluidRow(column(10, offset = 1,
                                   h2("Boxplot"),
                                   plotOutput('boxgraph.ttest'))),
                   
                   conditionalPanel(condition = "input.nonparm_ttest == 'no'",
                                    fluidRow(column(10, offset = 1,
                                                    h2("Q-Q plot"),
                                                    plotOutput('qqgraph.ttest'))),
                                    fluidRow(column(10, offset = 1,
                                                    h2("Sharpiro-Wilk's normality test"),
                                                    verbatimTextOutput('swtest.ttest')))
                   ),
                   
                   fluidRow(column(10, offset = 1,
                                   h2("Summary statistics"),
                                   verbatimTextOutput('summ.ttest'),
                                   conditionalPanel(condition = "input.sample_ttest === 'twoSamp' && input.nonparm_ttest === 'no'",
                                                    h2("Levene's homogeneity of variance test"),
                                                    verbatimTextOutput('levene.ttest')
                                   ),
                                   h2("Hypothesis testing"),
                                   verbatimTextOutput('out.ttest')
                   ))
                 
                 )
        ),
        tabPanel('One-Way ANOVA',
                 conditionalPanel(condition = "(input.var1_anova !== 'undefined' && input.var1_anova.length > 0 &&
                                  input.var2_anova !== 'undefined' && input.var2_anova.length > 0
                                  )",
                 
                   fluidRow(column(10, offset = 1,
                                   h2("Boxplot"),
                                   plotOutput('boxgraph.anova'))),
                   fluidRow(column(10, offset = 1,
                                   h2("Mean plot"),
                                   plotOutput('meangraph.anova'))),
                   conditionalPanel(condition = "input.nonparm_anova == 'no'",
                                    fluidRow(column(10, offset = 1,
                                                    h2("Q-Q plot"),
                                                    plotOutput('qqgraph.anova'))),
                                    fluidRow(column(10, offset = 1,
                                                    h2("Sharpiro-Wilk's normality test"),
                                                    verbatimTextOutput('swtest.anova'))),
                                    fluidRow(column(10, offset = 1,
                                                    h2("Residual plot"),
                                                    plotOutput('resgraph.anova'))),
                                    fluidRow(column(10, offset = 1,
                                                    h2("Bartlett's homogeneity test"),
                                                    verbatimTextOutput('bltest.anova')))
                   ),
                   fluidRow(column(10, offset = 1,
                                   h2("Summary statistics"),
                                   verbatimTextOutput('summ.anova'),
                                   h2("One-way ANOVA"),
                                   verbatimTextOutput('out.anova'),
                                   h2("Post-hoc pairwise comparisons"),
                                   verbatimTextOutput('posthoc.anova')
                   ))
                 
                 )
        ),
        tabPanel('Two-Way ANOVA',
                 conditionalPanel(condition = "(input.var1_tanova !== 'undefined' && input.var1_tanova.length > 0 &&
                                  input.var2_tanova !== 'undefined' && input.var2_tanova.length > 0 &&
                                  input.var3_tanova !== 'undefined' && input.var3_tanova.length > 0
                                  )",
                                  
                   h2("Q-Q and density plot"),
                   plotOutput('qqdplot.tanova'),
                   h2("Sharpiro-Wilk's normality test"),
                   verbatimTextOutput('swtest.tanova'),
                   h2("Residual plot"),
                   plotOutput('resplot.tanova'),
                   h2("Interaction plot"),
                   plotOutput('interplot.tanova'),
                   h2("Box plot"),
                   plotOutput('boxplot.tanova'),
                   h2("Summary"),
                   verbatimTextOutput('out.tanova'),
                   h2("Post-hoc pairwise comparisons"),
                   verbatimTextOutput('postout.tanova')
                 
                 )
        ),
        tabPanel('Tumor Growth Analysis',
                 conditionalPanel(condition = "(input.var1_tmg !== 'undefined' && input.var1_tmg.length > 0 &&
                                  input.var2_tmg !== 'undefined' && input.var2_tmg.length > 0 &&
                                  input.var3_tmg !== 'undefined' && input.var3_tmg.length > 0 &&
                                  input.var4_tmg !== 'undefined' && input.var4_tmg.length > 0
                                  )",
                                  
                   h2("Spider plot"),
                   plotOutput('spiderplot.tmg'),
                   h2("Spider plot by group"),
                   plotOutput('gspiderplot.tmg'),
                   h2("Mean growth plot"),
                   plotOutput('meanplot.tmg'),
                   h2("Linear mixed-effects model"),
                   h3("Summary"),
                   verbatimTextOutput('out.tmg'),
                   h3("Q-Q and density plot"),
                   plotOutput('qqdplot.tmg'),
                   h2("Sharpiro-Wilk's normality test"),
                   verbatimTextOutput('swtest.tmg'),
                   h3("Residual plot"),
                   plotOutput('resplot.tmg')
                   
                 )
        ),
        tabPanel('Categorical Test',
                 conditionalPanel(condition = "(input.var1_cat !== 'undefined' && input.var1_cat.length > 0 &&
                                  ((input.var2_cat !== 'undefined' && input.var2_cat.length > 0) ||
                                  (input.time_cat !== 'undefined' && input.time_cat.length > 0))
                                  )",
                                  
                   h2("Levels"),
                   verbatimTextOutput('summ.cat'),
                   h2("Contingency table"),
                   verbatimTextOutput('cont.cat'),
                   conditionalPanel(condition = "input.test_cat == 'chisq'",
                                    h2("Expected contingency table"),
                                    verbatimTextOutput('econt.cat')
                   ),
                   h2("Hypothesis testing"),
                   verbatimTextOutput('out.cat')
                 
                 )
        ),
        tabPanel('Correlation',
                 conditionalPanel(condition = "input.var_cor !== 'undefined' && input.var_cor.length > 1",
                                  h2("Q-Q plot"),
                                  plotOutput('qqgraph.cor'),
                                  h2("Sharpiro-Wilk's normality test"),
                                  verbatimTextOutput('swtest.cor'),
                                  
                                  h2("Correlation plot"),
                                  plotOutput('plot.cor'),
                                  h2("Transformation"),
                                  tableOutput('trans.cor'),
                                  h2("Summary"),
                                  verbatimTextOutput('out.cor')
                 )
        ),
        tabPanel('Regression',
                 conditionalPanel(condition = "input.indep_reg !== 'undefined' && input.indep_reg.length > 0",
                                  
                                  conditionalPanel(condition = "input.model_reg === 'linear'",
                                                   h2("Q-Q plot"),
                                                   plotOutput('qqgraph.reg'),
                                                   h2("Sharpiro-Wilk's normality test"),
                                                   verbatimTextOutput('swtest.reg'),
                                                   conditionalPanel(condition = "input.indep_reg !== 'undefined' && input.indep_reg.length === 1",
                                                    h2("Regression plot"),
                                                    plotOutput('regression.reg')
                                                   )
                                   ),
                                  h2("Diagnostic plots"),
                                  plotOutput('plot.reg'),
                                  h2("Summary"),
                                  verbatimTextOutput('out.reg')
                                  
                 )
        ),
        tabPanel('Survival Curve',
                 conditionalPanel(condition = "input.t_surv !== 'undefined' && input.t_surv.length > 0 &&
                                  input.s_surv !== 'undefined' && input.s_surv.length > 0",
                                  
                   conditionalPanel(condition = "input.group_surv == 'yes'",
                                    h2("Summary of the group variable"),
                                    verbatimTextOutput('summ.g.surv')
                   ),
                   h2("Kaplan-Meier curve"),
                   plotOutput('plot.surv'),
                   conditionalPanel(condition = "input.group_surv == 'no'",
                                    h2("Survival rate and Median survival"),
                                    verbatimTextOutput('out.surv')
                   ),
                   conditionalPanel(condition = "input.group_surv == 'yes'",
                                    h2("Survival rate, Median survival, and Hazard ratio"),
                                    verbatimTextOutput('out2.surv')
                   )
                 )
        )
      )
    )
  )
#)

#shinyServer(
server = 
  function(input, output, session) {
    
    #reactive({
    #  #options(warn=-1) or 0
    #  options(warn=as.numeric(input$warn))
    #})
      
    data <- reactive({
      study.data <- NULL
      if(input$rname.data == ""){
        inFile <- input$file1 
        if (is.null(inFile)){return(NULL)} 
        study.data <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                               quote=input$quote)
      }else{
        tname = strsplit(input$rname.data,"::")[[1]]
        tpos.pkg = which(rdlist$results[,1]==tname[1])
        tpos.dat = which(rdlist$results[,3]==tname[2])
        if(length(tpos.pkg)>0 & length(tpos.dat)>0){
          study.data <- eval(parse(text=input$rname.data))
        }
      }
      study.data
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var1_var",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "group.var",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var1_ttest",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "pair2.ttest",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var2.ttest",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "group2.ttest",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var1_anova",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var2_anova",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var1_tanova",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var2_tanova",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var3_tanova",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var1_tmg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var2_tmg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var3_tmg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var4_tmg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var1_cat",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var2_cat",
        choices=c("",names(data())))
      
    }) 
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "time_cat",
        choices=c("",names(data())))
      
    }) 
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "var_cor",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "log.cor",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "sqr.cor",
        choices=c("",names(data())))
      
    })
    
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "res.reg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "logit.reg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "ct.reg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "cs.reg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "indep_reg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "indepcat_reg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "scale.reg",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "t_surv",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "s_surv",
        choices=c("",names(data())))
      
    })
    
    # Updata value user could select
    observe({
      updateSelectInput(
        session,
        "g.surv",
        choices=c("",names(data())))
      
    })
    
    # Output a t distribution plot
    output$tplot <- renderPlot({
      # Display the Student's t distributions with various
      # degrees of freedom and compare to the normal distribution
      
      x <- seq(input$range[1], input$range[2], length=100)
      hx <- dnorm(x)  
      labels <- c('t distribution', 'normal distribution')
      
      plot(x, hx, type="l", lty=2, xlab="x value", col = "black",
           ylab="Density", main="t Distributions")
      
      
      lines(x, dt(x,input$df), lwd=2, col="red")
      
      
      legend("topright", inset=.05, title="Distributions",
             labels, lwd=2, lty=c(1, 2), col=c("red", "black"))
    })
    
    # Output a data table for the upload tab page
    output$contents <- renderTable({
      inFile <- input$file1 
      if (is.null(inFile))
        return(NULL)
      read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
      
    })
    
    #####
    # Data View
    
    # list all R data
    output$rdata.data <- DT::renderDataTable(
      #print(dim(df))
      DT::datatable(rdlist$results[,c(1,3,4)], options = list(pageLength = 10),editable=FALSE)
    )
    
    # display 10 rows initially
    output$tab.data <- DT::renderDataTable(
      #print(dim(df))
      DT::datatable(data(), options = list(pageLength = 10),editable=FALSE)
    )
    output$summ.data <- renderPrint({Hmisc::describe(data())})
    
    output$desc.data <- renderText(
      extract_help(input$rname.data, input$file1, to="html")
      #extract_help("ggplot2", "diamonds", to="html")
    )
    
    #####
    # Variable View
    
    # Output a histogram for the variables user chose
    output$hisgraph.var <- renderPlot({
      if(input$isgroup_var=="no"){
        var1 <- data()[,input$var1_var]
        if (is.null(var1)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        p1 <- hist(new.var1, breaks = input$bins.var)
        plot(p1, col=rgb(0,0,1,1/4))
      }else{
        var1 <- data.frame(v=data()[,input$var1_var],g=data()[,input$group.var])
        if (is.null(var1$v)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          tpos <- which(!is.na(var1$v) & var1$v>0)
          new.var1 <- var1[tpos,]
          new.var1$v <- log(new.var1$v)
        }else{
          new.var1 <- var1
          new.var1$v <- sqrt(new.var1$v)
        }
        
        new.var1$g <- factor(new.var1$g)
        wdata <- new.var1
        
        nlev <- length(levels(wdata$g))
        
        # change fill and outline color manually 
        ggplot(wdata, aes(x = v)) +
          geom_histogram(aes(color = g, fill = g), 
                         position = "identity", bins = input$bins.var, alpha = 0.4) +
          scale_color_manual(values = 1:nlev) + #c("#00AFBB", "#E7B800")) +
          scale_fill_manual(values = 1:nlev) #c("#00AFBB", "#E7B800"))
      }
    }) 
    
    # Output a density plot for the variables user chose
    output$dengraph.var <- renderPlot({
      if(input$isgroup_var=="no"){
        var1 <- data()[,input$var1_var]
        if (is.null(var1)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        p1 <- density(new.var1, na.rm=T)
        plot(p1)
      }else{
        var1 <- data.frame(v=data()[,input$var1_var],g=data()[,input$group.var])
        if (is.null(var1$v)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          tpos <- which(!is.na(var1$v) & var1$v>0)
          new.var1 <- var1[tpos,]
          new.var1$v <- log(new.var1$v)
        }else{
          new.var1 <- var1
          new.var1$v <- sqrt(new.var1$v)
        }
        
        new.var1$g <- factor(new.var1$g)
        wdata <- new.var1
        
        nlev <- length(levels(wdata$g))
        
        ggplot(wdata, aes(x = v, colour = g)) +
          geom_density()      
      }
    })    
    
    # Output a box plot for the variables user chose
    output$boxgraph.var <- renderPlot({
      if(input$isgroup_var=="no"){
        var1 <- na.omit(data()[,input$var1_var])
        if (is.null(var1)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- as.numeric(var1)
        }else if(input$transform.var == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        #boxplot(new.var1)
        
        wdata = data.frame(val=new.var1)
        
        ggplot(wdata, aes(x="",y=val)) +
          geom_boxplot() +
          xlab("") +
          ylab(input$var1_var) +
          geom_dotplot(binaxis='y', stackdir='center') + #,binwidth = 3) +
          #geom_hline(yintercept = 75,colour="red") +
          theme_classic()
        
      }else{
        var1 <- na.omit(data.frame(v=data()[,input$var1_var],g=data()[,input$group.var]))
        if (is.null(var1$v)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          tpos <- which(!is.na(var1$v) & var1$v>0)
          new.var1 <- var1[tpos,]
          new.var1$v <- log(new.var1$v)
        }else{
          new.var1 <- var1
          new.var1$v <- sqrt(new.var1$v)
        }
        
        new.var1$g <- factor(new.var1$g)
        wdata <- new.var1
        
        nlev <- length(levels(wdata$g))
        
        #bp <- ggplot(wdata, aes(x=g, y=v, fill=g)) + 
        #  geom_boxplot()#+
        #labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")
        #bp + theme_classic()
        
        ggplot(wdata, aes(x=g,y=v,fill=g)) +
          geom_boxplot() + #position=position_dodge(.1)) +
          xlab(input$group.var) +
          ylab(input$var1_var) +
          geom_dotplot(binaxis='y', stackdir='center') + #,binwidth = 3)+ 
          guides(fill=guide_legend(title=input$group.var)) +
          #geom_hline(yintercept = 75,colour="red") +
          theme_classic()
        
      }
    })    
    
    # Output a qq plot for the variables user chose
    output$qqgraph.var <- renderPlot({
      if(input$isgroup_var=="no"){
        var1 <- data()[,input$var1_var]
        if (is.null(var1)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        car::qqPlot(new.var1)
      }else{
        var1 <- data.frame(v=data()[,input$var1_var],g=data()[,input$group.var])
        if (is.null(var1$v)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          tpos <- which(!is.na(var1$v) & var1$v>0)
          new.var1 <- var1[tpos,]
          new.var1$v <- log(new.var1$v)
        }else{
          new.var1 <- var1
          new.var1$v <- sqrt(new.var1$v)
        }
        
        new.var1$g <- factor(new.var1$g)
        wdata <- new.var1
        
        nlev <- length(levels(wdata$g))
        
        gg <- ggplot(data = wdata, mapping = aes(sample = v, color = g, fill = g)) +
          stat_qq_band(alpha=0.5) +
          stat_qq_line() +
          stat_qq_point() +
          facet_wrap(~ g) +
          labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
        gg
      }
    })   
    
    # sw normality test
    output$swtest.var <- renderPrint({
      if(input$isgroup_var=="no"){
        var1 <- data()[,input$var1_var]
        if (is.null(var1)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        shapiro.test(new.var1)
      }else{
        var1 <- data.frame(v=data()[,input$var1_var],g=data()[,input$group.var])
        if (is.null(var1$v)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          tpos <- which(!is.na(var1$v) & var1$v>0)
          new.var1 <- var1[tpos,]
          new.var1$v <- log(new.var1$v)
        }else{
          new.var1 <- var1
          new.var1$v <- sqrt(new.var1$v)
        }
        
        new.var1$g <- factor(new.var1$g)
        wdata <- new.var1
        
        tapply(wdata$v, wdata$g, shapiro.test)
        
      }
    })
    
    # summary
    output$summ.var <- renderPrint({
      if(input$isgroup_var=="no"){
        var1 <- data()[,input$var1_var]
        if (is.null(var1)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        psych::describe(new.var1)
      }else{
        var1 <- data.frame(v=data()[,input$var1_var],g=data()[,input$group.var])
        if (is.null(var1$v)){return(NULL)}
        if(input$transform.var == 'none'){
          new.var1 <- var1
        }else if(input$transform.var == 'logtrans'){
          tpos <- which(!is.na(var1$v) & var1$v>0)
          new.var1 <- var1[tpos,]
          new.var1$v <- log(new.var1$v)
        }else{
          new.var1 <- var1
          new.var1$v <- sqrt(new.var1$v)
        }
        
        new.var1$g <- factor(new.var1$g)
        wdata <- new.var1
        
        tapply(wdata$v, wdata$g, psych::describe)
        
      }
    })
    
    ########
    # T-test
    
    # Output a box plot for the variables user chose
    output$boxgraph.ttest <- renderPlot({
      if(input$sample_ttest == "oneSamp"){
        var1 <- na.omit(data()[,input$var1_ttest])
        if (is.null(var1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        #boxplot(new.var1)
        #abline(h=input$null.ttest,lty=2,col="blue")
        
        wdata = data.frame(val=new.var1)
        
        ggplot(wdata, aes(x="",y=val)) +
          geom_boxplot() +
          xlab("") +
          ylab(input$var1_ttest) +
          geom_dotplot(binaxis='y', stackdir='center') + #,binwidth = 3) +
          geom_hline(yintercept = input$null.ttest,colour="red") +
          theme_classic()
        
      }else if(input$sample_ttest == "pairSamp"){
        var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$pair2.ttest])
        if (is.null(var1$v1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          new.var1 <- var1[tpos,]
          new.var1$v1 <- log(new.var1$v1)
          new.var1$v2 <- log(new.var1$v2)
        }else{
          new.var1 <- var1
          new.var1$v1 <- sqrt(new.var1$v1)
          new.var1$v2 <- sqrt(new.var1$v2)
        }
        
        tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$pair2.ttest,dim(new.var1)[1]))
        wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg),paired=c(1:dim(new.var1)[1],1:dim(new.var1)[1]))
        
        nlev <- length(levels(wdata$g))
        
        #bp <- ggplot(wdata, aes(x=g, y=v, fill=g)) + 
        #  geom_boxplot()#+
        #labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")
        #bp + theme_classic()
        
        ggplot(wdata,aes(x=g,y=v,fill=g)) +
          geom_boxplot() +
          xlab("") +
          ylab("") +
          geom_line(aes(group=paired), position = position_dodge(0.2)) +
          geom_point(aes(fill=g,group=paired), position = position_dodge(0.2)) +
          theme_classic()
        
      }else{
        if(input$isgroup_ttest == "no"){
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          #if (is.null(var1$v1)){return(NULL)}
          #if(input$transform.ttest == 'none'){
          #  new.var1 <- var1
          #}else if(input$transform.ttest == 'logtrans'){
          #  tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          #  new.var1 <- var1[tpos,]
          #  new.var1$v1 <- log(new.var1$v1)
          #  new.var1$v2 <- log(new.var1$v2)
          #}else{
          #  new.var1 <- var1
          #  new.var1$v1 <- sqrt(new.var1$v1)
          #  new.var1$v2 <- sqrt(new.var1$v2)
          #}
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          
          #nlev <- length(levels(wdata$g))
          
          #bp <- ggplot(wdata, aes(x=g, y=v, fill=g)) + 
          #  geom_boxplot()#+
          #labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")
          #bp + theme_classic()
          
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          var1 <- na.omit(data()[,input$var1_ttest])
          var2 <- na.omit(data()[,input$var2.ttest])
          if (is.null(var1) | is.null(var2)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
            new.var2 <- var2
          }else if(input$transform.ttest == 'logtrans'){
            tpos1 <- which(!is.na(var1) & var1>0)
            tpos2 <- which(!is.na(var2) & var2>0)
            new.var1 <- log(var1[tpos1])
            new.var2 <- log(var2[tpos2])
          }else{
            tpos1 <- which(!is.na(var1) & var1>=0)
            tpos2 <- which(!is.na(var2) & var2>=0)
            new.var1 <- sqrt(var1[tpos1])
            new.var2 <- sqrt(var2[tpos2])
          }
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          tg <- c(rep("g1",length(new.var1)),rep("g2",length(new.var2)))
          wdata <- data.frame(v=c(new.var1,new.var2),g=factor(tg))
          
          ggplot(wdata, aes(x=g,y=v,fill=g)) +
            geom_boxplot() + #position=position_dodge(.1)) +
            xlab("Group") +
            ylab(input$var1_ttest) +
            geom_dotplot(binaxis='y', stackdir='center') + #,binwidth = input$binwidth.ttest)+ 
            guides(fill=guide_legend(title="Group")) +
            theme_classic()
        }else{
          var1 <- data.frame(v=data()[,input$var1_ttest],g=data()[,input$group2.ttest])
          if (is.null(var1$v)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
          }else if(input$transform.ttest == 'logtrans'){
            tpos <- which(!is.na(var1$v) & var1$v>0)
            new.var1 <- var1[tpos,]
            new.var1$v <- log(new.var1$v)
          }else{
            new.var1 <- var1
            new.var1$v <- sqrt(new.var1$v)
          }
          
          new.var1$g <- factor(new.var1$g)
          wdata <- new.var1
          
          nlev <- length(levels(wdata$g))
          
          #bp <- ggplot(wdata, aes(x=g, y=v, fill=g)) + 
          #  geom_boxplot()#+
          #labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")
          #bp + theme_classic()
          
          ggplot(wdata, aes(x=g,y=v,fill=g)) +
            geom_boxplot() + #position=position_dodge(.1)) +
            xlab(input$group2.ttest) +
            ylab(input$var1_ttest) +
            geom_dotplot(binaxis='y', stackdir='center') + #,binwidth = input$binwidth.ttest)+ 
            guides(fill=guide_legend(title=input$group2.ttest)) +
            theme_classic()
          
        }
      }
    })    
    
    
    # qq plot: t-test
    output$qqgraph.ttest <- renderPlot({
      if(input$sample_ttest == "oneSamp"){
        var1 <- na.omit(data()[,input$var1_ttest])
        if (is.null(var1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        #boxplot(new.var1)
        #abline(h=input$null.ttest,lty=2,col="blue")
        
        wdata = data.frame(val=new.var1)

        car::qqPlot(wdata$val,xlab="",ylab=input$var1_ttest)
        
      }else if(input$sample_ttest == "pairSamp"){
        var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$pair2.ttest])
        if (is.null(var1$v1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          new.var1 <- var1[tpos,]
          new.var1$v1 <- log(new.var1$v1)
          new.var1$v2 <- log(new.var1$v2)
        }else{
          new.var1 <- var1
          new.var1$v1 <- sqrt(new.var1$v1)
          new.var1$v2 <- sqrt(new.var1$v2)
        }
        
        #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$pair2.ttest,dim(new.var1)[1]))
        #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg),paired=c(1:dim(new.var1)[1],1:dim(new.var1)[1]))
        
        diff.wdata <- na.omit(new.var1$v2 - new.var1$v1)
        
        #nlev <- length(levels(wdata$g))
        
        #bp <- ggplot(wdata, aes(x=g, y=v, fill=g)) + 
        #  geom_boxplot()#+
        #labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")
        #bp + theme_classic()
        
        car::qqPlot(diff.wdata,xlab="",ylab=paste0(input$pair2.ttest,"-",input$var1_ttest))
        
      }else{
        if(input$isgroup_ttest == "no"){
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          #if (is.null(var1$v1)){return(NULL)}
          #if(input$transform.ttest == 'none'){
          #  new.var1 <- var1
          #}else if(input$transform.ttest == 'logtrans'){
          #  tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          #  new.var1 <- var1[tpos,]
          #  new.var1$v1 <- log(new.var1$v1)
          #  new.var1$v2 <- log(new.var1$v2)
          #}else{
          #  new.var1 <- var1
          #  new.var1$v1 <- sqrt(new.var1$v1)
          #  new.var1$v2 <- sqrt(new.var1$v2)
          #}
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          
          #nlev <- length(levels(wdata$g))
          
          #bp <- ggplot(wdata, aes(x=g, y=v, fill=g)) + 
          #  geom_boxplot()#+
          #labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")
          #bp + theme_classic()
          
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          var1 <- na.omit(data()[,input$var1_ttest])
          var2 <- na.omit(data()[,input$var2.ttest])
          if (is.null(var1) | is.null(var2)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
            new.var2 <- var2
          }else if(input$transform.ttest == 'logtrans'){
            tpos1 <- which(!is.na(var1) & var1>0)
            tpos2 <- which(!is.na(var2) & var2>0)
            new.var1 <- log(var1[tpos1])
            new.var2 <- log(var2[tpos2])
          }else{
            tpos1 <- which(!is.na(var1) & var1>=0)
            tpos2 <- which(!is.na(var2) & var2>=0)
            new.var1 <- sqrt(var1[tpos1])
            new.var2 <- sqrt(var2[tpos2])
          }
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          tg <- c(rep("g1",length(new.var1)),rep("g2",length(new.var2)))
          wdata <- data.frame(v=c(new.var1,new.var2),g=factor(tg))
          
          gg <- ggplot(data = wdata, mapping = aes(sample = v, color = g, fill = g)) +
            stat_qq_band(alpha=0.5) +
            stat_qq_line() +
            stat_qq_point() +
            facet_wrap(~ g) +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          gg
          
        }else{
          var1 <- data.frame(v=data()[,input$var1_ttest],g=data()[,input$group2.ttest])
          if (is.null(var1$v)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
          }else if(input$transform.ttest == 'logtrans'){
            tpos <- which(!is.na(var1$v) & var1$v>0)
            new.var1 <- var1[tpos,]
            new.var1$v <- log(new.var1$v)
          }else{
            new.var1 <- var1
            new.var1$v <- sqrt(new.var1$v)
          }
          
          new.var1$g <- factor(new.var1$g)
          wdata <- new.var1
          
          nlev <- length(levels(wdata$g))
          
          #bp <- ggplot(wdata, aes(x=g, y=v, fill=g)) + 
          #  geom_boxplot()#+
          #labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")
          #bp + theme_classic()
 
          gg <- ggplot(data = wdata, mapping = aes(sample = v, color = g, fill = g)) +
            stat_qq_band(alpha=0.5) +
            stat_qq_line() +
            stat_qq_point() +
            facet_wrap(~ g) +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          gg
          
        }
      }
 
    })
    
    # shapiro test: t-test
    output$swtest.ttest <- renderPrint({
      if(input$sample_ttest == "oneSamp"){
        var1 <- na.omit(data()[,input$var1_ttest])
        if (is.null(var1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        #boxplot(new.var1)
        #abline(h=input$null.ttest,lty=2,col="blue")
        
        wdata = data.frame(val=new.var1)
        
        shapiro.test(wdata$val)
      }else if(input$sample_ttest == "pairSamp"){
        var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$pair2.ttest])
        if (is.null(var1$v1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          new.var1 <- var1[tpos,]
          new.var1$v1 <- log(new.var1$v1)
          new.var1$v2 <- log(new.var1$v2)
        }else{
          new.var1 <- var1
          new.var1$v1 <- sqrt(new.var1$v1)
          new.var1$v2 <- sqrt(new.var1$v2)
        }
        
        #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$pair2.ttest,dim(new.var1)[1]))
        #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg),paired=c(1:dim(new.var1)[1],1:dim(new.var1)[1]))
        
        diff.wdata <- na.omit(new.var1$v2 - new.var1$v1)

        shapiro.test(diff.wdata)
      }else{
        if(input$isgroup_ttest == "no"){
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          #if (is.null(var1$v1)){return(NULL)}
          #if(input$transform.ttest == 'none'){
          #  new.var1 <- var1
          #}else if(input$transform.ttest == 'logtrans'){
          #  tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          #  new.var1 <- var1[tpos,]
          #  new.var1$v1 <- log(new.var1$v1)
          #  new.var1$v2 <- log(new.var1$v2)
          #}else{
          #  new.var1 <- var1
          #  new.var1$v1 <- sqrt(new.var1$v1)
          #  new.var1$v2 <- sqrt(new.var1$v2)
          #}
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          
          #nlev <- length(levels(wdata$g))
          
          #bp <- ggplot(wdata, aes(x=g, y=v, fill=g)) + 
          #  geom_boxplot()#+
          #labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")
          #bp + theme_classic()
          
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          var1 <- na.omit(data()[,input$var1_ttest])
          var2 <- na.omit(data()[,input$var2.ttest])
          if (is.null(var1) | is.null(var2)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
            new.var2 <- var2
          }else if(input$transform.ttest == 'logtrans'){
            tpos1 <- which(!is.na(var1) & var1>0)
            tpos2 <- which(!is.na(var2) & var2>0)
            new.var1 <- log(var1[tpos1])
            new.var2 <- log(var2[tpos2])
          }else{
            tpos1 <- which(!is.na(var1) & var1>=0)
            tpos2 <- which(!is.na(var2) & var2>=0)
            new.var1 <- sqrt(var1[tpos1])
            new.var2 <- sqrt(var2[tpos2])
          }
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          tg <- c(rep("g1",length(new.var1)),rep("g2",length(new.var2)))
          wdata <- data.frame(v=c(new.var1,new.var2),g=factor(tg))

          tapply(wdata$v, wdata$g, shapiro.test)
        }else{
          var1 <- data.frame(v=data()[,input$var1_ttest],g=data()[,input$group2.ttest])
          if (is.null(var1$v)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
          }else if(input$transform.ttest == 'logtrans'){
            tpos <- which(!is.na(var1$v) & var1$v>0)
            new.var1 <- var1[tpos,]
            new.var1$v <- log(new.var1$v)
          }else{
            new.var1 <- var1
            new.var1$v <- sqrt(new.var1$v)
          }
          
          new.var1$g <- factor(new.var1$g)
          wdata <- new.var1

          tapply(wdata$v, wdata$g, shapiro.test)
        }
      }
    })
    
    
    # summary
    output$summ.ttest <- renderPrint({
      if(input$sample_ttest == "oneSamp"){
        var1 <- data()[,input$var1_ttest]
        if (is.null(var1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        psych::describe(new.var1)
      }else if(input$sample_ttest == "pairSamp"){
        var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$pair2.ttest])
        if (is.null(var1$v1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          new.var1 <- var1[tpos,]
          new.var1$v1 <- log(new.var1$v1)
          new.var1$v2 <- log(new.var1$v2)
        }else{
          new.var1 <- var1
          new.var1$v1 <- sqrt(new.var1$v1)
          new.var1$v2 <- sqrt(new.var1$v2)
        }
        
        tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$pair2.ttest,dim(new.var1)[1]))
        wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
        
        nlev <- length(levels(wdata$g))
        
        tapply(wdata$v, wdata$g, psych::describe)
      }else{
        if(input$isgroup_ttest == "no"){
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          #if (is.null(var1$v1)){return(NULL)}
          #if(input$transform.ttest == 'none'){
          #  new.var1 <- var1
          #}else if(input$transform.ttest == 'logtrans'){
          #  tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          #  new.var1 <- var1[tpos,]
          #  new.var1$v1 <- log(new.var1$v1)
          #  new.var1$v2 <- log(new.var1$v2)
          #}else{
          #  new.var1 <- var1
          #  new.var1$v1 <- sqrt(new.var1$v1)
          #  new.var1$v2 <- sqrt(new.var1$v2)
          #}
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          
          #nlev <- length(levels(wdata$g))
          
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          var1 <- na.omit(data()[,input$var1_ttest])
          var2 <- na.omit(data()[,input$var2.ttest])
          if (is.null(var1) | is.null(var2)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
            new.var2 <- var2
          }else if(input$transform.ttest == 'logtrans'){
            tpos1 <- which(!is.na(var1) & var1>0)
            tpos2 <- which(!is.na(var2) & var2>0)
            new.var1 <- log(var1[tpos1])
            new.var2 <- log(var2[tpos2])
          }else{
            tpos1 <- which(!is.na(var1) & var1>=0)
            tpos2 <- which(!is.na(var2) & var2>=0)
            new.var1 <- sqrt(var1[tpos1])
            new.var2 <- sqrt(var2[tpos2])
          }
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          tg <- c(rep("g1",length(new.var1)),rep("g2",length(new.var2)))
          wdata <- data.frame(v=c(new.var1,new.var2),g=factor(tg))
          
          tapply(wdata$v, wdata$g, psych::describe)
        }else{
          var1 <- data.frame(v=data()[,input$var1_ttest],g=data()[,input$group2.ttest])
          if (is.null(var1$v)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
          }else if(input$transform.ttest == 'logtrans'){
            tpos <- which(!is.na(var1$v) & var1$v>0)
            new.var1 <- var1[tpos,]
            new.var1$v <- log(new.var1$v)
          }else{
            new.var1 <- var1
            new.var1$v <- sqrt(new.var1$v)
          }
          
          new.var1$g <- factor(new.var1$g)
          wdata <- new.var1
          
          nlev <- length(levels(wdata$g))
          
          tapply(wdata$v, wdata$g, psych::describe)
        }
      }
    })    
    
    # equal variance test
    output$levene.ttest <- renderPrint({
      if(input$sample_ttest != "twoSamp"){
        return(NULL)
      }else{
        if(input$isgroup_ttest == "no"){
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          #if (is.null(var1$v1)){return(NULL)}
          #if(input$transform.ttest == 'none'){
          #  new.var1 <- var1
          #}else if(input$transform.ttest == 'logtrans'){
          #  tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          #  new.var1 <- var1[tpos,]
          #  new.var1$v1 <- log(new.var1$v1)
          #  new.var1$v2 <- log(new.var1$v2)
          #}else{
          #  new.var1 <- var1
          #  new.var1$v1 <- sqrt(new.var1$v1)
          #  new.var1$v2 <- sqrt(new.var1$v2)
          #}
          
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          var1 <- na.omit(data()[,input$var1_ttest])
          var2 <- na.omit(data()[,input$var2.ttest])
          if (is.null(var1) | is.null(var2)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
            new.var2 <- var2
          }else if(input$transform.ttest == 'logtrans'){
            tpos1 <- which(!is.na(var1) & var1>0)
            tpos2 <- which(!is.na(var2) & var2>0)
            new.var1 <- log(var1[tpos1])
            new.var2 <- log(var2[tpos2])
          }else{
            tpos1 <- which(!is.na(var1) & var1>=0)
            tpos2 <- which(!is.na(var2) & var2>=0)
            new.var1 <- sqrt(var1[tpos1])
            new.var2 <- sqrt(var2[tpos2])
          }
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          tg <- c(rep("g1",length(new.var1)),rep("g2",length(new.var2)))
          wdata <- data.frame(v=c(new.var1,new.var2),g=factor(tg))
          
          car::leveneTest(v~g,data=wdata)
        }else{
          var1 <- data.frame(v=data()[,input$var1_ttest],g=data()[,input$group2.ttest])
          if (is.null(var1$v)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
          }else if(input$transform.ttest == 'logtrans'){
            tpos <- which(!is.na(var1$v) & var1$v>0)
            new.var1 <- var1[tpos,]
            new.var1$v <- log(new.var1$v)
          }else{
            new.var1 <- var1
            new.var1$v <- sqrt(new.var1$v)
          }
          
          new.var1$g <- factor(new.var1$g)
          wdata <- new.var1
          
          car::leveneTest(v~g,data=wdata)
        }
      }
    })    
    
    # t-test
    output$out.ttest <- renderPrint({
      
      conf <- input$conf.ttest
      
      if(input$sample_ttest == "oneSamp"){
        var1 <- data()[,input$var1_ttest]
        if (is.null(var1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          new.var1 <- log(var1[which(!is.na(var1) & var1>0)])
        }else{
          new.var1 <- sqrt(var1)
        }
        if(input$nonparm_ttest=="no"){
          t.test(new.var1, alternative = input$tail.ttest, mu = input$null.ttest, conf.level = conf)
        }else{
          wilcox.test(new.var1, alternative = input$tail.ttest, mu = input$null.ttest, conf.level = conf)
        }
      }else if(input$sample_ttest == "pairSamp"){
        var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$pair2.ttest])
        if (is.null(var1$v1)){return(NULL)}
        if(input$transform.ttest == 'none'){
          new.var1 <- var1
        }else if(input$transform.ttest == 'logtrans'){
          tpos <- which(!is.na(var1$v1) & var1$v1>0 & !is.na(var1$v2) & var1$v2>0)
          new.var1 <- var1[tpos,]
          new.var1$v1 <- log(new.var1$v1)
          new.var1$v2 <- log(new.var1$v2)
        }else{
          new.var1 <- var1
          new.var1$v1 <- sqrt(new.var1$v1)
          new.var1$v2 <- sqrt(new.var1$v2)
        }
        
        tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$pair2.ttest,dim(new.var1)[1]))
        wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
        
        nlev <- length(levels(wdata$g))
        
        #t.test(v~g, data = wdata, alternative = input$tail.ttest, paired = TRUE, conf.level = conf)
        if(input$nonparm_ttest=="no"){
          t.test(new.var1$v1,new.var1$v2, alternative = input$tail.ttest, paired = TRUE, conf.level = conf)
        }else{
          wilcox.test(new.var1$v1,new.var1$v2, alternative = input$tail.ttest, paired = TRUE, conf.level = conf)
        }
      }else{
        ve <- ifelse(input$varequal.ttest == 'yes', TRUE, FALSE)
        
        if(input$isgroup_ttest == "no"){
          #var1 <- data.frame(v1=data()[,input$var1_ttest],v2=data()[,input$var2.ttest])
          var1 <- na.omit(data()[,input$var1_ttest])
          var2 <- na.omit(data()[,input$var2.ttest])
          if (is.null(var1) | is.null(var2)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
            new.var2 <- var2
          }else if(input$transform.ttest == 'logtrans'){
            tpos1 <- which(!is.na(var1) & var1>0)
            tpos2 <- which(!is.na(var2) & var2>0)
            new.var1 <- log(var1[tpos1])
            new.var2 <- log(var2[tpos2])
          }else{
            tpos1 <- which(!is.na(var1) & var1>=0)
            tpos2 <- which(!is.na(var2) & var2>=0)
            new.var1 <- sqrt(var1[tpos1])
            new.var2 <- sqrt(var2[tpos2])
          }
          
          #tg <- c(rep(input$var1_ttest,dim(new.var1)[1]),rep(input$var2.ttest,dim(new.var1)[1]))
          #wdata <- data.frame(v=c(new.var1$v1,new.var1$v2),g=factor(tg))
          
          #nlev <- length(levels(wdata$g))
          
          #t.test(v~g, data = wdata, alternative = input$tail.ttest, var.equal = ve, conf.level = conf)
          if(input$nonparm_ttest=="no"){
            t.test(new.var1,new.var2, alternative = input$tail.ttest, var.equal = ve, conf.level = conf)
          }else{
            wilcox.test(new.var1,new.var2, alternative = input$tail.ttest, var.equal = ve, conf.level = conf)
          }
        }else{
          var1 <- data.frame(v=data()[,input$var1_ttest],g=data()[,input$group2.ttest])
          if (is.null(var1$v)){return(NULL)}
          if(input$transform.ttest == 'none'){
            new.var1 <- var1
          }else if(input$transform.ttest == 'logtrans'){
            tpos <- which(!is.na(var1$v) & var1$v>0)
            new.var1 <- var1[tpos,]
            new.var1$v <- log(new.var1$v)
          }else{
            new.var1 <- var1
            new.var1$v <- sqrt(new.var1$v)
          }
          
          new.var1$g <- factor(new.var1$g) #,levels=c(input$var1_ttest,input$var2.ttest))
          wdata <- new.var1
          
          nlev <- length(levels(wdata$g))
          
          if(input$nonparm_ttest=="no"){
            t.test(v~g, data = wdata, alternative = input$tail.ttest, var.equal = ve, conf.level = conf)
          }else{
            wilcox.test(v~g, data = wdata, alternative = input$tail.ttest, var.equal = ve, conf.level = conf)
          }
        }
      }
    })
    
    #############
    # One-Way ANOVA
    
    # ggpubr::ggline(aa,x="gr",y="aa",add=c("mean_se","jitter"),error.plot='pointrange',color="indianred)
    
    # normality
    # red.res = aov(v~G,data=wdata)
    # summary(red.res)
    # res = residuals(red.res)
    # library(car)
    # qqPlot(res)
    # shapiro.test(res)
    
    # equal variance
    # plot(red.res,1)
    # bartlett.test(v~g,data=wdata)
    
    # post-hoc analysis
    # pairwise.t.test(wdata$v,wdata$g,p.adj="holm")
    
    # uneuqal variance
    # oneway.test(v~g,data=wdata)
    # pairwise.t.test(wdata$v,wdata$g,p.adj="holm",pool.sd=FALSE)
    
    # box plot
    output$boxgraph.anova <- renderPlot({
      #var1 <- data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova])
      var1 <- na.omit(data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova]))
      if (is.null(var1$v)){return(NULL)}
      if(input$transform.anova == 'none'){
        new.var1 <- var1
      }else if(input$transform.anova == 'logtrans'){
        tpos <- which(!is.na(var1$v) & var1$v>0)
        new.var1 <- var1[tpos,]
        new.var1$v <- log(new.var1$v)
      }else{
        new.var1 <- var1
        new.var1$v <- sqrt(new.var1$v)
      }
      
      new.var1$g <- factor(new.var1$g)
      wdata <- new.var1
      
      nlev <- length(levels(wdata$g))
      
      bp <- ggplot(wdata, aes(x=g, y=v, fill=g)) + 
        geom_boxplot()+
        labs(
          #title="Plot of length  per dose",
          y=input$var1_anova
          ,x=input$var2_anova
        )
      bp + theme_classic()
    })
    
    # mean plot
    output$meangraph.anova <- renderPlot({
      #var1 <- data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova])
      var1 <- na.omit(data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova]))
      if (is.null(var1$v)){return(NULL)}
      if(input$transform.anova == 'none'){
        new.var1 <- var1
      }else if(input$transform.anova == 'logtrans'){
        tpos <- which(!is.na(var1$v) & var1$v>0)
        new.var1 <- var1[tpos,]
        new.var1$v <- log(new.var1$v)
      }else{
        new.var1 <- var1
        new.var1$v <- sqrt(new.var1$v)
      }
      
      new.var1$g <- factor(new.var1$g)
      wdata <- new.var1
      
      nlev <- length(levels(wdata$g))
      
      ggpubr::ggline(wdata,x="g",y="v",xlab=input$var2_anova,ylab=input$var1_anova,add=c("mean_se","jitter"),error.plot='pointrange',color="indianred")
    })
    
    # one-way anova
    output$out.anova <- renderPrint({
      #var1 <- data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova])
      var1 <- na.omit(data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova]))
      if (is.null(var1$v)){return(NULL)}
      if(input$transform.anova == 'none'){
        new.var1 <- var1
      }else if(input$transform.anova == 'logtrans'){
        tpos <- which(!is.na(var1$v) & var1$v>0)
        new.var1 <- var1[tpos,]
        new.var1$v <- log(new.var1$v)
      }else{
        new.var1 <- var1
        new.var1$v <- sqrt(new.var1$v)
      }
      
      new.var1$g <- factor(new.var1$g)
      wdata <- new.var1
      
      nlev <- length(levels(wdata$g))
      
      if(input$nonparm_anova == "yes"){
        kruskal.test(v~g,data=wdata)
      }else{
        if(input$varequal.anova == "yes"){
          summary(aov(v~g,data=wdata))
        }else{
          oneway.test(v~g,data=wdata)
        }
      }
    })
    
    # qq plot anova
    output$qqgraph.anova <- renderPlot({
      #var1 <- data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova])
      var1 <- na.omit(data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova]))
      if (is.null(var1$v)){return(NULL)}
      if(input$transform.anova == 'none'){
        new.var1 <- var1
      }else if(input$transform.anova == 'logtrans'){
        tpos <- which(!is.na(var1$v) & var1$v>0)
        new.var1 <- var1[tpos,]
        new.var1$v <- log(new.var1$v)
      }else{
        new.var1 <- var1
        new.var1$v <- sqrt(new.var1$v)
      }
      
      new.var1$g <- factor(new.var1$g)
      wdata <- new.var1
      
      res <- residuals(aov(v~g,data=wdata))
      car::qqPlot(res,xlab="",ylab="")
    })
    
    # shapiro test
    output$swtest.anova <- renderPrint({
      #var1 <- data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova])
      var1 <- na.omit(data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova]))
      if (is.null(var1$v)){return(NULL)}
      if(input$transform.anova == 'none'){
        new.var1 <- var1
      }else if(input$transform.anova == 'logtrans'){
        tpos <- which(!is.na(var1$v) & var1$v>0)
        new.var1 <- var1[tpos,]
        new.var1$v <- log(new.var1$v)
      }else{
        new.var1 <- var1
        new.var1$v <- sqrt(new.var1$v)
      }
      
      new.var1$g <- factor(new.var1$g)
      wdata <- new.var1
      
      res <- residuals(aov(v~g,data=wdata))
      shapiro.test(res)
    })
    
    # residual plot anova
    output$resgraph.anova <- renderPlot({
      #var1 <- data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova])
      var1 <- na.omit(data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova]))
      if (is.null(var1$v)){return(NULL)}
      if(input$transform.anova == 'none'){
        new.var1 <- var1
      }else if(input$transform.anova == 'logtrans'){
        tpos <- which(!is.na(var1$v) & var1$v>0)
        new.var1 <- var1[tpos,]
        new.var1$v <- log(new.var1$v)
      }else{
        new.var1 <- var1
        new.var1$v <- sqrt(new.var1$v)
      }
      
      new.var1$g <- factor(new.var1$g)
      wdata <- new.var1
      red.res <- aov(v~g,data=wdata)
      
      plot(red.res,1)
    })
    
    # bartlett equal variance test
    output$bltest.anova <- renderPrint({
      #var1 <- data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova])
      var1 <- na.omit(data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova]))
      if (is.null(var1$v)){return(NULL)}
      if(input$transform.anova == 'none'){
        new.var1 <- var1
      }else if(input$transform.anova == 'logtrans'){
        tpos <- which(!is.na(var1$v) & var1$v>0)
        new.var1 <- var1[tpos,]
        new.var1$v <- log(new.var1$v)
      }else{
        new.var1 <- var1
        new.var1$v <- sqrt(new.var1$v)
      }
      
      new.var1$g <- factor(new.var1$g)
      wdata <- new.var1
      
      bartlett.test(v~g,data=wdata)
    })
    
    # summary
    output$summ.anova <- renderPrint({
      #var1 <- data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova])
      var1 <- na.omit(data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova]))
      if (is.null(var1$v)){return(NULL)}
      if(input$transform.anova == 'none'){
        new.var1 <- var1
      }else if(input$transform.anova == 'logtrans'){
        tpos <- which(!is.na(var1$v) & var1$v>0)
        new.var1 <- var1[tpos,]
        new.var1$v <- log(new.var1$v)
      }else{
        new.var1 <- var1
        new.var1$v <- sqrt(new.var1$v)
      }
      
      new.var1$g <- factor(new.var1$g)
      wdata <- new.var1
      
      tapply(wdata$v, wdata$g, psych::describe)
    })
    
    # post-hoc analysis
    output$posthoc.anova <- renderPrint({
      #var1 <- data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova])
      var1 <- na.omit(data.frame(v=data()[,input$var1_anova],g=data()[,input$var2_anova]))
      if (is.null(var1$v)){return(NULL)}
      if(input$transform.anova == 'none'){
        new.var1 <- var1
      }else if(input$transform.anova == 'logtrans'){
        tpos <- which(!is.na(var1$v) & var1$v>0)
        new.var1 <- var1[tpos,]
        new.var1$v <- log(new.var1$v)
      }else{
        new.var1 <- var1
        new.var1$v <- sqrt(new.var1$v)
      }
      
      new.var1$g <- factor(new.var1$g)
      wdata <- new.var1
      
      nlev <- length(levels(wdata$g))
      
      if(input$nonparm_anova == "yes"){
        pairwise.t.test(wdata$v,wdata$g,p.adj=input$posthoc.anova)
      }else{
        if(input$varequal.anova == "yes"){
          pairwise.t.test(wdata$v,wdata$g,p.adj=input$posthoc.anova)
        }else{
          pairwise.t.test(wdata$v,wdata$g,p.adj=input$posthoc.anova,pool.sd=FALSE)
        }
      }
    })
    
    ###############
    # Two-Way ANOVA
    
    # interaction plot  
    output$interplot.tanova <- renderPlot({
      var1 <- data()[,c(input$var1_tanova,input$var2_tanova,input$var3_tanova)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tanova == 'none'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tanova == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- log(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- sqrt(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- rank(new.var1[,input$var1_tanova])
      }
      
      new.var1[,input$var2_tanova] <- factor(new.var1[,input$var2_tanova])
      new.var1[,input$var3_tanova] <- factor(new.var1[,input$var3_tanova])
      wdata <- data.frame(response=new.var1[,1],fac1=new.var1[,2],group=new.var1[,3])
      
      interaction.plot(wdata$fac1,wdata$group,wdata$response
                       ,xlab=input$xlab.tanova,ylab=input$ylab.tanova,trace.label=input$tlab.tanova
                       ,col=c(2:(length(levels(wdata$group)))),lwd=2
      )
    })
    
    # box plot
    output$boxplot.tanova <- renderPlot({
      var1 <- data()[,c(input$var1_tanova,input$var2_tanova,input$var3_tanova)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tanova == 'none'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tanova == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- log(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- sqrt(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- rank(new.var1[,input$var1_tanova])
      }
      
      new.var1[,input$var2_tanova] <- factor(new.var1[,input$var2_tanova])
      new.var1[,input$var3_tanova] <- factor(new.var1[,input$var3_tanova])
      wdata <- data.frame(response=new.var1[,1],fac1=new.var1[,2],group=new.var1[,3])
      
      ggplot(wdata, aes(x=fac1, y=response, fill=group)) +
        geom_boxplot(position=position_dodge(input$pdodge.tanova)) +
        xlab(input$xlab.tanova) +
        ylab(input$ylab.tanova) +
        geom_dotplot(binaxis='y', stackdir='center',
                     position=position_dodge(input$pdodge.tanova)) +
        guides(fill=guide_legend(title=input$tlab.tanova))
    })
    
    # two-way anova
    output$out.tanova <- renderPrint({
      var1 <- data()[,c(input$var1_tanova,input$var2_tanova,input$var3_tanova)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tanova == 'none'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tanova == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- log(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- sqrt(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- rank(new.var1[,input$var1_tanova])
      }
      
      new.var1[,input$var2_tanova] <- factor(new.var1[,input$var2_tanova])
      new.var1[,input$var3_tanova] <- factor(new.var1[,input$var3_tanova])
      #wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      
      tform = as.formula(paste0(input$var1_tanova,"~",input$var2_tanova,"*",input$var3_tanova))
      afit = aov(tform,data=new.var1)
      summary(afit)
    })
    
    # two-way anova: pairwise post-hoc
    output$postout.tanova <- renderPrint({
      var1 <- data()[,c(input$var1_tanova,input$var2_tanova,input$var3_tanova)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tanova == 'none'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tanova == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- log(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- sqrt(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- rank(new.var1[,input$var1_tanova])
      }
      
      new.var1[,input$var2_tanova] <- factor(new.var1[,input$var2_tanova])
      new.var1[,input$var3_tanova] <- factor(new.var1[,input$var3_tanova])
      #wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      
      paircomp(input$var1_tanova,input$var2_tanova,input$var3_tanova,new.var1,input$posthoc.tanova)
    })
    
    # two-way anova: qq and density plot
    output$qqdplot.tanova <- renderPlot({
      var1 <- data()[,c(input$var1_tanova,input$var2_tanova,input$var3_tanova)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tanova == 'none'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tanova == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- log(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- sqrt(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- rank(new.var1[,input$var1_tanova])
      }
      
      new.var1[,input$var2_tanova] <- factor(new.var1[,input$var2_tanova])
      new.var1[,input$var3_tanova] <- factor(new.var1[,input$var3_tanova])
      #wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      
      tform = as.formula(paste0(input$var1_tanova,"~",input$var2_tanova,"*",input$var3_tanova))
      afit = aov(tform,data=new.var1)
      
      par(mfrow=c(1,2))
      car::qqPlot(afit$residuals)
      plot(density(afit$residuals),main="")
    })
    
    # two-way anova: shapiro test
    output$swtest.tanova <- renderPrint({
      var1 <- data()[,c(input$var1_tanova,input$var2_tanova,input$var3_tanova)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tanova == 'none'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tanova == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- log(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- sqrt(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- rank(new.var1[,input$var1_tanova])
      }
      
      new.var1[,input$var2_tanova] <- factor(new.var1[,input$var2_tanova])
      new.var1[,input$var3_tanova] <- factor(new.var1[,input$var3_tanova])
      #wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      
      tform = as.formula(paste0(input$var1_tanova,"~",input$var2_tanova,"*",input$var3_tanova))
      afit = aov(tform,data=new.var1)
      
      res <- afit$residuals
      shapiro.test(res)
    })
    
    # two-way anova: residual plot
    output$resplot.tanova <- renderPlot({
      var1 <- data()[,c(input$var1_tanova,input$var2_tanova,input$var3_tanova)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tanova == 'none'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tanova == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- log(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]) & var1[,input$var1_tanova]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- sqrt(new.var1[,input$var1_tanova])
      }else if(input$transform.tanova == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var1_tanova]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var1_tanova] <- rank(new.var1[,input$var1_tanova])
      }
      
      new.var1[,input$var2_tanova] <- factor(new.var1[,input$var2_tanova])
      new.var1[,input$var3_tanova] <- factor(new.var1[,input$var3_tanova])
      #wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      
      tform = as.formula(paste0(input$var1_tanova,"~",input$var2_tanova,"*",input$var3_tanova))
      afit = aov(tform,data=new.var1)
      
      plot(afit,1)
    })
    
    #######################
    # Tumor Growth Analysis
    
    # spider plot  
    output$spiderplot.tmg <- renderPlot({
      var1 <- data()[,c(input$var1_tmg,input$var2_tmg,input$var3_tmg,input$var4_tmg)]
      #if (is.null(var1$v)){return(NULL)}
      tpos <- which(!is.na(var1[,input$var2_tmg]))
      new.var1 <- var1[tpos,]
      
      new.var1[,input$var4_tmg] <- factor(new.var1[,input$var4_tmg])
      new.var1[,input$var3_tmg] <- as.numeric(as.character(new.var1[,input$var3_tmg]))
      wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      
      p=ggplot(wdata, aes(x=time, y=size, group=id)) +
        theme_bw(base_size=14) +
        theme(axis.title.x = element_text(face="bold"), axis.text.x = element_text(face="bold")) +
        theme(axis.title.y = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
        theme(plot.title = element_text(size=18, hjust=0.5)) +
        xlab(input$xlab.tmg) +
        ylab(input$ylab.tmg) +
        geom_line(aes(color=trt)) +
        geom_point(aes(shape=trt, color=trt), show.legend=TRUE)
      p
    })
    
    # individual spider plot  
    output$gspiderplot.tmg <- renderPlot({
      var1 <- data()[,c(input$var1_tmg,input$var2_tmg,input$var3_tmg,input$var4_tmg)]
      #if (is.null(var1$v)){return(NULL)}
      tpos <- which(!is.na(var1[,input$var2_tmg]))
      new.var1 <- var1[tpos,]
      
      new.var1[,input$var4_tmg] <- factor(new.var1[,input$var4_tmg])
      new.var1[,input$var3_tmg] <- as.numeric(as.character(new.var1[,input$var3_tmg]))
      wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      
      p=ggplot(wdata, aes(x=time, y=size, group=id)) +
        theme_bw(base_size=14) +
        theme(axis.title.x = element_text(face="bold"), axis.text.x = element_text(face="bold")) +
        theme(axis.title.y = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
        theme(plot.title = element_text(size=18, hjust=0.5)) +
        xlab(input$xlab.tmg) +
        ylab(input$ylab.tmg) +
        geom_line(aes(color=trt)) +
        geom_point(aes(shape=trt, color=trt), show.legend=TRUE) +
        facet_grid(~trt)
      p
    })
    
    # mean plot
    output$meanplot.tmg <- renderPlot({
      var1 <- data()[,c(input$var1_tmg,input$var2_tmg,input$var3_tmg,input$var4_tmg)]
      #if (is.null(var1$v)){return(NULL)}
      tpos <- which(!is.na(var1[,input$var2_tmg]))
      new.var1 <- var1[tpos,]
      
      new.var1[,input$var4_tmg] <- factor(new.var1[,input$var4_tmg])
      new.var1[,input$var3_tmg] <- as.numeric(as.character(new.var1[,input$var3_tmg]))
      wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      
      tlevels <- levels(wdata$trt)
      
      tdata <- summarySEM(wdata, measurevar="size", groupvars=c("trt","time"))
      
      pd <- position_dodge(input$pdodge.tmg) # move them .05 to the left and right
      
      if(input$bar.tmg=="msd"){
        p=ggplot(tdata, aes(x=time, y=size, colour=trt)) + 
          xlab("Time (in days)") +
          ylab("Tumor size") +
          geom_errorbar(aes(ymin=size-sd, ymax=size+sd), width=.3, position=pd) +
          geom_line(position=pd,size=.3) +
          geom_point(position=pd)
      }else if(input$bar.tmg=="mse"){
        p=ggplot(tdata, aes(x=time, y=size, colour=trt)) + 
          xlab("Time (in days)") +
          ylab("Tumor size") +
          geom_errorbar(aes(ymin=size-se, ymax=size+se), width=.3, position=pd) +
          geom_line(position=pd,size=.3) +
          geom_point(position=pd)
      }else if(input$bar.tmg=="mci"){
        p=ggplot(tdata, aes(x=time, y=size, colour=trt)) + 
          xlab("Time (in days)") +
          ylab("Tumor size") +
          geom_errorbar(aes(ymin=size-ci, ymax=size+ci), width=.3, position=pd) +
          geom_line(position=pd,size=.3) +
          geom_point(position=pd)
      }else if(input$bar.tmg=="mdci"){
        p=ggplot(tdata, aes(x=time, y=median, colour=trt)) + 
          xlab("Time (in days)") +
          ylab("Tumor size") +
          geom_errorbar(aes(ymin=mlow, ymax=mup), width=.3, position=pd) +
          geom_line(position=pd,size=.3) +
          geom_point(position=pd)
      }
      p
    })
    
    # mixed-effects model
    output$out.tmg <- renderPrint({
      var1 <- data()[,c(input$var1_tmg,input$var2_tmg,input$var3_tmg,input$var4_tmg)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tmg == 'none'){
        tpos <- which(!is.na(var1[,input$var2_tmg]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tmg == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]) & var1[,input$var2_tmg]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- log(new.var1[,input$var2_tmg])
      }else if(input$transform.tmg == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]) & var1[,input$var2_tmg]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- sqrt(new.var1[,input$var2_tmg])
      }else if(input$transform.tmg == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- rank(new.var1[,input$var2_tmg])
      }
      
      new.var1[,input$var4_tmg] <- factor(new.var1[,input$var4_tmg])
      new.var1[,input$var3_tmg] <- as.numeric(as.character(new.var1[,input$var3_tmg]))
      #wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      new.var1$ID <- new.var1[,input$var1_tmg] 
      
      tform = as.formula(paste0(input$var2_tmg,"~",input$var3_tmg,"*",input$var4_tmg))
      lmfit = lme(tform,random=~1|ID,data=new.var1,control=list(msVerbose=FALSE))
      anova(lmfit)
    })
    
    # mixed-effects model: qq and density plot
    output$qqdplot.tmg <- renderPlot({
      var1 <- data()[,c(input$var1_tmg,input$var2_tmg,input$var3_tmg,input$var4_tmg)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tmg == 'none'){
        tpos <- which(!is.na(var1[,input$var2_tmg]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tmg == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]) & var1[,input$var2_tmg]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- log(new.var1[,input$var2_tmg])
      }else if(input$transform.tmg == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]) & var1[,input$var2_tmg]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- sqrt(new.var1[,input$var2_tmg])
      }else if(input$transform.tmg == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- rank(new.var1[,input$var2_tmg])
      }
      
      new.var1[,input$var4_tmg] <- factor(new.var1[,input$var4_tmg])
      new.var1[,input$var3_tmg] <- as.numeric(as.character(new.var1[,input$var3_tmg]))
      #wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      new.var1$ID <- new.var1[,input$var1_tmg] 
      
      tform = as.formula(paste0(input$var2_tmg,"~",input$var3_tmg,"*",input$var4_tmg))
      lmfit = lme(tform,random=~1|ID,data=new.var1,control=list(msVerbose=FALSE))
      
      par(mfrow=c(1,2))
      car::qqPlot(lmfit$residuals)
      plot(density(lmfit$residuals),main="")
    })
    
    # mixed-effects model: shapiro test
    output$swtest.tmg <- renderPrint({
      var1 <- data()[,c(input$var1_tmg,input$var2_tmg,input$var3_tmg,input$var4_tmg)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tmg == 'none'){
        tpos <- which(!is.na(var1[,input$var2_tmg]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tmg == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]) & var1[,input$var2_tmg]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- log(new.var1[,input$var2_tmg])
      }else if(input$transform.tmg == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]) & var1[,input$var2_tmg]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- sqrt(new.var1[,input$var2_tmg])
      }else if(input$transform.tmg == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- rank(new.var1[,input$var2_tmg])
      }
      
      new.var1[,input$var4_tmg] <- factor(new.var1[,input$var4_tmg])
      new.var1[,input$var3_tmg] <- as.numeric(as.character(new.var1[,input$var3_tmg]))
      #wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      new.var1$ID <- new.var1[,input$var1_tmg] 
      
      tform = as.formula(paste0(input$var2_tmg,"~",input$var3_tmg,"*",input$var4_tmg))
      lmfit = lme(tform,random=~1|ID,data=new.var1,control=list(msVerbose=FALSE))
      
      res <- lmfit$residuals
      shapiro.test(res)
    })        
    
    # mixed-effects model: residual plot
    output$resplot.tmg <- renderPlot({
      var1 <- data()[,c(input$var1_tmg,input$var2_tmg,input$var3_tmg,input$var4_tmg)]
      #if (is.null(var1$v)){return(NULL)}
      if(input$transform.tmg == 'none'){
        tpos <- which(!is.na(var1[,input$var2_tmg]))
        new.var1 <- var1[tpos,]
      }else if(input$transform.tmg == 'logtrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]) & var1[,input$var2_tmg]>0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- log(new.var1[,input$var2_tmg])
      }else if(input$transform.tmg == 'srtrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]) & var1[,input$var2_tmg]>=0)
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- sqrt(new.var1[,input$var2_tmg])
      }else if(input$transform.tmg == 'ranktrans'){
        tpos <- which(!is.na(var1[,input$var2_tmg]))
        new.var1 <- var1[tpos,]
        new.var1[,input$var2_tmg] <- rank(new.var1[,input$var2_tmg])
      }
      
      new.var1[,input$var4_tmg] <- factor(new.var1[,input$var4_tmg])
      new.var1[,input$var3_tmg] <- as.numeric(as.character(new.var1[,input$var3_tmg]))
      #wdata <- data.frame(id=new.var1[,1],size=new.var1[,2],time=new.var1[,3],trt=new.var1[,4])
      new.var1$ID <- new.var1[,input$var1_tmg] 
      
      tform = as.formula(paste0(input$var2_tmg,"~",input$var3_tmg,"*",input$var4_tmg))
      lmfit = lme(tform,random=~1|ID,data=new.var1,control=list(msVerbose=FALSE))
      
      plot(lmfit)
    })
    
    #################
    # Categorical Test
    
    
    # Summary
    output$summ.cat <- renderPrint({
      if(input$test_cat!="trend"){
        var1 <- data()[,c(input$var1_cat,input$var2_cat)]
        var1[,input$var1_cat] = factor(var1[,input$var1_cat])
        var1[,input$var2_cat] = factor(var1[,input$var2_cat])
      }else{
        var1 <- data()[,c(input$var1_cat,input$time_cat)]
        var1[,input$var1_cat] = factor(var1[,input$var1_cat])
        var1[,input$time_cat] = factor(var1[,input$time_cat])
      }
      alevel = list()
      alevel[[names(var1)[1]]] = levels(var1[,1])
      alevel[[names(var1)[2]]] = levels(var1[,2])
      alevel
    })
    
    # Contingency table
    output$cont.cat <- renderPrint({
      if(input$test_cat!="trend"){
        var1 <- data()[,c(input$var1_cat,input$var2_cat)]
        var1[,input$var1_cat] = factor(var1[,input$var1_cat])
        var1[,input$var2_cat] = factor(var1[,input$var2_cat])
        table(var1[,input$var1_cat],var1[,input$var2_cat])
      }else{
        var1 <- data()[,c(input$var1_cat,input$time_cat)]
        var1[,input$var1_cat] = factor(var1[,input$var1_cat])
        var1[,input$time_cat] = factor(var1[,input$time_cat])
        table(var1[,input$var1_cat],var1[,input$time_cat])
      }
    })
    
    # Expected contingency table
    output$econt.cat <- renderPrint({
      if(input$test_cat!="trend"){
        var1 <- data()[,c(input$var1_cat,input$var2_cat)]
        var1[,input$var1_cat] = factor(var1[,input$var1_cat])
        var1[,input$var2_cat] = factor(var1[,input$var2_cat])
        ttab = table(var1[,input$var1_cat],var1[,input$var2_cat])
      }else{
        var1 <- data()[,c(input$var1_cat,input$time_cat)]
        var1[,input$var1_cat] = factor(var1[,input$var1_cat])
        var1[,input$time_cat] = factor(var1[,input$time_cat])
        ttab = table(var1[,input$var1_cat],var1[,input$time_cat])
      }
      if(input$test_cat == "chisq"){
        as.matrix(chisq.test(ttab)$expected)
      }else if(input$test_cat == "fisher"){
        as.matrix(c(0))
      }else if(input$test_cat == "mcnemar"){
        as.matrix(c(0))
      }else if(input$test_cat == "trend"){
        as.matrix(c(0))
      }
    })
    
    # Testing out
    output$out.cat <- renderPrint({
      if(input$test_cat!="trend"){
        var1 <- data()[,c(input$var1_cat,input$var2_cat)]
        var1[,input$var1_cat] = factor(var1[,input$var1_cat])
        var1[,input$var2_cat] = factor(var1[,input$var2_cat])
        ttab = table(var1[,input$var1_cat],var1[,input$var2_cat])
      }else{
        var1 <- data()[,c(input$var1_cat,input$time_cat)]
        var1[,input$var1_cat] = factor(var1[,input$var1_cat])
        var1[,input$time_cat] = factor(var1[,input$time_cat])
        ttab = table(var1[,input$var1_cat],var1[,input$time_cat])
      }
      if(input$test_cat == "chisq"){
        chisq.test(ttab)
      }else if(input$test_cat == "fisher"){
        fisher.test(ttab)
      }else if(input$test_cat == "mcnemar"){
        mcnemar.test(ttab)
      }else if(input$test_cat == "trend"){
        prop.trend.test(ttab[1,],apply(ttab,2,sum))
      }
    })    
    
    
    #############
    # Correlation
    
    # cor plot
    output$plot.cor <- renderPlot({
      var1 <- (data()[,input$var_cor])
      #if (is.null(var1$v)){return(NULL)}
      new.var1 <- var1
      if(length(input$log.cor)>0){
        ttvar = intersect(input$log.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tdd = var1[,ttvar[i]]
            tpos1 <- which(!is.na(tdd) & tdd>0)
            tpos2 <- which(is.na(tdd) | tdd<=0)
            tdd2 = log(tdd)
            tdd2[tpos2] = NA
            new.var1[,ttvar[i]] = tdd2
          }
        }
      }
      if(length(input$sqr.cor)>0){
        ttvar = intersect(input$sqr.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tdd = var1[,ttvar[i]]
            tpos1 <- which(!is.na(tdd) & tdd>=0)
            tpos2 <- which(is.na(tdd) | tdd<0)
            tdd2 = sqrt(tdd)
            tdd2[tpos2] = NA
            new.var1[,ttvar[i]] = tdd2
          }
        }
      }
      
      #pairs(new.var1, lower.panel=panel.lm, upper.panel=panel.cor, cex.labels=1,font.labels=2,method=input$cmet_cor)
      #pairs(new.var1, method=input$cmet_cor, lower.panel=panel.smooth, upper.panel=panel.cor, cex.labels=1,font.labels=2)
      
      if(input$cmet_cor=="pearson"){
        pairs(new.var1, lower.panel=panel.lm, upper.panel=panel.cor, cex.labels=1,font.labels=2,method=input$cmet_cor)
      }else{
        pairs(new.var1, method=input$cmet_cor, lower.panel=panel.smooth, upper.panel=panel.cor, cex.labels=1,font.labels=2)
      }
      
    })

    # Cor: qq plot regression
    output$qqgraph.cor <- renderPlot({
      var1 <- (data()[,input$var_cor])
      #if (is.null(var1$v)){return(NULL)}
      new.var1 <- var1
      if(length(input$log.cor)>0){
        ttvar = intersect(input$log.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tdd = var1[,ttvar[i]]
            tpos1 <- which(!is.na(tdd) & tdd>0)
            tpos2 <- which(is.na(tdd) | tdd<=0)
            tdd2 = log(tdd)
            tdd2[tpos2] = NA
            new.var1[,ttvar[i]] = tdd2
          }
        }
      }
      if(length(input$sqr.cor)>0){
        ttvar = intersect(input$sqr.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tdd = var1[,ttvar[i]]
            tpos1 <- which(!is.na(tdd) & tdd>=0)
            tpos2 <- which(is.na(tdd) | tdd<0)
            tdd2 = sqrt(tdd)
            tdd2[tpos2] = NA
            new.var1[,ttvar[i]] = tdd2
          }
        }
      }
      
      num.var = dim(new.var1)[2]
      nam.var = colnames(new.var1)
      num.col = 3
      num.row = ceiling(num.var/3)
      
      if(num.var<=3){
        par(mfrow=c(1,num.var))
      }else{
        par(mfrow=c(num.row,num.col))
      }
      
      for(i in 1:num.var){
        car::qqPlot(na.omit(new.var1[,i]),xlab="",ylab="",main=nam.var[i])
      }
    })
    
    # Cor: shapiro test
    output$swtest.cor <- renderPrint({
      var1 <- (data()[,input$var_cor])
      #if (is.null(var1$v)){return(NULL)}
      new.var1 <- var1
      if(length(input$log.cor)>0){
        ttvar = intersect(input$log.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tdd = var1[,ttvar[i]]
            tpos1 <- which(!is.na(tdd) & tdd>0)
            tpos2 <- which(is.na(tdd) | tdd<=0)
            tdd2 = log(tdd)
            tdd2[tpos2] = NA
            new.var1[,ttvar[i]] = tdd2
          }
        }
      }
      if(length(input$sqr.cor)>0){
        ttvar = intersect(input$sqr.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tdd = var1[,ttvar[i]]
            tpos1 <- which(!is.na(tdd) & tdd>=0)
            tpos2 <- which(is.na(tdd) | tdd<0)
            tdd2 = sqrt(tdd)
            tdd2[tpos2] = NA
            new.var1[,ttvar[i]] = tdd2
          }
        }
      }
 
      apply(new.var1,2,function(vv){shapiro.test(na.omit(vv))})

    })
    
    
    # correlation analysis
    output$out.cor <- renderPrint({
      var1 <- (data()[,input$var_cor])
      #if (is.null(var1$v)){return(NULL)}
      new.var1 <- var1
      if(length(input$log.cor)>0){
        ttvar = intersect(input$log.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tdd = var1[,ttvar[i]]
            tpos1 <- which(!is.na(tdd) & tdd>0)
            tpos2 <- which(is.na(tdd) | tdd<=0)
            tdd2 = log(tdd)
            tdd2[tpos2] = NA
            new.var1[,ttvar[i]] = tdd2
          }
        }
      }
      if(length(input$sqr.cor)>0){
        ttvar = intersect(input$sqr.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tdd = var1[,ttvar[i]]
            tpos1 <- which(!is.na(tdd) & tdd>=0)
            tpos2 <- which(is.na(tdd) | tdd<0)
            tdd2 = sqrt(tdd)
            tdd2[tpos2] = NA
            new.var1[,ttvar[i]] = tdd2
          }
        }
      }
      
      cmethod = "Pearson's correlation"
      if(input$cmet_cor=="spearman"){
        cmethod = "Spearman's correlation"
      }else if(input$cmet_cor=="kendall"){
        cmethod = "Kendall's correlation"
      }
      
      list(method=cmethod,coefficient=cor.mest(new.var1,method=input$cmet_cor),p.value=cor.mp(new.var1,method=input$cmet_cor))
    })
    
    # transformation summary
    output$trans.cor <- renderTable({
      ttext = matrix("No",length(input$var_cor),3)
      ttext[,1] = input$var_cor
      if(length(input$log.cor)>0){
        ttvar = intersect(input$log.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tpos = which(input$var_cor==ttvar[i])
            ttext[tpos,2] = "Yes"
          }
        }
      }
      if(length(input$sqr.cor)>0){
        ttvar = intersect(input$sqr.cor,input$var_cor)
        tlen = length(ttvar)
        if(tlen>0){
          for(i in 1:tlen){
            tpos = which(input$var_cor==ttvar[i])
            ttext[tpos,3] = "Yes"
          }
        }
      }
      colnames(ttext) = c("Variable","log-transformation","square root-transformation")
      ttext
    })
    
    ###################
    # Regression
    
    
    # qq plot regression
    output$qqgraph.reg <- renderPlot({
      if(input$model_reg == "cox"){
        var1 <- data()[,c(input$ct.reg,input$cs.reg,input$indep_reg)]
        
        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        tform = paste0("Surv(",input$ct.reg,",",input$cs.reg,")")
      }else{
        if(input$model_reg == "linear"){
          resvar = input$res.reg
        }else if(input$model_reg == "logistic"){
          resvar = input$logit.reg
        }
        var1 <- data()[,c(resvar,input$indep_reg)]
        
        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        if(input$restrans_reg == "logtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = log(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }else if(input$restrans_reg == "srtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = sqrt(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }
        if(input$model_reg == "logistic"){
          new.var1[,resvar] = factor(new.var1[,resvar])
        }
        tform = paste0(resvar)
      }
      
      if(length(input$scale.reg)>0){
        ttvar = intersect(input$scale.reg,input$indep_reg)
        if(length(ttvar)>0)
          new.var1[,ttvar] = scale(new.var1[,ttvar])
      }
      
      if(length(input$indep_reg)>0){
        for(i in 1:length(input$indep_reg)){
          if(i==1){
            tform = paste0(tform,"~",input$indep_reg[i])
          }else{
            tform = paste0(tform,"+",input$indep_reg[i])
          }
        }
      }else{
        tform = paste0(tform,"~1")
      }
      
      if(input$model_reg == "cox"){
        #tfit = coxph(as.formula(tform),data=new.var1)
        #tlen = length(input$indep_reg)
        #trow = tlen%/%3
        #trem = tlen%%3
        #if(trem>0){
        #  trow = trow+1
        #}
        #par(mfrow=c(trow,3))
        #plot(cox.zph(tfit))
      }else if(input$model_reg == "linear"){
        tfit = lm(as.formula(tform),data=new.var1)
        par(mfrow=c(2,2))
        res <- residuals(tfit)
        car::qqPlot(res,xlab="",ylab="",main="Residual")
        plot(density(res,na.rm=T),main="")
        car::qqPlot(new.var1[,resvar],xlab="",ylab="",main=resvar)
        plot(density(new.var1[,resvar],na.rm=T),main="")
      }else if(input$model_reg == "logistic"){
        #tfit = glm(as.formula(tform),data=new.var1,family="binomial")
        #par(mfrow=c(2,2))
        #plot(tfit)
      }
    })
    
    # regression plot regression
    output$regression.reg <- renderPlot({
      if(input$model_reg == "cox"){
        var1 <- data()[,c(input$ct.reg,input$cs.reg,input$indep_reg)]
        
        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        tform = paste0("Surv(",input$ct.reg,",",input$cs.reg,")")
      }else{
        if(input$model_reg == "linear"){
          resvar = input$res.reg
        }else if(input$model_reg == "logistic"){
          resvar = input$logit.reg
        }
        var1 <- data()[,c(resvar,input$indep_reg)]
        
        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        if(input$restrans_reg == "logtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = log(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }else if(input$restrans_reg == "srtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = sqrt(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }
        if(input$model_reg == "logistic"){
          new.var1[,resvar] = factor(new.var1[,resvar])
        }
        tform = paste0(resvar)
      }
      
      if(length(input$scale.reg)>0){
        ttvar = intersect(input$scale.reg,input$indep_reg)
        if(length(ttvar)>0)
          new.var1[,ttvar] = scale(new.var1[,ttvar])
      }
      
      if(length(input$indep_reg)>0){
        for(i in 1:length(input$indep_reg)){
          if(i==1){
            tform = paste0(tform,"~",input$indep_reg[i])
          }else{
            tform = paste0(tform,"+",input$indep_reg[i])
          }
        }
      }else{
        tform = paste0(tform,"~1")
      }
      
      if(input$model_reg == "cox"){
        #tfit = coxph(as.formula(tform),data=new.var1)
        #tlen = length(input$indep_reg)
        #trow = tlen%/%3
        #trem = tlen%%3
        #if(trem>0){
        #  trow = trow+1
        #}
        #par(mfrow=c(trow,3))
        #plot(cox.zph(tfit))
      }else if(input$model_reg == "linear" & length(input$indep_reg)==1){
        tfit = lm(as.formula(tform),data=new.var1)
        plot(new.var1[,input$indep_reg[1]],new.var1[,resvar],xlab=input$indep_reg[1],ylab=resvar)
        abline(tfit,col="red",lwd=2)
      }else if(input$model_reg == "logistic"){
        #tfit = glm(as.formula(tform),data=new.var1,family="binomial")
        #par(mfrow=c(2,2))
        #plot(tfit)
      }
    })
    
    
    # shapiro test
    output$swtest.reg <- renderPrint({
      if(input$model_reg == "cox"){
        var1 <- data()[,c(input$ct.reg,input$cs.reg,input$indep_reg)]
        
        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        tform = paste0("Surv(",input$ct.reg,",",input$cs.reg,")")
      }else{
        if(input$model_reg == "linear"){
          resvar = input$res.reg
        }else if(input$model_reg == "logistic"){
          resvar = input$logit.reg
        }
        var1 <- data()[,c(resvar,input$indep_reg)]
        
        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        if(input$restrans_reg == "logtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = log(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }else if(input$restrans_reg == "srtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = sqrt(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }
        if(input$model_reg == "logistic"){
          new.var1[,resvar] = factor(new.var1[,resvar])
        }
        tform = paste0(resvar)
      }
      
      if(length(input$scale.reg)>0){
        ttvar = intersect(input$scale.reg,input$indep_reg)
        if(length(ttvar)>0)
          new.var1[,ttvar] = scale(new.var1[,ttvar])
      }
      
      if(length(input$indep_reg)>0){
        for(i in 1:length(input$indep_reg)){
          if(i==1){
            tform = paste0(tform,"~",input$indep_reg[i])
          }else{
            tform = paste0(tform,"+",input$indep_reg[i])
          }
        }
      }else{
        tform = paste0(tform,"~1")
      }
      
      if(input$model_reg == "cox"){
        #tfit = coxph(as.formula(tform),data=new.var1)
        #tlen = length(input$indep_reg)
        #trow = tlen%/%3
        #trem = tlen%%3
        #if(trem>0){
        #  trow = trow+1
        #}
        #par(mfrow=c(trow,3))
        #plot(cox.zph(tfit))
      }else if(input$model_reg == "linear"){
        tfit = lm(as.formula(tform),data=new.var1)
        res <- residuals(tfit)
        #shapiro.test(res)
        sres = data.frame(residual=res,response=new.var1[,resvar])
        apply(sres,2,function(vv){shapiro.test(na.omit(vv))})
      }else if(input$model_reg == "logistic"){
        #tfit = glm(as.formula(tform),data=new.var1,family="binomial")
        #par(mfrow=c(2,2))
        #plot(tfit)
      }
    })
    
    # reg plot
    output$plot.reg <- renderPlot({
      if(input$model_reg == "cox"){
        var1 <- data()[,c(input$ct.reg,input$cs.reg,input$indep_reg)]
        
        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        tform = paste0("Surv(",input$ct.reg,",",input$cs.reg,")")
      }else{
        if(input$model_reg == "linear"){
          resvar = input$res.reg
        }else if(input$model_reg == "logistic"){
          resvar = input$logit.reg
        }
        var1 <- data()[,c(resvar,input$indep_reg)]
        
        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        if(input$restrans_reg == "logtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = log(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }else if(input$restrans_reg == "srtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = sqrt(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }
        if(input$model_reg == "logistic"){
          new.var1[,resvar] = factor(new.var1[,resvar])
        }
        tform = paste0(resvar)
      }
      
      if(length(input$scale.reg)>0){
        ttvar = intersect(input$scale.reg,input$indep_reg)
        if(length(ttvar)>0)
          new.var1[,ttvar] = scale(new.var1[,ttvar])
      }
      
      if(length(input$indep_reg)>0){
        for(i in 1:length(input$indep_reg)){
          if(i==1){
            tform = paste0(tform,"~",input$indep_reg[i])
          }else{
            tform = paste0(tform,"+",input$indep_reg[i])
          }
        }
      }else{
        tform = paste0(tform,"~1")
      }
      
      if(input$model_reg == "cox"){
        tfit = coxph(as.formula(tform),data=new.var1)
        tlen = length(input$indep_reg)
        trow = tlen%/%3
        trem = tlen%%3
        if(trem>0){
          trow = trow+1
        }
        par(mfrow=c(trow,3))
        plot(cox.zph(tfit))
      }else if(input$model_reg == "linear"){
        tfit = lm(as.formula(tform),data=new.var1)
        par(mfrow=c(2,2))
        plot(tfit)
      }else if(input$model_reg == "logistic"){
        tfit = glm(as.formula(tform),data=new.var1,family="binomial")
        par(mfrow=c(2,2))
        plot(tfit)
      }  
    })
    
    # Regression analysis
    output$out.reg <- renderPrint({
      if(input$model_reg == "cox"){
        var1 <- data()[,c(input$ct.reg,input$cs.reg,input$indep_reg)]
        
        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        tform = paste0("Surv(",input$ct.reg,",",input$cs.reg,")")
      }else{
        if(input$model_reg == "linear"){
          resvar = input$res.reg
        }else if(input$model_reg == "logistic"){
          resvar = input$logit.reg
        }
        var1 <- data()[,c(resvar,input$indep_reg)]

        if(length(input$indepcat_reg)>0){
          for(i in 1:length(input$indepcat_reg)){
            var1[,input$indepcat_reg[i]] = factor(var1[,input$indepcat_reg[i]])
          }
        }
        
        new.var1 <- var1
        if(input$restrans_reg == "logtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = log(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }else if(input$restrans_reg == "srtrans"){
          tdd = var1[,resvar]
          tpos1 <- which(!is.na(tdd) & tdd>0)
          tpos2 <- which(is.na(tdd) | tdd<=0)
          tdd2 = sqrt(tdd)
          tdd2[tpos2] = NA
          new.var1[,resvar] = tdd2
        }
        if(input$model_reg == "logistic"){
          new.var1[,resvar] = factor(new.var1[,resvar])
        }
        tform = paste0(resvar)
      }
      
      if(length(input$scale.reg)>0){
        ttvar = intersect(input$scale.reg,input$indep_reg)
        if(length(ttvar)>0)
          new.var1[,ttvar] = scale(new.var1[,ttvar])
      }
      
      if(length(input$indep_reg)>0){
        for(i in 1:length(input$indep_reg)){
          if(i==1){
            tform = paste0(tform,"~",input$indep_reg[i])
          }else{
            tform = paste0(tform,"+",input$indep_reg[i])
          }
        }
      }else{
        tform = paste0(tform,"~1")
      }
      
      out.reg = c()
      if(input$model_reg == "cox"){
        tfit = coxph(as.formula(tform),data=new.var1)
        tsumm = cbind(summary(tfit)$conf[,c(1)],summary(tfit)$conf[,c(3)],summary(tfit)$conf[,c(4)],summary(tfit)$coef[,5])
        #print(tsumm)
        colnames(tsumm) = c("HR","95% CI (lower)","95% CI (upper)","p")
        rownames(tsumm) = rownames(summary(tfit)$coef)
        out.reg = list(Output=summary(tfit),Summary=tsumm,PHassumption=cox.zph(tfit))
      }else if(input$model_reg == "linear"){
        tfit = lm(as.formula(tform),data=new.var1)
        tsumm = cbind(exp(summary(tfit)$coef[,1]),exp(confint(tfit)),summary(tfit)$coef[,4])
        colnames(tsumm) = c("Estimate","95% CI (lower)","95% CI (upper)","p")
        rownames(tsumm) = rownames(summary(tfit)$coef)
        out.reg = list(Output=summary(tfit),Summary=tsumm)
      }else if(input$model_reg == "logistic"){
        tfit = glm(as.formula(tform),data=new.var1,family="binomial")
        tsumm = cbind(exp(summary(tfit)$coef[,1]),exp(confint(tfit)),summary(tfit)$coef[,4])
        colnames(tsumm) = c("OR","95% CI (lower)","95% CI (upper)","p")
        rownames(tsumm) = rownames(summary(tfit)$coef)
        out.reg = list(Levels=levels(new.var1[,1]),Output=summary(tfit),Summary=tsumm)
      }  
      
      out.reg
    })
    
    ###################
    # Survival Curve
    
    # summary of group variable
    output$summ.g.surv <- renderPrint({
      if(input$group_surv == "yes"){
        var1 <- as.character(data()[,input$g.surv])
        if (is.null(var1)){return(NULL)}
        #psych::describe(var1)
        tout = list()
        tout$Unique_values = unique(var1)
        tout$Summary_table = table(var1)
        tout
      }
    })
    
    # median and survival rate calculation
    # group E N rate at 1, median, HR, p-value
    output$out.surv <- renderPrint({
      sout = c()
      if(input$group_surv == "no"){
        tdata = data.frame(t=data()[,input$t_surv],s=data()[,input$s_surv])
        if(input$unit1.surv != input$unit2.surv){
          if(input$unit1.surv == "months"){
            tmp.t = tdata$t * ((365*3+366)/(12*4))
          }else if(input$unit1.surv == "years"){
            tmp.t = tdata$t * (365*3+366)/4
          }else{
            tmp.t = tdata$t
          }
          if(input$unit2.surv == "days"){
            tdata$t = tmp.t
          }else if(input$unit2.surv == "months"){
            tdata$t = tmp.t/((365*3+366)/(12*4))
          }else if(input$unit2.surv == "years"){
            tdata$t = tmp.t/((365*3+366)/4)
          }
        }
        sfit = survfit(Surv(t, s) ~ 1, data = tdata)
        srate = summary(sfit, times = input$rate_surv)
        sout = list(Rate=srate,Median=sfit)
      }
      sout
    })
    
    # median and survival rate calculation
    # group E N rate at 1, median, HR, p-value
    output$out2.surv <- renderPrint({
      sout = c()
      if(input$group_surv == "yes"){
        tdata = data.frame(t=data()[,input$t_surv],s=data()[,input$s_surv],group=data()[,input$g.surv])
        if(input$unit1.surv != input$unit2.surv){
          if(input$unit1.surv == "months"){
            tmp.t = tdata$t * ((365*3+366)/(12*4))
          }else if(input$unit1.surv == "years"){
            tmp.t = tdata$t * (365*3+366)/4
          }else{
            tmp.t = tdata$t
          }
          if(input$unit2.surv == "days"){
            tdata$t = tmp.t
          }else if(input$unit2.surv == "months"){
            tdata$t = tmp.t/((365*3+366)/(12*4))
          }else if(input$unit2.surv == "years"){
            tdata$t = tmp.t/((365*3+366)/4)
          }
        }
        #print("####")
        #print(input$g.ref.surv)
        if(input$g.ref.surv==""){
          tdata$group = factor(tdata$group)
        }else{
          #print("#*#*#*#*")
          #print(input$g.ref.surv)
          tlevel = strsplit(input$g.ref.surv,"\\s+")[[1]]
          #print("******")
          #print(tlevel)
          tdata$group = factor(tdata$group,levels=tlevel)
        }
        sfit <- survfit(Surv(t, s) ~ group, data = tdata)
        srate = summary(sfit, times = input$rate_surv)
        scox <- coxph(Surv(t,s)~group,data=tdata)
        cout <- summary(scox)
        phout <- cox.zph(scox)
        tsumm = cbind(summary(scox)$conf[,c(1)],summary(scox)$conf[,c(3)],summary(scox)$conf[,c(4)],summary(scox)$coef[,5])
        #print(tsumm)
        colnames(tsumm) = c("HR","95% CI (lower)","95% CI (upper)","p")
        rownames(tsumm) = rownames(summary(scox)$coef)
        sout = list(Rate=srate,Median=sfit,Cox.output=cout,Summary=tsumm,PH.assmuption=phout)
      }
      sout
    })
    
    # survival curves
    output$plot.surv <- renderPlot({
      if(input$group_surv == "no"){
        tdata = data.frame(t=data()[,input$t_surv],s=data()[,input$s_surv])
        if(input$unit1.surv != input$unit2.surv){
          if(input$unit1.surv == "months"){
            tmp.t = tdata$t * ((365*3+366)/(12*4))
          }else if(input$unit1.surv == "years"){
            tmp.t = tdata$t * (365*3+366)/4
          }else{
            tmp.t = tdata$t
          }
          if(input$unit2.surv == "days"){
            tdata$t = tmp.t
          }else if(input$unit2.surv == "months"){
            tdata$t = tmp.t/((365*3+366)/(12*4))
          }else if(input$unit2.surv == "years"){
            tdata$t = tmp.t/((365*3+366)/4)
          }
        }
        sfit = survfit(Surv(t, s) ~ 1, data = tdata)
        plot_main <- 
          ggsurvplot(
            data = tdata, 
            fit = sfit,
            xlab = input$xlab.surv,
            ylab = input$ylab.surv,
            legend = "none",
            #xscale = 30.4,
            #break.x.by = 182.4, 
            risk.table = TRUE,
            risk.table.y.text = FALSE)
      }else{
        tdata = data.frame(t=data()[,input$t_surv],s=data()[,input$s_surv],group=data()[,input$g.surv])
        if(input$unit1.surv != input$unit2.surv){
          if(input$unit1.surv == "months"){
            tmp.t = tdata$t * ((365*3+366)/(12*4))
          }else if(input$unit1.surv == "years"){
            tmp.t = tdata$t * (365*3+366)/4
          }else{
            tmp.t = tdata$t
          }
          if(input$unit2.surv == "days"){
            tdata$t = tmp.t
          }else if(input$unit2.surv == "months"){
            tdata$t = tmp.t/((365*3+366)/(12*4))
          }else if(input$unit2.surv == "years"){
            tdata$t = tmp.t/((365*3+366)/4)
          }
        }
        #print("####")
        #print(input$g.ref.surv)
        if(input$g.ref.surv==""){
          tdata$group = factor(tdata$group)
        }else{
          #print("#*#*#*#*")
          #print(input$g.ref.surv)
          tlevel = strsplit(input$g.ref.surv,"\\s+")[[1]]
          #print("******")
          #print(tlevel)
          tdata$group = factor(tdata$group,levels=tlevel)
        }
        sfit <- survfit(Surv(t, s) ~ group, data = tdata)
        tpval <- ifelse(input$pval.surv == "yes",TRUE,FALSE)
        plot_main <-
          ggsurvplot(
            fit = sfit, 
            data = tdata,
            xlab = input$xlab.surv,
            ylab = input$ylab.surv,
            risk.table = TRUE,
            risk.table.y.text = FALSE,
            pval = tpval
          )
      }
      plot_main
    })
    
    ###################
    # Contingency Table
    
    # display 10 rows initially
    output$tab.ctab <- renderPrint({
      thead = input$header.ctab
      ttab = input$data.ctab
      #generate data frame
      data0 <- read.csv(textConnection(ttab),sep=" ",header=thead)
      #generate matrix
      as.matrix(data0)
    })
    
    # expected values
    output$etab.ctab <- renderPrint({
      thead = input$header.ctab
      ttab = input$data.ctab
      #generate data frame
      data0 <- read.csv(textConnection(ttab),sep=" ",header=thead)
      #generate matrix
      ttab2 = as.matrix(data0)
      if(input$type_ctab == "chisq"){
        tout = chisq.test(ttab2)$expected
      }else{
        tout = c(0)
      }
      #generate matrix
      as.matrix(tout)
    })
    
    # summary
    output$summ.ctab <- renderPrint({
      thead = input$header.ctab
      ttab = input$data.ctab
      #generate data frame
      data0 <- read.csv(textConnection(ttab),sep=" ",header=thead)
      #generate matrix
      Hmisc::describe(as.matrix(data0))
    })    
    
    # test outcome
    output$out.ctab <- renderPrint({
      thead = input$header.ctab
      ttab = input$data.ctab
      #generate data frame
      data0 <- read.csv(textConnection(ttab),sep=" ",header=thead)
      #generate matrix
      ttab2 = as.matrix(data0)
      if(input$type_ctab == "chisq"){
        tout = chisq.test(ttab2)
      }else if(input$type_ctab == "fisher"){
        tout = fisher.test(ttab2,conf.level=input$conf.ctab,alternative=input$tail.ctab)
      }else if(input$type_ctab == "mcnemar"){
        tout = mcnemar.test(ttab2)
      }else if(input$type_ctab == "trend"){
        tout = prop.trend.test(ttab2[1,],apply(ttab2,2,sum))
      }
      tout
    })
    
    ################################
    # Multiple Comparison Correction
    
    # display p-values
    output$tab.mcorrect <- renderPrint({
      ttab = input$data.mcorrect
      #generate data frame
      data0 <- read.csv(textConnection(ttab),sep=" ",header=FALSE)
      #generate matrix
      as.numeric(na.omit(as.numeric(c(unlist(data0)))))
    })
    
    # test outcome
    output$out.mcorrect <- renderPrint({
      ttab = input$data.mcorrect
      #generate data frame
      data0 <- read.csv(textConnection(ttab),sep=" ",header=FALSE)
      #generate matrix
      tpvals =  as.numeric(na.omit(as.numeric(c(unlist(data0)))))
      tout = tpvals
      tmethod = "None"
      tflag = 0
      if(input$holm_mcorrect=="yes"){
        tmp = p.adjust(tpvals,method="holm")
        tout = rbind(c(tout),(tmp))
        tmethod = c(tmethod,"Holm")
        tflag = tflag + 1
      }
      if(input$bonferroni_mcorrect=="yes"){
        tmp = p.adjust(tpvals,method="bonferroni")
        tout = rbind((tout),(tmp))
        tmethod = c(tmethod,"Bonferroni")
        tflag = tflag + 1
      }
      if(input$hommel_mcorrect=="yes"){
        tmp = p.adjust(tpvals,method="hommel")
        tout = rbind((tout),(tmp))
        tmethod = c(tmethod,"Hommel")
        tflag = tflag + 1
      }
      if(input$hochberg_mcorrect=="yes"){
        tmp = p.adjust(tpvals,method="hochberg")
        tout = rbind((tout),(tmp))
        tmethod = c(tmethod,"Hochberg")
        tflag = tflag + 1
      }
      if(input$BH_mcorrect=="yes"){
        tmp = p.adjust(tpvals,method="BH")
        tout = rbind((tout),(tmp))
        tmethod = c(tmethod,"BH")
        tflag = tflag + 1
      }
      if(input$BY_mcorrect=="yes"){
        tmp = p.adjust(tpvals,method="BY")
        tout = rbind((tout),(tmp))
        tmethod = c(tmethod,"BY")
        tflag = tflag + 1
      }
      if(tflag==0){
        tout = matrix(tpvals,1,length(tpvals))
      }
      tout2 = as.data.frame(tout)
      dimnames(tout2)[[2]] = paste0("p",1:length(tpvals))
      tout3 =cbind(Method=tmethod,tout2)
      tout3
    })
    
  }
#)

options(warn=-1)

# Run the application 
shinyApp(ui = ui, server = server)
