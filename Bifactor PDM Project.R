library(clusterGeneration)
library(psych)
library(lavaan)
library(matrixcalc)
library(semPlot)

########################################################################################

# Evaluating Model fit of Bifactor Models and Correlated-traits Model in Random Positive- Definite Covariance Matrices
### Completed for forum Spring 2019
### Since Model fit is overstated in observed data by bifactor models compared to other models, we conducted a small simulation
### study to examine bifactor model overfit where the data are generated with an underlaying CFA and bifactor network,
### with varying item size.
### Study results indicated that bifactor models, regardless of the true model structure, overfit the data

######################################################################################
#Practice generating data

#Generate Data - Just one Matrix 10 items

fun <- genPositiveDefMat("c-vine", dim=10)
fun2 <- genPositiveDefMat("c-vine", dim=10)
matrix10 <- replicate(10,fun)
data <- simplify2array(matrix10, higher=T)


###Generate Data Loop -- 50 matrcies
set.seed(7)
data <- list()
for(i in 1:1000){
  fun <- genPositiveDefMat("c-vine", dim=10)
  data <- c(data,list(fun))
}

########################################################################################
### Generate intial data with underlying CFA and Bifactor model structure 
### 10 items, 2 factors, 1,000 data sets generated and analyzed


### CFA LOOP
modelcfa <- paste("FA2~~FA1","FA1~~FA1", "FA2 ~~ FA2", sep = "\n")
FA1 <- "FA1 =~ 1*V1 + V2 + V3 + V4 + V5"
FA2 <- "FA2 =~ 1*V6 + V7 + V8 + V9 + V10"
modelcfa <- paste(modelcfa, FA1, FA2, sep = "\n")



for(i in 1:1000){
ind <- data[[i]]$Sigma
mdata <- matrix(ind, nrow = 10, ncol = 10,
                dimnames = list(paste("V", 1:10, sep=""), (paste("V", 1:10, sep="")))) 
fitcfa <- lavaan::sem(modelcfa, sample.cov =mdata, sample.nobs=1000)


if (!is.null(gof <- tryCatch((round(fitMeasures(fitcfa)[c("chisq", "df", "cfi", "tli", "rmsea","srmr")], digits = 3)), 
                error=function(e){cat("ERROR", conditionMessage(e), "\n")}))) {
  
  print(gof) }
}




### Bifactor loop
modelbi <- paste("FA1~~FA1", "FA2 ~~ FA2", sep = "\n")
g <- "g =~ 1*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10"
FA1 <- "FA1 =~ 1*V1 + V2 + V3 + V4 +V5"
FA2 <- "FA2 =~ 1*V6 + V7 + V8 + V9 +V10"
modelbi <- paste(modelbi,g, FA1, FA2, sep = "\n")
for(i in 1:1000){
  ind <- data[[i]]$Sigma
  bdata <- matrix(ind, nrow = 10, ncol = 10,
                  dimnames = list(paste("V", 1:10, sep=""), (paste("V", 1:10, sep="")))) 
  
  fitbi <- lavaan::sem(modelbi, sample.cov =bdata, sample.nobs=1000)
  
  if (!is.null(gof <- tryCatch((round(fitMeasures(fitbi)[c("chisq", "df", "cfi", "tli", "rmsea","srmr")], digits = 3)), 
                               error=function(e){cat("ERROR", conditionMessage(e), "\n")}))) {
    if(c((gof[3]>.8),(gof[4]>.8))){
      print(gof)}}
}



#####################################################################################
# Expand research for poster presentation Fall 2019 Odum institute at UNC

### This time, vary sample size (n = 100, 500, and 1,000) and number of items (i = 10 and 20)
### Study results found in poster. 

# Pirani S., & Bauer, D., (2019 November). Evaluating Bifactor and Correlated- Traits Models with Random Positive- 
#    Definite Covariance Matrices. Poster presented at Thomas M. Carsey Graduate Student Symposium 
#    in Chapel Hill, NC

######################################################################################


# For poster - 10 items, 50 random matrices, change sample size =  100, 500, 1000

set.seed(7)
data10 <- list()
for(i in 1:50){
  fun <- genPositiveDefMat("c-vine", dim=10)
  data10 <- c(data10,list(fun))
}


#CFA

modelcfa <- paste("FA2~~FA1","FA1~~FA1", "FA2 ~~ FA2", sep = "\n")
FA1 <- "FA1 =~ 1*V1 + V2 + V3 + V4 + V5"
FA2 <- "FA2 =~ 1*V6 + V7 + V8 + V9 + V10"
modelcfa <- paste(modelcfa, FA1, FA2, sep = "\n")

for(i in 1:50){
  ind <- data10[[i]]$Sigma
  mdata <- matrix(ind, nrow = 10, ncol = 10,
                  dimnames = list(paste("V", 1:10, sep=""), (paste("V", 1:10, sep="")))) 
  
  fitcfa <- lavaan::sem(modelcfa, sample.cov =mdata, sample.nobs=100)
  
  print(tryCatch((round(fitMeasures(fitcfa)[c("chisq", "df", "cfi", "tli", "rmsea","srmr")], digits = 3)), 
                 error=function(e){cat("ERROR", conditionMessage(e), "\n")}))
}

# Bifactor 

modelbi <- paste("FA1~~FA1", "FA2 ~~ FA2", sep = "\n")
g <- "g =~ 1*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10"
FA1 <- "FA1 =~ 1*V1 + V2 + V3 + V4 +V5"
FA2 <- "FA2 =~ 1*V6 + V7 + V8 + V9 +V10"
modelbi <- paste(modelbi,g, FA1, FA2, sep = "\n")

for(i in 1:50){
  ind <- data10[[i]]$Sigma
  bdata <- matrix(ind, nrow = 10, ncol = 10,
                  dimnames = list(paste("V", 1:10, sep=""), (paste("V", 1:10, sep="")))) 
  
  fitbi <- lavaan::sem(modelbi, sample.cov =bdata, sample.nobs=100)
  
  print(tryCatch((round(fitMeasures(fitbi)[c("chisq", "df", "cfi", "tli", "rmsea")], digits = 3)), 
                 error=function(e){cat("ERROR", conditionMessage(e), "\n")}))
  
}

##############################################################################################3


# For poster - 20 items (sample 100, 500, 1000) -- CHANGE SAMPLE SIZE ONLY

set.seed(7)
data20 <- list()
for(i in 1:50){
  fun <- genPositiveDefMat("c-vine", dim=10)
  data20 <- c(data20,list(fun))
}


#CFA

modelcfa <- paste("FA2~~FA1","FA1~~FA1", "FA2 ~~ FA2", sep = "\n")
FA1 <- "FA1 =~ 1*V1 + V2 + V3 + V4 + V5 + V6 +V7 + V8 + V9 +V10 "
FA2 <- "FA2 =~ 1*V11 +V12 +V13 +V14+ V15 +V16 +V17 +V18 +V19 +V20"
modelcfa <- paste(modelcfa, FA1, FA2, sep = "\n")

for(i in 1:50){
  ind <- data20[[i]]$Sigma
  mdata <- matrix(ind, nrow = 20, ncol = 20,
                  dimnames = list(paste("V", 1:20, sep=""), (paste("V", 1:20, sep="")))) 
  
  fitcfa <- lavaan::sem(modelcfa, sample.cov =mdata, sample.nobs=100)
  
  print(tryCatch((round(fitMeasures(fitcfa)[c("chisq", "df", "cfi", "tli", "rmsea","srmr")], digits = 3)), 
                 error=function(e){cat("ERROR", conditionMessage(e), "\n")}))
}

#Bifactor 

modelbi <- paste("FA1~~FA1", "FA2 ~~ FA2", sep = "\n")
g <- "g =~ 1*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
          V11 +V12 +V13 +V14+ V15 +V16 +V17 +V18 +V19 +V20"
FA1 <- "FA1 =~ 1*V1 + V2 + V3 + V4 + V5 + V6 +V7 + V8 + V9 +V10 "
FA2 <- "FA2 =~ 1*V11 +V12 +V13 +V14+ V15 +V16 +V17 +V18 +V19 +V20"
modelbi <- paste(modelbi,g, FA1, FA2, sep = "\n")

for(i in 1:50){
  ind <- data20[[i]]$Sigma
  bdata <- matrix(ind, nrow = 20, ncol = 20,
                  dimnames = list(paste("V", 1:20, sep=""), (paste("V", 1:20, sep="")))) 
  
  fitbi <- lavaan::sem(modelbi, sample.cov =bdata, sample.nobs=100)
  
  print(tryCatch((round(fitMeasures(fitbi)[c("chisq", "df", "cfi", "tli", "rmsea")], digits = 3)), 
                 error=function(e){cat("ERROR", conditionMessage(e), "\n")}))
  
}



##################################################################################



