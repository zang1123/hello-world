#Author: Dan Zang, PhD
#Affiliation: Integrated Laboratory Systems, Inc. 
#Contact: dzang@ils-inc.com
#Title: ARQSARModeling.R 
#Purpose: this code takes an input of PaDEL descriptors and uses the SVM model 
#to predict androgen receptor (AR) antagonism 
#Notes: random forest was used to determine the most relevant variables from 
#the PaDEL descriptors. 

#define output direcory

outDir<-"C:/ARQSARModeling"

# Load the data set containing active and inactive chemicals

ToxCastData<-read.table(file.path(outDir, "ToxCast-Descriptors-Active-Inactive.txt"), header=T, sep="\t", as.is=T)

# There are 1806 chemicals (206 actives and 1600 inactives)
# There are 1447 columns (1444 descriptors+"AUC.Antagonist"+"CODE"+ 
# "CASRN")

dim(ToxCastData)   # 1806   1447

# Use only the descriptors for data processing

ToxCastDescriptors<- ToxCastData[, c(1:1444)]

dim(ToxCastDescriptors)   # 1806   1444

# Load package caret for data processing

library(caret)

# Remove variables with low variance

ToxCastDescNonZero <- ToxCastDescriptors[,-nearZeroVar(ToxCastDescriptors)]

# Variables reduced to 1014 from original 1444

dim(ToxCastDescNonZero)     # 1806   1014
# Remove correlated variables with correlation coefficient > 0.9

ToxCastDescNonCorr<-ToxCastDescNonZero[,-findCorrelation(cor(ToxCastDescNonZero),.9)]


# Variables reduced to 509 from 1014

dim(ToxCastDescNonCorr)   # 1806  509
# Names of 509 descriptors after removing low variance and correlated variables

DescNames<-names(ToxCastDescNonCorr)


# The subset of ToxCast chemicals with Active label (206 chemicals)

ToxCastActive<-ToxCastData[1:206, DescNames]

dim(ToxCastActive)      #    206   509

# The subset of ToxCast chemicals with Inactive label (1600 chemicals)

ToxCastInactive<-ToxCastData[207:1806, DescNames]

dim(ToxCastInactive)      #    1600   509


# Sampling of Active chemicals. Training: 140; Test: 66.

ActiveTest<-sample(1:206, 66, replace = FALSE)

# Training set of active chemicals

ToxCastActiveTraining<-ToxCastActive[-ActiveTest,]

dim(ToxCastActiveTraining)           #    140   509

# Test set of active chemicals

ToxCastActiveTest<-ToxCastActive[ActiveTest,]

dim(ToxCastActiveTest)           #    66   509


# Sampling of Inactive chemicals. Training: 1100; Test: 500.

InactiveTest<-sample(1:1600, 500, replace = FALSE)

# Training set of inactive chemicals

ToxCastInactiveTraining<-ToxCastInactive[-InactiveTest,]

dim(ToxCastInactiveTraining)           #    1100   509

# Test set of inactive chemicals

ToxCastInactiveTest<-ToxCastInactive[InactiveTest,]

dim(ToxCastInactiveTest)           #    500   509


# Training set with 1240 chemicals (140 actives and 1100 inactives)

ToxCastTraining<-rbind(ToxCastActiveTraining, ToxCastInactiveTraining)

dim(ToxCastTraining)    # 1240   509

# Test set with 566 chemicals (66 actives and 500 inactives)

ToxCastTest<-rbind(ToxCastActiveTest, ToxCastInactiveTest)

dim(ToxCastTest)    # 566    509


# The labels for training and test sets

TrainingClass <- factor(c( rep("Active", 140), rep("Inactive", 1100)))

TestClass <- factor(c( rep("Active", 66), rep("Inactive", 500)))


# Set class weights 

wts <- 7/ table(TrainingClass)

# Variable importance ranking by randomForest

# Load package randomForest 

library(randomForest)

# Build a random forest model

ToxCast.rf<-randomForest(ToxCastTraining, TrainingClass, importance=TRUE, proximity=TRUE, classwt=wts)
# The importance of variables is ranked based on Active accuracy [1], Inactive accuracy # [2], MeanDecreaseAccuracy [3] and MeanDecreaseGini [4]. Here we use the rank from # MeanDecreaseAccuracy

# The 509 descriptors are sorted based on their importance

VarSorted<-sort(ToxCast.rf$importance[,3], decreasing = TRUE)

# The top 30 descriptors are employed to build classification models

VarTop30<-names(VarSorted[1:30])

# Training data consists of 1240 chemicals and 30 descriptors

traindata<- ToxCastTraining[,VarTop30]

dim(traindata)     # 1240   30

# Test data consists of 566 chemicals and 30 descriptors

testdata<- ToxCastTest[,VarTop30]

dim(testdata)       #  566  30

# Load package 1071 for data processing

library(e1071)

# Build a support vector machine model

SVMmodel <- svm(traindata, TrainingClass, cost = 90, gamma = 0.01, class.weights = wts)

# Predict the training set

PredTrain<-predict(SVMmodel, traindata)

table(PredTrain, TrainingClass)

# Predict the test set

PredTest<-predict(SVMmodel, testdata)

table(PredTest, TestClass)
