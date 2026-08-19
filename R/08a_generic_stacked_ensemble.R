##Generic Stacked ensemble Script

rm(list = ls())
library(data.table)
library(gbm)
library(xgboost)
library(caret)
library(mgcv)
library(dplyr)
library(glmnet)
library(matrixStats)
library(quadprog)
library(gtools)
library(Cubist)
library(nnet)
library(randomForest)
library(kernlab)
library(monomvn)
library(gam)
library(rJava) # ..gz file in kp folder and need to ask chatgpt how install it properly
library(bartMachine)
library(h2o) # ..gz file in kp folder
library(kknn)
library(earth)
##library(gpls)
library(brnn)
library(bst)
##library(monmlp)

##library(mboost)
library(plyr)
##library(import)

library(ranger)
library(e1071)

#install.packages("e1071")
##library(partDSA)
##BiocManager::install(c("fastICA"))
##BiocManager::install(c("gpls"))

#~~~~~~~~~~~~~#
# i. Setup ####
#~~~~~~~~~~~~~#
# set seed #
set.seed(12345)

##Set working directory
setwd("C:/Users/naekaphirat/Desktop/KP")

#set output directory
model_date = format(Sys.Date(), "%Y_%m_%d")
outputdir <-  paste0('stacked_ensemble/',model_date, '/')
dir.create(outputdir, showWarnings = F, recursive = T)

library(data.table)
#Load data
mydata <- data.table(read.csv('outlier_data2.csv', stringsAsFactors = F))

#exclude outliers
mydata <- mydata[mydata$is_outlier == 0,]

#remove previous covariates
mydata <- mydata[, -c(15:19)]

#Set minimum sample size
mydata <- mydata[mydata$sample_size >10,
                 ]
#specify child models to include (for example...)

#child_models <- c('ridge', 'ranger', 'gam', 'knn', 'gbm_h2o')
#child_models <- c('gam', 'xgboost', 'enet', 'ridge', 'lasso', 'rf', 'nnet', 'cubist', 
#                  'pcaNNet', 'knn', 'gamLoess', 'gbm_h2o',
#                  'blasso', 'gbm', 'bagEarth', 'svmRadial', 'kknn', 'gaussprPoly',
#                  'lmStepAIC', 'svmLinear', 'svmPoly', 'glmnet_h2o', 'brnn', 'bstSm',
#                  'gamSpline', 'BstLm')

#Not working properly methodology
#child_models <- c('xgbtree', , 'bartMachine', gaussprLinear', 'pcr', 'ranger')



#Waiting list
#  'icr' 'monmlp'

child_models <- c('ridge', 'gam', 'gamLoess', 'gaussprPoly') #

#child_models <- c('monmlp')

#specify the stacker you want to use out of CWM (constrained weighted mean, from quadratic programming), 
# RWM (weighted mean based on R-sqr) GBM, GLM, nnet
stacker <- 'RWM'

#specify the family you are modelling (currently can use binomial or gaussian)
family <- 'binomial'

#specify transformation of the data to do - 'logit', 'log' or NULL
#nb if using logit this is only compatable for xgboost, random forest, cubist and neural networks
#child models with binomial data as these are using the probability, GAM and penalised regressions are using n and d
transformation <- 'log'

#Centre scale covariates?
centre_scale <- FALSE

#Include year in your models?
include_year <-  TRUE

#load covariates - either centre-scaled or standard
covs <- read.csv('Gram_covs_1990to2021.csv', stringsAsFactors = F)


#specify holdout method, currently can use random or country
holdout_method <- 'random'

#specify covariates you want to include in the model

covs_to_include <- c('year_id', 'pigs_pc', 'age_std_hiv_prev',
                     'sev_scalar_agestd_hiv', 'physicians_pc',
                     'smok_daily_prev_agstd_both', 'location_id') #

#specify what your columns are
p <- 'val'          #the proportion of your indicator successes
n <- 'resistant'            #the number of your indicator successes
d <- 'sample_size'  #the denominator (sample size)
w <- NULL            #the weights to use

#Specify which years you are modelling for
min_year <- 1990
max_year <- 2025

#rename some colums to avoid confusion
colnames(mydata)[colnames(mydata)==d] <- 'd' 
if(!is.null(p)) {colnames(mydata)[colnames(mydata)==p] <- 'p'} 
if(!is.null(n)) {colnames(mydata)[colnames(mydata)==n] <- 'n'} 

#if you dont have n but have p and d
if(is.null(n) &!is.null(p)&!is.null(d)){mydata$n <- mydata$p*mydata$d}

#if you dont have p but have n and d
if(is.null(p) &!is.null(n)&!is.null(d)){mydata$p <- mydata$n/mydata$d}

#if you havent specified a weights column set to 1
if(is.null(w)){
  mydata$w <- 1
} else {
  colnames(mydata)[colnames(mydata)==w] <- 'w' 
}

#perform transformations as specified
if(is.null(transformation)){
} else if(transformation == 'log'){
  if(family == 'binomial'){
    mydata$n <- log(mydata$n)
    mydata$d <- log(mydata$d)
    mydata$p <- log(mydata$p)
  } else if(family == 'gaussian')
    mydata$n <- log(mydata$n)
} else if(transformation == 'logit'){  #ln(p/1-p)
  if(family == 'binomial'){
    mydata$p <- log(mydata$p/(1-mydata$p))
  } else if(family == 'gaussian'){
    message('should not be using logit transformation with gaussian data')
  }
}

#restrict covs to those included
covs <- covs[colnames(covs) %in% covs_to_include | colnames(covs)=='location_id' | colnames(covs) =='year_id']
covs <- data.table(covs)
covs <- na.omit(covs)
head(covs)

# transform covs
#centre scale the covariates if desired
if(centre_scale == TRUE){
  covs <- data.frame(covs)
  covs[colnames(covs) %in% covs_to_include] <- data.frame(scale(covs[colnames(covs) %in% covs_to_include]))
  covs$year <- scale(covs$year_id)
  covs <-  data.table(covs)
}

#merge covs onto data
mydata <- merge(mydata, covs, by = c('location_id', 'year_id'))
mydata <- data.table(mydata)

## remove NAs
if(family == 'binomial'){
  mydata    <- na.omit(mydata, c('n', 'd', 'p', names(covs)))
}

if(family == 'gaussian'){
  mydata    <- na.omit(mydata, c('n', names(covs)))
}

## shuffle the data into five random folds
if(holdout_method == 'random'){
  mydata <- mydata[sample(nrow(mydata)),]
  mydata[,fold_id := cut(seq(1,nrow(mydata)),breaks=5,labels=FALSE)]
}

if(holdout_method == 'country'){
  country <- unique(mydata[, country])
  country <- country[sample(length(country))]
  fold_id <- cut(seq(1,length(country)),breaks=5,labels=FALSE)
  folds <- as.data.table(cbind(country, fold_id))
  
  mydata <- merge(mydata, folds, by = c('country'))
  mydata$fold_id <- as.numeric(mydata$fold_id)
  rm(country, fold_id)
}

# Limit to required years
mydata <- mydata[mydata$year_id >= min_year & mydata$year_id <= max_year,]
covs <- covs[covs$year_id >= min_year & covs$year_id <= max_year,]

## add a row id column
mydata[, a_rowid := seq(1:nrow(mydata))]

if(include_year == TRUE){
  covs_to_include <- c('year', covs_to_include)
}
mydata$year <- mydata$year_id
covs$year <- covs$year_id

##plot( mydata$maternal_educ_yrs_pc, mydata$p, pch=19)
##abline(lm(mydata$p~mydata$maternal_educ_yrs_pc), col="red") # regression line (y~x)
##lines(lowess(mydata$maternal_educ_yrs_pc,mydata$p), col="blue") # lowess line (x,y)

##cor(mydata$maternal_educ_yrs_pc,mydata$p, method = c("pearson", "kendall", "spearman"))
##cor.test(mydata$maternal_educ_yrs_pc,mydata$p, method=c("pearson", "kendall", "spearman"))

#~~~~~~~~~~~~~~~~~~~~~#
# Fit child models ####
#~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~#
# 1. XGBoost ####
#~~~~~~~~~~~~~~~#

if('xgboost' %in% child_models){
  dir.create(paste0(outputdir, '/xgboost'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  #tune the XGBoost
  # message("Model tuning xgboost")
  # #set the options for parameters to look at whilst tuning - can adjust these as required
  xg_grid <- expand.grid(nrounds = c(50, 100, 200),
                         max_depth = c(4, 6, 8, 10, 12),
                         eta = (3:8) / 100,
                         colsample_bytree = .5,
                         min_child_weight = 1,
                         subsample = 1,
                         gamma = 0)
  
  # Set cross validation options, default to 5 times repeated 5-fold cross validation
  # Selection function is "oneSE" to pick simplest model within one standard error of minimum
  # then imput this into the training model
  train_control <- trainControl(selectionFunction = "oneSE",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]))
  
  # Fit model
  xg_fit <- train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = F,
                  tuneGrid = xg_grid,
                  metric = "RMSE",
                  method = "xgbTree",
                  ##   objective = if(family == 'binomial'){"reg:logistic"}else if(family == 'gaussian'){"reg:linear"}else{message('Family of model not compatiable')},
                  weights = mydata$w)
  
  # Save model fit object
  saveRDS(xg_fit, paste0(outputdir, "/xgboost/xg_fit.RDS"))
  
  # Save the best parameters to csv file
  write.csv(xg_fit$bestTune, paste0(outputdir, 'xgboost/xgboost_best_tune_.csv'))
  xg_best_tune <- xg_fit$bestTune
  
  #set up final parameters based on the model tuning
  xg_grid_final <- expand.grid(nrounds = xg_best_tune$nrounds,
                               max_depth = xg_best_tune$max_depth,
                               eta = xg_best_tune$eta,
                               colsample_bytree = .5,
                               min_child_weight = 1,
                               subsample = 1,
                               gamma = 0)
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  #Fit final model
  message("Fitting xgboost on final tuned hyperparameters")
  xg_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = F,
                        tuneGrid = xg_grid_final,
                        metric = "RMSE",
                        method = "xgbTree",
                        ##              objective = if(family == 'binomial'){"reg:logistic"}else if(family == 'gaussian'){"reg:linear"}else{message('Family of model not compatiable')},
                        weights = mydata$w)
  
  # Plot the covariate importance of final model
  cov_plot <-
    ggplot(varImp(xg_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/xgboost/_covariate_importance.png'),
         plot = cov_plot)
  
  # Extract out of sample and in sample predictions
  mydata[, 'xgboost_cv_pred'   := arrange(xg_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'xgboost_full_pred' := predict(xg_fit_final, mydata)]
  #fill in the data
  #mydata[fold_id==i, 'xgboost_cv_pred' := predict(xg_fit_final, mydata[fold_id==i,],type = 'response')] 
  
  #save model fit
  xg_fit_final$model_name <- "xgboost"
  saveRDS(xg_fit_final, paste0(outputdir, '/xgboost/full_xgboost.RDS'))
  
  #predict out for all locations
  covs[, 'xgboost' := predict(xg_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, xg_best_tune, xg_fit, xg_fit_final, xg_grid, xg_grid_final)
}

#~~~~~~~~~~~#
# 2. GAM ####
#~~~~~~~~~~~#
if('gam' %in% child_models){
  dir.create(paste0(outputdir, '/gam/'), showWarnings = F)
  
  #If there are any binary covariates then remove them from the cov list and add as additional terms
  
  #set response variable
  if(family == 'binomial'){
    response <- cbind(sucesses = mydata$n, 
                      failures = mydata$d - mydata$n)
  } else if (family == 'gaussian'){
    response <- mydata$n
  }
  
  #build the GAM formula:
  #response ~ 1 + s(covariates, spline arguments)
  #TO DO: Look at tuning this model with different splines
  gam_formula <- paste0('response ~ 1+ s(', paste(covs_to_include, collapse = ", bs = 'ts', k = 3) + s("), ", bs = 'ts', k = 3)")
  gam_formula <- as.formula(gam_formula)
  
  # Fit full model
  #has some sort of parrallelisation inbuilt - set using this 
  full_gam = mgcv::gam(gam_formula, 
                       data = mydata, 
                       family = if(family =='binomial'){'quasibinomial'}else if(family == 'gaussian'){'gaussian'}, 
                       weights = mydata$w, 
                       control = list(nthreads = 2))
  full_gam$model_name = 'GAM'
  
  #predict using full model fit earlier
  mydata[,'gam_full_pred' := predict(full_gam, mydata, type = 'response')]
  
  #fit the model on the holdouts
  for(i in 1:5){
    if(family == 'binomial'){
      response <- cbind(successes = mydata$n[mydata$fold_id!=i], 
                        failures = mydata$d[mydata$fold_id!=i] - mydata$n[mydata$fold_id!=i])
    } else if (family == 'gaussian'){
      response <- mydata$n[mydata$fold_id!=i] 
    }      
    
    baby_gam = mgcv::gam(gam_formula, 
                         data = mydata[mydata$fold_id!=i], 
                         family = if(family =='binomial'){'quasibinomial'}else if(family == 'gaussian'){'gaussian'}, 
                         weights = mydata$weight[mydata$fold_id!=i], 
                         control = list(nthreads = 2))
    
    #fill in the data
    mydata[fold_id==i, 'gam_cv_pred' := predict(baby_gam, mydata[fold_id==i,],type = 'response')] 
    
  }
  
  #save full model fit
  saveRDS(full_gam, paste0(outputdir, '/gam/full_gam.RDS'))
  
  #predict out for all locations
  covs[,'gam' := predict(full_gam, covs, type = 'response')]
  
  #plot out GAM results to analyse
  pdf(paste0(outputdir, '/gam/plots.pdf'))
  gam.check(full_gam)
  dev.off()
  
  rm(baby_gam, full_gam, gam_formula, response)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Penalised regression (E-net/Ridge/Lasso) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#alpha 0 = Ridge, alpha 1 = Lasso, inbetween = e-net
if('enet' %in% child_models | 'ridge' %in% child_models | 'lasso' %in% child_models){
  
  dir.create(paste0(outputdir, '/glmnet'),showWarnings = F)
  
  #define the response to be modeled (2 variable matrix)
  if(family == 'binomial'){
    response <- cbind(failures   = mydata$d - mydata$n, 
                      successes = mydata$n)
  }else if(family == 'gaussian'){
    response <- mydata$n
  }
  
  #define variables to include (as a matrix)
  vars <- as.matrix(mydata[, covs_to_include, with = F])
  colnames(vars) <- covs_to_include
  
  #train model to select lambda and alpha,
  #no function to select alpha in glmnet so just run with a few options and compare MSE
  #let the function select its own ranges of lambda (can select manually of desired)
  #use 5 fold CV for this (?? could change to 10 ??)
  cv_lambda0 = cv.glmnet(x = vars , y= response, family = family, alpha = 0, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda0.25 = cv.glmnet(x = vars , y= response, family = family, alpha = 0.25, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda0.5 = cv.glmnet(x = vars , y= response, family = family, alpha = 0.5, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda0.75 = cv.glmnet(x = vars , y= response, family = family, alpha = 0.75, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda1 = cv.glmnet(x = vars , y= response, family = family, alpha = 1, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  
  #plot out the lambda and alpha options
  ##LOOK AT THESE PLOTS to select your prefered penalised regression model (cannot use multiple as they will be correlated)
  pdf(paste0(outputdir, '/glmnet/parameter_selection.pdf'))
  par(mfrow=c(3,2))
  plot(cv_lambda0)
  plot(cv_lambda0.25)
  plot(cv_lambda0.5)
  plot(cv_lambda0.75)
  plot(cv_lambda1)
  plot(log(cv_lambda0$lambda),cv_lambda0$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv_lambda0$name)
  points(log(cv_lambda0.25$lambda),cv_lambda0.25$cvm,pch=19,col="pink")
  points(log(cv_lambda0.5$lambda),cv_lambda0.5$cvm,pch=19,col="blue")
  points(log(cv_lambda0.75$lambda),cv_lambda0.75$cvm,pch=19,col="yellow")
  points(log(cv_lambda1$lambda),cv_lambda1$cvm,pch=19,col="green")
  legend("bottomright",legend=c("alpha= 1","alpha= .75", "alpha= .5", "alpha= .25","alpha 0"),pch=19,col=c("green","yellow","blue","pink","red"))
  dev.off()
  
  #fit the full model using selected lambda and alpha
  if('ridge' %in% child_models){full_ridge = glmnet(x = vars , y= response, family = family, alpha = 0, weights = mydata$w)}
  if('enet' %in% child_models){full_enet = glmnet(x = vars , y= response, family = family, alpha = 0.5, weights = mydata$w)}
  if('lasso' %in% child_models){full_lasso = glmnet(x = vars , y= response, family = family, alpha = 1, weights = mydata$w)}
  
  #predict full model results (requires matrix)
  # used 'response' which gives the same results as inverse logit (link)
  if('ridge' %in% child_models){mydata[,'ridge_full_pred' := predict(full_ridge,newx = vars, s = cv_lambda0$lambda.1se, type = 'response')]}
  if('lasso' %in% child_models){mydata[,'lasso_full_pred' := predict(full_lasso,newx = vars, s = cv_lambda1$lambda.1se, type = 'response')]}
  if('enet' %in% child_models){mydata[,'enet_full_pred' := predict(full_enet,newx = vars, s = cv_lambda0.5$lambda.1se, type = 'response')]}
  
  #fit the model on the holdouts
  for(i in 1:5){
    if(family == 'binomial'){
      response <- cbind(failures = mydata$d[mydata$fold_id!=i] - mydata$n[mydata$fold_id!=i], 
                        successes = mydata$n[mydata$fold_id!=i])
    } else if(family == 'gaussian'){
      response <- mydata$n[mydata$fold_id!=i] 
    }
    
    vars <- as.matrix(mydata[fold_id != i, covs_to_include, with = F])
    colnames(vars) <- covs_to_include
    
    if('ridge' %in% child_models){baby_ridge = glmnet(x = vars , y= response, family = family, lambda = cv_lambda0$lambda.1se, alpha = 0, weights = mydata$w[mydata$fold_id!=i])}
    if('lasso' %in% child_models){baby_lasso = glmnet(x = vars , y= response, family = family, lambda = cv_lambda1$lambda.1se, alpha = 1, weights = mydata$w[mydata$fold_id!=i])}
    if('enet' %in% child_models){baby_enet = glmnet(x = vars , y= response, family = family, lambda = cv_lambda0.5$lambda.1se, alpha = 0.5, weights = mydata$w[mydata$fold_id!=i])}
    
    new_vars <- as.matrix(mydata[fold_id == i, covs_to_include, with = F])
    
    #fill in the data
    if('ridge' %in% child_models){mydata[fold_id==i,'ridge_cv_pred' := predict(baby_ridge,newx = new_vars, s = cv_lambda0$lambda.1se, type = 'response')]}
    if('lasso' %in% child_models){mydata[fold_id==i,'lasso_cv_pred' := predict(baby_lasso,newx = new_vars, s = cv_lambda1$lambda.1se, type = 'response')]}
    if('enet' %in% child_models){mydata[fold_id==i,'enet_cv_pred' := predict(baby_enet,newx = new_vars, s = cv_lambda0.5$lambda.1se, type = 'response')]}
  }
  
  #save the model and relevent coefficients
  if('ridge' %in% child_models){saveRDS(cv_lambda0, paste0(outputdir, '/glmnet/full_ridge.rds'))}
  if('enet' %in% child_models){saveRDS(cv_lambda0.5, paste0(outputdir, '/glmnet/full_enet.rds'))}
  if('lasso' %in% child_models){saveRDS(cv_lambda1, paste0(outputdir, '/glmnet/full_lasso.rds'))}
  
  #predict out for all locations - can change the lambdas if required
  all_names <- names(covs) 
  new_covs <- as.matrix(covs)
  names(new_covs) <- all_names
  if('ridge' %in% child_models){covs[,'ridge' := predict(full_ridge,newx = new_covs[,rownames(full_ridge$beta)], s = cv_lambda0$lambda.1se, type = 'response')]}
  if('enet' %in% child_models){covs[,'enet' := predict(full_enet,newx = new_covs[,rownames(full_enet$beta)], s = cv_lambda0.5$lambda.1se, type = 'response')]}
  if('lasso' %in% child_models){covs[,'lasso' := predict(full_lasso,newx = new_covs[,rownames(full_lasso$beta)], s = cv_lambda1$lambda.1se, type = 'response')]}
  
  rm(cv_lambda1, cv_lambda0.5, cv_lambda0, cv_lambda0.25, cv_lambda0.75, full_lasso, full_enet, full_ridge, baby_lasso, baby_ridge, baby_enet, new_vars, response, vars, i, new_covs, all_names)    
}

#~~~~~~~~~~~~~~~~~~~~~#
# 4. Random forest ####
#~~~~~~~~~~~~~~~~~~~~~#
if('rf' %in% child_models){
  dir.create(paste0(outputdir, '/rf'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.mtry=c(1:length(covs_to_include)))
  
  # Fit model
  rf_fit <- train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = T,
                  tuneGrid = tunegrid,
                  metric = "RMSE",
                  method = "rf",
                  weights = mydata$w)
  
  # Save model fit object 
  saveRDS(rf_fit, paste0(outputdir, "/rf/rf_fit.RDS"))
  png(paste0(outputdir, '/rf/rf_fit.png'))
  plot(rf_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- rf_fit$bestTune$mtry
  write.csv(mtry_tune, paste0(outputdir, '/rf/rf_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  rf_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = T,
                        tuneGrid = tunegrid_final,
                        metric = "RMSE",
                        method = "rf",
                        importance=T,
                        weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'rf_cv_pred'   := arrange(rf_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'rf_full_pred' := predict(rf_fit_final, mydata)]
  
  #save model fit
  rf_fit_final$model_name <- "rf"
  saveRDS(rf_fit_final, paste0(outputdir, '/rf/full_rf.RDS'))
  
  cov_plot <-
    ggplot(varImp(rf_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/rf/rf_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'rf' := predict(rf_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, rf_fit, rf_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Neural networks ####
#~~~~~~~~~~~~~~~~~~~~~~~#
if('nnet' %in% child_models){
  dir.create(paste0(outputdir, '/nnet'), showWarnings = F)
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.decay = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001), .size = c(3,4, 5, 6, 7, 8, 9,10,11,12))
  
  # Fit model
  nn_fit <- train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = T,
                  tuneGrid = tunegrid,
                  metric = "RMSE",
                  method = "nnet",
                  linout = FALSE,
                  maxit = 1000,
                  weights = mydata$w)
  
  # Save model fit object 
  saveRDS(nn_fit, paste0(outputdir, "/nnet/nn_fit.RDS"))
  png(paste0(outputdir, '/nnet/nn_fit.png'))
  plot(nn_fit)  
  dev.off()
  
  # Save the best parameters to csv file
  write.csv(nn_fit$bestTune, paste0(outputdir, '/nnet/nnet_best_tune.csv'))
  
  #specify the parameters
  tunegrid_final <- expand.grid(.decay=nn_fit$bestTune$decay, .size=nn_fit$bestTune$size)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  nn_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = T,
                        tuneGrid = tunegrid_final,
                        metric = "RMSE",
                        method = "nnet",
                        linout = FALSE,
                        maxit = 1000,
                        weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'nnet_cv_pred'   := arrange(nn_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'nnet_full_pred' := predict(nn_fit_final, mydata)]
  
  #save model fit
  nn_fit_final$model_name <- "nn"
  saveRDS(nn_fit_final, paste0(outputdir, '/nnet/full_nn.RDS'))
  
  cov_plot <-
    ggplot(varImp(nn_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/nnet/nn_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'nnet' := predict(nn_fit_final, covs)]
  
  rm(nn_fit, nn_fit_final, train_control, train_control_final, tunegrid, tunegrid_final, cov_plot)
}

#~~~~~~~~~~~~~~~~~~~~#
# 6. Cubist model ####
#~~~~~~~~~~~~~~~~~~~~#
if('cubist' %in% child_models){
  dir.create(paste0(outputdir, '/cubist'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.committees = c(seq(1, 40, 5)), 
                          .neighbors  = c(0, 3, 6, 9))
  
  # Fit model
  cubist_fit <- train(form,
                      data = mydata,
                      trControl = train_control,
                      verbose = T,
                      tuneGrid = tunegrid,
                      metric = "RMSE",
                      method = "cubist",
                      control = Cubist::cubistControl(),
                      weights = mydata$w)
  
  # Save model fit object 
  saveRDS(cubist_fit, paste0(outputdir, "/cubist/cubist_fit.RDS"))
  png(paste0(outputdir, '/cubist/cubist_fit.png'))
  plot(cubist_fit)  
  dev.off()
  
  # Save the best parameters to csv file
  write.csv(cubist_fit$bestTune, paste0(outputdir, '/cubist/cubist_best_tune.csv'))
  
  #specify the parameters
  tunegrid_final <- expand.grid(.committees=cubist_fit$bestTune$committees, .neighbors=cubist_fit$bestTune$neighbors)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  cubist_fit_final <- train(form,
                            data = mydata,
                            trControl = train_control_final,
                            verbose = T,
                            tuneGrid = tunegrid_final,
                            metric = "RMSE",
                            method = "cubist",
                            weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'cubist_cv_pred'   := arrange(cubist_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'cubist_full_pred' := predict(cubist_fit_final, mydata)]
  
  #save model fit
  cubist_fit_final$model_name <- "cubist"
  saveRDS(cubist_fit_final, paste0(outputdir, '/cubist/full_cubist.RDS'))
  
  cov_plot <-
    ggplot(varImp(cubist_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/cubist/cubist_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'cubist' := predict(cubist_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, cubist_fit, cubist_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~#
# 7. pcaNNet model ####
#~~~~~~~~~~~~~~~~~~~~#
if('pcaNNet' %in% child_models){
  dir.create(paste0(outputdir, '/pcaNNet'), showWarnings = F)
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.decay = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001), .size = c(4, 5, 6, 7, 8, 9))
  
  # Fit model
  nn_fit <- train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = T,
                  tuneGrid = tunegrid,
                  metric = "RMSE",
                  method = "pcaNNet",
                  linout = FALSE,
                  maxit = 1000,
                  weights = mydata$w)
  
  # Save model fit object 
  saveRDS(nn_fit, paste0(outputdir, "/pcaNNet/nn_fit.RDS"))
  png(paste0(outputdir, '/pcaNNet/nn_fit.png'))
  plot(nn_fit)  
  dev.off()
  
  # Save the best parameters to csv file
  write.csv(nn_fit$bestTune, paste0(outputdir, '/pcaNNet/pcaNNet_best_tune.csv'))
  
  #specify the parameters
  tunegrid_final <- expand.grid(.decay=nn_fit$bestTune$decay, .size=nn_fit$bestTune$size)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  nn_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = T,
                        tuneGrid = tunegrid_final,
                        metric = "RMSE",
                        method = "pcaNNet",
                        linout = FALSE,
                        maxit = 1000,
                        weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'pcaNNet_cv_pred'   := arrange(nn_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'pcaNNet_full_pred' := predict(nn_fit_final, mydata)]
  
  #save model fit
  nn_fit_final$model_name <- "nn"
  saveRDS(nn_fit_final, paste0(outputdir, '/pcaNNet/full_nn.RDS'))
  
  cov_plot <-
    ggplot(varImp(nn_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/pcaNNet/nn_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'pcaNNet' := predict(nn_fit_final, covs)]
  
  rm(nn_fit, nn_fit_final, train_control, train_control_final, tunegrid, tunegrid_final, cov_plot)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 8. eXtreme Gradient Boosting model ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if('xgbTree' %in% child_models){
  dir.create(paste0(outputdir, '/xgbTree'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  
  #tune the xgbTree
  # message("Model tuning xgbTree")
  # #set the options for parameters to look at whilst tuning - can adjust these as required
  xg_grid <- expand.grid(nrounds = c(50, 100, 200),
                         max_depth = c(4, 6, 8, 10, 12),
                         eta = (3:8) / 100,
                         colsample_bytree = .5,
                         min_child_weight = 1,
                         subsample = 1,
                         gamma = 0)
  
  
  # Set cross validation options, default to 5 times repeated 5-fold cross validation
  # Selection function is "oneSE" to pick simplest model within one standard error of minimum
  # then imput this into the training model
  train_control <- trainControl(selectionFunction = "oneSE",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]))
  
  # Fit model
  xg_fit <- train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = F,
                  tuneGrid = xg_grid,
                  metric = "RMSE",
                  method = "xgbTree",
                  objective = if(family == 'binomial'){"reg:logistic"}else if(family == 'gaussian'){"reg:linear"}else{message('Family of model not compatiable')},
                  weights = mydata$w)
  
  # Save model fit object
  saveRDS(xg_fit, paste0(outputdir, "/xgbTree/xg_fit.RDS"))
  
  # Save the best parameters to csv file
  write.csv(xg_fit$bestTune, paste0(outputdir, 'xgbTree/xgbTree_best_tune_.csv'))
  xg_best_tune <- xg_fit$bestTune
  
  #set up final parameters based on the model tuning
  xg_grid_final <- expand.grid(nrounds = xg_best_tune$nrounds,
                               max_depth = xg_best_tune$max_depth,
                               eta = xg_best_tune$eta,
                               colsample_bytree = .5,
                               min_child_weight = 1,
                               subsample = 1,
                               gamma = 0)
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  #Fit final model
  message("Fitting xgbTree on final tuned hyperparameters")
  xg_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = F,
                        tuneGrid = xg_grid_final,
                        metric = "RMSE",
                        method = "xgbTree",
                        objective = if(family == 'binomial'){"reg:logistic"}else if(family == 'gaussian'){"reg:linear"}else{message('Family of model not compatiable')},
                        weights = mydata$w)
  
  # Plot the covariate importance of final model
  cov_plot <-
    ggplot(varImp(xg_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/xgbTree/_covariate_importance.png'),
         plot = cov_plot)
  
  # Extract out of sample and in sample predictions
  mydata[, 'xgbTree_cv_pred'   := arrange(xg_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'xgbTree_full_pred' := predict(xg_fit_final, mydata)]
  
  #save model fit
  xg_fit_final$model_name <- "xgbTree"
  saveRDS(xg_fit_final, paste0(outputdir, '/xgbTree/full_xgbTree.RDS'))
  
  #predict out for all locations
  covs[, 'xgbTree' := predict(xg_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, xg_best_tune, xg_fit, xg_fit_final, xg_grid, xg_grid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 9. Support Vector Machines with Radial Basis Function Kernel  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if('svmRadial' %in% child_models){
  dir.create(paste0(outputdir, '/svmRadial'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  ##tunegrid <- expand.grid(sigma = c(.01, .05, 0.2, 0.25), C = seq(0, 3, 0.1))
  tunegrid <- expand.grid(sigma = c(0.1,0.25,0.5), C = seq(0, 5, 0.5))
  
  ##tunegrid <- expand.grid(sigma = c(.01, .015, 0.2),C = c(0.75, 0.9, 1, 1.1, 1.25))
  
  # Fit model
  svmRadial_fit <- train(form,
                         data = mydata,
                         trControl = train_control,
                         verbose = T,
                         tuneGrid = tunegrid,
                         metric = "RMSE",
                         method = "svmRadial",
                         weights = mydata$w,
                         preProcess = c("center","scale"))
  
  # Save model fit object 
  saveRDS(svmRadial_fit, paste0(outputdir, "/svmRadial/svmRadial_fit.RDS"))
  png(paste0(outputdir, '/svmRadial/svmRadial_fit.png'))
  plot(svmRadial_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- svmRadial_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/svmRadial/svmRadial_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  svmRadial_fit_final <- train(form,
                               data = mydata,
                               trControl = train_control_final,
                               verbose = T,
                               tuneGrid = tunegrid_final,
                               metric = "RMSE",
                               method = "svmRadial",
                               importance=T,
                               weights = mydata$w,
                               preProcess = c("center","scale"))
  
  # Extract out of sample and in sample predictions
  mydata[, 'svmRadial_cv_pred'   := arrange(svmRadial_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'svmRadial_full_pred' := predict(svmRadial_fit_final, mydata)]
  
  #save model fit
  svmRadial_fit_final$model_name <- "svmRadial"
  saveRDS(svmRadial_fit_final, paste0(outputdir, '/svmRadial/full_svmRadial.RDS'))
  
  cov_plot <-
    ggplot(varImp(svmRadial_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/svmRadial/svmRadial_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'svmRadial' := predict(svmRadial_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, svmRadial_fit, svmRadial_fit_final, tunegrid)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 10. K-nearest neighbour model ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if('knn' %in% child_models){
  dir.create(paste0(outputdir, '/knn'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- data.frame(k = seq(1,25))
  
  # Fit model
  knn_fit <- train(form,
                   data = mydata,
                   trControl = train_control,
                   verbose = T,
                   tuneGrid = tunegrid,
                   metric = "RMSE",
                   method = "knn",
                   weights = mydata$w)
  
  # Save model fit object 
  saveRDS(knn_fit, paste0(outputdir, "/knn/knn_fit.RDS"))
  png(paste0(outputdir, '/knn/knn_fit.png'))
  plot(knn_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- knn_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/knn/knn_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  knn_fit_final <- train(form,
                         data = mydata,
                         trControl = train_control_final,
                         verbose = T,
                         tuneGrid = tunegrid_final,
                         metric = "RMSE",
                         method = "knn",
                         importance=T,
                         weights = mydata$w)
  
  
  # Extract out of sample and in sample predictions
  mydata[, 'knn_cv_pred'   := arrange(knn_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'knn_full_pred' := predict(knn_fit_final, mydata)]
  
  #save model fit
  knn_fit_final$model_name <- "knn"
  saveRDS(knn_fit_final, paste0(outputdir, '/knn/full_knn.RDS'))
  
  cov_plot <-
    ggplot(varImp(knn_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/knn/knn_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'knn' := predict(knn_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, knn_fit, knn_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 11. Generalized Additive Model using LOESS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if('gamLoess' %in% child_models){
  dir.create(paste0(outputdir, '/gamLoess'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(span = seq(0.15, 0.65, len = 10), degree = 1)
  
  # Fit model
  gamLoess_fit <- train(form,
                        data = mydata,
                        trControl = train_control,
                        verbose = T,
                        tuneGrid = tunegrid,
                        metric = "RMSE",
                        method = "gamLoess",
                        weights = mydata$w)
  
  # Save model fit object 
  saveRDS(gamLoess_fit, paste0(outputdir, "/gamLoess/gamLoess_fit.RDS"))
  png(paste0(outputdir, '/gamLoess/gamLoess_fit.png'))
  plot(gamLoess_fit)  
  dev.off()
  
  #specify the parameters
  ##mtry_tune <- gamLoess_fit$bestTune$mtry
  mtry_tune <- gamLoess_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/gamLoess/gamLoess_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  gamLoess_fit_final <- train(form,
                              data = mydata,
                              trControl = train_control_final,
                              verbose = T,
                              tuneGrid = tunegrid_final,
                              metric = "RMSE",
                              method = "gamLoess",
                              importance=T,
                              weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'gamLoess_cv_pred'   := arrange(gamLoess_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'gamLoess_full_pred' := predict(gamLoess_fit_final, mydata)]
  
  #save model fit
  gamLoess_fit_final$model_name <- "gamLoess"
  saveRDS(gamLoess_fit_final, paste0(outputdir, '/gamLoess/full_gamLoess.RDS'))
  
  cov_plot <-
    ggplot(varImp(gamLoess_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/gamLoess/gamLoess_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'gamLoess' := predict(gamLoess_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, gamLoess_fit, gamLoess_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 12. Bayesian Additive Regression Trees ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if('bartMachine' %in% child_models){
  dir.create(paste0(outputdir, '/bartMachine'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(num_trees = c(10, 15, 20, 100), k = 2, alpha = 0.95, beta = 2, nu = 3)
  
  # Fit model
  bartMachine_fit <- train(form,
                           data = mydata,
                           trControl = train_control,
                           verbose = T,
                           tuneGrid = tunegrid,
                           metric = "RMSE",
                           method = "bartMachine",
                           weights = mydata$w)
  
  # Save model fit object 
  saveRDS(bartMachine_fit, paste0(outputdir, "/bartMachine/bartMachine_fit.RDS"))
  png(paste0(outputdir, '/bartMachine/bartMachine_fit.png'))
  plot(bartMachine_fit)  
  dev.off()
  
  #specify the parameters
  ##mtry_tune <- bartMachine_fit$bestTune$mtry
  mtry_tune <- bartMachine_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/bartMachine/bartMachine_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  tunegrid_final 
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  bartMachine_fit_final <- train(form,
                                 data = mydata,
                                 trControl = train_control_final,
                                 verbose = T,
                                 tuneGrid = tunegrid_final,
                                 metric = "RMSE",
                                 method = "bartMachine",
                                 weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'bartMachine_cv_pred'   := arrange(bartMachine_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'bartMachine_full_pred' := predict(bartMachine_fit_final, mydata)]
  
  #save model fit
  bartMachine_fit_final$model_name <- "bartMachine"
  saveRDS(bartMachine_fit_final, paste0(outputdir, '/bartMachine/full_bartMachine.RDS'))
  
  cov_plot <-
    ggplot(varImp(bartMachine_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/bartMachine/bartMachine_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'bartMachine' := predict(bartMachine_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, bartMachine_fit, bartMachine_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 13. Gaussian Process ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if('gaussprLinear' %in% child_models){
  dir.create(paste0(outputdir, '/gaussprLinear'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  gaussprLinear_fit_final <- train(form,
                                   data = mydata,
                                   trControl = train_control_final,
                                   verbose = T,
                                   metric = "RMSE",
                                   method = "gaussprLinear",
                                   weights = mydata$w)
  
  # Save model fit object 
  saveRDS(gaussprLinear_fit_final, paste0(outputdir, "/gaussprLinear/gaussprLinear_fit.RDS"))
  dev.off()
  
  # Extract out of sample and in sample predictions
  mydata[, 'gaussprLinear_cv_pred'   := arrange(gaussprLinear_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'gaussprLinear_full_pred' := predict(gaussprLinear_fit_final, mydata)]
  
  #save model fit
  gaussprLinear_fit_final$model_name <- "gaussprLinear"
  saveRDS(gaussprLinear_fit_final, paste0(outputdir, '/gaussprLinear/full_gaussprLinear.RDS'))
  
  cov_plot <-
    ggplot(varImp(gaussprLinear_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/gaussprLinear/gaussprLinear_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'gaussprLinear' := predict(gaussprLinear_fit_final, covs)]
  
  ##rm(form, cov_plot, train_control_final, gaussprLinear_fit_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 14. Gaussian Process with Polynomial Kernel ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if('gaussprPoly' %in% child_models){
  dir.create(paste0(outputdir, '/gaussprPoly'), showWarnings = F)
  
  #set response variable
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if (family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  # Create model formula
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <-  expand.grid(degree = seq(1, min(5, 3)),scale = 10 ^((1:5) - 4))
  
  
  # Fit model
  gaussprPoly_fit <- train(form,
                           data = mydata,
                           trControl = train_control,
                           verbose = T,
                           tuneGrid = tunegrid,
                           metric = "RMSE",
                           method = "gaussprPoly",
                           weights = mydata$w)
  
  # Save model fit object 
  saveRDS(gaussprPoly_fit, paste0(outputdir, "/gaussprPoly/gaussprPoly_fit.RDS"))
  png(paste0(outputdir, '/gaussprPoly/gaussprPoly_fit.png'))
  plot(gaussprPoly_fit)  
  dev.off()
  
  #specify the parameters
  ##mtry_tune <- gaussprPoly_fit$bestTune$mtry
  mtry_tune <- gaussprPoly_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/gaussprPoly/gaussprPoly_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  gaussprPoly_fit_final <- train(form,
                                 data = mydata,
                                 trControl = train_control_final,
                                 verbose = T,
                                 tuneGrid = tunegrid_final,
                                 metric = "RMSE",
                                 method = "gaussprPoly",
                                 importance=T,
                                 weights = mydata$w)
  
  min<-min(gaussprPoly_fit_final$pred$pred)
  max<-max(gaussprPoly_fit_final$pred$pred)
  gaussprPoly_fit_final$pred$pred<- (gaussprPoly_fit_final$pred$pred-min)/(max-min)
  
  # Extract out of sample and in sample predictions
  mydata[, 'gaussprPoly_cv_pred'   := arrange(gaussprPoly_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'gaussprPoly_full_pred' := predict(gaussprPoly_fit_final, mydata)]
  
  #save model fit
  gaussprPoly_fit_final$model_name <- "gaussprPoly"
  saveRDS(gaussprPoly_fit_final, paste0(outputdir, '/gaussprPoly/full_gaussprPoly.RDS'))
  
  cov_plot <-
    ggplot(varImp(gaussprPoly_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/gaussprPoly/gaussprPoly_covariate_importance.png'),
         
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'gaussprPoly' := predict(gaussprPoly_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, gaussprPoly_fit, gaussprPoly_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 15. Gradient Boosting Machines ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if('gbm_h2o' %in% child_models){
  
  h2o.init()
  
  dir.create(paste0(outputdir, '/gbm_h2o'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(ntrees=150,max_depth = c(4,6,8,10),min_rows = 2,learn_rate = 0.01,col_sample_rate = 0.8)
  
  warnings()
  # Fit model
  gbm_h2o_fit <- train(form,
                       data = mydata,
                       trControl = train_control,
                       verbose = T,
                       tuneGrid = tunegrid,
                       metric = "RMSE",
                       method = "gbm_h2o",
                       weights = mydata$w)
  
  # Save model fit object 
  saveRDS(gbm_h2o_fit, paste0(outputdir, "/gbm_h2o/gbm_h2o_fit.RDS"))
  png(paste0(outputdir, '/gbm_h2o/gbm_h2o_fit.png'))
  plot(gbm_h2o_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- gbm_h2o_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/gbm_h2o/gbm_h2o_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  gbm_h2o_fit_final <- train(form,
                             data = mydata,
                             trControl = train_control_final,
                             verbose = T,
                             tuneGrid = tunegrid_final,
                             metric = "RMSE",
                             method = "gbm_h2o",
                             weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'gbm_h2o_cv_pred'   := arrange(gbm_h2o_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'gbm_h2o_full_pred' := predict(gbm_h2o_fit_final, mydata)]
  
  #save model fit
  gbm_h2o_fit_final$model_name <- "gbm_h2o"
  saveRDS(gbm_h2o_fit_final, paste0(outputdir, '/gbm_h2o/full_gbm_h2o.RDS'))
  
  cov_plot <-
    ggplot(varImp(gbm_h2o_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/gbm_h2o/gbm_h2o_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'gbm_h2o' := predict(gbm_h2o_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, gbm_h2o_fit, gbm_h2o_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 16. The Bayesian lasso ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if('blasso' %in% child_models){
  dir.create(paste0(outputdir, '/blasso'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <-  expand.grid(sparsity = seq(.3, .7, length = 10))
  
  # Fit model
  blasso_fit <- train(form,
                      data = mydata,
                      trControl = train_control,
                      tuneGrid = tunegrid,
                      metric = "RMSE",
                      method = "blasso",
                      weights = mydata$w)
  
  # Save model fit object 
  saveRDS(blasso_fit, paste0(outputdir, "/blasso/blasso_fit.RDS"))
  png(paste0(outputdir, '/blasso/blasso_fit.png'))
  plot(blasso_fit)  
  dev.off()
  
  #specify the parameters
  ##mtry_tune <- blasso_fit$bestTune$mtry
  mtry_tune <- blasso_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/blasso/blasso_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  blasso_fit_final <- train(form,
                            data = mydata,
                            trControl = train_control_final,
                            tuneGrid = tunegrid_final,
                            metric = "RMSE",
                            method = "blasso",
                            weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'blasso_cv_pred'   := arrange(blasso_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'blasso_full_pred' := predict(blasso_fit_final, mydata)]
  
  #save model fit
  blasso_fit_final$model_name <- "blasso"
  saveRDS(blasso_fit_final, paste0(outputdir, '/blasso/full_blasso.RDS'))
  
  cov_plot <-
    ggplot(varImp(blasso_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/blasso/blasso_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'blasso' := predict(blasso_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, blasso_fit, blasso_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~#
# 17. gbm model ####
#~~~~~~~~~~~~~~~~~~~~#
if('gbm' %in% child_models){
  dir.create(paste0(outputdir, '/gbm'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  #tunegrid <- expand.grid(interaction.depth = seq(1, 7, by = 2),n.trees = seq(100, 1000, by = 50),n.minobsinnode = 10, shrinkage = c(0.01, 0.1))
  
  tunegrid <- expand.grid(n.trees = (1:30)*50, interaction.depth = c(5, 9, 13, 17), shrinkage = 0.1,n.minobsinnode = 20)
  
  # Fit model
  gbm_fit <- train(form,
                   data = mydata,
                   trControl = train_control,
                   tuneGrid = tunegrid,
                   metric = "RMSE",
                   method = "gbm",
                   weights = mydata$w)
  
  # Save model fit object 
  saveRDS(gbm_fit, paste0(outputdir, "/gbm/gbm_fit.RDS"))
  png(paste0(outputdir, '/gbm/gbm_fit.png'))
  plot(gbm_fit)  
  dev.off()
  
  # Save the best parameters to csv file
  write.csv(gbm_fit$bestTune, paste0(outputdir, '/gbm/gbm_best_tune.csv'))
  
  #specify the parameters
  #tunegrid_final <- expand.grid(.committees=gbm_fit$bestTune$committees, .neighbors=gbm_fit$bestTune$neighbors)
  mtry_tune <- gbm_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/gbm/gbm_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  gbm_fit_final <- train(form,
                         data = mydata,
                         trControl = train_control_final,
                         tuneGrid = tunegrid_final,
                         metric = "RMSE",
                         method = "gbm",
                         weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'gbm_cv_pred'   := arrange(gbm_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'gbm_full_pred' := predict(gbm_fit_final, mydata)]
  
  #save model fit
  gbm_fit_final$model_name <- "gbm"
  saveRDS(gbm_fit_final, paste0(outputdir, '/gbm/full_gbm.RDS'))
  
  cov_plot <-
    ggplot(varImp(gbm_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/gbm/gbm_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'gbm' := predict(gbm_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, gbm_fit, gbm_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~#
# 18. k-Nearest Neighbors ####
#~~~~~~~~~~~~~~~~~~~~~#
if('kknn' %in% child_models){
  dir.create(paste0(outputdir, '/kknn'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tuneGrid <- expand.grid(kmax = seq(5, 100, length = 10),            # allows to test a range of k values
                          distance = 1:5,        # allows to test a range of distance values
                          kernel = c('gaussian',  # different weighting types in kknn
                                     'triangular',
                                     'rectangular',
                                     'epanechnikov',
                                     'optimal'))
  # Fit model
  kknn_fit <- train(form,
                    data = mydata,
                    trControl = train_control,
                    tuneGrid = tuneGrid,
                    metric = "RMSE",
                    method = "kknn",
                    weights = mydata$w)
  
  # Save model fit object 
  saveRDS(kknn_fit, paste0(outputdir, "/kknn/kknn_fit.RDS"))
  png(paste0(outputdir, '/kknn/kknn_fit.png'))
  plot(kknn_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- kknn_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/kknn/kknn_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  kknn_fit_final <- train(form,
                          data = mydata,
                          trControl = train_control_final,
                          tuneGrid = tunegrid_final,
                          metric = "RMSE",
                          method = "kknn",
                          ##                        importance=T,
                          weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'kknn_cv_pred'   := arrange(kknn_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'kknn_full_pred' := predict(kknn_fit_final, mydata)]
  
  #save model fit
  kknn_fit_final$model_name <- "kknn"
  saveRDS(kknn_fit_final, paste0(outputdir, '/kknn/full_kknn.RDS'))
  
  cov_plot <-
    ggplot(varImp(kknn_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/kknn/kknn_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'kknn' := predict(kknn_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, kknn_fit, kknn_fit_final, tuneGrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~#
# 19. Bagged MARS ####
#~~~~~~~~~~~~~~~~~~~~~#
if('bagEarth' %in% child_models){
  dir.create(paste0(outputdir, '/bagEarth'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <-expand.grid(degree = 1, nprune = (1:10) * 2)
  
  # Fit model
  bagEarth_fit <- train(form,
                        data = mydata,
                        trControl = train_control,
                        tuneGrid = tunegrid,
                        metric = "RMSE",
                        method = "bagEarth",
                        weights = mydata$w)
  
  # Save model fit object 
  saveRDS(bagEarth_fit, paste0(outputdir, "/bagEarth/bagEarth_fit.RDS"))
  png(paste0(outputdir, '/bagEarth/bagEarth_fit.png'))
  plot(bagEarth_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- bagEarth_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/bagEarth/bagEarth_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  bagEarth_fit_final <- train(form,
                              data = mydata,
                              trControl = train_control_final,
                              tuneGrid = tunegrid_final,
                              metric = "RMSE",
                              method = "bagEarth",
                              weights = mydata$w)
  
  bagEarth_fit_final
  # Extract out of sample and in sample predictions
  mydata[, 'bagEarth_cv_pred'   := arrange(bagEarth_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'bagEarth_full_pred' := predict(bagEarth_fit_final, mydata)]
  
  #save model fit
  bagEarth_fit_final$model_name <- "bagEarth"
  saveRDS(bagEarth_fit_final, paste0(outputdir, '/bagEarth/full_bagEarth.RDS'))
  
  #predict out for all locations
  covs[, 'bagEarth' := predict(bagEarth_fit_final, covs)]
  
  rm(form, train_control, train_control_final, bagEarth_fit, bagEarth_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~#
# 20. Principal component regression #### 
#~~~~~~~~~~~~~~~~~~~~~#
if('pcr' %in% child_models){
  dir.create(paste0(outputdir, '/pcr'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <-expand.grid(ncomp = 1:10)
  
  # Fit model
  pcr_fit <- train(form,
                   data = mydata,
                   trControl = train_control,
                   tuneGrid = tunegrid,
                   metric = "RMSE",
                   method = "pcr",
                   weights = mydata$w)
  
  # Save model fit object 
  saveRDS(pcr_fit, paste0(outputdir, "/pcr/pcr_fit.RDS"))
  png(paste0(outputdir, '/pcr/pcr_fit.png'))
  plot(pcr_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- pcr_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/pcr/pcr_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  pcr_fit_final <- train(form,
                         data = mydata,
                         trControl = train_control_final,
                         tuneGrid = tunegrid_final,
                         metric = "RMSE",
                         method = "pcr",
                         ##                        importance=T,
                         weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'pcr_cv_pred'   := arrange(pcr_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'pcr_full_pred' := predict(pcr_fit_final, mydata)]
  
  #save model fit
  pcr_fit_final$model_name <- "pcr"
  saveRDS(pcr_fit_final, paste0(outputdir, '/pcr/full_pcr.RDS'))
  
  cov_plot <-
    ggplot(varImp(pcr_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/pcr/pcr_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'pcr' := predict(pcr_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, pcr_fit, pcr_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~#
# 21. Linear Regression with Stepwise Selection####
#~~~~~~~~~~~~~~~~~~~~~#
if('lmStepAIC' %in% child_models){
  dir.create(paste0(outputdir, '/lmStepAIC'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  # Fit model
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  lmStepAIC_fit_final <- train(form,
                               data = mydata,
                               trControl = train_control_final,
                               metric = "RMSE",
                               method = "lmStepAIC",
                               weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'lmStepAIC_cv_pred'   := arrange(lmStepAIC_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'lmStepAIC_full_pred' := predict(lmStepAIC_fit_final, mydata)]
  
  #save model fit
  lmStepAIC_fit_final$model_name <- "lmStepAIC"
  saveRDS(lmStepAIC_fit_final, paste0(outputdir, '/lmStepAIC/full_lmStepAIC.RDS'))
  
  cov_plot <-
    ggplot(varImp(lmStepAIC_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/lmStepAIC/lmStepAIC_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'lmStepAIC' := predict(lmStepAIC_fit_final, covs)]
  
  rm(form, cov_plot, train_control_final,  lmStepAIC_fit_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 22. Support Vector Machines with Linear Kernel  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if('svmLinear' %in% child_models){
  dir.create(paste0(outputdir, '/svmLinear'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(C = c(0.1, 1, 10, 100, 1000))
  
  # Fit model
  svmLinear_fit <- train(form,
                         data = mydata,
                         trControl = train_control,
                         verbose = T,
                         tuneGrid = tunegrid,
                         metric = "RMSE",
                         method = "svmLinear",
                         weights = mydata$w,
                         preProcess = c("center","scale"))
  
  # Save model fit object 
  saveRDS(svmLinear_fit, paste0(outputdir, "/svmLinear/svmLinear_fit.RDS"))
  png(paste0(outputdir, '/svmLinear/svmLinear_fit.png'))
  plot(svmLinear_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- svmLinear_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/svmLinear/svmLinear_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  svmLinear_fit_final <- train(form,
                               data = mydata,
                               trControl = train_control_final,
                               verbose = T,
                               tuneGrid = tunegrid_final,
                               metric = "RMSE",
                               method = "svmLinear",
                               weights = mydata$w,
                               preProcess = c("center","scale"))
  
  # Extract out of sample and in sample predictions
  mydata[, 'svmLinear_cv_pred'   := arrange(svmLinear_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'svmLinear_full_pred' := predict(svmLinear_fit_final, mydata)]
  
  #save model fit
  svmLinear_fit_final$model_name <- "svmLinear"
  saveRDS(svmLinear_fit_final, paste0(outputdir, '/svmLinear/full_svmLinear.RDS'))
  
  cov_plot <-
    ggplot(varImp(svmLinear_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/svmLinear/svmLinear_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'svmLinear' := predict(svmLinear_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, svmLinear_fit, svmLinear_fit_final, tunegrid)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 23. Support Vector Machines with Polynomial Kernel  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if('svmPoly' %in% child_models){
  dir.create(paste0(outputdir, '/svmPoly'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  C <- c(0.1,1,10,100)
  degree <- c(1,2,3)
  scale <- 1
  
  tunegrid <- expand.grid(C=C,degree=degree,scale=scale)
  
  # Fit model
  svmPoly_fit <- train(form,
                       data = mydata,
                       trControl = train_control,
                       verbose = T,
                       tuneGrid = tunegrid,
                       metric = "RMSE",
                       method = "svmPoly",
                       weights = mydata$w,
                       preProcess = c("center","scale"))
  
  # Save model fit object 
  saveRDS(svmPoly_fit, paste0(outputdir, "/svmPoly/svmPoly_fit.RDS"))
  png(paste0(outputdir, '/svmPoly/svmPoly_fit.png'))
  plot(svmPoly_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- svmPoly_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/svmPoly/svmPoly_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  svmPoly_fit_final <- train(form,
                             data = mydata,
                             trControl = train_control_final,
                             verbose = T,
                             tuneGrid = tunegrid_final,
                             metric = "RMSE",
                             method = "svmPoly",
                             weights = mydata$w,
                             preProcess = c("center","scale"))
  
  # Extract out of sample and in sample predictions
  mydata[, 'svmPoly_cv_pred'   := arrange(svmPoly_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'svmPoly_full_pred' := predict(svmPoly_fit_final, mydata)]
  
  #save model fit
  svmPoly_fit_final$model_name <- "svmPoly"
  saveRDS(svmPoly_fit_final, paste0(outputdir, '/svmPoly/full_svmPoly.RDS'))
  
  cov_plot <-
    ggplot(varImp(svmPoly_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/svmPoly/svmPoly_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'svmPoly' := predict(svmPoly_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, svmPoly_fit, svmPoly_fit_final, tunegrid)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 24. glmnet ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

if('glmnet_h2o' %in% child_models){
  
  h2o.init()
  
  dir.create(paste0(outputdir, '/glmnet_h2o'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(alpha = 0:1, lambda = seq(0.0001, 1, length = 10))
  
  # Fit model
  glmnet_h2o_fit <- train(form,
                          data = mydata,
                          trControl = train_control,
                          tuneGrid = tunegrid,
                          metric = "RMSE",
                          method = "glmnet_h2o",
                          weights = mydata$w)
  
  # Save model fit object 
  saveRDS(glmnet_h2o_fit, paste0(outputdir, "/glmnet_h2o/glmnet_h2o_fit.RDS"))
  png(paste0(outputdir, '/glmnet_h2o/glmnet_h2o_fit.png'))
  plot(glmnet_h2o_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- glmnet_h2o_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/glmnet_h2o/glmnet_h2o_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  glmnet_h2o_fit_final <- train(form,
                                data = mydata,
                                trControl = train_control_final,
                                tuneGrid = tunegrid_final,
                                metric = "RMSE",
                                method = "glmnet_h2o",
                                weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'glmnet_h2o_cv_pred'   := arrange(glmnet_h2o_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'glmnet_h2o_full_pred' := predict(glmnet_h2o_fit_final, mydata)]
  
  #save model fit
  glmnet_h2o_fit_final$model_name <- "glmnet_h2o"
  saveRDS(glmnet_h2o_fit_final, paste0(outputdir, '/glmnet_h2o/full_glmnet_h2o.RDS'))
  
  #predict out for all locations
  covs[, 'glmnet_h2o' := predict(glmnet_h2o_fit_final, covs)]
  
  rm(form, train_control, train_control_final, glmnet_h2o_fit, glmnet_h2o_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 25. Bayesian Regularized Neural Networks ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if('brnn' %in% child_models){
  
  dir.create(paste0(outputdir, '/brnn'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(neurons=c(1,2,3,4,5))
  
  warnings()
  # Fit model
  brnn_fit <- train(form,
                    data = mydata,
                    trControl = train_control,
                    tuneGrid = tunegrid,
                    metric = "RMSE",
                    method = "brnn",
                    weights = mydata$w)
  
  # Save model fit object 
  saveRDS(brnn_fit, paste0(outputdir, "/brnn/brnn_fit.RDS"))
  png(paste0(outputdir, '/brnn/brnn_fit.png'))
  plot(brnn_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- brnn_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/brnn/brnn_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  brnn_fit_final <- train(form,
                          data = mydata,
                          trControl = train_control_final,
                          tuneGrid = tunegrid_final,
                          metric = "RMSE",
                          method = "brnn",
                          weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'brnn_cv_pred'   := arrange(brnn_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'brnn_full_pred' := predict(brnn_fit_final, mydata)]
  
  #save model fit
  brnn_fit_final$model_name <- "brnn"
  saveRDS(brnn_fit_final, paste0(outputdir, '/brnn/full_brnn.RDS'))
  
  cov_plot <-
    ggplot(varImp(brnn_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/brnn/brnn_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'brnn' := predict(brnn_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, brnn_fit, brnn_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 26. Boosted Linear Model ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if('bstSm' %in% child_models){
  
  dir.create(paste0(outputdir, '/bstSm'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(mstop = 1000, nu = c(0.1,0.5,1,1.5))
  
  # Fit model
  bstSm_fit <- train(form,
                     data = mydata,
                     trControl = train_control,
                     tuneGrid = tunegrid,
                     metric = "RMSE",
                     method = "bstSm",
                     weights = mydata$w)
  
  # Save model fit object 
  saveRDS(bstSm_fit, paste0(outputdir, "/bstSm/bstSm_fit.RDS"))
  png(paste0(outputdir, '/bstSm/bstSm_fit.png'))
  plot(bstSm_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- bstSm_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/bstSm/bstSm_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  bstSm_fit_final <- train(form,
                           data = mydata,
                           trControl = train_control_final,
                           tuneGrid = tunegrid_final,
                           metric = "RMSE",
                           method = "bstSm",
                           weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'bstSm_cv_pred'   := arrange(bstSm_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'bstSm_full_pred' := predict(bstSm_fit_final, mydata)]
  
  #save model fit
  bstSm_fit_final$model_name <- "bstSm"
  saveRDS(bstSm_fit_final, paste0(outputdir, '/bstSm/full_bstSm.RDS'))
  
  cov_plot <-
    ggplot(varImp(bstSm_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/bstSm/bstSm_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'bstSm' := predict(bstSm_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, bstSm_fit, bstSm_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~#
# 27. icr ####
#~~~~~~~~~~~~~~~~#

if('icr' %in% child_models){
  
  dir.create(paste0(outputdir, '/icr'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid<-expand.grid(n.comp = c(2,5,10,15))
  
  # Fit model
  icr_fit <- train(form,
                   data = mydata,
                   trControl = train_control,
                   tuneGrid = tunegrid,
                   metric = "RMSE",
                   method = "icr",
                   weights = mydata$w)
  
  # Save model fit object 
  saveRDS(icr_fit, paste0(outputdir, "/icr/icr_fit.RDS"))
  png(paste0(outputdir, '/icr/icr_fit.png'))
  plot(icr_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- icr_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/icr/icr_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  icr_fit_final <- train(form,
                         data = mydata,
                         trControl = train_control_final,
                         tuneGrid = tunegrid_final,
                         metric = "RMSE",
                         method = "icr",
                         weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'icr_cv_pred'   := arrange(icr_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'icr_full_pred' := predict(icr_fit_final, mydata)]
  
  #save model fit
  icr_fit_final$model_name <- "icr"
  saveRDS(icr_fit_final, paste0(outputdir, '/icr/full_icr.RDS'))
  
  cov_plot <-
    ggplot(varImp(icr_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/icr/icr_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'icr' := predict(icr_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, icr_fit, icr_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~#
# 28. monmlp ####
#~~~~~~~~~~~~~~~~#

if('monmlp' %in% child_models){
  
  dir.create(paste0(outputdir, '/monmlp'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid<-expand.grid(hidden1=3, n.ensemble = c(5,10,15,20))
  
  # Fit model
  monmlp_fit <- train(form,
                      data = mydata,
                      trControl = train_control,
                      tuneGrid = tunegrid,
                      metric = "RMSE",
                      method = "monmlp",
                      weights = mydata$w)
  
  # Save model fit object 
  saveRDS(monmlp_fit, paste0(outputdir, "/monmlp/monmlp_fit.RDS"))
  png(paste0(outputdir, '/monmlp/monmlp_fit.png'))
  plot(monmlp_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- monmlp_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/monmlp/monmlp_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  monmlp_fit_final <- train(form,
                            data = mydata,
                            trControl = train_control_final,
                            tuneGrid = tunegrid_final,
                            metric = "RMSE",
                            method = "monmlp",
                            weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'monmlp_cv_pred'   := arrange(monmlp_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'monmlp_full_pred' := predict(monmlp_fit_final, mydata)]
  
  #save model fit
  monmlp_fit_final$model_name <- "monmlp"
  saveRDS(monmlp_fit_final, paste0(outputdir, '/monmlp/full_monmlp.RDS'))
  
  cov_plot <-
    ggplot(varImp(monmlp_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/monmlp/monmlp_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'monmlp' := predict(monmlp_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, monmlp_fit, monmlp_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~#
# 29. gamSpline ####
#~~~~~~~~~~~~~~~~#

if('gamSpline' %in% child_models){
  
  dir.create(paste0(outputdir, '/gamSpline'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid<-expand.grid(df = c(2,3,4,5,6,7,8,9,10))
  
  # Fit model
  gamSpline_fit <- train(form,
                         data = mydata,
                         trControl = train_control,
                         tuneGrid = tunegrid,
                         metric = "RMSE",
                         method = "gamSpline",
                         weights = mydata$w)
  
  # Save model fit object 
  saveRDS(gamSpline_fit, paste0(outputdir, "/gamSpline/gamSpline_fit.RDS"))
  png(paste0(outputdir, '/gamSpline/gamSpline_fit.png'))
  plot(gamSpline_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- gamSpline_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/gamSpline/gamSpline_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  gamSpline_fit_final <- train(form,
                               data = mydata,
                               trControl = train_control_final,
                               tuneGrid = tunegrid_final,
                               metric = "RMSE",
                               method = "gamSpline",
                               weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'gamSpline_cv_pred'   := arrange(gamSpline_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'gamSpline_full_pred' := predict(gamSpline_fit_final, mydata)]
  
  #save model fit
  gamSpline_fit_final$model_name <- "gamSpline"
  saveRDS(gamSpline_fit_final, paste0(outputdir, '/gamSpline/full_gamSpline.RDS'))
  
  cov_plot <-
    ggplot(varImp(gamSpline_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/gamSpline/gamSpline_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'gamSpline' := predict(gamSpline_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, gamSpline_fit, gamSpline_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~#
# 30. gamboost ####
#~~~~~~~~~~~~~~~~#

if('gamboost' %in% child_models){
  
  dir.create(paste0(outputdir, '/gamboost'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid<-expand.grid(mstop = c(50, 100, 150),prune = c("yes", "no"))
  
  # Fit model
  gamboost_fit <- train(form,
                        data = mydata,
                        trControl = train_control,
                        tuneGrid = tunegrid,
                        metric = "RMSE",
                        method = "gamboost",
                        weights = mydata$w)
  
  # Save model fit object 
  saveRDS(gamboost_fit, paste0(outputdir, "/gamboost/gamboost_fit.RDS"))
  png(paste0(outputdir, '/gamboost/gamboost_fit.png'))
  plot(gamboost_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- gamboost_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/gamboost/gamboost_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  gamboost_fit_final <- train(form,
                              data = mydata,
                              trControl = train_control_final,
                              tuneGrid = tunegrid_final,
                              metric = "RMSE",
                              method = "gamboost",
                              weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'gamboost_cv_pred'   := arrange(gamboost_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'gamboost_full_pred' := predict(gamboost_fit_final, mydata)]
  
  #save model fit
  gamboost_fit_final$model_name <- "gamboost"
  saveRDS(gamboost_fit_final, paste0(outputdir, '/gamboost/full_gamboost.RDS'))
  
  cov_plot <-
    ggplot(varImp(gamboost_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/gamboost/gamboost_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'gamboost' := predict(gamboost_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, gamboost_fit, gamboost_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~#
# 31. BstLm ####
#~~~~~~~~~~~~~~~~#

if('BstLm' %in% child_models){
  
  dir.create(paste0(outputdir, '/BstLm'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid<-expand.grid(mstop = c(50, 100, 150),nu = c(0.01,0.1))
  
  # Fit model
  BstLm_fit <- train(form,
                     data = mydata,
                     trControl = train_control,
                     tuneGrid = tunegrid,
                     metric = "RMSE",
                     method = "BstLm",
                     weights = mydata$w)
  
  # Save model fit object 
  saveRDS(BstLm_fit, paste0(outputdir, "/BstLm/BstLm_fit.RDS"))
  png(paste0(outputdir, '/BstLm/BstLm_fit.png'))
  plot(BstLm_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- BstLm_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/BstLm/BstLm_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  BstLm_fit_final <- train(form,
                           data = mydata,
                           trControl = train_control_final,
                           tuneGrid = tunegrid_final,
                           metric = "RMSE",
                           method = "BstLm",
                           weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'BstLm_cv_pred'   := arrange(BstLm_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'BstLm_full_pred' := predict(BstLm_fit_final, mydata)]
  
  #save model fit
  BstLm_fit_final$model_name <- "BstLm"
  saveRDS(BstLm_fit_final, paste0(outputdir, '/BstLm/full_BstLm.RDS'))
  
  cov_plot <-
    ggplot(varImp(BstLm_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/BstLm/BstLm_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'BstLm' := predict(BstLm_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, BstLm_fit, BstLm_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~#
# 13. Gaussian Process####
#~~~~~~~~~~~~~~~~~~~~~#
if('gaussprLinear' %in% child_models){
  dir.create(paste0(outputdir, '/gaussprLinear'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  # Fit model
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  gaussprLinear_fit_final <- train(form,
                                   data = mydata,
                                   trControl = train_control_final,
                                   metric = "RMSE",
                                   method = "gaussprLinear",
                                   weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'gaussprLinear_cv_pred'   := arrange(gaussprLinear_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'gaussprLinear_full_pred' := predict(gaussprLinear_fit_final, mydata)]
  
  #save model fit
  gaussprLinear_fit_final$model_name <- "gaussprLinear"
  saveRDS(gaussprLinear_fit_final, paste0(outputdir, '/gaussprLinear/full_gaussprLinear.RDS'))
  
  cov_plot <-
    ggplot(varImp(gaussprLinear_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/gaussprLinear/gaussprLinear_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'gaussprLinear' := predict(gaussprLinear_fit_final, covs)]
  
  rm(form, cov_plot, train_control_final,  gaussprLinear_fit_final)
}

#~~~~~~~~~~~~~~~~#
# 32. partDSA ####
#~~~~~~~~~~~~~~~~#

if('partDSA' %in% child_models){
  
  dir.create(paste0(outputdir, '/partDSA'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid<-expand.grid(cut.off.growth=c(10,15,20,25,30,35,40,45),MPD=c(0.01,0.05,0.1,0.15,0.2,0.25,0.5))
  
  # Fit model
  partDSA_fit <- train(form,
                       data = mydata,
                       trControl = train_control,
                       tuneGrid = tunegrid,
                       metric = "RMSE",
                       method = "partDSA",
                       weights = mydata$w)
  
  # Save model fit object 
  saveRDS(partDSA_fit, paste0(outputdir, "/partDSA/partDSA_fit.RDS"))
  png(paste0(outputdir, '/partDSA/partDSA_fit.png'))
  plot(partDSA_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- partDSA_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/partDSA/partDSA_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  partDSA_fit_final <- train(form,
                             data = mydata,
                             trControl = train_control_final,
                             tuneGrid = tunegrid_final,
                             metric = "RMSE",
                             method = "partDSA",
                             weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'partDSA_cv_pred'   := arrange(partDSA_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'partDSA_full_pred' := predict(partDSA_fit_final, mydata)]
  
  #save model fit
  partDSA_fit_final$model_name <- "partDSA"
  saveRDS(partDSA_fit_final, paste0(outputdir, '/partDSA/full_partDSA.RDS'))
  
  cov_plot <-
    ggplot(varImp(partDSA_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/partDSA/partDSA_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'partDSA' := predict(partDSA_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, partDSA_fit, partDSA_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~#
# 33. ranger ####
#~~~~~~~~~~~~~~~~#

if('ranger' %in% child_models){
  
  dir.create(paste0(outputdir, '/ranger'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  #tunegrid<-expand.grid(.mtry=c(1:length(covs_to_include)),.splitrule = c("variance","gini"),.min.node.size = c(5,10,15))
  tunegrid <- expand.grid(
    .mtry = seq_along(covs_to_include),
    .splitrule = if (family == "binomial") c("gini") else c("variance"),
    .min.node.size = c(5, 10, 15) #(5, 10, 15)
  )
  # Fit model
  ranger_fit <- train(form,
                      data = mydata,
                      trControl = train_control,
                      tuneGrid = tunegrid,
                      metric = "RMSE",
                      method = "ranger",
                      weights = mydata$w)
  #  num.trees = 1000)
  
  # Save model fit object 
  saveRDS(ranger_fit, paste0(outputdir, "/ranger/ranger_fit.RDS"))
  png(paste0(outputdir, '/ranger/ranger_fit.png'))
  plot(ranger_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- ranger_fit$bestTune
  write.csv(mtry_tune, paste0(outputdir, '/ranger/ranger_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  ranger_fit_final <- train(form,
                            data = mydata,
                            trControl = train_control_final,
                            tuneGrid = tunegrid_final,
                            metric = "RMSE",
                            method = "ranger",
                            weights = mydata$w)
  #num.trees = 1500
  #)
  
  # Extract out of sample and in sample predictions
  mydata[, 'ranger_cv_pred'   := arrange(ranger_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'ranger_full_pred' := predict(ranger_fit_final, mydata)]
  
  #save model fit
  ranger_fit_final$model_name <- "ranger"
  saveRDS(ranger_fit_final, paste0(outputdir, '/ranger/full_ranger.RDS'))
  
  # cov_plot <-
  #  ggplot(varImp(ranger_fit_final, scale = FALSE)) +
  #   labs(x = "Covariate", y = "Relative Importance") +
  #   theme_bw()
  #  ggsave(filename = paste0(outputdir, '/ranger/ranger_covariate_importance.png'),
  #       plot = cov_plot)
  
  #predict out for all locations
  covs[, 'ranger' := predict(ranger_fit_final, covs)]
  
  rm(form, train_control, train_control_final, ranger_fit, ranger_fit_final, tunegrid, tunegrid_final)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check the correlation of the stackers ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#if any of the preds are highly correlated the remove one of the correlated models

for(i in 1:length(child_models)){
  for(j in 1:length(child_models)){
    if(i==j){
    } else{
      if(cor(mydata[,get(paste0(child_models[i], '_cv_pred'))],mydata[,get(paste0(child_models[j], '_cv_pred'))])^2>0.8){message(paste0(child_models[i],  ' and ', child_models[j], ' correlated, remove one'))}
    }
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Print out correlations between data and predictions ####
# this is to aid selection of child models               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for(i in 1:length(child_models)){
  if(family == 'binomial'){
    message(paste0(child_models[i], ' correlation: ', round(cor(mydata$p, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)))
  }
  if(family == 'gaussian'){
    message(paste0(child_models[i], ' correlation: ', round(cor(mydata$n, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)))
  }
}

#save the temporary files
if(centre_scale == TRUE){
  mydata$year <-  NULL
  covs$year <- NULL
}
write.csv(mydata, paste0(outputdir, '/fitted_child_models.csv'), row.names = F)
write.csv(covs, paste0(outputdir, '/child_model_preds.csv'), row.names = F)

mydata <- read.csv(paste0(outputdir, '/fitted_child_models.csv'), stringsAsFactors =  F)
covs <- read.csv(paste0(outputdir, '/child_model_preds.csv'), stringsAsFactors = F)
mydata <- data.table(mydata)
covs <- data.table(covs)

if('gaussprPoly' %in% child_models){
  min<-min(covs$gaussprPoly)
  max<-max(covs$gaussprPoly)
  covs$gaussprPoly<- (covs$gaussprPoly-min)/(max-min)
}

if('gamLoess' %in% child_models){
  min<-min(covs$gamLoess)
  max<-max(covs$gamLoess)
  covs$gamLoess<- (covs$gamLoess-min)/(max-min)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Combined the child model estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# what methods to use to create the stackers
# use either a GPR or constrained weighted mean or other model
# ensure beta coefficients are contrained to sum to 1

# specify child models to use (if need to remove some due to correlation), want max 5ish, remove others from dataframe
# child_models <- child_models[child_models!='xgboost']
# child_models <- child_models[child_models!='gam']
# child_models <- child_models[child_models!='lasso']
# child_models <- child_models[child_models!='ridge']
# child_models <- child_models[child_models!='enet']
# child_models <- child_models[child_models!='rf']
# child_models <- child_models[child_models!='nnet']
# child_models <- child_models[child_models!='cubist']

#remove unwanted child models from the data frame to calculate child model weights
stackers <- data.table(mydata)
##stackers[, (colnames(stackers)[grep('xgboost', colnames(stackers))]) := NULL]
# stackers[, (colnames(stackers)[grep('gam', colnames(stackers))]) := NULL]
# stackers[, (colnames(stackers)[grep('lasso', colnames(stackers))]) := NULL]
# stackers[, (colnames(stackers)[grep('ridge', colnames(stackers))]) := NULL]
# stackers[, (colnames(stackers)[grep('enet', colnames(stackers))]) := NULL]
# stackers[, (colnames(stackers)[grep('rf', colnames(stackers))]) := NULL]
# stackers[, (colnames(stackers)[grep('nnet', colnames(stackers))]) := NULL]
# stackers[, (colnames(stackers)[grep('cubist', colnames(stackers))]) := NULL]

stackers <- data.frame(stackers)
X <- as.matrix(stackers[colnames(stackers)[(grep('cv_pred', colnames(stackers)))]])

Y = if(family == 'binomial'){stackers$p}else if(family == 'gaussian'){stackers$n}

#stack the predictions for all locations
C <- data.frame(covs)
C <- as.matrix(C[c(child_models)])

#Select which stacker you are using and stack the estimates

if(stacker == 'CWM'){
  # The following code is from http://zoonek.free.fr/blosxom/R/2012-06-01_Optimization.html
  # Calculate coefficients (child stacker weights)
  # Coefficients must sum to 1 and >=0
  
  s <- solve.QP( 
    t(X) %*% X,t(Y) %*% X, 
    cbind(  # One constraint per COLUMN
      matrix(1, nr=length(child_models), nc=1),
      diag(length(child_models)),
      -diag(length(child_models))
    ),
    c(1, 
      rep(0.000001, length(child_models)),
      rep(-1, length(child_models))), 
    meq = 0 # Only the first constraint is an equality if meq = 1 so set to 0, the others are >=
  )
  
  #calculate weighted stackers
  #for the fits with data
  mydata$stacked_preds <- rowWeightedMeans(X, w = s$solution)   # or can use: mydata$crossprod <- crossprod(t(X), s$solution)
  
  #Calculate the stacked predictions
  covs$cv_custom_stage_1 <- rowWeightedMeans(C, w = s$solution)
}

if(stacker == 'RWM'){
  r2 <- rep(NA, length(child_models))
  for(i in 1:length(child_models)){
    if(family == 'binomial'){
      r2[i] <- round(cor(mydata$p, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)
    }
    if(family == 'gaussian'){
      r2[i] <- round(cor(mydata$n, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)
    }
  }
  
  total <-  sum(r2)
  weights <- r2/total
  mydata$stacked_preds <- rowWeightedMeans(X, w = weights)   
  covs$cv_custom_stage_1 <- rowWeightedMeans(C, w = weights)
}

if(stacker == 'GBM'){
  if(family == 'gaussian'){
    form <- as.formula(paste0('ddd_per_1000 ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  #if(family == 'binomial'){
  #  form <- as.formula(paste0('logit(p) ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  #}
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  model_gbm<- train(form, data = mydata, method='gbm',trControl = train_control_final, tuneLength=3)
  
  ##mydata[, 'stacked_preds'   := inv.logit(predict(model_gbm, mydata))]
  covs <- data.frame(covs)
  colnames(covs)[colnames(covs) %in% child_models] <- colnames(stackers)[(grep('cv_pred', colnames(stackers)))]
  covs$cv_custom_stage_1 <- inv.logit(predict(model_gbm, covs[colnames(covs)[(grep('cv_pred', colnames(covs)))]])) 
  covs <- data.table(covs)
  
  #save the covariate importance plot
  jpeg(paste0(outputdir, 'stacker_cov_importance.jpeg'),
       height = 10, width = 10, units = 'cm', res = 150)
  ggplot(varImp(model_gbm, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  
  dev.off()
}

if(stacker == 'GLM'){
  if(family == 'gaussian'){
    form <- as.formula(paste0('ddd_per_1000 ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  model_glm<- 
    train(form, data = mydata, method='glm', trControl=train_control_final, tuneLength=3)
  
  mydata[, 'stacked_preds'   := predict(model_glm, mydata)]
  
  covs <- data.frame(covs)
  colnames(covs)[colnames(covs) %in% child_models] <- colnames(stackers)[(grep('cv_pred', colnames(stackers)))]
  covs$cv_custom_stage_1 <- predict(model_glm, covs[colnames(covs)[(grep('cv_pred', colnames(covs)))]]) 
  covs <- data.table(covs)
}

if(stacker == 'nnet'){
  if(family == 'gaussian'){
    form <- as.formula(paste0('ddd_per_1000 ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  nnetGrid <-  expand.grid(size = seq(from = 1, to = 10, by = 1),
                           decay = seq(from = 0.1, to = 0.5, by = 0.1))
  model_nnet0<- 
    train(form, data = mydata, method='nnet', trControl=train_control_final, tuneGrid = nnetGrid)
  
  #specify the parameters
  tunegrid_final <- expand.grid(.decay=model_nnet0$bestTune$decay, .size=model_nnet0$bestTune$size)
  tunegrid_final
  #fit the final model
  model_nnet <- train(form,
                      data = mydata,
                      trControl = train_control_final,
                      tuneGrid = tunegrid_final,
                      method = "nnet",
                      tuneLength=3)
  
  model_nnet$pred
  
  mydata[, 'stacked_preds'   := arrange(model_nnet$obs, rowIndex)[,"pred"]]
  mydata[, 'stacked_preds'   := predict(model_nnet, mydata)]
  
  head(mydata)
  
  covs <- data.frame(covs)
  colnames(covs)[colnames(covs) %in% child_models] <- colnames(stackers)[(grep('cv_pred', colnames(stackers)))]
  covs$cv_custom_stage_1 <- predict(model_nnet, covs[colnames(covs)[(grep('cv_pred', colnames(covs)))]]) 
  covs <- data.table(covs)
}

#clean up the covs dataset
stg1 <-  covs[, .(location_id, year_id,
                  age_group_id = rep(22, length(covs$location_id)),
                  sex_id = rep(3, length(covs$location_id)),
                  cv_custom_stage_1 = cv_custom_stage_1)]

#remove covariates from the dataset
#mydata[, colnames(mydata)[grep('^cv_', colnames(mydata))] := NULL]

#check that the estimates are within expected range
#max(stg1$cv_custom_stage_1, na.rm = T)

#save prediction
write.csv(mydata, paste0(outputdir, '/fitted_stackers.csv'), row.names = F)
write.csv(stg1, paste0(outputdir, '/custom_stage1_df.csv'), row.names = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calculate metrics for all models ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#get the RMSE and R2 of each model into a data frame and save it
child_model_metrics <- data.frame(child_models)
child_model_metrics$r_sq <- NA

mydata <- data.frame(mydata)
for(i in 1:length(child_models)){
  cm <- paste0(child_models, '_cv_pred')[i]
  pred <- mydata[c(cm)]
  pred <-  unlist(pred)
  child_model_metrics$rmse[child_model_metrics$child_models == child_models[i]] <- round(RMSE(pred, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n}),4)
  child_model_metrics$r_sq[child_model_metrics$child_models == child_models[i]] <- round(cor(pred, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n})^2,2)
}

child_model_metrics$child_models <- as.character(child_model_metrics$child_models)
child_model_metrics <- rbind(child_model_metrics, c('Stackers', NA, NA) )
child_model_metrics$rmse[child_model_metrics$child_models=='Stackers'] <- round(RMSE(mydata$stacked_preds, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n}),4)
child_model_metrics$r_sq[child_model_metrics$child_models=='Stackers'] <- round(cor(mydata$stacked_preds, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n})^2,2)

write.csv(child_model_metrics, paste0(outputdir, '/national_stacker_metrics.csv'), row.names = F)

#~~~~~~~~~~~~~~~~~~~~~#
# Plot the results ####
#~~~~~~~~~~~~~~~~~~~~~#
library(foreign)

#locs <- read.dbf("Z:/AMR/Shapefiles/GBD_2025/GBD2025_analysis_final.dbf")

library(sf)
library(sp)

setwd("C:/Users/naekaphirat/Desktop/KP/shp")
locs <- st_read("GBD2023_analysis_final_loc_set_22.shp")
setwd("C:/Users/naekaphirat/Desktop/KP")

covs <- data.table(covs)

# Do any need reverse transforming?
# covs[, gam := inv.logit(gam)]
# covs[, enet := inv.logit(enet)]
# covs[, rf := inv.logit(rf)]
# covs[, nnet := exp(nnet)]
# covs[, cubist := exp(cubist)]
# covs[, xgboost := exp(xgboost)]
#covs[, cv_custom_stage_1 := exp(cv_custom_stage_1)]

#merge on the data points
mydata <- data.frame(mydata)
input <- mydata[c('location_id', 'year_id', 'n')]
covs <- merge(covs, input, by = c('location_id', 'year_id'), all.x = T, all.y = T, allow.cartesian = T)

colnames(covs) <-  gsub('_cv_pred', '', colnames(covs))
#reshape long
covs <- melt(covs, id.vars = c('location_id', 'year_id','n'),
             measure.vars = c(child_models, 'cv_custom_stage_1'))

#merge on locations
covs <- merge(covs, locs, by.x = 'location_id', by.y = 'loc_id', all.x = T, all.y = F)
covs <- covs[covs$level == 3,]

# Plot out the stg1s data and check
pdf(paste0(outputdir, 'stacker_results.pdf'),
    height = 8.3, width = 11.7)

#plot out a page for each region
for(i in 1:length(unique(covs$region_id))){
  subset <- covs[covs$region_id == unique(covs$region_id)[i],]
  print(
    ggplot(subset)+
      geom_line(aes(x=year_id, y = value, group = variable, colour = variable))+
      geom_point(aes(x=year_id, y = (n)))+
      facet_wrap(~loc_name, nrow = ceiling(sqrt(length(unique(subset$location_id)))))+
      ylim(0,10)
    
  )
}
dev.off()

