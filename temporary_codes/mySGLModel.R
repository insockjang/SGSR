require(glmnet)
mySGLModel <- setRefClass(Class = "mySGLModel",
                          contains="PredictiveModel",
                          fields="model",
                          methods = list(
                            initialize = function(...){
                              return(.self)
                            },
                            
                            rawModel = function(){
                              return(.self$model)
                            },
                            
                            customTrain = function(featureData, responseData, alpha = 1, lambda = lambda, nfolds = 5,iterations = 20,groups = NA, ...){
                              # groups should be list 
                              # each group must contains genes (each group = gene set)
                              penalty<-rep(1,ncol(featureData))
                              path<-c()                               
                              mse<-list()
                              k1 <- 0
                              while(k1 <= iterations){
                                k1<-k1+1
                                MSE<-c()
                                set.seed(2)
                                fit<-cv.glmnet(featureData,responseData,alpha=1,lambda = lambdas,nfolds=5,penalty.factor = penalty)
                                a<-min(fit$cvm)
                                M<-foreach(kkk = 1:length(groups)) %dopar% {
                                  group <- groups[[kkk]]
                                  b<-match(group,colnames(featureData))
                                  if(length(which(is.na(b)==0))>0){
                                    penalty2<-penalty
                                    penalty2[b[which(is.na(b)==0)]]<-0
                                    set.seed(2)
                                    fit<-cv.glmnet(featureData,responseData,alpha=1,lambda = lambdas,nfolds=5,penalty.factor = penalty2)
                                    return(min(fit$cvm))
                                  }else{
                                    return(a)
                                  }     
                                }
                                MSE<-do.call("c",M)                                 
                                if(min(MSE)<a){
                                  group1<-groups[[which.min(MSE)]]
                                  path<-rbind(path,c(which.min(MSE),min(MSE)))
                                  b<-match(group1,colnames(featureData))
                                  penalty[b[which(is.na(b)==0)]]<-0
                                  mse[[k1]]<-MSE
                                }else{
                                  break
                                }
                              }
                              
                              .self$model <- cv.glmnet(featureData,responseData, alpha = alpha, lambda = lambda, nfolds = nfolds,penalty.factor = penalty,...) 
                              optParam <- c(.self$model$cvm[which.min(.self$model$cvm)],alpha,.self$model$lambda[which.min(.self$model$cvm)])
                              names(optParam)<-c("MSE","alpha","lambdaOpt")
                              .self$model$optParam <- optParam                               
                              .self$model$optPath <- path                               
                              
                            },
                            
                            customPredict = function(featureData){
                              predictedResponse <- predict(.self$model, featureData, s="lambda.min")
                              return(predictedResponse)
                            },
                            
                            getCoefficients = function(){
                              return(coef(.self$model,s = "lambda.min"))
                            }
                            
                          )
                          
)
