#' Load results from a caret train object saved in an .Rda file.
#' 
#' Loads an Rda file in the enclosing function environment.
#' Allows for assignment of the Rda file to a variable, without polluting the global environment.
#' 
#' @param dataset A character string specifying the dataset associated with the desired results object.
#' @param par A character string specifying the parameter or model associated with the desired results object.
#' @return The 'results' slot from a caret train object.
#' @examples
#' results <- loadResults("BreastCancer", "nnet")
#' @export
loadResults <- function(dataset, par){
load(sprintf("results_%s_%s.Rda", dataset, par))
results
}

#' Programmatically pass arguments to a script for training caret models.
#' 
#' @param script A character string specifying the name of the script you wish to run.
#' @param dataset A character string specifying the name of the dataset you are training the models on.
#' @param par Character string vector specifying the parameters of interest. Examples of parameters include different values of hyperparameters or different caret methods.
#' @return Response from Domino indicating if runs are succesfully triggered. 
#' @examples
#' models <- c("nnet", "rf", "gbm")
#' launcher(script = "variable_model.R", dataset = "BreastCancer", par = models)
#' @export
launcher <- function(script, dataset, par){
  lapply(par, function(m) domino.run(script, dataset, m))
}

#' Compare results from caret models trained on multiple instances on Domino.
#' 
#' Assumes results are downloaded as .Rda files to a single machine. Loads Rda and files and rbinds the results objects into a single data.frame for the sake of comparison.
#' 
#' @param dataset A character string specifying the name of the dataset you are training the models on.
#' @param par Character string vector specifying the parameters of interest. Examples of parameters include different values of hyperparameters or different caret methods.
#' @param metric A character string specifying the metric you wish to rank caret results by. 
#' @return data.frame containing results from caret train objects, ranked in descending order according to the chosen accuracy metric.
#' @examples
#' models <- c("nnet", "rf", "gbm")
#' compare(dataset = "BreastCancer", par = models, metric = "ROC")
#' @export
compare <- function(dataset, par, metric){
  result <- lapply(par, loadResults, dataset=dataset)
  result <- do.call(rbind, result)
  result[order(-result[[metric]]), ]
}


load_models <- function(dataset, method){
  load(sprintf("model_obj_%s_%s.Rda", dataset, method))
  model
}

#' Create a list of caret train objects trained on multiple Domino instances.
#' 
#' Assumes caret train objects are downloaded as .Rda files to a single machine. Loads Rda files and returns all caret models in a single list.
#' 
#' @param dataset A character string specifying the name of the dataset you are training the models on.
#' @param models Character string vector specifying the names of the models you wish to load.
#' @return a list of class "ensemble" containing caret train objects.
#' @examples
#' models <- c("nnet", "rf", "gbm")
#' modelList(dataset = "BreastCancer", models = models)
#' @export
modelList <- function(dataset, models){
  ens <- lapply(models, load_models, dataset=dataset)
  class(ens) <- "ensemble"
  ens
}

#' Predict method for objects of class "ensemble"
#' 
#' 
#' @param model_list An object of class "ensemble" containing a list of caret train objects
#' @param newdata A data.frame containing the features you wish to predict on. 
#' @param type A character string specifying the response; for now only supports "prob"
#' @param outcome An integer indicating the column index to select from each model's prediction matrix. 
#' @return a matrix of predictions, with each column vector corresponding to a different model.
#' @examples
#' model_list <- modelList(dataset="BreastCancer", models=models)
#' pred_mat <- predict(model_list, newdata=testing, type="prob", outcome=2)
#' @export
predict.ensemble <- function(model_list, newdata, type, outcome){
  
  pred_fact <- function(newdata, type, outcome){
    function(model){
      out <- predict(model, newdata=newdata, type=type)
      out <- out[,outcome]
      out
    }
  }
  
  
  predict_filt <- pred_fact(newdata=newdata, type = "prob", outcome=2)
  out <- sapply(model_list, predict_filt)
  if(class(out)=="numeric"){
    names(out) <- sapply(model_list, function(x) x$method)}else
    {colnames(out) <- sapply(model_list, function(x) x$method) }
  out
}


featureImp <- function(dataset, method){
  load(sprintf("model_obj_%s_%s.Rda", dataset, method))
  print(varImp(model))
}


featuRe <- function(dataset, models){
print(lapply(models, featureImp, dataset=dataset))
}





#used for API
live_blender <- function(dataset, models, Cl.thickness, Cell.size, Cell.shape, Marg.adhesion, Epith.c.size, Bare.nuclei, Bl.cromatin, Normal.nucleoli, Mitoses){
  df <- data.frame(Cl.thickness = Cl.thickness, Cell.size = Cell.size, Cell.shape = Cell.shape, Marg.adhesion = Marg.adhesion, Epith.c.size = Epith.c.size, Bare.nuclei = Bare.nuclei, Bl.cromatin =Bl.cromatin, Normal.nucleoli =Normal.nucleoli, Mitoses = Mitoses)
  df <- data.frame(lapply(df, as.factor))
  model_list <- modelList(dataset=dataset, models=models)
  pred_vec <- predict(model_list, newdata=df, type="prob", outcome=2)
  #prediction
  mean(pred_vec)
  }

#' Stack a set of models using a regularized logistic regression
#' 
#' 
#' @param models character vector specifying which models to stack.
#' @param train_mat A matrix containing the scored training data created from calling predict.ensemble on the training data, along with a model list. 
#' @param pred_mat A matrix containing the scored test data created from calling predict.ensemble on the testing data, along with a model list. 
#' @param cost An integer specifying the cost parameter of the regularized logistic regression. See package LiblineaR for details.
#' @param type An integer specifying the type of regularization. See package LiblineaR for details.
#' @return a numeric vector of probabilties indicating the probability of the target class.
#' @examples
#' pred_mat <- predict(model_list, newdata=testing, type="prob", outcome=2)
#' train_mat <- predict(model_list, newdata=training, type="prob", outcome=2)
#' log_ensemble <- stackerReg(models=models, train_mat=train_mat, pred_mat=pred_mat, cost=1, type=6)
#' @export
stackerReg <- function(models, train_mat, pred_mat, cost, type){
  l1m <- LiblineaR::LiblineaR(data = train_mat[ ,models], target = training$Class, type = type, cost = cost, epsilon = 0.01,
                   svr_eps = NULL, bias = TRUE, wi = NULL, cross = 0, verbose = FALSE)
  log_ensemble <- predict(l1m, newx=pred_mat[ ,models], proba=TRUE)#combine model predictions on test data, using logit model
  return(log_ensemble$probabilities[,2])
}


topModel <- function(n, ordered_acc, pred_mat, reference){
  if(n>1){
  top_mean <- apply(pred_mat[,ordered_acc$model[1:n]], 1, mean)
  pred <- ROCR::prediction(top_mean, reference)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  auc.tmp <- ROCR::performance(pred, "auc")
  auc <- as.numeric(auc.tmp@y.values)
  cbind(as.integer(n), auc, paste(ordered_acc$model[1:n], collapse=", "))  
  }else{
    top_mean <-pred_mat[,ordered_acc$model[1]]
    pred <- ROCR::prediction(top_mean, reference)
    perf <- ROCR::performance(pred, "tpr", "fpr")
    auc.tmp <- ROCR::performance(pred, "auc")
    auc <- as.numeric(auc.tmp@y.values)
    cbind(as.integer(n), auc, paste(ordered_acc$model[1]))}
}

topComb <- function(ordered_acc, pred_mat, reference){
  top_ens <- sapply(1:length(models), topModel, ordered_acc= ordered_acc, pred_mat=pred_mat, reference=reference)
  top_m <- data.frame(t(top_ens))
  colnames(top_m) <- c("number_models", "AUC", "model_names")
  top_m
 }

#' Calculate the error under the receiver operator curve and plot the ROC, given a vector of predicitons and a reference vector.
#' 
#' @param pred The name of the prediction vector.
#' @param reference The name of the reference vector.
#' @return Area under the receiver operator curve. Also plots the ROC as a side effect. 
#' @examples
#' AUC <- aucFun(pred, reference)
aucFun <- function(pred, reference) {
  pred <- ROCR::prediction(pred, reference)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  plot(perf)
  auc.tmp <- ROCR::performance(pred, "auc")
  auc <- as.numeric(auc.tmp@y.values)
  print(auc)}




