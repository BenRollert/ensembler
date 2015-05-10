#' Load 'results' object from an .Rda file; can assign output to a variable.
#' 
#' @param dataset A character string specifying the dataset associated with the desired results object.
#' @param par A character string specifying the parameter or model associated with the desired results object.
#' @return The 'results' slot from a caret train object.
#' @examples
#' results <- loadResults("BreastCancer", "nnet")
loadResults <- function(dataset, par){
load(sprintf("results_%s_%s.Rda", dataset, par))
results
}


launcher <- function(script, dataset, par){
  lapply(par, function(m) domino.run(script, dataset, m))
}

compare <- function(dataset, par, metric){
  result <- lapply(par, loadResults, dataset=dataset)
  result <- do.call(rbind, result)
  result[order(-result[[metric]]), ]
}


load_models <- function(dataset, method){
  load(sprintf("model_obj_%s_%s.Rda", dataset, method))
  model
}


modelList <- function(dataset, models){
  ens <- lapply(models, load_models, dataset=dataset)
  class(ens) <- "ensemble"
  ens
}


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


stackerReg <- function(models, train_mat, pred_mat, cost, type){
  l1m <- LiblineaR(data = train_mat[ ,models], target = training$Class, type = type, cost = cost, epsilon = 0.01,
                   svr_eps = NULL, bias = TRUE, wi = NULL, cross = 0, verbose = FALSE)
  log_ensemble <- predict(l1m, newx=pred_mat[ ,models], proba=TRUE)#combine model predictions on test data, using logit model
  return(log_ensemble$probabilities[,2])
}


topModel <- function(n, ordered_acc, pred_mat, reference){
  if(n>1){
  top_mean <- apply(pred_mat[,ordered_acc$model[1:n]], 1, mean)
  pred <- prediction(top_mean, reference)
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred, "auc")
  auc <- as.numeric(auc.tmp@y.values)
  cbind(as.integer(n), auc, paste(ordered_acc$model[1:n], collapse=", "))  
  }else{
    top_mean <-pred_mat[,ordered_acc$model[1]]
    pred <- prediction(top_mean, reference)
    perf <- performance(pred, "tpr", "fpr")
    auc.tmp <- performance(pred, "auc")
    auc <- as.numeric(auc.tmp@y.values)
    cbind(as.integer(n), auc, paste(ordered_acc$model[1]))}
}

topComb <- function(ordered_acc, pred_mat, reference){
  top_ens <- sapply(1:length(models), topModel, ordered_acc= ordered_acc, pred_mat=pred_mat, reference=reference)
  top_m <- data.table(t(top_ens))
  setnames(top_m, c("number_models", "AUC", "model_names"))
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
  pred <- prediction(pred, reference)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  auc.tmp <- performance(pred, "auc")
  auc <- as.numeric(auc.tmp@y.values)
  print(auc)}




