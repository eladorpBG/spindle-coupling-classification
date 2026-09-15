library(ggpubr)
library(rstatix)
library(survey)
library(reshape2)
library(caret)
library(randomForest)
library(zoo)
library(pROC)
library(ROCR)
library(xgboost)
library(fastDummies)
library(SHAPforxgboost)
library(ggplot2)
library(PRROC)
library(shapviz)
library(splitstackshape)




set.seed(123)

## Load data
adults_data_2025 <- read.csv("~/spindles/adults_data_2025.csv")
spndl_data <- adults_data_2025
spndl_data$coupling_label <- as.factor(spndl_data$coupling_label)
spndl_data <- subset(spndl_data, !grepl("O", channel))
# spndl_data$channel <- ifelse(grepl("F3", spndl_data$channel), "F3", spndl_data$channel)
# spndl_data$channel <- ifelse(grepl("C3", spndl_data$channel), "C3", spndl_data$channel)
# spndl_data$channel <- ifelse(grepl("F4", spndl_data$channel), "F4", spndl_data$channel)
# spndl_data$channel <- ifelse(grepl("C4", spndl_data$channel), "C4", spndl_data$channel)
spndl_data <- tibble::rowid_to_column(spndl_data, "ID")

# ages <- ages_data
# merged_data <- merge(spndl_data, ages, by = "patient_id")
# merged_data <- merged_data[order(match(merged_data$ID , spndl_data$ID )), ]
# 
# spndl_data <- merged_data

## Divide to train and test
spndl_data$refrPeriod <- na.aggregate(spndl_data$refrPeriod, by = spndl_data$coupling_label)
spndl_data <- subset(spndl_data, select= -numCycles)

# set the seed to make your partition reproducible
#train <- spndl_data[sample(spndl_data$ID, size = smp_size),]
# trainIndex <- createDataPartition(spndl_data$coupling_label, p = 0.66, list = FALSE, times = 1)

stratified_list <- stratified(spndl_data, c('coupling_label', 'patient_id'), size = 0.66, bothSets=TRUE)


# Split the dataset into training and testing sets
train <- as.data.frame(stratified_list[[1]])
test <- as.data.frame(stratified_list[[2]])

mean(table(spndl_data$patient_id[spndl_data$coupling_label == 1]))
prop.table(table(spndl_data$coupling_label))
prop.table(table(train$coupling_label))
prop.table(table(test$coupling_label))

## Classification model - XGBoost

train_mat <- as.matrix(train[,sapply(train, is.numeric)])
train_mat <- train_mat[, 6:ncol(train_mat)]
# train_categ_encoded <- dummy_cols(train$channel)
# train_mat <- cbind(train_mat, as.matrix(train_categ_encoded[,2:ncol(train_categ_encoded)]))

test_mat <- as.matrix(test[,sapply(test, is.numeric)])
test_mat <- test_mat[, 6:ncol(test_mat)]
# test_categ_encoded <- dummy_cols(test$channel)
# test_mat <- cbind(test_mat, as.matrix(test_categ_encoded[,2:ncol(test_categ_encoded)]))

xgb_train <- xgb.DMatrix(data = train_mat, label = as.numeric(train$coupling_label) -1)
xgb_test <- xgb.DMatrix(data = test_mat, label = as.numeric(test$coupling_label) -1)
imbalance_weight <- sum(train$coupling_label == 0) / sum(train$coupling_label == 1)


## Random Tuning

train_xgb <- function(params) {
  
  # Perform cross-validation with xgb.cv
  cv_results <- xgb.cv(
    params = params,         # Parameters for the model
    data = xgb_train,             # Training data (xgb.DMatrix)
    nrounds = params$nrounds, # Number of boosting rounds
    nfold = 5,                    # Number of cross-validation folds
    verbose = 1,                  # Verbosity (1 for detailed output)
    scale_pos_weight = imbalance_weight,  # Handle class imbalance
    early_stopping_rounds = 10,   # Stop if no improvement for 10 rounds
    maximize = TRUE,              # Maximize the evaluation metric
    metrics = "auc",              # Evaluation metric
    booster = "dart"
  )
  
  # Output cross-validation results
  print(cv_results)
  # model <- xgb.train(
  #   params = params,
  #   data = xgb_train,
  #   nrounds = params$nrounds,
  #   watchlist = list(validation1 = xgb_train),
  #   early_stopping_rounds = 10,
  #   verbose = 0,
  #   scale_pos_weight = imbalance_weight
  # )
  eval_log <- as.data.frame(cv_results$evaluation_log)
  auc <- eval_log[cv_results$best_iteration, "test_auc_mean"]
  # auc <- auc[,]
  iter <- cv_results$best_iteration
  list(as.numeric(auc), as.numeric(iter))
}

random_search <- 50
results <- data.frame()

for (i in 1:random_search) {
  print("######################################")
  print(paste("Iteration num", i))
  print("######################################")
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = runif(1, 0.01, 0.3),
    max_depth = sample(3:10, 1),
    min_child_weight = sample(1:5, 1),
    subsample = runif(1, 0.6, 1),
    colsample_bytree = runif(1, 0.6, 1),
    nrounds = sample(c(500, 1000, 2500), 1)
  )
  score <- train_xgb(params)
  results <- rbind(results, cbind(as.data.frame(params), as.numeric(score[1]), as.numeric(score[2])))
}

write.csv(as.data.frame(results), "tuning_results_adults.csv")

tuning_results_adults <- read_csv("~/spindles/tuning_results_adults.csv")
results <- tuning_results_adults

max_score_idx = which.max(results$`as.numeric(score[1])`)
best_params <- results[max_score_idx,-(ncol(results)-1)]
best_params_list <- as.list(best_params[,-ncol(best_params)])
# 
# ## Train model 
# xgb_model <- xgb.train(
#   booster = "dart",
#   params = best_params_list,
#   data = xgb_train,
#   nrounds = best_params$`as.numeric(score[2])`,
#   verbose = 1,
#   scale_pos_weight = imbalance_weight
# )
# xgb_model
# 
# # Training performance
# xgb_preds <- predict(xgb_model, train_mat, reshape = TRUE )
# xgb_preds_df <- as.data.frame(xgb_preds)
# colnames(xgb_preds_df) <- "1"
# xgb_preds_df$"0" <- 1 - xgb_preds_df$"1"
# xgb_preds_df$PredictedClass <- apply(xgb_preds_df, 1, function(y) colnames(xgb_preds_df)[which.max(y)])
# confusion_matrix_train <- confusionMatrix(as.factor(xgb_preds_df$PredictedClass), train$coupling_label)
# print(confusion_matrix_train)
# 
# xgb_prediction <- prediction(xgb_preds, train$coupling_label);
# xgb_auc <- performance(xgb_prediction, measure = "auc")@y.values[[1]]
# print(xgb_auc)
# roc_curve <- performance(xgb_prediction, measure="tpr", x.measure="fpr")
# pdf()
# plot(roc_curve)
# dev.off()
# 
# # Test performance
# xgb_preds <- predict(xgb_model, test_mat, reshape = TRUE )
# xgb_preds_df <- as.data.frame(xgb_preds)
# colnames(xgb_preds_df) <- "1"
# xgb_preds_df$"0" <- 1 - xgb_preds_df$"1"
# xgb_preds_df$PredictedClass <- apply(xgb_preds_df, 1, function(y) colnames(xgb_preds_df)[which.max(y)])
# confusion_matrix_test <- confusionMatrix(as.factor(xgb_preds_df$PredictedClass), test$coupling_label)
# print(confusion_matrix_test)
# 
# xgb_prediction <- prediction(xgb_preds, test$coupling_label);
# xgb_auc <- performance(xgb_prediction, measure = "auc")@y.values[[1]]
# print(xgb_auc)
# roc_curve <- performance(xgb_prediction, measure="tpr", x.measure="fpr")
# png()
# plot(roc_curve)
# dev.off()
# 
# # PR curve
# 
# # Precision-recall curve and AUC
# pr_result <- pr.curve(scores.class0 = xgb_preds, weights.class0 = (as.numeric(test$coupling_label) - 1), curve = TRUE)
# 
# # Display AUC of the PR curve
# print(paste("AUC-PR:", pr_result$auc.integral))
# 
# # Plot the precision-recall curve
# png()
# plot(pr_result, main = "Precision-Recall Curve")
# dev.off()
# 
# 
# importance_matrix = xgb.importance(colnames(xgb_train), model = xgb_model)
# importance_matrix
# pdf()
# xgb.plot.importance(importance_matrix)
# dev.off()
# 
# ## Shap analysis
# 
# 
# # Compute SHAP values
# shap_long <- shap.prep(xgb_model = xgb_model,
#                        X_train = test_mat)
# 
# png("shap_summary.png")
# shap.plot.summary(shap_long)
# dev.off()
# 
# for (feature_name in colnames(test_mat)) {
#   png(paste0(feature_name, ".png"))
#   plot <- shap.plot.dependence(data_long = shap_long, x = feature_name, color_feature ="ages", alpha = 0.5, size0 = 0.3, dilute=10)
#   print(plot)
#   dev.off()
# }
# 
# png(paste0("fano_zoom", ".png"))
# plot <- shap.plot.dependence(data_long = shap_long, x = "fano",
#                              color_feature ="ages",
#                              alpha = 0.5, size0 = 0.3, dilute=10) +
#   scale_x_continuous(limits = c(0, 0.5))
# print(plot)
# dev.off() 
# 
# ## Get missclassified
# 
# data_inds <- data_widxs
# data_inds <- subset(data_inds, !grepl("O", channel))
# data_inds <- tibble::rowid_to_column(data_inds, "ID")
# 
# xgb_preds_df$id <- test$ID
# sorted_votes <-  xgb_preds_df[order(xgb_preds_df$`0`),] # Coupled is first 
# sorted_votes_inds <- as.numeric(sorted_votes$id)
# test_by_votes <- test[order(match(test$ID , sorted_votes$id)), ]
# test_by_votes$coupled_vote <- sorted_votes$`1`
# test_by_votes$coupled_prediction <- sorted_votes$PredictedClass
# missclass <- test_by_votes[test_by_votes$coupling_label != test_by_votes$coupled_prediction,]
# filtered_df <- data_inds %>% filter(ID %in% missclass$ID)
# filtered_df <- filtered_df[order(match(filtered_df$ID , missclass$ID )), ]
# 
# 
# missclass$spindle_idx <- filtered_df$spindle_idx
# # Extract first 5 rows
# first_ten <- head(missclass, 10)
# 
# # Extract last 5 rows
# last_ten <- tail(missclass, 10)
# 
# # Combine them into a new dataframe
# top_miss <- rbind(first_ten, last_ten)
# 
# # Print the new dataframe
# write.csv(top_miss, "top_miss2.csv")
# write.csv(data_inds, "data_inds.csv")
# 
# fixed_two <- fixed2
# 
# ## Force plots
# top_miss
# miss_mat <- as.matrix(top_miss[,sapply(top_miss, is.numeric)])
# miss_mat <- miss_mat[, 8:ncol(miss_mat)-2]
# shp <- shapviz(xgb_model, X_pred = miss_mat, X = miss_mat)
