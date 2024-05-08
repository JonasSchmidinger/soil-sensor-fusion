# Download this R-script to run it in your local environment 
# Should work in RStudio, but not tested in other environments
# This script was "cleaned" for publication i.e. added comments, removed unnecessary code and made it more readable
# I checked it and it should run but maybe I missed something, so please let me know if you encounter any issues
# If there are any bugs or issues, please contact Jonas.Schmidinger@uni-osnabrueck.de I will then share the messy undocumented script that was used to generate the results in the paper


# Check if packages are installed if not install them
library(randomForest)
library(data.table)
library(dplyr)
library(Metrics)
library(stats)
library(Cubist)
library(foreach)
library(parallel)
library(doParallel)
library(e1071)
library(xgboost)
library(penalized)
library(future)



# Determine parameters manually for local environment testing:
# cores  #Determine number of cores to be used for parallel computing
# debug_run # Set to TRUE if you want to run condensed script for debugging purposes, else run full functions as in paper
# soil_property # If debug_run = F, determine which soil property to predict (1 = SOM, 2 = pH, 3 = P, 4 = moisture, 5 = K, 6 = Mg)
# nested_folds # If debug_run = F, set number of folds for nested cross validation

# If debug_run = T, the script will run a condensed version of the code for debugging purposes
debug_run <- T

# E.g. to run for SOM prediction:
# debug_run <- F 
# cores <- 128 
# nested_folds <- 5 
# debug_run <- F
# soil_property <- 1



#############
#'* Creating a couple of functions used later on *
# *r2_general* for calculating the R^2 value, can be negative if results are worse than mean(y_observed)
# Input a vector with predictions as *preds* and a vector with actual values as *actual*
# Output is R^2 value
r2_general <- function(preds, actual) {
  return(1 - sum((preds - actual)^2) / sum((actual - mean(actual))^2))
}

# *sample_cv* for data splitting in CV
# Input dataframe to be split as *df_input* and the number of folds as *folds_input*
# Output is a vector of the row-length of the *df_input* with numbers 1 to *folds_input* in random order, based on these indices the data can be split
sample_cv <- function(df_input, folds_input) {
  desired_length <- nrow(df_input)
  values <- 1:folds_input
  vector_fill <- rep(values, times = floor(desired_length / length(values)))
  if (length(vector_fill) != desired_length) {
    rest_vector <- sample(1:folds_input, desired_length - length(vector_fill))
    vector_fill <- append(vector_fill, rest_vector)
  }
  permuted_vec <- sample(vector_fill)
  return(permuted_vec)
}

# *delete_commas* for better storage of the PSS combination names in a dataframe
# Input the dataframe *input_df* that contains Covariates_ID column about the PSS combinations
# Output is the same dataframe with corrected names for sensors
delete_commas <- function(input_df) {
  for (i in 1:length(input_df$covariates_ID)) {
    cnter <- 0
    index_charac <- nchar(input_df$covariates_ID[i]) - cnter
    while (substr(input_df$covariates_ID[i], index_charac, index_charac) == ",") {
      cnter <- cnter + 1
      index_charac <- nchar(input_df$covariates_ID[i]) - cnter
    }
    replacement_character <- substr(input_df$covariates_ID[i], 1, nchar(input_df$covariates_ID[i]) - cnter)
    input_df$covariates_ID[i] <- replacement_character
  }
  return(input_df)
}

# *PCA_selection* for selecting the number of components for a sensor if present in the dataframe
# Input the training data and test (or validation) fold as *training_input* and *testing_input*
# Input the hyperparameters i.e. number of components to be used for the specific run as *param_grid_input*
# Input sensor name, i.e. "NIR" in our case, for which the number of desired components are selected
# Output is a list with the training data and test (or validation) fold with selected number of components for a sensor (i.e. NIR in our case) 
PCA_selection <- function(training_input, testing_input, param_grid_input, sensor_name) {
  if (any(grepl(sensor_name, names(training_input)) | grepl(paste0(".*", sensor_name, ".*"), names(training_input)))) {
    my_vector <- names(training_input)
    matching_entries <- my_vector[grep(sensor_name, my_vector)]
    anti_join_entries <- matching_entries[1:param_grid_input[1, paste0(sensor_name, "_PCA")]]
    sensor_anti_join <- setdiff(matching_entries, anti_join_entries)
    training_input <- training_input %>% dplyr::select(-one_of(sensor_anti_join))
    testing_input <- testing_input %>% dplyr::select(-one_of(sensor_anti_join))
  }
  
  return(list(training_input, testing_input))
}

# *sensor_in_grid* for evaluating whether a sensor i.e., in our case "NIR", is part of a PSS combination that is evaluated. If not, the hyperparameter grid excludes search of principle components
# Input dataframe consisting of the predictors and soil properties as *input_data*
# Input dataframe consisting of hyperparemeters to be tuned in the grid as *grid_input*
# Input name of PSS, i.e. "NIR" in our case, as *sensor_name*

sensor_in_grid <- function(input_data,grid_input,sensor_name){
  if (any(grepl(sensor_name, names(input_data)) | grepl(paste0(".*", sensor_name, ".*"), names(input_data)))) {
    grid_output <- grid_input
  } else {
    grid_output <- grid_input[grid_input$NIR_PCA == 10,]
    
  }
  return(grid_output)
}
#

#'* Organize data*
# Read dataframe containing the target soil properties and all predictors
Joint_Table <- readRDS(file = "Output/Joint_Table.df")
Joint_Table


# Group all predictors into a vector that contains specific string in their name i.e. all NIR and sent (Sentinal-2) predictors
# Will be used to match predictors to the sensor
NIR_names <- names(Joint_Table)[grep("NIR", names(Joint_Table))]
Sentinel_names <- names(Joint_Table)[grep("sent", names(Joint_Table))]
Sentinel_names
# Define list in which for each sensors, the belonging predictor IDs are associated
groups <- list(
  XRF = c("Al_XRF", "Fe_XRF", "Mn_XRF", "Ca_XRF", "P_XRF", "Si_XRF", "Zn_XRF", "Pb_XRF"),
  ISE = "pH_Veris",
  `γ` = c("K40_Gamma", "U238_Gamma", "Cs137_Gamma", "TC_Gamma"),
  CSMoist = "Moisture_Veris",
  EC = "EC_rm",
  NIR = NIR_names,
  RS = c(Sentinel_names)
)

PSS_groups <- list(
  XRF = c("Al_XRF", "Fe_XRF", "Mn_XRF", "Ca_XRF", "P_XRF", "Si_XRF", "Zn_XRF", "Pb_XRF"),
  ISE = "pH_Veris",
  `γ` = c("K40_Gamma", "U238_Gamma", "Cs137_Gamma", "TC_Gamma"),
  CSMoist = "Moisture_Veris",
  EC = "EC_rm",
  NIR = NIR_names
)

# Generate all possible combinations of PSSs, combinations_chosen contains dataframe with all possible combinations of PSSs
colnames_vector <- as.vector(names(PSS_groups))

combinations_chosen <- data.table::rbindlist(
  sapply(c(1, 2, 3, 4, 5), function(i) as.data.frame(t(combn(colnames_vector, i)))),
  fill = TRUE
)

# Generate all possible combinations of PSSs with inclusion of RS
combinations_chosen_RS <- cbind(RS_column = "RS", combinations_chosen)

#############
#'* Run main body of the script *
# *combination_performance* for running the main code by preparing data and storing modeling results
# Input the dataframe containing all the combinations of PSSs as *combinations_without_RS* and those combinations also containing RS as *combinations_with_RS*
# Input a vector of soil attributes that we want to predict as *soil_attributes*
# Input the dataframe containing the target soil properties and all predictors as *input*
# Input the list of sensors with their associated names as *group_list*
# Input the number of folds for outer cross-validation as *folds_for_CV*
# Input the number of cores to be used for parallel computing as *num_cores*
# Input the hyperparameters for Cubist as *param_grid_Cubist*
# Input the hyperparameters for Random Forest as *param_grid_RF*
# Input the hyperparameters for Support Vector Regression as *param_grid_SVR*
# Input the hyperparameters for Ridge Regression as *param_grid_ridge*
# Input the hyperparameters for XGBoost as *param_grid_XGB*
# Input the number of folds for nested cross-validation as *folds_for_nCV*
# Output is a dataframe containing evaluation results for all combinations and models 

combination_performance <- function(combinations_without_RS, combinations_with_RS, soil_attributes, input, group_list, folds_for_CV, num_cores, param_grid_Cubist, param_grid_RF, param_grid_SVR, param_grid_ridge, param_grid_XGB, folds_for_nCV) {
  # Prepare cluster for parallel computing
  active_Rstudio <- search()
  package_vector <- gsub("^package:", "", grep("^package:", active_Rstudio, value = TRUE))
  num_cores <- num_cores
  cl <- makeClusterPSOCK(num_cores, revtunnel = TRUE, verbose = TRUE)
  registerDoParallel(cl)
  
  # Count how many PSS were used in run, appears as "cov_number" in finale dataframe
  covariates_used <- 5 - (apply(combinations_without_RS, 1, function(row) sum(is.na(row))))
  covariates_used_add <- covariates_used
  
  # Differentiate in binary loop if run is with or without RS
  # Allows potential parallelization of the loop but is not used in this script
  identifier_fusion <- c("without RS", "with RS")
  
  # Store combinations used
  cov_ID_result.df <- data.frame()
  
  # Store final results within this dataframe
  result_df <- data.frame()
  
  # Loop over all identified soil attributes (e.g. pH, SOM etc.), where j is the individual soil attribute 
  for (j in soil_attributes) {
    # Loop over the vector which defines whether it is run with RS or without RS
    for (fi in identifier_fusion) {
      print(fi)
      
      # Select combinations based on whether RS is included or not
      if (fi == "without RS") {
        combinations <- combinations_without_RS
      } else {
        combinations <- combinations_with_RS
        covariates_used <- covariates_used + 1
      }
      
      # Prepare the foreach loop in which the main computing is done. Loop runs in parallel depending on the number of selected cores
      # Function named *modelling_for_all_samples* does the modeling and evaluation of the models
      inbetween_result.df <- foreach(i = 1:nrow(combinations), .combine = rbind, .packages = package_vector, .export = ls(globalenv())) %dopar% {
        modelling_for_all_samples(combinations, covariates_used, group_list, j, input, fi, i, folds_for_CV, covariates_used_add, param_grid_Cubist, param_grid_RF, param_grid_SVR, param_grid_ridge, param_grid_XGB, folds_for_nCV)
      }
    
      result_df <- rbind(inbetween_result.df, result_df)
    }
    
    # In the following the results from the modelling are aggregated and labeled
    combinations_together <- rbind(combinations_without_RS, combinations_without_RS)
    empty_charac <- as.data.frame(replace(combinations_together, is.na(combinations_together), ""))
    covariates_ID <- apply(empty_charac, 1, paste, collapse = ",")
    cov_ID_result.df <- cbind(result_df, covariates_ID)
  }
  result_df <- cbind(result_df, covariates_ID)
  result_df <- delete_commas(result_df)
  stopCluster(cl)
  return(result_df)
}


# *modelling_for_all_samples* does the predictive modeling and returns the evaluation results for a soil attribute
# Inputs for *modelling_for_all_samples* are passed to function within the body of *combination_performance*

modelling_for_all_samples <- function(combinations, covariates_used, group_list, j, input, fi, i, folds_for_CV, covariates_used_add, param_grid_Cubist, param_grid_RF, param_grid_SVR, param_grid_ridge, param_grid_XGB, folds_for_nCV) {
  sink(NULL)
  worker_id <- Sys.getpid()
  
  # Identify the sensor-acronyms for the specific combinations
  identifier_cov <- as.character(combinations[i, 1:covariates_used[[i]]])
  
  fill_vec <- c()
  for (l in identifier_cov) { 
    
    # Change the sensor-acronym to the corresponding predictor names in the dataset
    fill_vec <- append(fill_vec, group_list[[l]])
  }
  identifier <- append(j, fill_vec) 
  
  # Create dataframe with selected predictions and target soil property
  selected_columns.df <- input[identifier] 
  
  # Prepare outer folds 
  set.seed(150)
  indices_split <- sample_cv(input, folds_for_CV) 
  
  # Prepare vectors in which predictions and observations of outer folds are stored
  Cubist_y_pred_vector <- c()
  y_obs_vector <- c()
  RF_y_pred_vector <- c()
  SVR_y_pred_vector <- c()
  ridge_y_pred_vector <- c()
  XGB_y_pred_vector <- c()
  Meta_Model_y_pred_vector <- c()
  
  # Start outer loop 
  for (fold in 1:folds_for_CV) {
    # Allocate outer folds i.e., define training and test folds
    train_data <- selected_columns.df[indices_split != fold, ]
    test_data <- selected_columns.df[indices_split == fold, ]
    
    # Prepare vectors to store performances of inner loop (ncv = nested cross-validation)
    Cubist_ncv_perform <- c()
    RF_ncv_perform <- c()
    SVR_ncv_perform <- c()
    ridge_ncv_perform <- c()
    XGB_ncv_perform <- c()
    
    # Build lists in which to store predictions. Later used for model stacking
    Cubist_list <- list()
    RF_list <- list()
    SVR_list <- list()
    ridge_list <- list()
    XGB_list <- list()
    
    # Check whether NIR is part of PSS combination. If so grid will be reduced excluding the component search for NIR
    param_grid_Cubist_loop<- sensor_in_grid(train_data,param_grid_Cubist,"NIR")
    
    param_grid_RF_loop<- sensor_in_grid(train_data,param_grid_RF,"NIR")
    
    # Identify the length of the grid search for the nested/inner loop for RF and cubist
    total_iterations <- nrow(param_grid_Cubist_loop)
    
    for (hyperparameter in 1:total_iterations) {

      iteration_start_time <- Sys.time()
      
      # Prepare inner folds for (Cubist & RF)
      set.seed(150)
      indices_split_ncv <- sample_cv(train_data, folds_for_nCV)
      
      # Prepare vectors in which predictions and observations for inner folds are stored (Cubist & RF)
      Cubist_y_pred_vector_ncv <- c()
      y_obs_vector_ncv <- c()
      RF_y_pred_vector_ncv <- c()
      
      # Start inner loop (for Cubist & RF)
      for (fold_ncv in 1:folds_for_nCV) {
        
        # Allocate inner folds i.e., define inner training and validation folds
        train_data_ncv <- train_data[indices_split_ncv != fold_ncv, ]
        validation_data_ncv <- train_data[indices_split_ncv == fold_ncv, ]
        
        # Selecting NIR components from grid if NIR present in PSS combination
        train_val <- PCA_selection(train_data_ncv, validation_data_ncv, param_grid_Cubist_loop[hyperparameter, ], "NIR")
        train_data_ncv <- train_val[[1]]
        validation_data_ncv <- train_val[[2]]
        
        # Fit Cubist with hyperparameter combinations from grid
        cubist_model_cnv <- cubist(base::subset(train_data_ncv, select = -eval(parse(text = j))), train_data_ncv[[j]],
                                   control = cubistControl(rules = param_grid_Cubist_loop$rules[hyperparameter],
                                                           sample = param_grid_Cubist_loop$sample[hyperparameter]),
                                   committees = 80
        )

        # Fit RF with hyperparameter combinations from grid
        RF_model_cnv <- randomForest(base::subset(train_data_ncv, select = -eval(parse(text = j))), train_data_ncv[[j]],
                                     nodesize = param_grid_RF_loop$nodesize[hyperparameter],
                                     ntree = 500,
                                     samplesize = nrow(train_data_ncv)*param_grid_RF_loop$samplesize[hyperparameter],
                                     mtry = ceiling((ncol(train_data_ncv) - 1) / (param_grid_RF_loop$mtry[hyperparameter]))
        )
        
        # Use RF and Cubist to predict validation fold
        Cubist_y_pred_ncv <- predict(cubist_model_cnv,
                                     newdata = base::subset(validation_data_ncv, select = -eval(parse(text = j))),
                                     neighbors = param_grid_Cubist_loop$neigh[hyperparameter]
        )
        RF_y_pred_ncv <- predict(RF_model_cnv, newdata = base::subset(validation_data_ncv, select = -eval(parse(text = j))))
        
        # Append the predictions and the validation-observations to the pre-defined vector
        # Will be needed to calculate the RMSE by comparing predictions with observations from validation-folds
        Cubist_y_pred_vector_ncv <- append(Cubist_y_pred_vector_ncv, Cubist_y_pred_ncv)
        RF_y_pred_vector_ncv <- append(RF_y_pred_vector_ncv, RF_y_pred_ncv)
        y_obs_vector_ncv <- append(y_obs_vector_ncv, validation_data_ncv[[j]])
      }
      
      # Store predictions in list
      # Will be needed for the model stacking later on, as it uses the predictions of the inner loop as feature
      Cubist_list[[hyperparameter]] <- Cubist_y_pred_vector_ncv
      RF_list[[hyperparameter]] <- RF_y_pred_vector_ncv
      
      # Calculate RMSE associated to hyperparameter combination and append to pre-defined vector
      # Will be needed for hyperparameter selection later on
      Cubist_rmse_value_ncv <- Metrics::rmse(y_obs_vector_ncv, Cubist_y_pred_vector_ncv)
      Cubist_ncv_perform <- append(Cubist_ncv_perform, Cubist_rmse_value_ncv)      
      RF_rmse_value_ncv <- Metrics::rmse(y_obs_vector_ncv, RF_y_pred_vector_ncv)
      RF_ncv_perform <- append(RF_ncv_perform, RF_rmse_value_ncv)
      
      iteration_end_time <- Sys.time()
      iteration_time <- iteration_end_time - iteration_start_time
      if (hyperparameter %% 100 == 0) {
        #cat("Worker:", worker_id, "Iteration", hyperparameter, "of", total_iterations, "finished at", format(Sys.time(), "%H:%M:%S"), "ETA:", format(iteration_time * (total_iterations - hyperparameter), units = "auto"), "\n")
      }
    }
    
    
    ##########################################################################################################################################
    
    #cat("Worker:", worker_id, "Trainig Ridge started at", format(Sys.time(), "%H:%M:%S"), "\n")
    
    # The same pipeline as for RF and Cubist will now be repeated for the other models (i.e. Ridge, SVR, XGBoost), we refrain from describing each step 
    param_grid_ridge_loop<- sensor_in_grid(train_data,param_grid_ridge,"NIR")
    
    total_iterations <- nrow(param_grid_ridge_loop)
    for (hyperparameter in 1:total_iterations) {
      iteration_start_time <- Sys.time()
      set.seed(150)
      indices_split_ncv <- sample_cv(train_data, folds_for_nCV)
      y_obs_vector_ncv <- c()
      ridge_y_pred_vector_ncv <- c()
      
      
      for (fold_ncv in 1:folds_for_nCV) {
        train_data_ncv <- train_data[indices_split_ncv != fold_ncv, ]
        validation_data_ncv <- train_data[indices_split_ncv == fold_ncv, ]
        
        train_val <- PCA_selection(train_data_ncv, validation_data_ncv, param_grid_ridge_loop[hyperparameter, ], "NIR")
        train_data_ncv <- train_val[[1]]
        validation_data_ncv <- train_val[[2]]
        
        ridge_model_cnv <- penalized(train_data_ncv[[j]], as.matrix(base::subset(train_data_ncv, select = -eval(parse(text = j)))),
                                     lambda2 = param_grid_ridge_loop$lambda[hyperparameter]
        )

        ridge_y_pred_ncv_matrix <- predict(ridge_model_cnv, penalized = as.matrix(base::subset(validation_data_ncv, select = -eval(parse(text = j)))))
        ridge_y_pred_ncv <- ridge_y_pred_ncv_matrix[, 1] # We only want to extract the predictions from the matrix
        
        ridge_y_pred_vector_ncv <- append(ridge_y_pred_vector_ncv, ridge_y_pred_ncv)
        y_obs_vector_ncv <- append(y_obs_vector_ncv, validation_data_ncv[[j]])
      }
      
      ridge_list[[hyperparameter]] <- ridge_y_pred_vector_ncv
      ridge_rmse_value_ncv <- Metrics::rmse(y_obs_vector_ncv, ridge_y_pred_vector_ncv)
      ridge_ncv_perform <- append(ridge_ncv_perform, ridge_rmse_value_ncv)
      
      iteration_end_time <- Sys.time()
      iteration_time <- iteration_end_time - iteration_start_time
      if (hyperparameter %% 100 == 0) {
        #cat("Worker:", worker_id, "Iteration", hyperparameter, "of", total_iterations, "finished at", format(Sys.time(), "%H:%M:%S"), "ETA:", format(iteration_time * (total_iterations - hyperparameter), units = "auto"), "\n")
      }
    }
    
    
    #cat("Worker:", worker_id, "Trainig SVR started at", format(Sys.time(), "%H:%M:%S"), "\n")
    param_grid_SVR_loop<- sensor_in_grid(train_data,param_grid_SVR,"NIR")
    
    total_iterations <- nrow(param_grid_SVR_loop)
    
    for (hyperparameter in 1:total_iterations) {
      iteration_start_time <- Sys.time()
      
      set.seed(150)
      indices_split_ncv <- sample_cv(train_data, folds_for_nCV)
      y_obs_vector_ncv <- c()
      SVR_y_pred_vector_ncv <- c()
      
      for (fold_ncv in 1:folds_for_nCV) {
        train_data_ncv <- train_data[indices_split_ncv != fold_ncv, ]
        validation_data_ncv <- train_data[indices_split_ncv == fold_ncv, ]
        
        train_val <- PCA_selection(train_data_ncv, validation_data_ncv, param_grid_SVR_loop[hyperparameter, ], "NIR")
        train_data_ncv <- train_val[[1]]
        validation_data_ncv <- train_val[[2]]
        
        SVR_model_cnv <- svm(
          y = train_data_ncv[[j]], x = base::subset(train_data_ncv, select = -eval(parse(text = j))),
          kernel = param_grid_SVR_loop$kernel[hyperparameter],
          cost = param_grid_SVR_loop$cost[hyperparameter],
          gamma = param_grid_SVR_loop$gamma[hyperparameter]
        )
        

        SVR_y_pred_ncv <- predict(SVR_model_cnv, newdata = base::subset(validation_data_ncv, select = -eval(parse(text = j))))
        
        SVR_y_pred_vector_ncv <- append(SVR_y_pred_vector_ncv, SVR_y_pred_ncv)
        y_obs_vector_ncv <- append(y_obs_vector_ncv, validation_data_ncv[[j]])
      }
      
      SVR_list[[hyperparameter]] <- SVR_y_pred_vector_ncv
      SVR_rmse_value_ncv <- Metrics::rmse(y_obs_vector_ncv, SVR_y_pred_vector_ncv)
      SVR_ncv_perform <- append(SVR_ncv_perform, SVR_rmse_value_ncv)
      
      iteration_end_time <- Sys.time()
      iteration_time <- iteration_end_time - iteration_start_time
      if (hyperparameter %% 100 == 0) {
        #cat("Worker:", worker_id, "Iteration", hyperparameter, "of", total_iterations, "finished at", format(Sys.time(), "%H:%M:%S"), "ETA:", format(iteration_time * (total_iterations - hyperparameter), units = "auto"), "\n")
      }
    }
    #cat("Worker:", worker_id, "Trainig XGBoost started at", format(Sys.time(), "%H:%M:%S"), "\n")
    param_grid_XGB_loop<- sensor_in_grid(train_data,param_grid_XGB,"NIR")
    
    total_iterations <- nrow(param_grid_XGB_loop)
    for (hyperparameter in 1:total_iterations) {
      iteration_start_time <- Sys.time()
      set.seed(150)
      indices_split_ncv <- sample_cv(train_data, folds_for_nCV)
      y_obs_vector_ncv <- c()
      XGB_y_pred_vector_ncv <- c()
      
      for (fold_ncv in 1:folds_for_nCV) {
        train_data_ncv <- train_data[indices_split_ncv != fold_ncv, ]
        validation_data_ncv <- train_data[indices_split_ncv == fold_ncv, ]
        
        train_val <- PCA_selection(train_data_ncv, validation_data_ncv, param_grid_XGB_loop[hyperparameter, ], "NIR")
        train_data_ncv <- train_val[[1]]
        validation_data_ncv <- train_val[[2]]
        
        XGB_model_cnv <- xgboost(
          data = as.matrix(base::subset(train_data_ncv, select = -eval(parse(text = j)))),
          label = train_data_ncv[, j],
          nrounds = 1000,
          min_child_weight = 3,
          max_depth = param_grid_XGB_loop$max_depth[hyperparameter],
          eta = param_grid_XGB_loop$eta[hyperparameter],
          gamma = param_grid_XGB_loop$gamma[hyperparameter],
          subsample = param_grid_XGB_loop$subsample[hyperparameter],
          colsample_bytree = param_grid_XGB_loop$colsample_bytree[hyperparameter],
          objective = "reg:squarederror",
          nthread = 1,
          verbose = F
        )

        XGB_y_pred_ncv <- predict(XGB_model_cnv, as.matrix(base::subset(validation_data_ncv, select = -eval(parse(text = j)))))
        
        XGB_y_pred_vector_ncv <- append(XGB_y_pred_vector_ncv, XGB_y_pred_ncv)
        y_obs_vector_ncv <- append(y_obs_vector_ncv, validation_data_ncv[[j]])
      }
      
      XGB_list[[hyperparameter]] <- XGB_y_pred_vector_ncv
      XGB_rmse_value_ncv <- Metrics::rmse(y_obs_vector_ncv, XGB_y_pred_vector_ncv)
      XGB_ncv_perform <- append(XGB_ncv_perform, XGB_rmse_value_ncv)
      
      iteration_end_time <- Sys.time()
      iteration_time <- iteration_end_time - iteration_start_time
      if (hyperparameter %% 100 == 0) {
        #cat("Worker:", worker_id, "Iteration", hyperparameter, "of", total_iterations, "finished at", format(Sys.time(), "%H:%M:%S"), "ETA:", format(iteration_time * (total_iterations - hyperparameter), units = "auto"), "\n")
      }
    }
    
    
    # Determine which hyperparemter combination had lowest RMSE, store their position i.e. which position in the grid dataframe this lowest performance has been achieved 
    # Will be used to identify the hyperparameters 
    lowest_Cubist_rmse <- which.min(Cubist_ncv_perform)
    lowest_RF_rmse <- which.min(RF_ncv_perform)
    lowest_SVR_rmse <- which.min(SVR_ncv_perform)
    lowest_ridge_rmse <- which.min(ridge_ncv_perform)
    lowest_XGB_rmse <- which.min(XGB_ncv_perform)
    
    # Store which hyperparameter combination achieved the best results
    # Will be used to fit train models with thiese hyperparameters in outer loop
    Cubist_selected_hyperparameter <- param_grid_Cubist_loop[lowest_Cubist_rmse, ]
    RF_selected_hyperparameter <- param_grid_RF_loop[lowest_RF_rmse, ]
    SVR_selected_hyperparameter <- param_grid_SVR_loop[lowest_SVR_rmse, ]
    ridge_selected_hyperparameter <- param_grid_ridge_loop[lowest_ridge_rmse, ]
    XGB_selected_hyperparameter <- param_grid_XGB_loop[lowest_XGB_rmse, ]
    
    # Store predictions associated to best hyperparameter combination 
    # Will be feed to meta-model/model stacking as training data
    Meta_Model_training_data <- data.frame(
      "Predictions" = y_obs_vector_ncv,
      "Cubist_variable" = Cubist_list[[lowest_Cubist_rmse]],
      "RF_variable" = RF_list[[lowest_RF_rmse]],
      "SVR_variable" = SVR_list[[lowest_SVR_rmse]],
      "ridge_variable" = ridge_list[[lowest_ridge_rmse]],
      "XGB_variable" = XGB_list[[lowest_XGB_rmse]]
    )
    
    # Store observations of test data from outer loop in pre defined vector
    # This will be used to calculate RMSE and R^2 
    y_obs_vector <- append(y_obs_vector, test_data[[j]])
    
    # Fit models with best hyperparameters from inner fold
    # Fit Cubist
    Cubist_data <- PCA_selection(train_data, test_data, Cubist_selected_hyperparameter, "NIR")
    cubist_model <- cubist(base::subset(Cubist_data[[1]], select = -eval(parse(text = j))), Cubist_data[[1]][[j]],
                           control = cubistControl(rules = Cubist_selected_hyperparameter$rules[1],
                                                   sample = Cubist_selected_hyperparameter$sample[1]),
                           committees = 80
    )
    
    # Fit RF
    RF_data <- PCA_selection(train_data, test_data, RF_selected_hyperparameter, "NIR")
    RF_model <- randomForest(base::subset(RF_data[[1]], select = -eval(parse(text = j))), RF_data[[1]][[j]],
                             nodesize = RF_selected_hyperparameter$nodesize[1],
                             ntree = 500,
                             samplesize = nrow(RF_data)*RF_selected_hyperparameter$samplesize[1],
                             mtry = ceiling((ncol(RF_data[[1]]) - 1) / (RF_selected_hyperparameter$mtry[1]))
    )
    
    # Fit SVR
    SVR_data <- PCA_selection(train_data, test_data, SVR_selected_hyperparameter, "NIR")
    svr_model <- svm(
      x = base::subset(SVR_data[[1]], select = -eval(parse(text = j))), y = SVR_data[[1]][[j]],
      kernel = SVR_selected_hyperparameter$kernel[1],
      cost = SVR_selected_hyperparameter$cost[1],
      gamma = SVR_selected_hyperparameter$gamma[1]
    )
    
    # Fit RR
    ridge_data <- PCA_selection(train_data, test_data, ridge_selected_hyperparameter, "NIR")
    ridge_model <- penalized(ridge_data[[1]][[j]], as.matrix(base::subset(ridge_data[[1]], select = -eval(parse(text = j)))),
                             lambda2 = ridge_selected_hyperparameter$lambda[1]
    )
    
    # Fit XGBoost
    XGB_data <- PCA_selection(train_data, test_data, XGB_selected_hyperparameter, "NIR")
    XGB_model <- xgboost(
      data = as.matrix(base::subset(XGB_data[[1]], select = -eval(parse(text = j)))),
      label = XGB_data[[1]][, j],
      nrounds = 1000,
      min_child_weight = 3,
      max_depth = XGB_selected_hyperparameter$max_depth[1],
      eta = XGB_selected_hyperparameter$eta[1],
      gamma = XGB_selected_hyperparameter$gamma[1],
      subsample = XGB_selected_hyperparameter$subsample[1],
      colsample_bytree = XGB_selected_hyperparameter$colsample_bytree[1],
      objective = "reg:squarederror",
      nthread = 1,
      verbose = F
    )
    
    # Make predictions on the test set with each base-model
    Cubist_y_pred <- predict(cubist_model,
                             newdata = base::subset(Cubist_data[[2]], select = -eval(parse(text = j))),
                             Cubist_selected_hyperparameter$neigh[1]
    )
    RF_y_pred <- predict(RF_model, newdata = base::subset(RF_data[[2]], select = -eval(parse(text = j))))
    SVR_y_pred <- predict(svr_model, newdata = base::subset(SVR_data[[2]], select = -eval(parse(text = j))))
    ridge_y_pred_matrix <- predict(ridge_model, penalized = as.matrix(base::subset(ridge_data[[2]], select = -eval(parse(text = j)))))
    ridge_y_pred <- ridge_y_pred_matrix[, 1] # For RR extra step to extract predictions
    XGB_y_pred <- predict(XGB_model, newdata = as.matrix(base::subset(XGB_data[[2]], select = -eval(parse(text = j)))))
    
    # Store predictions of outer loop in pre-defined vectors 
    # Will be used for model evaluation comparison to test set observations and for the meta model/model stacking, as predictors
    Cubist_y_pred_vector <- append(Cubist_y_pred_vector, Cubist_y_pred)
    RF_y_pred_vector <- append(RF_y_pred_vector, RF_y_pred)
    SVR_y_pred_vector <- append(SVR_y_pred_vector, SVR_y_pred)
    ridge_y_pred_vector <- append(ridge_y_pred_vector, ridge_y_pred)
    XGB_y_pred_vector <- append(XGB_y_pred_vector, XGB_y_pred)
    
    # Fit Meta Model (i.e. MLR) for model stacking, predictors are the inner loop predictions of the base-models
    Meta_Model <- lm(Predictions ~ Cubist_variable + RF_variable + SVR_variable + ridge_variable + XGB_variable, data = Meta_Model_training_data)
    
    # Make predictions with the meta model/model stacking in which the fitted meta model is used to calibrate the predictions from the outer loop
    Meta_Model_y_pred <- predict(Meta_Model, newdata = data.frame(
      "Cubist_variable" = Cubist_y_pred,
      "RF_variable" = RF_y_pred,
      "SVR_variable" = SVR_y_pred,
      "ridge_variable" = ridge_y_pred,
      "XGB_variable" = XGB_y_pred
    ))
    
    #Store predictons for model evaluation
    Meta_Model_y_pred_vector <- append(Meta_Model_y_pred_vector, Meta_Model_y_pred)
  }
  
  # Calculate RMSE and R^2 model performances for all outer loop base-models and the model stacking 
  Cubist_R2_value <- r2_general(preds = Cubist_y_pred_vector, actual = y_obs_vector)
  rmse_value_Cubist <- Metrics::rmse(y_obs_vector, Cubist_y_pred_vector)
  R2_value_RF <- r2_general(preds = RF_y_pred_vector, actual = y_obs_vector)
  rmse_value_RF <- Metrics::rmse(y_obs_vector, RF_y_pred_vector)
  R2_value_SVR <- r2_general(preds = SVR_y_pred_vector, actual = y_obs_vector)
  rmse_value_SVR <- Metrics::rmse(y_obs_vector, SVR_y_pred_vector)
  R2_value_ridge <- r2_general(preds = ridge_y_pred_vector, actual = y_obs_vector)
  rmse_value_ridge <- Metrics::rmse(y_obs_vector, ridge_y_pred_vector)
  R2_value_Meta <- r2_general(preds = Meta_Model_y_pred_vector, actual = y_obs_vector)
  rmse_value_Meta <- Metrics::rmse(y_obs_vector, Meta_Model_y_pred_vector)
  R2_value_XGB <- r2_general(preds = XGB_y_pred_vector, actual = y_obs_vector)
  rmse_value_XGB <- Metrics::rmse(y_obs_vector, XGB_y_pred_vector)
  
  # Store evaluation values (RMSE and R^2) in dataframe, as well as which soil property (*Soil_attribute = j*), which PSSs have been used (*cov_number = covariates_used_add[[i]]*) and if RS has been used (*fusion = fi*)
  df.temp <- data.frame(
    R2_MM = R2_value_Meta, MM_rmse = rmse_value_Meta,
    R2_Cubist = Cubist_R2_value, rmse_Cubist = rmse_value_Cubist,
    R2_RF = R2_value_RF, RF_rmse = rmse_value_RF,
    R2_ridge = R2_value_ridge, ridge_rmse = rmse_value_ridge,
    R2_SVR = R2_value_SVR, SVR_rmse = rmse_value_SVR,
    R2_XGB = R2_value_XGB, XGB_rmse = rmse_value_XGB,
    Soil_attribute = j, cov_number = covariates_used_add[[i]], fusion = fi
  )
  
  # Return results for aggregation in *combination_performance*
  return(df.temp)
}

# Set hyperparameter grids for all base-models

hyperparameter_grid_RF <- expand.grid(
  NIR_PCA = c(5,10,15,20),
  nodesize = c(3, 6, 9, 12),
  mtry = c(1.25, 1.5, 2, 3, 4, 5),
  samplesize = c(0.6, 0.8, 1)
)

hyperparameter_grid_Cubist <- expand.grid(
  NIR_PCA = c(5,10,15,20),
  rules = c(4, 6,12,24),
  neigh = c(0, 1, 2, 4, 6, 8),
  sample = c(60,80,0)
)

hyperparameter_grid_SVR <- expand.grid(
  NIR_PCA = c(5,10,15,20),
  kernel = c("linear", "radial"),
  cost = c(0.1,1,10,100),
  gamma = c(0.01,0.1,1,10)
)

hyperparameter_grid_ridge <- expand.grid(
  NIR_PCA = c(5,10,15,20),
  lambda = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
)


hyperparameter_grid_xgb <- expand.grid(
  NIR_PCA = c(5,10,15,20),
  max_depth = c(4, 6, 8, 10),
  eta = c(0.005, 0.01, 0.05,0.1),
  gamma = c(0.5,1,2),
  colsample_bytree= c(0.5, 0.75, 1.0),
  subsample = c(0.6, 0.8, 1)
)


# Only for debugging purpose i.e. if debug_run = TRUE
# Set hyperparameter grids for all base-models with one parameter for fast run
if (debug_run) {
  hyperparameter_grid_RF <- expand.grid(
    NIR_PCA = c(10),
    nodesize = c(3),
    mtry = c(1.5),
    samplesize = c(0.8)
  )
  
  hyperparameter_grid_Cubist <- expand.grid(
    NIR_PCA = c(10),
    rules = c(4),
    neigh = c(0),
    sample = c(80)
    
  )
  
  hyperparameter_grid_SVR <- expand.grid(
    NIR_PCA = c(10),
    kernel = c("linear"),
    cost = c(1),
    gamma = c(0.1)
  )
  
  hyperparameter_grid_ridge <- expand.grid(
    NIR_PCA = c(10),
    lambda = c(0.001)
  )
  
  hyperparameter_grid_xgb <- expand.grid(
    NIR_PCA = c(10),
    max_depth = c(9),
    eta = c(0.3),
    gamma = c(0.1),
    lambda = c(0.9),
    subsample = c(0.9),
    colsample_bytree= c(0.9)
  )
}









# Run script i.e. main functions
# Here done separately for each soil property to run them in parallel on different nodes
# We used 6*128 cores to run script i.e. one node (128 cores) was used for modeling of one soil property 
# Each run had an estimated running time of ~ 34 hours with a node

if (debug_run == FALSE) {
  if (soil_property == 1) {
    # Property 1, SOM coded as Humus
    cubist_result_combination_SOM <- combination_performance(combinations_chosen, combinations_chosen_RS, c("Humus"), Joint_Table, groups, 5, cores, hyperparameter_grid_Cubist, hyperparameter_grid_RF, hyperparameter_grid_SVR, hyperparameter_grid_ridge, hyperparameter_grid_xgb, nested_folds)
    saveRDS(cubist_result_combination_SOM, file = "Output/cubist_result_combination_SOM")

  }
  
  if (soil_property == 2) {
    # Property 2, pH coded as pH_lab
    cubist_result_combination_pH_lab <- combination_performance(combinations_chosen, combinations_chosen_RS, c("pH_lab"), Joint_Table, groups, 5, cores, hyperparameter_grid_Cubist, hyperparameter_grid_RF, hyperparameter_grid_SVR, hyperparameter_grid_ridge, hyperparameter_grid_xgb, nested_folds)
    cubist_result_combination_pH_lab
    saveRDS(cubist_result_combination_pH_lab, file = "Output/cubist_result_combination_pH_lab")

  }
  
  if (soil_property == 3) {
    # Property 3, P coded as P_lab
    cubist_result_combination_P_lab <- combination_performance(combinations_chosen, combinations_chosen_RS, c("P_lab"), Joint_Table, groups, 5, cores, hyperparameter_grid_Cubist, hyperparameter_grid_RF, hyperparameter_grid_SVR, hyperparameter_grid_ridge, hyperparameter_grid_xgb, nested_folds)
    cubist_result_combination_P_lab
    saveRDS(cubist_result_combination_P_lab, file = "Output/cubist_result_combination_P_lab")

  }
  
  if (soil_property == 4) {
    # Property 4, moisture coded as Moisture
    cubist_result_combination_Moisture <- combination_performance(combinations_chosen, combinations_chosen_RS, c("Moisture"), Joint_Table, groups, 5, cores, hyperparameter_grid_Cubist, hyperparameter_grid_RF, hyperparameter_grid_SVR, hyperparameter_grid_ridge, hyperparameter_grid_xgb, nested_folds)
    saveRDS(cubist_result_combination_Moisture, file = "Output/cubist_result_combination_Moisture")

  }
  
  if (soil_property == 5) {
    # Property 5, K coded as K_lab
    cubist_result_combination_K_lab <- combination_performance(combinations_chosen, combinations_chosen_RS, c("K_lab"), Joint_Table, groups, 5, cores, hyperparameter_grid_Cubist, hyperparameter_grid_RF, hyperparameter_grid_SVR, hyperparameter_grid_ridge, hyperparameter_grid_xgb, nested_folds)
    cubist_result_combination_K_lab
    saveRDS(cubist_result_combination_K_lab, file = "Output/cubist_result_combination_K_lab")

  }
  
  if (soil_property == 6) {
    # Property 6, Mg coded as Mg_lab
    cubist_result_combination_Mg_lab <- combination_performance(combinations_chosen, combinations_chosen_RS, c("Mg_lab"), Joint_Table, groups, 5, cores, hyperparameter_grid_Cubist, hyperparameter_grid_RF, hyperparameter_grid_SVR, hyperparameter_grid_ridge, hyperparameter_grid_xgb, nested_folds)
    cubist_result_combination_Mg_lab
    saveRDS(cubist_result_combination_Mg_lab, file = "Output/cubist_result_combination_Mg_lab")

  }
}

# Run the script for debugging purposes with debug_run set to TRUE
# Here only condensed grid and K=2, L=2 are used for CV and nested CV for the prediction of SOC, denoted as "Humus"
# It has an approximate running time of few minutes with one core
# Results are nonsensical but it allows to run whole script quickly
if (debug_run) {
  test <- combination_performance(combinations_chosen, combinations_chosen_RS, c("Humus"), Joint_Table, groups, 2, 2, hyperparameter_grid_Cubist, hyperparameter_grid_RF, hyperparameter_grid_SVR, hyperparameter_grid_ridge, hyperparameter_grid_xgb, 2)
}



