# Modular code chunks for stroke_biomarker_analysis.Rmd.
# These named chunks are loaded with knitr::read_chunk();
# knit the R Markdown entry point rather than source() this file directly.

## ---- machine-learning-development ----
# Verify required packages
ml_packages <- c(
  "glmnet",
  "e1071",
  "randomForest",
  "limma"
)

ml_package_status <- vapply(
  ml_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)

if (!all(ml_package_status)) {
  stop("One or more required machine learning packages are missing.")
}
# Verify expression matrix and clinical metadata
stopifnot(
  is.matrix(expr_gene),
  is.numeric(expr_gene),
  ncol(expr_gene) == 63,
  nrow(metadata) == 63,
  identical(colnames(expr_gene), metadata$Sample_ID),
  !anyNA(expr_gene),
  !anyNA(metadata$Group),
  sum(metadata$Group == "Control") == 24,
  sum(metadata$Group == "Stroke") == 39
)
# Set factor levels for classification
metadata$Group <- factor(
  metadata$Group,
  levels = c("Control", "Stroke")
)
# Machine learning settings
ml_seed <- 2026

number_outer_folds <- 5
number_outer_repeats <- 2
number_inner_folds <- 3

candidate_gene_counts <- c(10L, 20L)
# --------------------------------------------------
# Function 1: Create stratified cross-validation folds
# --------------------------------------------------
make_stratified_folds <- function(group, k, seed) {

  group <- factor(
    group,
    levels = c("Control", "Stroke")
  )

  if (any(table(group) < k)) {
    stop("Each class must contain at least k samples.")
  }

  set.seed(seed)

  fold_assignment <- integer(length(group))

  for (class_name in levels(group)) {

    class_indices <- which(group == class_name)

    fold_assignment[class_indices] <- sample(
      rep(
        seq_len(k),
        length.out = length(class_indices)
      )
    )

  }

  lapply(
    seq_len(k),
    function(fold_number) {
      which(fold_assignment == fold_number)
    }
  )

}
# --------------------------------------------------
# Function 2: Calculate ROC-AUC
# --------------------------------------------------
calculate_auc <- function(actual, predicted_probability) {

  actual <- as.character(actual)

  positive <- actual == "Stroke"

  n_positive <- sum(positive)
  n_negative <- sum(!positive)

  if (
    n_positive == 0 ||
    n_negative == 0 ||
    any(!is.finite(predicted_probability))
  ) {
    return(NA_real_)
  }

  probability_ranks <- rank(
    predicted_probability,
    ties.method = "average"
  )

  auc_value <- (
    sum(probability_ranks[positive]) -
      n_positive * (n_positive + 1) / 2
  ) / (n_positive * n_negative)

  as.numeric(auc_value)

}
# --------------------------------------------------
# Function 3: Select genes using training data only
# --------------------------------------------------
select_top_genes <- function(sample_indices) {

  training_metadata <- metadata[
    sample_indices,
    ,
    drop = FALSE
  ]

  training_expression <- expr_gene[
    ,
    sample_indices,
    drop = FALSE
  ]
  # Create design matrix using training samples only
  training_design <- model.matrix(
    ~ Group,
    data = training_metadata
  )
  # Fit Limma model using training samples only
  training_fit <- limma::lmFit(
    training_expression,
    training_design
  )

  training_fit <- limma::eBayes(
    training_fit,
    robust = TRUE
  )
  # Rank genes by differential expression p-value
  ranked_results <- limma::topTable(
    training_fit,
    coef = "GroupStroke",
    number = 20,
    adjust.method = "BH",
    sort.by = "P"
  )

  selected_genes <- rownames(ranked_results)

  if (
    length(selected_genes) < max(candidate_gene_counts)
  ) {
    stop("Insufficient genes available for feature selection.")
  }

  selected_genes

}
# --------------------------------------------------
# Function 4: Fit classifier and obtain probabilities
# --------------------------------------------------
predict_classifier_probability <- function(
  model_name,
  training_x,
  training_y,
  validation_x,
  parameters
) {

  training_y <- factor(
    as.character(training_y),
    levels = c("Control", "Stroke")
  )

  validation_x <- as.matrix(validation_x)
  training_x <- as.matrix(training_x)
  # Elastic Net Logistic Regression
  if (model_name == "Elastic Net") {

    fitted_model <- glmnet::glmnet(
      x = training_x,
      y = as.integer(training_y == "Stroke"),
      family = "binomial",
      alpha = parameters$alpha,
      lambda = parameters$lambda,
      standardize = TRUE
    )

    predicted_probability <- predict(
      fitted_model,
      newx = validation_x,
      s = parameters$lambda,
      type = "response"
    )

    return(as.numeric(predicted_probability[, 1]))

  }
  # Scale data using training data only
  training_center <- colMeans(training_x)

  training_scale <- apply(
    training_x,
    2,
    sd
  )

  training_scale[
    !is.finite(training_scale) |
      training_scale == 0
  ] <- 1

  scaled_training_x <- scale(
    training_x,
    center = training_center,
    scale = training_scale
  )

  scaled_validation_x <- scale(
    validation_x,
    center = training_center,
    scale = training_scale
  )
  # Support Vector Machine
  if (model_name == "SVM") {

    fitted_model <- e1071::svm(
      x = scaled_training_x,
      y = training_y,
      kernel = "radial",
      cost = parameters$cost,
      gamma = parameters$gamma,
      probability = TRUE,
      scale = FALSE
    )

    svm_prediction <- predict(
      fitted_model,
      newdata = scaled_validation_x,
      probability = TRUE
    )

    svm_probabilities <- attr(
      svm_prediction,
      "probabilities"
    )

    if (
      is.null(svm_probabilities) ||
      !"Stroke" %in% colnames(svm_probabilities)
    ) {
      stop("SVM probability predictions are unavailable.")
    }

    return(
      as.numeric(svm_probabilities[, "Stroke"])
    )

  }
  # Random Forest
  if (model_name == "Random Forest") {

    selected_mtry <- max(
      1,
      min(
        ncol(training_x),
        round(
          ncol(training_x) *
            parameters$mtry_fraction
        )
      )
    )

    fitted_model <- randomForest::randomForest(
      x = training_x,
      y = training_y,
      ntree = 300,
      mtry = selected_mtry,
      nodesize = parameters$nodesize,
      importance = FALSE
    )

    rf_probabilities <- predict(
      fitted_model,
      newdata = validation_x,
      type = "prob"
    )

    return(
      as.numeric(rf_probabilities[, "Stroke"])
    )

  }

  stop("Unknown machine learning model.")

}
# --------------------------------------------------
# Function 5: Tune model using inner cross-validation
# --------------------------------------------------
tune_classifier <- function(
  model_name,
  parameter_grid,
  inner_cache,
  seed_base
) {

  candidate_results <- vector(
    "list",
    nrow(parameter_grid)
  )

  for (candidate in seq_len(nrow(parameter_grid))) {

    parameters <- as.list(
      parameter_grid[candidate, , drop = FALSE]
    )

    inner_auc <- rep(
      NA_real_,
      length(inner_cache)
    )

    for (fold_index in seq_along(inner_cache)) {

      fold_data <- inner_cache[[fold_index]]

      selected_genes <- fold_data$ranked_genes[
        seq_len(parameters$n_features)
      ]

      training_x <- t(
        expr_gene[
          selected_genes,
          fold_data$training_indices,
          drop = FALSE
        ]
      )

      validation_x <- t(
        expr_gene[
          selected_genes,
          fold_data$validation_indices,
          drop = FALSE
        ]
      )

      training_y <- metadata$Group[
        fold_data$training_indices
      ]

      validation_y <- metadata$Group[
        fold_data$validation_indices
      ]
      # Fit and predict within this inner fold
      inner_auc[fold_index] <- tryCatch({
        set.seed(
          seed_base +
            candidate * 100 +
            fold_index
        )
        probability <- predict_classifier_probability(
          model_name = model_name,
          training_x = training_x,
          training_y = training_y,
          validation_x = validation_x,
          parameters = parameters
        )
        calculate_auc(
          validation_y,
          probability
        )
      }, error = function(e) {
        NA_real_
      })

    }
    # Require successful results from every inner fold
    mean_auc <- if (all(is.finite(inner_auc))) {
      mean(inner_auc)
    } else {
      NA_real_
    }

    sd_auc <- if (all(is.finite(inner_auc))) {
      sd(inner_auc)
    } else {
      NA_real_
    }

    candidate_results[[candidate]] <- data.frame(
      Model = model_name,
      parameter_grid[candidate, , drop = FALSE],
      Mean_Inner_AUC = mean_auc,
      SD_Inner_AUC = sd_auc,
      Successful_Inner_Folds = sum(
        is.finite(inner_auc)
      ),
      stringsAsFactors = FALSE
    )

  }

  tuning_scores <- do.call(
    rbind,
    candidate_results
  )

  valid_candidates <- which(
    is.finite(tuning_scores$Mean_Inner_AUC)
  )

  if (length(valid_candidates) == 0) {
    stop(
      paste(
        "No valid tuning configuration for",
        model_name
      )
    )
  }
  # Select the configuration with the highest mean inner AUC
  best_index <- valid_candidates[
    which.max(
      tuning_scores$Mean_Inner_AUC[valid_candidates]
    )
  ]

  best_parameters <- as.list(
    tuning_scores[
      best_index,
      names(parameter_grid),
      drop = FALSE
    ]
  )

  list(
    best_parameters = best_parameters,
    best_inner_auc = tuning_scores$Mean_Inner_AUC[
      best_index
    ],
    tuning_scores = tuning_scores
  )

}
# --------------------------------------------------
# Define candidate parameter grids
# --------------------------------------------------
elastic_net_grid <- expand.grid(
  n_features = candidate_gene_counts,
  alpha = c(0.5, 1.0),
  lambda = c(0.001, 0.01, 0.1),
  KEEP.OUT.ATTRS = FALSE
)

svm_grid <- expand.grid(
  n_features = candidate_gene_counts,
  cost = c(0.1, 1, 10),
  gamma = c(0.01, 0.1),
  KEEP.OUT.ATTRS = FALSE
)

random_forest_grid <- expand.grid(
  n_features = candidate_gene_counts,
  mtry_fraction = c(0.2, 0.5),
  nodesize = c(1, 3),
  KEEP.OUT.ATTRS = FALSE
)
# Store parameter grids
model_grids <- list(
  "Elastic Net" = elastic_net_grid,
  "SVM" = svm_grid,
  "Random Forest" = random_forest_grid
)
# --------------------------------------------------
# Nested repeated cross-validation
# --------------------------------------------------
outer_prediction_results <- list()

tuning_results_list <- list()

selected_features_list <- list()

prediction_index <- 1
tuning_index <- 1
feature_index <- 1

number_samples <- nrow(metadata)

for (repeat_id in seq_len(number_outer_repeats)) {

  cat(
    "\nOuter cross-validation repeat:",
    repeat_id,
    "\n"
  )
  # Create outer stratified folds
  outer_folds <- make_stratified_folds(
    group = metadata$Group,
    k = number_outer_folds,
    seed = ml_seed + repeat_id
  )

  for (outer_fold in seq_len(number_outer_folds)) {

    cat(
      "\nRepeat",
      repeat_id,
      "- Outer fold",
      outer_fold,
      "\n"
    )
    # Split outer training and validation samples
    outer_validation_indices <- outer_folds[[outer_fold]]

    outer_training_indices <- setdiff(
      seq_len(number_samples),
      outer_validation_indices
    )
    # Create inner stratified folds within outer training data
    inner_folds <- make_stratified_folds(
      group = metadata$Group[
        outer_training_indices
      ],
      k = number_inner_folds,
      seed = ml_seed +
        repeat_id * 100 +
        outer_fold
    )
    # Precalculate gene rankings using each inner training fold
    inner_cache <- lapply(
      seq_along(inner_folds),
      function(inner_fold) {
        inner_validation_local <- inner_folds[[inner_fold]]
        inner_training_local <- setdiff(
          seq_along(outer_training_indices),
          inner_validation_local
        )
        inner_training_indices <- outer_training_indices[
          inner_training_local
        ]
        inner_validation_indices <- outer_training_indices[
          inner_validation_local
        ]
        list(
          training_indices = inner_training_indices,
          validation_indices = inner_validation_indices,
          ranked_genes = select_top_genes(
            inner_training_indices
          )
        )
      }
    )
    # Tune and evaluate each model
    for (model_name in names(model_grids)) {

      cat(
        "Tuning model:",
        model_name,
        "\n"
      )

      tuning_output <- tune_classifier(
        model_name = model_name,
        parameter_grid = model_grids[[model_name]],
        inner_cache = inner_cache,
        seed_base = ml_seed +
          repeat_id * 10000 +
          outer_fold * 1000
      )
# Store tuning results in one consistent schema across all models
tuning_table <- tuning_output$tuning_scores
parameter_names <- c("alpha", "lambda", "cost", "gamma", "mtry_fraction", "nodesize")
parameter_columns <- lapply(parameter_names, function(nm) {
  if (nm %in% names(tuning_table)) as.numeric(tuning_table[[nm]]) else rep(NA_real_, nrow(tuning_table))
})
names(parameter_columns) <- c("Alpha", "Lambda", "Cost", "Gamma", "Mtry_Fraction", "Nodesize")
tuning_results_list[[tuning_index]] <- data.frame(
  Model = as.character(tuning_table$Model), Repeat = repeat_id,
  Outer_Fold = outer_fold, N_Features = as.integer(tuning_table$n_features),
  parameter_columns,
  Mean_Inner_AUC = as.numeric(tuning_table$Mean_Inner_AUC),
  SD_Inner_AUC = as.numeric(tuning_table$SD_Inner_AUC),
  Successful_Inner_Folds = as.integer(tuning_table$Successful_Inner_Folds),
  stringsAsFactors = FALSE
)
tuning_index <- tuning_index + 1

      best_parameters <- tuning_output$best_parameters
      # Perform gene selection again using the complete outer training set
      outer_ranked_genes <- select_top_genes(
        outer_training_indices
      )

      final_genes <- outer_ranked_genes[
        seq_len(best_parameters$n_features)
      ]
      # Construct outer training and validation expression matrices
      outer_training_x <- t(
        expr_gene[
          final_genes,
          outer_training_indices,
          drop = FALSE
        ]
      )

      outer_validation_x <- t(
        expr_gene[
          final_genes,
          outer_validation_indices,
          drop = FALSE
        ]
      )

      outer_training_y <- metadata$Group[
        outer_training_indices
      ]

      outer_validation_y <- metadata$Group[
        outer_validation_indices
      ]
      # Fit the selected model and predict outer validation samples
      set.seed(
        ml_seed +
          repeat_id * 100000 +
          outer_fold * 1000 +
          match(model_name, names(model_grids))
      )

      outer_probabilities <- predict_classifier_probability(
        model_name = model_name,
        training_x = outer_training_x,
        training_y = outer_training_y,
        validation_x = outer_validation_x,
        parameters = best_parameters
      )

      if (
        length(outer_probabilities) !=
          length(outer_validation_indices) ||
        any(!is.finite(outer_probabilities))
      ) {

        stop(
          paste(
            "Invalid outer predictions generated by",
            model_name
          )
        )

      }
      # Save out-of-fold predictions
      outer_prediction_results[[prediction_index]] <- data.frame(
        Repeat = repeat_id,
        Outer_Fold = outer_fold,
        Sample_ID = metadata$Sample_ID[
          outer_validation_indices
        ],
        Actual = as.character(
          outer_validation_y
        ),
        Model = model_name,
        Probability_Stroke = as.numeric(
          outer_probabilities
        ),
        Predicted_Class = ifelse(
          outer_probabilities >= 0.5,
          "Stroke",
          "Control"
        ),
        stringsAsFactors = FALSE
      )

      prediction_index <- prediction_index + 1
      # Save selected genes for reproducibility
      selected_features_list[[feature_index]] <- data.frame(
        Repeat = repeat_id,
        Outer_Fold = outer_fold,
        Model = model_name,
        Feature_Rank = seq_along(final_genes),
        Gene_Symbol = final_genes,
        stringsAsFactors = FALSE
      )

      feature_index <- feature_index + 1

      cat(
        "Best inner AUC:",
        round(tuning_output$best_inner_auc, 4),
        "\n"
      )

    }

  }

}
# --------------------------------------------------
# Combine and export results
# --------------------------------------------------
ml_outer_predictions <- do.call(
  rbind,
  outer_prediction_results
)

ml_tuning_results <- do.call(
  rbind,
  tuning_results_list
)

ml_selected_features <- do.call(
  rbind,
  selected_features_list
)

rownames(ml_outer_predictions) <- NULL
rownames(ml_tuning_results) <- NULL
rownames(ml_selected_features) <- NULL
# Verify prediction structure
expected_prediction_rows <- (
  number_samples *
    number_outer_repeats *
    length(model_grids)
)

stopifnot(
  nrow(ml_outer_predictions) ==
    expected_prediction_rows,
  !anyNA(ml_outer_predictions$Probability_Stroke),
  all(
    ml_outer_predictions$Probability_Stroke >= 0 &
      ml_outer_predictions$Probability_Stroke <= 1
  ),
  length(unique(ml_outer_predictions$Model)) == 3
)
# Define output files
ml_prediction_file <- file.path(
  data_dir,
  "GSE16561_ml_outer_cv_predictions.csv"
)

ml_tuning_file <- file.path(
  data_dir,
  "GSE16561_ml_tuning_results.csv"
)

ml_features_file <- file.path(
  data_dir,
  "GSE16561_ml_selected_features.csv"
)

ml_settings_file <- file.path(
  data_dir,
  "GSE16561_ml_model_settings.csv"
)
# Export out-of-fold predictions
write.csv(
  ml_outer_predictions,
  ml_prediction_file,
  row.names = FALSE
)
# Export tuning results
write.csv(
  ml_tuning_results,
  ml_tuning_file,
  row.names = FALSE
)
# Export selected genes
write.csv(
  ml_selected_features,
  ml_features_file,
  row.names = FALSE
)
# Record model development settings
ml_model_settings <- data.frame(
  Setting = c(
    "Dataset",
    "Number of samples",
    "Control samples",
    "Stroke samples",
    "Number of genes available",
    "Outer cross-validation folds",
    "Outer cross-validation repeats",
    "Inner cross-validation folds",
    "Candidate gene counts",
    "Gene selection method",
    "Model selection metric",
    "Elastic Net",
    "Support Vector Machine",
    "Random Forest"
  ),
  Value = c(
    "GSE16561",
    nrow(metadata),
    sum(metadata$Group == "Control"),
    sum(metadata$Group == "Stroke"),
    nrow(expr_gene),
    number_outer_folds,
    number_outer_repeats,
    number_inner_folds,
    paste(candidate_gene_counts, collapse = ", "),
    "Limma training-fold differential expression ranking",
    "Mean inner-fold ROC-AUC",
    "Binomial Elastic Net Logistic Regression",
    "Radial Basis Function Kernel",
    "300 trees"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  ml_model_settings,
  ml_settings_file,
  row.names = FALSE
)
# Save complete development results for later evaluation
ml_results_rds <- file.path(
  data_dir,
  "GSE16561_ml_development_results.rds"
)

saveRDS(
  list(
    predictions = ml_outer_predictions,
    tuning_results = ml_tuning_results,
    selected_features = ml_selected_features,
    model_settings = ml_model_settings
  ),
  ml_results_rds
)
# Verify all output files
ml_output_files <- c(
  ml_prediction_file,
  ml_tuning_file,
  ml_features_file,
  ml_settings_file,
  ml_results_rds
)

verify_output_files(ml_output_files)
# Display model development summary
cat("\nMachine Learning Development Summary\n")

cat("------------------------------------\n")

cat(
  "Samples:",
  number_samples,
  "\n"
)

cat(
  "Outer validation predictions:",
  nrow(ml_outer_predictions),
  "\n"
)

cat(
  "Models developed:",
  length(model_grids),
  "\n"
)

cat(
  "Repeated outer folds:",
  number_outer_repeats * number_outer_folds,
  "\n"
)

cat("\nPrediction counts by model:\n")

print(
  table(ml_outer_predictions$Model)
)

cat("\nSelected model features:\n")

print(
  aggregate(
    Feature_Rank ~ Model,
    data = ml_selected_features,
    FUN = max
  )
)

cat("\nAll machine learning output files verified.\n")

## ---- machine-learning-evaluation ----
# --------------------------------------------------
# Load and validate outer-fold predictions
# --------------------------------------------------
ml_eval_predictions <- read.csv(
  ml_prediction_file,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Repeat",
  "Sample_ID",
  "Actual",
  "Model",
  "Probability_Stroke",
  "Predicted_Class"
)

stopifnot(
  all(required_columns %in% names(ml_eval_predictions)),
  !anyNA(ml_eval_predictions$Probability_Stroke),
  all(
    ml_eval_predictions$Actual %in% c("Control", "Stroke")
  ),
  all(
    ml_eval_predictions$Predicted_Class %in% c("Control", "Stroke")
  )
)
# Separate predictions by model and repeat
evaluation_groups <- split(
  ml_eval_predictions,
  list(
    ml_eval_predictions$Model,
    ml_eval_predictions$Repeat
  ),
  drop = TRUE
)
# --------------------------------------------------
# Calculate performance metrics
# --------------------------------------------------
calculate_performance <- function(prediction_data) {
  actual <- factor(prediction_data$Actual, levels = c("Control", "Stroke"))
  predicted <- factor(prediction_data$Predicted_Class, levels = c("Control", "Stroke"))
  cm <- table(actual, predicted)
  TN <- unname(cm["Control", "Control"])
  FP <- unname(cm["Control", "Stroke"])
  FN <- unname(cm["Stroke", "Control"])
  TP <- unname(cm["Stroke", "Stroke"])

  accuracy <- (TP + TN) / sum(cm)
  sensitivity <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  specificity <- if (TN + FP > 0) TN / (TN + FP) else NA_real_
  precision <- if (TP + FP > 0) TP / (TP + FP) else 0
  f1_score <- if (is.finite(sensitivity) && precision + sensitivity > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else 0

  data.frame(
    Model = prediction_data$Model[1], Repeat = prediction_data$Repeat[1],
    Samples = nrow(prediction_data), TP = TP, TN = TN, FP = FP, FN = FN,
    Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity,
    Precision = precision, F1_Score = f1_score,
    ROC_AUC = calculate_auc(actual, prediction_data$Probability_Stroke)
  )
}
# Calculate metrics for every model and repeat
performance_by_repeat <- do.call(
  rbind,
  lapply(
    evaluation_groups,
    calculate_performance
  )
)

rownames(performance_by_repeat) <- NULL
# --------------------------------------------------
# Summarise performance across the two repeats
# --------------------------------------------------
metric_columns <- c(
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "Precision",
  "F1_Score",
  "ROC_AUC"
)

model_groups <- split(
  performance_by_repeat,
  performance_by_repeat$Model
)

performance_summary <- do.call(
  rbind,
  lapply(model_groups, function(model_data) {
    summary_row <- data.frame(
      Model = model_data$Model[1]
    )
    for (metric in metric_columns) {
      summary_row[[paste0(metric, "_Mean")]] <-
        mean(model_data[[metric]], na.rm = TRUE)