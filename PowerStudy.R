# Enhanced Power Study Implementation - Supporting Multiple Datasets
# Modified to use VEER1 dataset from cancerdata package with sample removal

library(ranger)
library(vita)
library(golubEsets)
library(cancerdata)  # Added for VEER1 dataset
library(parallel)
library(doParallel)
library(foreach)
library(Biobase)     # Added for ExpressionSet handling
library(ggplot2)     # For plotting
library(doRNG)

# ────────────────────────────────────────────────────────────────────────────
# 1. Enhanced function to load datasets (leukemia + VEER1)
# ────────────────────────────────────────────────────────────────────────────
load_real_dataset <- function(dataset_name, samples_to_remove = NULL) {
  if (dataset_name == "leukemia") {
    data(Golub_Merge)
    
    # Extract expression data and phenotype
    expr_data <- exprs(Golub_Merge)
    expr_data <- t(expr_data)  # transpose to samples x genes
    
    # Get outcome variable (ALL vs AML)
    outcome_data <- pData(Golub_Merge)$ALL.AML
    
    # Remove any samples with missing outcome
    valid_samples <- !is.na(outcome_data)
    expr_data     <- expr_data[valid_samples, ]
    outcome_data  <- outcome_data[valid_samples]
    
    # Convert to factor
    outcome_data <- as.factor(outcome_data)
    
    cat("Leukemia dataset loaded:\n")
    cat("  Samples:", nrow(expr_data), "\n")
    cat("  Genes:",   ncol(expr_data), "\n")
    cat("  Outcome classes (pre‐simulation):\n")
    print(table(outcome_data))
    
    return(list(X = expr_data, y = outcome_data, n = nrow(expr_data), p = ncol(expr_data)))
    
  } else if (dataset_name == "veer1") {
    cat("Loading VEER1 dataset from cancerdata package...\n")
    
    # Load the VEER1 dataset
    data(VEER1)
    
    # Extract expression data and phenotype information
    expr_data <- exprs(VEER1)  # genes x samples
    expr_data <- t(expr_data)  # transpose to samples x genes
    
    # Get phenotype data
    pheno_data <- pData(VEER1)
    
    cat("Available phenotype columns:\n")
    print(colnames(pheno_data))
    cat("\nFirst few rows of phenotype data:\n")
    print(head(pheno_data))
    
    # Extract survival outcome
    # Common column names for survival in breast cancer datasets:
    # Look for columns containing survival, time, event, status, etc.
    survival_cols <- grep("surv|time|event|status|outcome|class", 
                          colnames(pheno_data), ignore.case = TRUE)
    
    if (length(survival_cols) > 0) {
      cat("\nPotential survival columns found:\n")
      for (i in survival_cols) {
        col_name <- colnames(pheno_data)[i]
        cat("  ", col_name, ":", class(pheno_data[, i]), "\n")
        if (is.factor(pheno_data[, i]) || is.character(pheno_data[, i])) {
          print(table(pheno_data[, i], useNA = "always"))
        } else {
          print(summary(pheno_data[, i]))
        }
        cat("\n")
      }
      
      # Use the first survival-related column
      outcome_data <- pheno_data[, survival_cols[1]]
    } else {
      # If no clear survival column, create a binary outcome based on sample names or use a placeholder
      cat("No clear survival column found. Checking sample names for survival information...\n")
      sample_names <- rownames(pheno_data)
      cat("Sample names pattern:\n")
      print(head(sample_names, 10))
      
      # Look for survival patterns in sample names (>5yr, <5yr, etc.)
      if (any(grepl("5.*yr|survival", sample_names, ignore.case = TRUE))) {
        outcome_data <- ifelse(grepl(">5.*yr|long.*survival", sample_names, ignore.case = TRUE), 1,
                               ifelse(grepl("<5.*yr|short.*survival", sample_names, ignore.case = TRUE), 0, NA))
        outcome_data <- outcome_data[!is.na(outcome_data)]
      } else {
        # Create a random binary outcome for demonstration
        cat("Creating random binary outcome for demonstration purposes...\n")
        set.seed(123)
        outcome_data <- rbinom(nrow(expr_data), 1, 0.5)
      }
    }
    
    # Remove specified samples if provided
    if (!is.null(samples_to_remove)) {
      cat("Removing samples:", samples_to_remove, "\n")
      
      # Ensure sample indices are valid
      valid_indices <- samples_to_remove[samples_to_remove <= nrow(expr_data)]
      if (length(valid_indices) != length(samples_to_remove)) {
        cat("Warning: Some sample indices were out of range and ignored\n")
      }
      
      if (length(valid_indices) > 0) {
        # Remove samples
        expr_data <- expr_data[-valid_indices, ]
        outcome_data <- outcome_data[-valid_indices]
        
        cat("Removed", length(valid_indices), "samples\n")
      }
    }
    
    # Convert to factor if not already
    if (!is.factor(outcome_data)) {
      outcome_data <- as.factor(outcome_data)
    }
    
    # Ensure we have valid data
    if (length(outcome_data) != nrow(expr_data)) {
      min_length <- min(length(outcome_data), nrow(expr_data))
      expr_data <- expr_data[1:min_length, ]
      outcome_data <- outcome_data[1:min_length]
      cat("Adjusted data to match dimensions\n")
    }
    
    # Remove genes with too many missing values
    gene_missing_prop <- colMeans(is.na(expr_data))
    keep_genes <- gene_missing_prop < 0.5  # Keep genes with <50% missing values
    expr_data <- expr_data[, keep_genes]
    
    # Impute remaining missing values with column means
    for (j in 1:ncol(expr_data)) {
      missing_idx <- is.na(expr_data[, j])
      if (any(missing_idx)) {
        expr_data[missing_idx, j] <- mean(expr_data[, j], na.rm = TRUE)
      }
    }
    
    cat("VEER1 dataset loaded:\n")
    cat("  Samples:", nrow(expr_data), "\n")
    cat("  Genes:",   ncol(expr_data), "\n")
    cat("  Outcome classes (pre‐simulation):\n")
    print(table(outcome_data))
    
    return(list(X = expr_data, y = outcome_data, 
                n = nrow(expr_data), p = ncol(expr_data)))
  }
  
  stop("Dataset not recognized: ", dataset_name, 
       ". Available options: 'leukemia', 'veer1'")
}

# ────────────────────────────────────────────────────────────────────────────
# 2. Outcome generation function (unchanged)
# ────────────────────────────────────────────────────────────────────────────
generate_outcome_logit <- function(X, effect_vars, effect_sizes) {
  n <- nrow(X)
  
  # Initialize linear predictor
  linear_pred <- rep(0, n)
  
  # Vectorized matrix multiplication (much faster)
  if (length(effect_vars) > 0 && length(effect_sizes) > 0) {
    # Handle missing values
    X_clean <- X[, effect_vars, drop = FALSE]
    X_clean[is.na(X_clean)] <- 0
    
    # Matrix multiplication for efficiency
    linear_pred <- as.vector(X_clean %*% effect_sizes)
  } else {
    linear_pred <- rep(0, n)
  }
  
  # Convert to probabilities using logistic function
  prob <- plogis(linear_pred)
  
  # Generate binary outcome
  outcome <- rbinom(n, 1, prob)
  return(as.factor(outcome))
}

# ────────────────────────────────────────────────────────────────────────────
# 3. P-value calculation (unchanged)
# ────────────────────────────────────────────────────────────────────────────
calculate_pvalues_corrected <- function(importance_values) {
  # Create null distribution M = M1 ∪ M2 ∪ M3
  M1 <- importance_values[importance_values < 0]
  M2 <- importance_values[importance_values == 0]
  M3 <- -M1  # Mirror negative values
  
  # Combine to form null distribution
  M <- c(M1, M2, M3)
  
  if (length(M) == 0) {
    # No null distribution available
    return(rep(0.5, length(importance_values)))
  }
  
  # Calculate p-values as proportion of null distribution >= observed value
  pvalues <- sapply(importance_values, function(x) {
    mean(M >= x)
  })
  return(pvalues)
}

# ────────────────────────────────────────────────────────────────────────────
# 4. Helper functions for different correlation structures (unchanged)
# ────────────────────────────────────────────────────────────────────────────

# Function to create block correlation structure
create_block_correlation <- function(X, effect_vars, block_size = 10, within_block_cor = 0.7, debug = FALSE) {
  X_new <- X
  n <- nrow(X)
  
  # Validate correlation parameter
  if (abs(within_block_cor) >= 1) {
    stop("within_block_cor must be between -1 and 1 (exclusive)")
  }
  
  # Create correlated blocks among effect variables
  n_blocks <- ceiling(length(effect_vars) / block_size)
  
  for (block_idx in 1:n_blocks) {
    start_idx <- (block_idx - 1) * block_size + 1
    end_idx <- min(block_idx * block_size, length(effect_vars))
    block_vars <- effect_vars[start_idx:end_idx]
    
    if (length(block_vars) > 1) {
      # Create correlated structure within this block
      base_var <- X[, block_vars[1]]  # Use first variable as base
      base_var <- scale(base_var)[, 1]  # Standardize for correlation calculation
      
      for (i in 2:length(block_vars)) {
        var_idx <- block_vars[i]
        # Create correlation with base variable using Cholesky-like approach
        noise <- rnorm(n, 0, 1)
        noise <- scale(noise)[, 1]  # Standardize noise
        
        # Create desired correlation: X_new = ρ * base + sqrt(1-ρ²) * noise
        rho <- within_block_cor
        X_new[, var_idx] <- rho * base_var + sqrt(1 - rho^2) * noise
      }
    }
  }
  
  if (debug) {
    cat("Created", n_blocks, "correlated blocks of size ~", block_size, "\n")
    cat("Target within-block correlation:", within_block_cor, "\n")
    
    # Verify actual correlations in first block
    if (length(effect_vars) >= 2) {
      first_block_end <- min(block_size, length(effect_vars))
      first_block_vars <- effect_vars[1:first_block_end]
      if (length(first_block_vars) > 1) {
        actual_cor <- cor(X_new[, first_block_vars[1]], X_new[, first_block_vars[2]])
        cat("Actual correlation in first block:", round(actual_cor, 3), "\n")
      }
    }
  }
  
  return(X_new)
}

# Function to induce specific correlation patterns
induce_correlation_patterns <- function(X, effect_vars, pairwise_cor = 0.8, effect_null_cor = 0.5, debug = FALSE) {
  X_new <- X
  n <- nrow(X)
  
  # Validate correlation parameters
  if (abs(pairwise_cor) >= 1) {
    stop("pairwise_cor must be between -1 and 1 (exclusive)")
  }
  if (abs(effect_null_cor) >= 1) {
    stop("effect_null_cor must be between -1 and 1 (exclusive)")
  }
  
  # Pattern 1: High correlation between adjacent effect variables
  pairs_created <- 0
  for (i in 1:(length(effect_vars) - 1)) {
    if (i %% 2 == 1) {  # Every other pair
      var1_idx <- effect_vars[i]
      var2_idx <- effect_vars[i + 1]
      
      # Create desired correlation between these two variables
      base_var <- X[, var1_idx]
      base_var <- scale(base_var)[, 1]  # Standardize
      
      noise <- rnorm(n, 0, 1)
      noise <- scale(noise)[, 1]  # Standardize noise
      
      # Create exact correlation: X_new = ρ * base + sqrt(1-ρ²) * noise
      rho <- pairwise_cor
      X_new[, var2_idx] <- rho * base_var + sqrt(1 - rho^2) * noise
      pairs_created <- pairs_created + 1
    }
  }
  
  # Pattern 2: Create some correlation with null variables
  nulls_correlated <- 0
  if (length(effect_vars) < ncol(X)) {
    null_vars <- setdiff(1:ncol(X), effect_vars)
    n_correlated_nulls <- min(20, length(null_vars))
    
    for (i in 1:n_correlated_nulls) {
      effect_var <- effect_vars[i %% length(effect_vars) + 1]
      null_var <- null_vars[i]
      
      # Create moderate correlation between effect and null variable
      base_var <- X[, effect_var]
      base_var <- scale(base_var)[, 1]  # Standardize
      
      noise <- rnorm(n, 0, 1)
      noise <- scale(noise)[, 1]  # Standardize noise
      
      # Create exact correlation
      rho <- effect_null_cor
      X_new[, null_var] <- rho * base_var + sqrt(1 - rho^2) * noise
      nulls_correlated <- nulls_correlated + 1
    }
  }
  
  if (debug) {
    cat("Induced correlation patterns:\n")
    cat("  Effect variable pairs created:", pairs_created, "(target cor =", pairwise_cor, ")\n")
    cat("  Effect-null correlations created:", nulls_correlated, "(target cor =", effect_null_cor, ")\n")
    
    # Verify actual correlations for first pair
    if (pairs_created > 0 && length(effect_vars) >= 2) {
      actual_cor <- cor(X_new[, effect_vars[1]], X_new[, effect_vars[2]])
      cat("  Actual correlation in first pair:", round(actual_cor, 3), "\n")
    }
    
    # Verify effect-null correlation
    if (nulls_correlated > 0) {
      null_vars <- setdiff(1:ncol(X), effect_vars)
      if (length(null_vars) > 0) {
        actual_effect_null_cor <- cor(X_new[, effect_vars[1]], X_new[, null_vars[1]])
        cat("  Actual effect-null correlation:", round(actual_effect_null_cor, 3), "\n")
      }
    }
  }
  
  return(X_new)
}

# Function to print correlation diagnostics
print_correlation_diagnostics <- function(X, effect_vars, n_effect_vars) {
  cat("\n=== CORRELATION DIAGNOSTICS ===\n")
  
  # Check correlations among first few effect variables
  if (n_effect_vars >= 2) {
    effect_subset <- X[, effect_vars[1:min(10, n_effect_vars)]]
    effect_cor_matrix <- cor(effect_subset)
    
    # Average correlation among effect variables
    effect_cors <- effect_cor_matrix[upper.tri(effect_cor_matrix)]
    cat("Effect variables correlations:\n")
    cat("  Mean absolute correlation:", round(mean(abs(effect_cors)), 3), "\n")
    cat("  Range:", round(range(effect_cors), 3), "\n")
  }
  
  # Check correlations between effect and null variables
  if (ncol(X) > n_effect_vars) {
    null_start <- max(n_effect_vars + 1, 100)
    null_end <- min(ncol(X), null_start + 9)
    
    if (null_start <= ncol(X)) {
      effect_subset <- X[, effect_vars[1:min(5, n_effect_vars)]]
      null_subset <- X[, null_start:null_end]
      
      cross_cors <- cor(effect_subset, null_subset)
      cat("Effect-null correlations:\n")
      cat("  Mean absolute correlation:", round(mean(abs(cross_cors)), 3), "\n")
      cat("  Range:", round(range(cross_cors), 3), "\n")
    }
  }
  cat("===============================\n\n")
}

# ────────────────────────────────────────────────────────────────────────────
# 5. Enhanced power analysis with correlation control (UPDATED)
# ────────────────────────────────────────────────────────────────────────────
power_analysis_with_correlation_control <- function(dataset_name, 
                                                    samples_to_remove = NULL,
                                                    n_replications  = 1000,
                                                    effect_sizes    = c(-1, 1, -2, 2, -3, 3, -4, 4),
                                                    n_trees         = 5000,
                                                    alpha           = 0.05,
                                                    mtry_value      = 500,
                                                    n_cores         = NULL,
                                                    scale_data      = TRUE,
                                                    correlation_method = "none",
                                                    # Correlation control parameters
                                                    block_size      = 10,
                                                    within_block_cor = 0.7,
                                                    pairwise_cor    = 0.8,
                                                    effect_null_cor = 0.5,
                                                    debug           = FALSE,
                                                    seed            = 2025) {
  set.seed(seed)
  
  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  
  # Load dataset
  dataset <- load_real_dataset(dataset_name, samples_to_remove)
  X_raw <- dataset$X
  n     <- dataset$n
  p     <- dataset$p
  
  # Data scaling
  if (scale_data) {
    gene_sds <- apply(X_raw, 2, sd, na.rm = TRUE)
    gene_sds[gene_sds == 0] <- 1
    X_scaled <- scale(X_raw, center = TRUE, scale = gene_sds)
    cat("  Data scaling: Standard normal (mean=0, sd=1)\n")
  } else {
    X_scaled <- X_raw
    cat("  Data scaling: None (using raw log-expression values)\n")
  }
  
  # Set up effect parameters (but NOT the actual variables yet)
  n_effect_sizes <- length(effect_sizes)
  vars_per_effect <- 10
  n_effect_vars   <- min(n_effect_sizes * vars_per_effect, p)
  
  cat("\nPower analysis setup:\n")
  cat("  Dataset:", dataset_name, "\n")
  cat("  Dataset dimensions:", n, "x", p, "\n")
  if (!is.null(samples_to_remove)) {
    cat("  Removed samples:", paste(samples_to_remove, collapse = ", "), "\n")
  }
  cat("  Number of effect variables:", n_effect_vars, "\n")
  cat("  Effect-size levels:", unique(effect_sizes), "\n")
  cat("  mtry =", mtry_value, "\n")
  cat("  n_trees =", n_trees, "\n")
  cat("  Alpha level =", alpha, "\n")
  cat("  Correlation method:", correlation_method, "\n")
  
  # Print correlation parameters if relevant
  if (correlation_method %in% c("block", "induced")) {
    cat("  Correlation parameters:\n")
    if (correlation_method == "block") {
      cat("    Block size:", block_size, "\n")
      cat("    Within-block correlation:", within_block_cor, "\n")
    } else if (correlation_method == "induced") {
      cat("    Pairwise correlation:", pairwise_cor, "\n")
      cat("    Effect-null correlation:", effect_null_cor, "\n")
    }
  }
  cat("\n")
  
  if (correlation_method == "none") {
    # Original approach: permute each column to remove all correlations
    cat("  Removing all correlations via permutation\n")
    X_final <- apply(X_scaled, 2, sample)
    
  } else if (correlation_method == "preserve") {
    # Preserve original correlations by not permuting
    cat("  Preserving original correlations (no permutation)\n")
    X_final <- X_scaled
    
  } else if (correlation_method == "block") {
    # First permute, then create block correlation structure
    cat("  Creating block correlation structure\n")
    X_permuted <- apply(X_scaled, 2, sample)
    # Note: We'll apply block correlation to randomly selected effect vars in each replication
    X_final <- X_permuted
    
  } else if (correlation_method == "induced") {
    # First permute, then induce specific correlation patterns
    cat("  Inducing custom correlation patterns\n")
    X_permuted <- apply(X_scaled, 2, sample)
    X_final <- X_permuted
    
  } else {
    stop("Unknown correlation_method: ", correlation_method, 
         ". Use 'none', 'preserve', 'block', or 'induced'")
  }
  
  # Set up parallel processing
  cl <- makeCluster(n_cores, type = "PSOCK")
  registerDoParallel(cl)
  registerDoRNG(seed = seed + 1)
  
  bad_class_balance <- 0L
  
  # Run replications in parallel
  results_list <- foreach(rep_id = 1:n_replications,
                          .combine  = 'rbind',
                          .packages = c('ranger', 'vita'),
                          .export   = c('generate_outcome_logit', 'calculate_pvalues_corrected',
                                        'create_block_correlation', 'induce_correlation_patterns'),
                          .inorder  = FALSE) %dopar% {
                            
                            tryCatch({
                              effect_vars <- sample(1:p, n_effect_vars, replace = FALSE)
                              effect_size_vector <- rep(effect_sizes, each = vars_per_effect, length.out = n_effect_vars)
                              
                              # Apply correlation structure if needed (for block/induced methods)
                              X_rep <- X_final  # Start with the base transformed data
                              
                              if (correlation_method == "block") {
                                X_rep <- create_block_correlation(X_rep, effect_vars, 
                                                                  block_size = block_size,
                                                                  within_block_cor = within_block_cor,
                                                                  debug = FALSE)
                              } else if (correlation_method == "induced") {
                                X_rep <- induce_correlation_patterns(X_rep, effect_vars,
                                                                     pairwise_cor = pairwise_cor,
                                                                     effect_null_cor = effect_null_cor,
                                                                     debug = FALSE)
                              }
                              
                              # Generate outcome using the selected effect variables
                              y <- generate_outcome_logit(X_rep, effect_vars, effect_size_vector)
                              
                              # Check class balance
                              taby <- table(y)
                              if (length(taby) < 2) {
                                bad_class_balance <<- bad_class_balance + 1L
                                warning("Outcome has only one class in replication ", rep_id)
                                return(NULL)
                              }
                              
                              # Combine into data.frame for RF
                              data_combined <- data.frame(y = y, X_rep)
                              
                              # Fit three forests
                              rf_air <- ranger(y ~ ., data = data_combined,
                                               num.trees     = n_trees,
                                               mtry          = mtry_value,
                                               importance    = "impurity_corrected",
                                               min.node.size = 1,
                                               num.threads   = 1)
                              
                              rf_perm <- ranger(y ~ ., data = data_combined,
                                                num.trees     = n_trees,
                                                mtry          = mtry_value,
                                                importance    = "permutation",
                                                min.node.size = 1,
                                                num.threads   = 1)
                              
                              rf_holdout <- holdoutRF(y ~ ., data = data_combined,
                                                      num.trees     = n_trees,
                                                      mtry          = mtry_value,
                                                      min.node.size = 1)
                              
                              # Extract importance values
                              imp_air     <- rf_air$variable.importance
                              imp_perm    <- rf_perm$variable.importance
                              imp_holdout <- rf_holdout$variable.importance
                              
                              # Calculate p-values
                              pvals_air     <- calculate_pvalues_corrected(imp_air)
                              pvals_perm    <- calculate_pvalues_corrected(imp_perm)
                              pvals_holdout <- calculate_pvalues_corrected(imp_holdout)
                              
                              # Create full effect-size vector based on THIS replication's effect vars
                              effect_size_full <- numeric(p)
                              effect_size_full[effect_vars] <- effect_size_vector
                              
                              # Return results
                              data.frame(
                                replication       = rep_id,
                                variable          = 1:p,
                                effect_size       = effect_size_full,
                                pvalue_air        = pvals_air,
                                pvalue_perm       = pvals_perm,
                                pvalue_holdout    = pvals_holdout,
                                rejected_air      = pvals_air <= alpha,
                                rejected_perm     = pvals_perm <= alpha,
                                rejected_holdout  = pvals_holdout <= alpha
                              )
                              
                            }, error = function(e) {
                              cat("Error in replication", rep_id, ":", e$message, "\n")
                              return(NULL)
                            })
                          }
  
  # Clean up
  stopCluster(cl)
  
  if (bad_class_balance > 0) {
    cat("\n### WARNING: ", bad_class_balance,
        " replications produced a single-class outcome.\n")
  }
  
  # Remove NULL replications
  results_list <- results_list[!sapply(results_list, is.null), ]
  return(results_list)
}


# ────────────────────────────────────────────────────────────────────────────
# 6. Summarize function (unchanged)
# ────────────────────────────────────────────────────────────────────────────
summarize_power_results_corrected <- function(results) {
  # Type I error (variables with effect_size == 0)
  type1_data <- results[results$effect_size == 0, ]
  
  type1_error <- data.frame(
    AIR         = mean(type1_data$rejected_air),
    Permutation = mean(type1_data$rejected_perm),
    Holdout     = mean(type1_data$rejected_holdout)
  )
  
  # Power by effect size – including 0 (null variables)
  results$abs_effect_size <- abs(results$effect_size)
  power_summary <- aggregate(
    cbind(rejected_air, rejected_perm, rejected_holdout) ~ abs_effect_size,
    data = results,
    FUN  = mean
  )
  names(power_summary)[1] <- "effect_size"
  
  cat("\nType I Error (≈ 0.05):\n")
  print(round(type1_error, 4))
  
  cat("\nPower by Effect Size (including 0 for null):\n")
  print(round(power_summary, 3))
  
  return(list(type1_error = type1_error, power = power_summary))
}

# ────────────────────────────────────────────────────────────────────────────
# 7. Plot function (unchanged)
# ────────────────────────────────────────────────────────────────────────────
plot_power_results_final <- function(summary_results, dataset_name = "") {
  power_data <- summary_results$power
  
  # Reshape data for plotting
  power_long <- data.frame(
    effect_size = rep(power_data$effect_size, 3),
    method      = rep(c("AIR", "Permutation", "Holdout"), each = nrow(power_data)),
    power       = c(power_data$rejected_air,
                    power_data$rejected_perm,
                    power_data$rejected_holdout)
  )
  power_long$facet_label <- paste("Power Analysis Results")
  
  # Determine y-axis limit
  max_power <- max(power_long$power, na.rm = TRUE)
  y_max     <- ceiling(max_power / 0.05) * 0.05
  if (y_max < 0.25) y_max <- 0.25
  
  p <- ggplot(power_long, aes(x = effect_size, y = power, color = method, shape = method)) +
    geom_point(size = 3) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 0.05, color = "red", linewidth = 1) +
    facet_wrap(~facet_label) +
    labs(
      x = "Effect Size",
      y = "Power (Proportion of Rejections)",
      color = NULL,
      shape = NULL
    ) +
    scale_y_continuous(breaks = seq(0, y_max, by = 0.05), limits = c(0, y_max)) +
    theme_light(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey80", colour = "black", linewidth = 1),
      strip.text       = element_text(face = "bold", colour = "black", size = rel(1.2)),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(linewidth = 0.5),
      panel.grid.minor = element_line(linewidth = 0.25),
      axis.line        = element_line(linewidth = 1, colour = "black"),
      axis.ticks       = element_line(linewidth = 1, colour = "black"),
      axis.text        = element_text(face = "bold", size = rel(1)),
      axis.title       = element_text(face = "bold", size = rel(1.1)),
      axis.title.y     = element_text(margin = margin(r = 15)),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(face = "bold", size = rel(1)),
      aspect.ratio     = 1
    )
  
  print(p)
  return(p)
}

# ────────────────────────────────────────────────────────────────────────────
# 8. Enhanced wrapper function with all correlation parameters (updated for VEER1)
# ────────────────────────────────────────────────────────────────────────────
run_correlation_study <- function(dataset_name    = "veer1",
                                  samples_to_remove = c(54, 32),
                                  n_replications  = 100,
                                  correlation_method = "preserve",
                                  # Correlation control parameters
                                  block_size      = 10,
                                  within_block_cor = 0.7,
                                  pairwise_cor    = 0.8,
                                  effect_null_cor = 0.5,
                                  n_cores         = NULL,
                                  mtry_value      = 500,
                                  scale_data      = TRUE,
                                  debug           = FALSE) {
  
  cat("Running correlation study for", dataset_name, "\n")
  cat("  Correlation method:", correlation_method, "\n")
  if (!is.null(samples_to_remove)) {
    cat("  Removing samples:", paste(samples_to_remove, collapse = ", "), "\n")
  }
  cat("  Replications:", n_replications, "\n")
  cat("  Data scaling:", ifelse(scale_data, "Standard normal", "Raw values"), "\n")
  cat("  Debug mode:", debug, "\n\n")
  
  results <- power_analysis_with_correlation_control(
    dataset_name       = dataset_name,
    samples_to_remove  = samples_to_remove,
    n_replications     = n_replications,
    correlation_method = correlation_method,
    block_size         = block_size,
    within_block_cor   = within_block_cor,
    pairwise_cor       = pairwise_cor,
    effect_null_cor    = effect_null_cor,
    mtry_value         = mtry_value,
    scale_data         = scale_data,
    debug              = debug
  )
  
  # Summarize and plot
  summary <- summarize_power_results_corrected(results)
  plot_power_results_final(summary)
  return(list(raw_results = results, summary = summary))
}

# ────────────────────────────────────────────────────────────────────────────
# 9. Enhanced debug function for single replication (UPDATED)
# ────────────────────────────────────────────────────────────────────────────
debug_single_replication <- function(dataset_name = "veer1", 
                                     samples_to_remove = c(54, 32), 
                                     scale_data = TRUE,
                                     correlation_method = "none") {
  cat("\n=== DEBUG: Single Replication Analysis ===\n")
  
  dataset <- load_real_dataset(dataset_name, samples_to_remove)
  X <- dataset$X
  n <- dataset$n
  p <- dataset$p
  
  # Data scaling based on scale_data parameter
  if (scale_data) {
    # Scale each gene to standard normal (mean=0, sd=1)
    gene_sds <- apply(X, 2, sd, na.rm = TRUE)
    gene_sds[gene_sds == 0] <- 1  # Handle zero variance genes
    X_scaled <- scale(X, center = TRUE, scale = gene_sds)
    cat("Using standard normal scaling (mean=0, sd=1)\n")
  } else {
    # Use raw data without scaling
    X_scaled <- X
    cat("Using raw log-expression values (no scaling)\n")
  }
  
  # Apply correlation method
  set.seed(999)
  if (correlation_method == "none") {
    X_final <- apply(X_scaled, 2, sample)
    cat("Applied permutation to remove correlations\n")
  } else if (correlation_method == "preserve") {
    X_final <- X_scaled
    cat("Preserved original correlations\n")
  } else if (correlation_method == "block") {
    X_permuted <- apply(X_scaled, 2, sample)
    X_final <- X_permuted  # Will apply block correlation after selecting effect vars
  } else if (correlation_method == "induced") {
    X_permuted <- apply(X_scaled, 2, sample)
    X_final <- X_permuted  # Will apply induced correlation after selecting effect vars
  }
  
  # RANDOMLY SELECT effect variables
  n_effect_vars <- min(80, p)
  set.seed(123)
  effect_vars <- sample(1:p, n_effect_vars, replace = FALSE)
  effect_sizes <- rep(c(-1, 1, -2, 2, -3, 3, -4, 4), each = 10, length.out = n_effect_vars)
  
  cat("\nRandomly selected effect variables (first 20):", effect_vars[1:20], "\n")
  
  # Apply correlation structures if needed
  if (correlation_method == "block") {
    X_final <- create_block_correlation(X_final, effect_vars, debug = TRUE)
  } else if (correlation_method == "induced") {
    X_final <- induce_correlation_patterns(X_final, effect_vars, debug = TRUE)
  }
  
  # Generate outcome
  y <- generate_outcome_logit(X_final, effect_vars, effect_sizes)
  
  # Check linear predictor
  linear_pred <- rep(0, n)
  for (i in seq_along(effect_vars)) {
    x_raw <- X_final[, effect_vars[i]]
    x_raw[is.na(x_raw)] <- 0
    linear_pred <- linear_pred + effect_sizes[i] * x_raw
  }
  cat("\nLinear predictor stats:\n")
  cat("  Range:", round(range(linear_pred), 2), "\n")
  cat("  Mean:", round(mean(linear_pred), 2), "\n")
  cat("  SD:", round(sd(linear_pred), 2), "\n")
  
  prob <- plogis(linear_pred)
  cat("\nProbability stats:\n")
  cat("  Range:", round(range(prob), 2), "\n")
  cat("  Class distribution:", table(y), "\n")
  cat("  Class proportions:", round(prop.table(table(y)), 2), "\n")
  
  # Print correlation diagnostics
  print_correlation_diagnostics(X_final, effect_vars, n_effect_vars)
  
  # Test variable importance detection
  cat("\n=== Testing importance detection ===\n")
  data_combined <- data.frame(y = y, X_final)
  rf_test <- ranger(y ~ ., data = data_combined,
                    num.trees = 1000,
                    mtry = 500,
                    importance = "impurity_corrected",
                    min.node.size = 1)
  
  imp <- rf_test$variable.importance
  cat("Mean importance for effect vars:", mean(imp[effect_vars]), "\n")
  cat("Mean importance for null vars:", mean(imp[-effect_vars]), "\n")
  cat("Top 10 variables by importance:", which(rank(-imp) <= 10), "\n")
  cat("How many are true effect vars:", sum(which(rank(-imp) <= 10) %in% effect_vars), "\n")
}

# ────────────────────────────────────────────────────────────────────────────
# 10. VEER1 dataset exploration function
# ────────────────────────────────────────────────────────────────────────────
explore_veer1_dataset <- function() {
  cat("=== VEER1 DATASET EXPLORATION ===\n")
  
  # Load the dataset
  data(VEER1)
  
  # Basic information
  cat("Dataset class:", class(VEER1), "\n")
  cat("Expression data dimensions:", dim(exprs(VEER1)), "\n")
  cat("Feature data dimensions:", dim(fData(VEER1)), "\n")
  cat("Phenotype data dimensions:", dim(pData(VEER1)), "\n\n")
  
  # Expression data summary
  expr_data <- exprs(VEER1)
  cat("Expression data summary:\n")
  cat("  Range:", round(range(expr_data, na.rm = TRUE), 2), "\n")
  cat("  Mean:", round(mean(expr_data, na.rm = TRUE), 2), "\n")
  cat("  Missing values:", sum(is.na(expr_data)), "\n\n")
  
  # Phenotype data
  pheno_data <- pData(VEER1)
  cat("Phenotype columns:\n")
  for (i in 1:ncol(pheno_data)) {
    col_name <- colnames(pheno_data)[i]
    col_class <- class(pheno_data[, i])
    cat("  ", col_name, "(", col_class, ")\n")
    
    if (is.factor(pheno_data[, i]) || is.character(pheno_data[, i])) {
      if (length(unique(pheno_data[, i])) <= 10) {
        print(table(pheno_data[, i], useNA = "always"))
      } else {
        cat("    Too many unique values to display\n")
      }
    } else {
      print(summary(pheno_data[, i]))
    }
    cat("\n")
  }
  
  # Feature data
  feature_data <- fData(VEER1)
  cat("Feature data columns:\n")
  print(colnames(feature_data))
  cat("\nFirst few feature entries:\n")
  print(head(feature_data))
  
  return(invisible(VEER1))
}

# ────────────────────────────────────────────────────────────────────────────
# 11. Example usage and testing (updated for VEER1)
# ────────────────────────────────────────────────────────────────────────────

# Explore the VEER1 dataset first
cat("Exploring VEER1 dataset structure...\n")
explore_veer1_dataset()

# Test 1: No correlation with VEER1 dataset and sample removal
results_no_corr_veer1 <- run_correlation_study(
  dataset_name = "veer1",
  samples_to_remove = c(54, 32),
  correlation_method = "none",
  n_replications = 1000,
  debug = FALSE
)

# Test 2: Preserve original correlations with VEER1
results_preserve_corr_veer1 <- run_correlation_study(
  dataset_name = "leukemia",
  samples_to_remove = c(54, 32),
  correlation_method = "preserve",
  n_replications = 200,
  debug = FALSE
)

# Test 3: Block correlation structure with custom parameters using VEER1
results_block_corr_veer1 <- run_correlation_study(
  dataset_name = "veer1",
  samples_to_remove = c(54, 32),
  correlation_method = "block",
  block_size = 10,
  within_block_cor = 0.8,
  n_replications = 200,
  debug = FALSE
)

# Test 4: Induced correlation patterns with custom parameters using VEER1
results_induced_corr_veer1 <- run_correlation_study(
  dataset_name = "veer1",
  samples_to_remove = c(54, 32),
  correlation_method = "induced",
  pairwise_cor = 0.3,
  effect_null_cor = 0.3,
  n_replications = 200, 
  debug = TRUE
)

# Test 5: Compare multiple correlation levels for block method using VEER1
correlation_levels <- c(0.3, 0.5, 0.7, 0.9)
block_results_veer1 <- list()
for (corr_level in correlation_levels) {
  cat("\n=== Testing VEER1 block correlation =", corr_level, "===\n")
  block_results_veer1[[paste0("corr_", corr_level)]] <- run_correlation_study(
    dataset_name = "veer1",
    samples_to_remove = c(54, 32),
    correlation_method = "block",
    within_block_cor = corr_level,
    n_replications = 30,
    debug = TRUE
  )
}

# Test 6: Debug single replication with different correlation methods using VEER1
debug_single_replication("veer1", samples_to_remove = c(54, 32), correlation_method = "none")
debug_single_replication("veer1", samples_to_remove = c(54, 32), correlation_method = "preserve")
debug_single_replication("veer1", samples_to_remove = c(54, 32), correlation_method = "block")

# Test 7: Compare VEER1 with and without sample removal
cat("\n=== Comparing VEER1 with and without sample removal ===\n")
results_with_removal <- run_correlation_study(
  dataset_name = "veer1",
  samples_to_remove = c(54, 32),
  correlation_method = "none",
  n_replications = 200,
  debug = TRUE
)

results_without_removal <- run_correlation_study(
  dataset_name = "veer1",
  samples_to_remove = NULL,
  correlation_method = "preserve",
  n_replications = 30,
  debug = TRUE
)

# Test 8: Leukemia dataset for comparison (original functionality)
results_leukemia <- run_correlation_study(
  dataset_name = "leukemia",
  samples_to_remove = NULL,
  correlation_method = "preserve",
  n_replications = 200,
  debug = FALSE
)

# Test 3: Block correlation structure with custom parameters using VEER1
results_block_corr_veer1 <- run_correlation_study(
  dataset_name = "leukemia",
  samples_to_remove = NULL,
  correlation_method = "block",
  block_size = 10,
  within_block_cor = 0.1,
  n_replications = 200,
  debug = FALSE
)

# Test 4: Induced correlation patterns with custom parameters using VEER1
results_induced_corr_veer1 <- run_correlation_study(
  dataset_name = "leukemia",
  samples_to_remove = NULL,
  correlation_method = "induced",
  pairwise_cor = 0.9,
  effect_null_cor = 0.9,
  n_replications = 200, 
  debug = FALSE
)
