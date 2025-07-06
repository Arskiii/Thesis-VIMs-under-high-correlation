# ═══════════════════════════════════════════════════════════════════════════════
# CONSOLIDATED RANDOM FOREST VARIABLE IMPORTANCE SIMULATION
# ═══════════════════════════════════════════════════════════════════════════════

library(ranger)
library(MASS)
library(Matrix)
library(ggplot2)

# ─── HELPER FUNCTIONS ──────────────────────────────────────────────────────────

# Single replicate: fit four forests
run_one <- function(df, formula, ntree = 50, is_classification = TRUE, node_size = 1) {
  target <- all.vars(formula)[1]
  
  # (a) Gini-impurity importance
  rf_imp <- ranger(
    formula, data = df,
    num.trees = ntree,
    min.node.size = node_size,
    importance = "impurity",
    respect.unordered.factors = "partition"
  )
  
  # (b) Permutation (OOB) importance
  rf_perm <- ranger(
    formula, data = df,
    num.trees = ntree,
    min.node.size = node_size,
    importance = "permutation",
    respect.unordered.factors = "partition"
  )
  
  # (c) Holdout importance using vita::holdoutRF() - UPDATED
  rf_holdout <- holdoutRF(
    formula, data = df,
    num.trees = ntree,
    min.node.size = node_size,
    respect.unordered.factors = "partition"
  )
  
  # (d) AIR importance (aka "impurity_corrected" in ranger >=0.14)
  rf_air <- ranger(
    formula, data = df,
    num.trees = ntree,
    min.node.size = node_size,
    importance = "impurity_corrected",
    respect.unordered.factors = "partition"
  )
  
  list(
    impurity = rf_imp$variable.importance,
    permutation = rf_perm$variable.importance,
    holdout = rf_holdout$variable.importance,
    AIR = rf_air$variable.importance
  )
}

# Generate correlated data with mixed variable types
make_correlated_mixed <- function(n, p_vec, ord_levels, nom_levels, rho) {
  p_tot <- length(p_vec) + length(ord_levels) + length(nom_levels) + 1
  Sigma <- matrix(rho, p_tot, p_tot)
  diag(Sigma) <- 1
  Z <- mvrnorm(n, mu = rep(0, p_tot), Sigma = Sigma)
  U <- pnorm(Z)
  
  cols <- list()
  col <- 1
  
  # Binary variables
  for (i in seq_along(p_vec)) {
    cols[[paste0("B", p_vec[i])]] <- as.integer(U[, col] > (1 - p_vec[i]))
    col <- col + 1
  }
  
  # Ordered factors
  for (i in seq_along(ord_levels)) {
    k <- ord_levels[i]
    breaks <- seq(0, 1, length.out = k + 1)
    cols[[paste0("O", k)]] <- factor(findInterval(U[, col], breaks, rightmost.closed = TRUE),
                                     levels = 1:k, ordered = TRUE)
    col <- col + 1
  }
  
  # Nominal factors
  for (i in seq_along(nom_levels)) {
    k <- nom_levels[i]
    breaks <- seq(0, 1, length.out = k + 1)
    cols[[paste0("N", k)]] <- factor(findInterval(U[, col], breaks, rightmost.closed = TRUE),
                                     levels = 1:k)
    col <- col + 1
  }
  
  # Continuous
  cols[["C"]] <- Z[, col]
  
  as.data.frame(cols)
}

# ─── MAIN SIMULATION FUNCTION ──────────────────────────────────────────────────

#' Comprehensive Random Forest Variable Importance Simulation
#'
#' @param n_reps Number of simulation replicates
#' @param correlation_method Character: "none", "uniform", or "blockwise"
#' @param rho Numeric: For "none" use NULL, for "uniform" use single value, for "blockwise" use vector of 4 values
#' @param n_sample Sample size per replicate
#' @param causal_variable Character: Variable name that should be causal (e.g., "X8") or NULL for no causal effect
#' @param beta Numeric: Effect size for causal variable (log odds ratio for classification) or NULL
#' @param method Character: "maf" (minor allele frequency), "categorical", or "mixed"
#' @param is_classification Logical: TRUE for classification, FALSE for regression
#' @param ntree Number of trees in random forest
#' @param node_size Minimum node size
#' @param seed Random seed for reproducibility
#'
#' @return List containing simulation results and plots
rf_importance_simulation <- function(
    n_reps = 1000,
    correlation_method = "none",  # "none", "uniform", "blockwise"
    rho = NULL,                   # NULL, single value, or vector of 4 values
    n_sample = 100,
    causal_variable = NULL,       # NULL or variable name like "X8"
    beta = NULL,                  # NULL or effect size
    method = "maf",               # "maf", "categorical", "mixed"
    is_classification = TRUE,
    ntree = 50,
    node_size = NULL,
    seed = 2025
) {
  
  set.seed(seed)
  
  # Set default node size
  if (is.null(node_size)) {
    node_size <- if (is_classification) 1 else 5
  }
  
  # Define variable sequences based on method
  if (method == "maf") {
    var_seq <- seq(0.05, 0.50, by = 0.05)
    var_names <- paste0("X", seq_along(var_seq))
  } else if (method == "categorical") {
    var_seq <- c(2, 3, 4, 5, 6, 7, 8, 10, 20, 30)
    var_names <- paste0("X", seq_along(var_seq))
  } else if (method == "mixed") {
    p_vec <- c(0.05, 0.1, 0.2, 0.5)
    ord_levels <- c(5, 10)
    nom_levels <- c(5, 8, 10)
  } else {
    stop("Method must be 'maf', 'categorical', or 'mixed'")
  }
  
  # Setup correlation structure
  if (correlation_method == "none") {
    if (method != "mixed") {
      # Independent variables
      correlation_matrix <- NULL
    }
  } else if (correlation_method == "uniform") {
    if (is.null(rho)) stop("rho must be specified for uniform correlation")
    if (method != "mixed") {
      k <- length(var_seq)
      correlation_matrix <- matrix(rho, k, k)
      diag(correlation_matrix) <- 1
    }
  } else if (correlation_method == "blockwise") {
    if (is.null(rho) || length(rho) != 4) {
      stop("rho must be a vector of 4 values for blockwise correlation")
    }
    # Build block-diagonal correlation matrix
    block_sizes <- rep(10, 4)
    Sigma_blocks <- lapply(seq_along(rho), function(i) {
      m <- block_sizes[i]
      rho_val <- rho[i]
      M <- matrix(rho_val, nrow = m, ncol = m)
      diag(M) <- 1
      M
    })
    correlation_matrix <- as.matrix(bdiag(Sigma_blocks))
    var_seq <- rep(0.1, 40)  # 40 variables for blockwise
    var_names <- paste0("X", seq_along(var_seq), "_ρ", rep(rho, each = 10))
  }
  
  # Storage for results
  res_list <- vector("list", n_reps)
  pb <- txtProgressBar(min = 1, max = n_reps, style = 3)
  
  # Main simulation loop
  for (i in seq_len(n_reps)) {
    
    if (method == "maf") {
      if (correlation_method == "none") {
        # Independent binary variables
        Xmat <- sapply(var_seq, function(p) rbinom(n_sample, 1, p))
      } else if (correlation_method == "uniform") {
        # Correlated binary variables via latent MVN
        Z <- mvrnorm(n_sample, rep(0, length(var_seq)), correlation_matrix)
        Xmat <- sapply(seq_along(var_seq), function(j) {
          as.integer(Z[, j] < qnorm(var_seq[j]))
        })
      } else if (correlation_method == "blockwise") {
        # Blockwise correlated binary variables
        Z <- mvrnorm(n = n_sample, mu = rep(0, length(var_seq)), Sigma = correlation_matrix)
        Xmat <- ifelse(Z > 0, 1L, 0L)
      }
      colnames(Xmat) <- var_names
      
    } else if (method == "categorical") {
      if (correlation_method == "none") {
        # Independent categorical variables
        Xmat <- sapply(var_seq, function(k) sample(1:k, size = n_sample, replace = TRUE))
      } else if (correlation_method == "uniform") {
        # Correlated categorical variables via latent MVN
        Z <- mvrnorm(n_sample, rep(0, length(var_seq)), correlation_matrix)
        U <- pnorm(Z)
        Xmat <- sapply(seq_along(var_seq), function(j) {
          k <- var_seq[j]
          breaks <- seq(0, 1, length.out = k + 1)
          findInterval(U[, j], breaks, rightmost.closed = TRUE)
        })
      }
      colnames(Xmat) <- var_names
      
    } else if (method == "mixed") {
      if (correlation_method == "none") {
        # Independent mixed variables
        Bmat <- sapply(p_vec, function(p) rbinom(n_sample, 1, p))
        colnames(Bmat) <- paste0("B", p_vec)
        
        O5 <- factor(rep(1:5, length.out = n_sample), levels = 1:5, ordered = TRUE)
        O10 <- factor(rep(1:10, length.out = n_sample), levels = 1:10, ordered = TRUE)
        N5 <- factor(rep(1:5, length.out = n_sample))
        N8 <- factor(rep(1:8, length.out = n_sample))
        N10 <- factor(rep(1:8, length.out = n_sample))
        C <- rnorm(n_sample)
        
        Xmat <- data.frame(Bmat, O5, O10, N5, N8, N10, C)
      } else {
        # Correlated mixed variables
        if (is.null(rho)) rho <- 0.9
        Xmat <- make_correlated_mixed(n_sample, p_vec, ord_levels, nom_levels, rho)
      }
    }
    
    # Generate outcome variable
    if (is_classification) {
      if (is.null(causal_variable) || is.null(beta)) {
        # No causal effect
        Y <- factor(rbinom(n_sample, 1, 0.5), labels = c("class0", "class1"))
      } else {
        # Causal effect - handle different data structures
        if (method == "mixed") {
          # For mixed method, Xmat is a data.frame
          if (!causal_variable %in% names(Xmat)) {
            stop(paste("Causal variable", causal_variable, "not found in mixed variables.",
                       "Available variables:", paste(names(Xmat), collapse = ", ")))
          }
          X_causal <- Xmat[[causal_variable]]
        } else {
          # For other methods, Xmat is a matrix
          if (is.matrix(Xmat)) {
            if (!causal_variable %in% colnames(Xmat)) {
              stop(paste("Causal variable", causal_variable, "not found in matrix columns.",
                         "Available variables:", paste(colnames(Xmat), collapse = ", ")))
            }
            X_causal <- Xmat[, causal_variable]
          } else {
            # Handle case where Xmat might be a data.frame
            if (!causal_variable %in% names(Xmat)) {
              stop(paste("Causal variable", causal_variable, "not found.",
                         "Available variables:", paste(names(Xmat), collapse = ", ")))
            }
            X_causal <- Xmat[[causal_variable]]
          }
        }
        intercept <- - beta * mean(X_causal)
        p <- plogis(intercept + beta * X_causal)
        Y <- factor(rbinom(n_sample, 1, p), labels = c("class0", "class1"))
      }
    } else {
      # Regression
      if (is.null(causal_variable) || is.null(beta)) {
        # No causal effect
        Y <- rnorm(n_sample)
      } else {
        # Causal effect - handle different data structures
        if (method == "mixed") {
          # For mixed method, Xmat is a data.frame
          if (!causal_variable %in% names(Xmat)) {
            stop(paste("Causal variable", causal_variable, "not found in mixed variables.",
                       "Available variables:", paste(names(Xmat), collapse = ", ")))
          }
          X_causal <- Xmat[[causal_variable]]
        } else {
          # For other methods, Xmat is a matrix
          if (is.matrix(Xmat)) {
            if (!causal_variable %in% colnames(Xmat)) {
              stop(paste("Causal variable", causal_variable, "not found in matrix columns.",
                         "Available variables:", paste(colnames(Xmat), collapse = ", ")))
            }
            X_causal <- Xmat[, causal_variable]
          } else {
            # Handle case where Xmat might be a data.frame
            if (!causal_variable %in% names(Xmat)) {
              stop(paste("Causal variable", causal_variable, "not found.",
                         "Available variables:", paste(names(Xmat), collapse = ", ")))
            }
            X_causal <- Xmat[[causal_variable]]
          }
        }
        intercept <- 0  # You can adjust this if needed
        Y <- intercept + beta * X_causal + rnorm(n_sample)
      }
    }
    
    # Create data frame
    df <- data.frame(Y, Xmat)
    
    # Run simulation for this replicate
    res_list[[i]] <- run_one(df, Y ~ ., ntree = ntree, 
                             is_classification = is_classification, 
                             node_size = node_size)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("All", n_reps, "reps done!\n")
  
  # ─── PROCESS RESULTS ────────────────────────────────────────────────────────
  
  # Flatten results into data frame
  all_imp <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    rep_list <- res_list[[i]]
    do.call(rbind, lapply(names(rep_list), function(meth) {
      data.frame(
        method = meth,
        variable = names(rep_list[[meth]]),
        imp = as.numeric(rep_list[[meth]]),
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # Add variable properties and create factor variables for plotting
  if (correlation_method == "blockwise") {
    all_imp$corr <- as.numeric(sub("^.*_ρ(.+)$", "\\1", all_imp$variable))
    
    # Special handling for causal variable case
    if (!is.null(causal_variable)) {
      # Determine which block the causal variable is in
      causal_block_rho <- as.numeric(sub("^.*_ρ(.+)$", "\\1", causal_variable))
      causal_block_vars <- paste0("X", ((which(rho == causal_block_rho) - 1) * 10 + 1):
                                    (which(rho == causal_block_rho) * 10), "_ρ", causal_block_rho)
      
      # Create three categories
      all_imp$var_category <- ifelse(all_imp$variable == causal_variable, "Causal Variable",
                                     ifelse(all_imp$variable %in% causal_block_vars, "Correlated Neighbors", 
                                            "Other Variables"))
      all_imp$var_category <- factor(all_imp$var_category, 
                                     levels = c("Causal Variable", "Correlated Neighbors", "Other Variables"))
    } else {
      # Original blockwise plotting by correlation level
      all_imp$corr_f <- factor(sprintf("%.2f", all_imp$corr), levels = sprintf("%.2f", rho))
    }
  } else if (method == "maf") {
    # Create property mapping for MAF variables
    maf_mapping <- setNames(var_seq, var_names)
    all_imp$maf <- maf_mapping[all_imp$variable]
    all_imp$maf_f <- factor(sprintf("%.2f", all_imp$maf), 
                            levels = sprintf("%.2f", var_seq))
  } else if (method == "categorical") {
    # Create property mapping for categorical variables
    cat_mapping <- setNames(var_seq, var_names)
    all_imp$n_cats <- cat_mapping[all_imp$variable]
    all_imp$n_cats_f <- factor(all_imp$n_cats, levels = var_seq)
  } else if (method == "mixed") {
    # Extract variable type and property - keep original variable names as factor
    var_mapping <- c(
      setNames(p_vec, paste0("B", p_vec)),
      O5 = 5, O10 = 10, N5 = 5, N8 = 8, N10 = 10, C = 1
    )
    all_imp$var_prop <- var_mapping[all_imp$variable]
    # Keep variable as factor for x-axis ordering
    all_imp$variable <- factor(all_imp$variable, levels = names(var_mapping))
  }
  
  # Fix method factor levels
  all_imp$method <- factor(all_imp$method,
                           levels = c("impurity", "AIR", "permutation", "holdout"))
  
  # ─── HELPER FUNCTIONS FOR PLOTTING ─────────────────────────────────────────
  
  # Function to create consistent plot styling like Utils.R
  create_importance_plot <- function(data, x_var, x_label, title_suffix = "") {
    ggplot(data, aes_string(x = x_var, y = "imp")) +
      geom_boxplot(outlier.size = 0.3, size = 1.2, fatten = 1.2) +
      geom_hline(yintercept = 0, color = "red", size = 1) +
      facet_wrap(~method, ncol = 2, scales = "free_y") +
      labs(
        x = x_label,
        y = "Variable Importance Measure"
      ) +
      theme_light(base_size = 12) +
      theme(
        strip.background = element_rect(fill = "grey80", colour = "black", size = 1),
        strip.text       = element_text(face = "bold", colour = "black", size = rel(1)),
        panel.border     = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_line(size = 0.25),
        axis.line        = element_line(size = 1, colour = "black"),
        axis.ticks       = element_line(size = 1, colour = "black"),
        axis.text        = element_text(face = "bold", size = rel(1)),
        axis.title       = element_text(face = "bold", size = rel(1.1)),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        aspect.ratio     = 1
      )
  }
  
  # ─── CREATE PLOTS ───────────────────────────────────────────────────────────
  
  if (correlation_method == "blockwise") {
    if (!is.null(causal_variable)) {
      # Special plot for causal variable case: group by variable category
      plot_obj <- create_importance_plot(all_imp, "var_category", "Variable Category")
    } else {
      # Original blockwise plot by correlation level
      plot_obj <- create_importance_plot(all_imp, "corr_f", "Block Correlation")
    }
  } else if (method == "maf") {
    plot_obj <- create_importance_plot(all_imp, "maf_f", "Minor Allele Frequency")
  } else if (method == "categorical") {
    plot_obj <- create_importance_plot(all_imp, "n_cats_f", "Number of Categories")
  } else if (method == "mixed") {
    plot_obj <- create_importance_plot(all_imp, "variable", "Variable Type")
  }
  
  # Return results
  return(list(
    results = all_imp,
    plot = plot_obj,
    parameters = list(
      n_reps = n_reps,
      correlation_method = correlation_method,
      rho = rho,
      n_sample = n_sample,
      causal_variable = causal_variable,
      beta = beta,
      method = method,
      is_classification = is_classification,
      ntree = ntree,
      node_size = node_size,
      seed = seed
    )
  ))
}

# ─── EXAMPLE USAGE ──────────────────────────────────────────────────────────────

# Example 1: Basic MAF study with no correlation
result1 <- rf_importance_simulation(
  n_reps = 100,
  correlation_method = "none",
  method = "maf",
  is_classification = FALSE
)
print(result1$plot)

# Example 2: Categorical variables with uniform correlation
# result2 <- rf_importance_simulation(
#   n_reps = 1000,
#   correlation_method = "uniform",
#   rho = 0.9,
#   method = "categorical",
#   is_classification = TRUE
# )
# print(result2$plot)

# Example 3: Mixed variables with causal effect
# result3 <- rf_importance_simulation(
#   n_reps = 1000,
#   correlation_method = "uniform",
#   rho = 0.5,
#   method = "mixed",
#   causal_variable = "B0.1",
#   beta = log(2),
#   is_classification = TRUE
# )
# print(result3$plot)

# Example 4a: Blockwise correlation study without causal effect
result4a <- rf_importance_simulation(
  n_reps = 100,
  correlation_method = "blockwise",
  rho = c(0.1, 0.3, 0.5, 0.8),
  method = "maf",
  is_classification = TRUE
)
print(result4a$plot)

# Example 4b: Blockwise correlation study with causal effect (shows leakage)
result4b <- rf_importance_simulation(
  n_reps = 100,
  correlation_method = "blockwise",
  rho = c(0.1, 0.3, 0.5, 0.8),
  method = "maf",
  causal_variable = "X33_ρ0.8",  # Variable 25 is in the 0.5 correlation block (variables 21-30)
  beta = 2,
  is_classification = FALSE
)
print(result4b$plot)  # This will show: Causal Variable | Correlated Neighbors | Other Variables

# Example 5: Regression with MAF
result5 <- rf_importance_simulation(
   n_reps = 1000,
   correlation_method = "none",
   method = "maf",
   is_classification = FALSE
)
print(result5$plot)
