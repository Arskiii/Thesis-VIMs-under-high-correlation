This repo contains the code for my Econometrics and Operational Research thesis: "Actual Impurity Reduction under Strong Variable Correlation: Accuracy and Power Analysis."

**Abstract:**

**I used ChatGPT for generating small parts of the code, and structuring and commenting the code**

## Function of the Files

### NullCase.R 
This code includes: 
- retrieving the impurity, AIR, permutation, and holdout importance from the Ranger package (Nembrini et al. 2018),
- the creation of the synthetic data,
- generating the variable importance measure scores on the synthetic data,
- plotting the results

### PowerStudy.R
This code includes:
- retrieving the impurity, AIR, permutation, holdout importance from the Ranger package,
- retrieving the leukemia dataset from the golubEsets package (Golub et al. 1999) and breast cancer dataset from the cancerdata package (Van de Vijver et al. 2002),
  - To install these packages, first install "BiocManager" in R (install.packages("BiocManager")). Then install the golubEsets and cancerdata packages in the following way:
      - BiocManager::install("golubEsets")
      - BiocManager::install("cancerdata")
- process the downloaded datasets (if samples need to be removed),
- generate and measure the power for the AIR, permutation, and holdout importance,
- plotting the results

## Code Files Structure

In the files you'll find a few big functions. You can run these functions and instantly go to the end of the file where you can provide the specifics of what scenario you want to simulate. You can find various examples in the code for the scenarios you can simulate.
Variables to design your simulation:

### NullCase.R
  - *n_reps* — number of repetitions
  - *correlation_method* — choose between "none" (independent), "uniform" or "blockwise" (the three different correlation scenarios described in my thesis)
  - *rho* — choose the correlation between 0-1
  - *method* — choose between "maf", "categorial" or "mixed" (the three null cases described in my thesis)
  - *causal_variable* — choose a variable that should be causal to the outcome such as "X33_ρ0.8" (in blockwise correlation with rho = c(0.1, 0.3, 0.5, 0.8)) or B0.1 (in mixed method)
  - *beta* — for y=alpha+beta*x choose your beta
  - *is_classification* — TRUE for classification, FALSE for regression

**PowerStudy.R**
- *dataset_name* — choose "leukemia" or "breast_cancer"
- *correlation_method* — choose "block" (blockwise correlation), "induced" (uniform correlation), "preserve" (native correlation), "none" (independent)
- *block_size* — choose a positive natural number, only applicable for blockwise correlation
- *within_block_cor* — choose within block correlation 0-1, only applicable for blockwise correlation
- *n_replications* — number of repetitions
- *debug* — TRUE for extra debugging output
- *samples_to_remove* — to remove samples from the dataset, only applicable for the leukemia dataset
- *pairwise_cor* — set correlation 0-1 between causal and null variables, only applicable for induced
- *effect_null_cor* — set correlation 0-1 between null variables, only applicable for induced
