# Pairwise_imputation
Code for reproducing the results in "A Pairwise Imputation Strategy for Retaining Predictive Features When Combining Multiple Datasets" by Wu et al.

## Reproducing the results in Section 3.1.1
Sim_sec31_simple.R contains the code for reproducing the simulation results in the first part of section 3.1.1, where the performance of using pairwise polynomial, merged polynomial, pairwise linear, and merged linear imputation model to impute the study-specific genes were compared, and all the genes in each study were predictive of the outcome with complete external validation set.

Sim_sec31_sparse_sig.R contains the the code for reproducing the simulation results in the second part of section 3.1.1, where additional "noise" genes were included.

## Reproducing the results in Section 3.1.2
Sim_sec32_MVN.R and Sim_sec32_sine_cosine.R contains the code reproducing the simulation results in section 3.1.2, where the validation set also had missing genes.

## Wilcoxon test on simulation results
Wilcoxon.R contains the code for performing and visualizing paired Wilcoxon test on the prediction RMSE between two methods.
