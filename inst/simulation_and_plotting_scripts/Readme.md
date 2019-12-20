# Simulation Pipelines

Example scripts for generating simulation replicates, creating and storing output files, and making plots of auRC and auPRC for method comparisons and class imbalance experiments.

Scripts:
1. compare-methods_auRC-auPRC(continuous).R
2. compare-methods_auRC-auPRC(discrete).R
3. make_auRC-auPRC_boxplots(method-comparison).R
4. generate_auRC-auPRC_replicates(imbalanced-data).R
5. make_auRC-auPRC_boxplots(imbalanced-hitmiss-nbds).R

Details:
1. Iteratively generates simulated continuous data sets with desired effects, runs several different feature selection methods, computes auRC and auPRC, stores auRC and auPRC for each method dataframe, and saves results.
2. Same as (1), but for GWAS data.
3. Uses data from (1) OR (2) and makes side-by-side boxplots of auRC and auPRC for each method and saves plot.
4. Generates 30 GWAS data sets for each level of imbalance in response (pct.imbalance=0.1, 0.2, 0.3, 0.4, 0.5), runs vwak-NPDR without balanced hit/miss neighborhoods (separate.hitmiss.nbds=F) OR runs vwak-NPDR with balanced hit/miss neighborhoods (separate.hitmiss.nbds=T), computes auRC and auPRC for each replicate data set, computes the best k for each attribute in each replicate data set, and saves auRC/auPRC and best-k matrix. Note: need to run twice: one with separate.hitmiss.nbds=T and the other with  separate.hitmiss.nbds=F.
5. Uses auRC/auPRC files from (4) to create box plots that compare balanced hit/miss vs imbalanced hit/miss NPDR for different levels of class imbalance (pct.imbalance=0.1, 0.2, 0.3, 0.4, 0.5). Saves plots.