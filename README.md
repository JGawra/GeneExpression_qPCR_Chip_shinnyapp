About: Gene Expression Analysis Shiny App
What is this?
This program is an interactive platform for gene expression data analysis. Aslo can be used for any different type of data that have sample ID in the 1st column and Treatment group in the 2nd column and 3rd, 4th and so on different measurment for numeric data. You can perform outlier detection, imputation, statistical testing (ANOVA, Kruskal-Wallis, post-hoc tests), PCA, and generate publication-ready plots and summaries.
You can access to the app here:

https://janan91.shinyapps.io/GeneExpression_qPCR_Chip_shinnyapp_V11/


Input File Format
Upload your data as a CSV or TSV file with the following structure:

First column: Sample ID
Second column: Treatment or group (categorical variable)
Remaining columns: Gene expression values (one column per gene, numeric)
Missing values are allowed and will be imputed as specified. Treatments should be named consistently.

Analyses Available
Outlier detection and replacement (per group or globally)
Imputation of missing values (mean, median, or zero)
Normality testing (Shapiro-Wilk)
Significant gene identification (ANOVA, Kruskal-Wallis, automatic based on normality, with optional FDR correction)
Post-hoc pairwise comparisons (Tukey, Dunn, t-test, Wilcoxon)
Principal Component Analysis (PCA) and biplots
Heatmaps (all genes or significant genes only)
Boxplots, violin plots, and custom qPCR-style plots
Outputs
Interactive plots in the browser
Downloadable plots (PNG, PDF for all genes)
Excel summary including normality, statistical test results, and significant gene lists
Tables of outliers, imputed values, and significant genes
Contact
Developed by Janan Gawra. For help, contact via LinkedIn.
https://www.linkedin.com/in/janangawra/
