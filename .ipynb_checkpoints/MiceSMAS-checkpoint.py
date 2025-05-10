#Here adjusted p_values is calculated by Applying Benjamini-Hochberg (BH) correction for p-values
#and no  log2FC,pvalue criteria is used
# directly MAS score is calculated and top genes are found.

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

# Load the results file
results_file_path = "gene_expression_p_values_with_ensembl_ids.csv"  # Update with your file path
results = pd.read_csv(results_file_path)

# Debug: Print column names and data types
print("Columns in the dataset:", results.columns)
print("Data types:", results.dtypes)

# Ensure numeric data for 'logFC' and 'p-value'
results['logFC'] = pd.to_numeric(results['logFC'], errors='coerce')
results['p-value'] = pd.to_numeric(results['p-value'], errors='coerce')

# Drop rows with missing or invalid values in 'logFC' or 'p-value'
results = results.dropna(subset=['logFC', 'p-value'])

# Step 1: Apply Benjamini-Hochberg (BH) correction for p-values
_, p_values_corrected, _, _ = multipletests(results['p-value'], method='fdr_bh')

# Step 2: Add the corrected p-values to the results DataFrame
results['adjusted p-value'] = p_values_corrected

# Step 3: Set the hyperparameters M and A
M = 1  # Example value for M
A = 1  # Example value for A

# Step 4: Calculate the Supervised Magnitude-Altitude Scoring (MAS)
results['log2FC_abs_M'] = np.abs(results['logFC']) ** M
results['log10_pBH_abs_A'] = np.abs(np.log10(results['adjusted p-value'])) ** A
results['MAS'] = results['log2FC_abs_M'] * results['log10_pBH_abs_A']

# Debug: Check for issues in MAS calculation
print("Top 5 rows of results after MAS calculation:")
print(results.head())

# Step 5: Sort by MAS and output the top genes
results_sorted = results.sort_values(by='MAS', ascending=False)
top_genes = results_sorted[['Row.names', 'ensembl_gene_id', 'logFC', 'p-value', 'adjusted p-value', 'MAS']].head(10)

# Save and print the top genes
top_genes.to_csv("top_genes_with_MAS.csv", index=False)
print("Top genes with MAS:", top_genes)
