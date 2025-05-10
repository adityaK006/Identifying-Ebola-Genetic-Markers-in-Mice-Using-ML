# Here adjusted p_values is not calculated by Applying Benjamini-Hochberg (BH) correction 
# as normalised p values are used to calculate SMAS
# the thresholds are log2FC>=-0.2 and p_values<=0.3
#venv = base

import pandas as pd
import numpy as np

results_file_path = "gene_expression_p_values_with_ensembl_ids.csv"  
results = pd.read_csv(results_file_path)

print("Columns in the dataset:", results.columns)
print("Data types:", results.dtypes)

results['logFC'] = pd.to_numeric(results['logFC'], errors='coerce')
results['p-value'] = pd.to_numeric(results['p-value'], errors='coerce')

results = results.dropna(subset=['logFC', 'p-value'])

logFC_threshold = -0.2  
pvalue_threshold = 0.3  

filtered_results = results[
    (results['logFC'].abs() >= logFC_threshold) & 
    (results['p-value'] <= pvalue_threshold)
]

print(f"Number of genes passing thresholds: {filtered_results.shape[0]}")

# Supervised Magnitude-Altitude Scoring (MAS)
M = 1  # Example value for M
A = 1  # Example value for A

filtered_results['log2FC_abs_M'] = np.abs(filtered_results['logFC']) ** M
filtered_results['log10_p_abs_A'] = np.abs(np.log10(filtered_results['p-value'])) ** A
filtered_results['MAS'] = filtered_results['log2FC_abs_M'] * filtered_results['log10_p_abs_A']


print("Top 5 rows of filtered results after MAS calculation:")
print(filtered_results.head())


filtered_results_sorted = filtered_results.sort_values(by='MAS', ascending=False)
top_genes = filtered_results_sorted[['Row.names', 'ensembl_gene_id', 'logFC', 'p-value', 'MAS']].head(10)


# top_genes.to_csv("top_filtered_genes_with_MAS4.csv", index=False)
# print("Top genes with MAS:", top_genes)
