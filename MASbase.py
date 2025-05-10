import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Step 1: Load the merged dataset
file_path = "ML base code/merged_dataset.csv"
data = pd.read_csv(file_path)

# Step 2: Extract relevant metadata columns
metadata_columns = ['ensembl_gene_id_day3', 'gene_symbol_day3', 'entrezgene_id_day3', 'uniprotswissprot_day3', 'description_day3']
data_day3_metadata = data[metadata_columns]

# Replace "day3" with "day6" for each element in the metadata_columns list
metadata_columns_day6 = [col.replace("day3", "day6") for col in metadata_columns]

# Extract Day 6 metadata columns
data_day6_metadata = data[metadata_columns_day6]

# Step 3: Extract expression columns for Day 3 and Day 6 samples
day3_expression = data.filter(regex="^F1_\d{3}$")  # Regex pattern to match Day 3 expression columns
day6_expression = data.filter(regex="^F1_\d{3}$|^F2_\d{3}$")  # Regex pattern to match Day 6 expression columns

# Check if data is present
print("Day 3 expression data shape:", day3_expression.shape)
print("Day 6 expression data shape:", day6_expression.shape)

# Step 4: Assign group labels (0 = Day 3, 1 = Day 6)
day3_labels = np.zeros(day3_expression.shape[1])
day6_labels = np.ones(day6_expression.shape[1])

# Combine expression data for both days and create labels
expression_data = pd.concat([day3_expression, day6_expression], axis=1)
labels = np.concatenate([day3_labels, day6_labels])

# Step 5: Perform differential expression analysis (t-test) with BH correction
p_values = []
logFC_values = []

# Iterate over the rows (genes) of the expression data
for gene in expression_data.index:
    day3_vals = day3_expression.loc[gene].dropna().values  # Remove NaN values
    day6_vals = day6_expression.loc[gene].dropna().values  # Remove NaN values

    # Skip if either of the values is empty or has zero variance
    if len(day3_vals) == 0 or len(day6_vals) == 0:
        continue  # Skip gene if there are missing values in either group

    # Perform variance check (skip if zero variance in either group)
    if np.var(day3_vals) == 0 or np.var(day6_vals) == 0:
        continue  # Skip genes with no variance
    
    # Perform t-test
    t_stat, p_val = stats.ttest_ind(day3_vals, day6_vals)
    
    # Calculate log fold change (logFC)
    logFC = np.log2(np.mean(day6_vals) / np.mean(day3_vals))
    
    p_values.append(p_val)
    logFC_values.append(logFC)

# Step 6: Apply BH correction (Benjamini-Hochberg) for multiple testing correction
_, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

# Step 7: Filter genes based on logFC > 1 and adjusted p-value < 0.05
selected_genes = pd.DataFrame({
    'Gene': data['gene_symbol_day3'].values,
    'logFC': logFC_values,
    'p-value': p_values,
    'adjusted p-value': p_values_corrected
})

# Check the distribution of logFC and adjusted p-value
print(selected_genes[['logFC','p-value', 'adjusted p-value']].describe())

## Step 7: Relaxed filter (logFC > 0.5 and adjusted p-value < 0.1)
selected_genes = selected_genes[(selected_genes['logFC'].abs() > -0.17) & (selected_genes['adjusted p-value'] < 0.89)]

# Step 8: Calculate Magnitude-Altitude Score (MAS)
M = 1  # Hyperparameter M
A = 1  # Hyperparameter A
selected_genes['MAS'] = np.abs(selected_genes['logFC']) * np.abs(np.log10(selected_genes['adjusted p-value'] + 1e-10)) ** A

# Step 9: Rank genes by MAS score
selected_genes = selected_genes.sort_values(by='MAS', ascending=False)

# Step 10: Output top d genes (you can set d to any value you want, e.g., top 10)
top_d_genes = selected_genes.head(10)

# Step 11: Save results to CSV (optional)
top_d_genes.to_csv("top_d_genes.csv", index=False)

# Output the top d genes with their MAS scores
print(top_d_genes)
