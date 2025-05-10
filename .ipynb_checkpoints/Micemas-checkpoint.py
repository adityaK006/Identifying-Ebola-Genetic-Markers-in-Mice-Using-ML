import pandas as pd
import numpy as np
from scipy import stats

# Load the datasets for Day 3 and Day 6
day3_file_path = "ML base code/Miceday03.xlsx"  # Update with your Day 3 file path
day6_file_path = "ML base code/Miceday06.xlsx"  # Update with your Day 6 file path

# Load Day 3 and Day 6 data
day3_data = pd.read_excel(day3_file_path)
day6_data = pd.read_excel(day6_file_path)

# Extract expression columns for Day 3 and Day 6
day3_expression = day3_data.filter(regex=r"^F1_\d{3}$")  # Day 3 columns (F1_XXX)
day6_expression = day6_data.filter(regex=r"^F1_\d{3}$|^F2_\d{3}$")  # Day 6 columns (F1_XXX, F2_XXX)

# Extract gene details (`ensembl_gene_id`) for both datasets
day3_gene_ids = day3_data[['Row.names', 'ensembl_gene_id']]  # Day 3 gene details
day6_gene_ids = day6_data[['Row.names', 'ensembl_gene_id']]  # Day 6 gene details

# Merge the gene details from Day 3 and Day 6 to ensure all genes are accounted for
merged_gene_ids = pd.concat([day3_gene_ids, day6_gene_ids]).drop_duplicates().reset_index(drop=True)

# Ensure log-transformation is applied to expression values
day3_expression_log = np.log1p(day3_expression)  # Apply log transformation to Day 3 data
day6_expression_log = np.log1p(day6_expression)  # Apply log transformation to Day 6 data

# Step 1: Differential Expression Analysis (Welch's t-test)
p_values = []
logFC_values = []
gene_details = []  # To store the gene details

# Step 2: Iterate through genes and compute t-test and log fold change
for gene in day3_expression_log.index:  # Ensure the index corresponds to gene identifiers
    # Fetch Day 3 and Day 6 values for the gene
    day3_vals = day3_expression_log.loc[gene].values
    day6_vals = day6_expression_log.loc[gene].values

    # Perform Welch's t-test
    t_stat, p_val = stats.ttest_ind(day3_vals, day6_vals, equal_var=False)

    # Calculate log fold change (logFC)
    logFC = np.log2(np.mean(day6_vals) / np.mean(day3_vals))

    # Append results
    p_values.append(p_val)
    logFC_values.append(logFC)

    # Fetch gene details from the merged dataset
    gene_row = merged_gene_ids[merged_gene_ids['Row.names'] == gene]
    if not gene_row.empty:
        gene_details.append(gene_row.iloc[0].to_dict())
    else:
        gene_details.append({'Row.names': gene, 'ensembl_gene_id': 'Unknown'})

# Step 3: Create a DataFrame with the calculated values and gene details
results = pd.DataFrame({
    'Row.names': [detail['Row.names'] for detail in gene_details],
    'ensembl_gene_id': [detail['ensembl_gene_id'] for detail in gene_details],
    'logFC': logFC_values,
    'p-value': p_values
})

# Step 4: Output the results (with gene details included)
print(results)

# Optionally, save the results to a CSV file
results.to_csv("gene_expression_p_values_with_ensembl_ids.csv", index=False)
