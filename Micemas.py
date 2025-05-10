import pandas as pd
import numpy as np
from scipy import stats

# Load the datasets for Day 3 and Day 6
day3_file_path = "ML base code/Miceday03.xlsx"  
day6_file_path = "ML base code/Miceday06.xlsx" 

day3_data = pd.read_excel(day3_file_path)
day6_data = pd.read_excel(day6_file_path)

day3_expression = day3_data.filter(regex=r"^F1_\d{3}$")  
day6_expression = day6_data.filter(regex=r"^F1_\d{3}$|^F2_\d{3}$")  

day3_gene_ids = day3_data[['Row.names', 'ensembl_gene_id']]  
day6_gene_ids = day6_data[['Row.names', 'ensembl_gene_id']]  

merged_gene_ids = pd.concat([day3_gene_ids, day6_gene_ids]).drop_duplicates().reset_index(drop=True)


day3_expression_log = np.log1p(day3_expression)  
day6_expression_log = np.log1p(day6_expression)  

#Differential Expression Analysis (Welch's t-test)
p_values = []
logFC_values = []
gene_details = []  

for gene in day3_expression_log.index:  
    
    day3_vals = day3_expression_log.loc[gene].values
    day6_vals = day6_expression_log.loc[gene].values

    #Perform Welch's t-test
    t_stat, p_val = stats.ttest_ind(day3_vals, day6_vals, equal_var=False)
    
    logFC = np.log2(np.mean(day6_vals) / np.mean(day3_vals))

    p_values.append(p_val)
    logFC_values.append(logFC)

    gene_row = merged_gene_ids[merged_gene_ids['Row.names'] == gene]
    if not gene_row.empty:
        gene_details.append(gene_row.iloc[0].to_dict())
    else:
        gene_details.append({'Row.names': gene, 'ensembl_gene_id': 'Unknown'})

results = pd.DataFrame({
    'Row.names': [detail['Row.names'] for detail in gene_details],
    'ensembl_gene_id': [detail['ensembl_gene_id'] for detail in gene_details],
    'logFC': logFC_values,
    'p-value': p_values
})

print(results)

results.to_csv("gene_expression_p_values_with_ensembl_ids.csv", index=False)
