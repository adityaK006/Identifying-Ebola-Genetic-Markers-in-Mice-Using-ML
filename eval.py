# import pandas as pd
# import numpy as np
# from sklearn.metrics import roc_curve, auc
# import matplotlib.pyplot as plt

# # Load the results file
# results_file_path = "gene_expression_p_values_with_ensembl_ids.csv"  # Update with your file path
# results = pd.read_csv(results_file_path)

# # Debug: Print column names and data types
# print("Columns in the dataset:", results.columns)
# print("Data types:", results.dtypes)

# # Ensure numeric data for 'logFC' and 'p-value'
# results['logFC'] = pd.to_numeric(results['logFC'], errors='coerce')
# results['p-value'] = pd.to_numeric(results['p-value'], errors='coerce')

# # Drop rows with missing or invalid values in 'logFC' or 'p-value'
# results = results.dropna(subset=['logFC', 'p-value'])

# # Step 1: Apply thresholds for filtering
# logFC_threshold = -0.2  # Filter for substantial expression differences
# pvalue_threshold = 0.3  # Significant p-value threshold

# # Filter the results based on thresholds
# filtered_results = results[
#     (results['logFC'].abs() >= logFC_threshold) & 
#     (results['p-value'] <= pvalue_threshold)
# ]

# # Debug: Check the number of genes passing the thresholds
# print(f"Number of genes passing thresholds: {filtered_results.shape[0]}")

# # Step 2: Calculate the Supervised Magnitude-Altitude Scoring (MAS)
# M = 1  # Example value for M
# A = 1  # Example value for A

# filtered_results['log2FC_abs_M'] = np.abs(filtered_results['logFC']) ** M
# filtered_results['log10_p_abs_A'] = np.abs(np.log10(filtered_results['p-value'])) ** A
# filtered_results['MAS'] = filtered_results['log2FC_abs_M'] * filtered_results['log10_p_abs_A']

# # Debug: Check for issues in MAS calculation
# print("Top 5 rows of filtered results after MAS calculation:")
# print(filtered_results.head())

# # Step 3: Sort by MAS and output the top genes
# filtered_results_sorted = filtered_results.sort_values(by='MAS', ascending=False)
# top_genes = filtered_results_sorted[['Row.names', 'ensembl_gene_id', 'logFC', 'p-value', 'MAS']].head(10)

# # Save and print the top genes
# # top_genes.to_csv("top_filtered_genes_with_MAS4.csv", index=False)
# # print("Top genes with MAS:", top_genes)

# # Step 4: ROC Curve and AUC
# # Now compute true_labels for only the filtered results
# true_labels = (filtered_results['p-value'] <= 0.05).astype(int)  # True labels: 1 for significant genes, 0 for others

# # Generate the ROC curve
# fpr, tpr, thresholds = roc_curve(true_labels, filtered_results['MAS'])
# roc_auc = auc(fpr, tpr)

# # Plot the ROC curve
# plt.figure()
# plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
# plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver Operating Characteristic')
# plt.legend(loc="lower right")
# plt.show()

# # Print the AUC score
# print(f"AUC (Area Under Curve) for MAS: {roc_auc:.2f}")



# import pandas as pd
# import numpy as np
# from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, classification_report
# import matplotlib.pyplot as plt

# # Load the results file
# results_file_path = "gene_expression_p_values_with_ensembl_ids.csv"  # Update with your file path
# results = pd.read_csv(results_file_path)

# # Debug: Print column names and data types
# print("Columns in the dataset:", results.columns)
# print("Data types:", results.dtypes)

# # Ensure numeric data for 'logFC' and 'p-value'
# results['logFC'] = pd.to_numeric(results['logFC'], errors='coerce')
# results['p-value'] = pd.to_numeric(results['p-value'], errors='coerce')

# # Drop rows with missing or invalid values in 'logFC' or 'p-value'
# results = results.dropna(subset=['logFC', 'p-value'])

# # Step 1: Apply thresholds for filtering
# logFC_threshold = -0.2  # Filter for substantial expression differences
# pvalue_threshold = 0.3  # Significant p-value threshold

# # Filter the results based on thresholds
# filtered_results = results[
#     (results['logFC'].abs() >= logFC_threshold) & 
#     (results['p-value'] <= pvalue_threshold)
# ]

# # Debug: Check the number of genes passing the thresholds
# print(f"Number of genes passing thresholds: {filtered_results.shape[0]}")

# # Step 2: Calculate the Supervised Magnitude-Altitude Scoring (MAS)
# M = 1  # Example value for M
# A = 1  # Example value for A

# filtered_results['log2FC_abs_M'] = np.abs(filtered_results['logFC']) ** M
# filtered_results['log10_p_abs_A'] = np.abs(np.log10(filtered_results['p-value'])) ** A
# filtered_results['MAS'] = filtered_results['log2FC_abs_M'] * filtered_results['log10_p_abs_A']

# # Debug: Check for issues in MAS calculation
# print("Top 5 rows of filtered results after MAS calculation:")
# print(filtered_results.head())

# # Step 3: Sort by MAS and output the top genes
# filtered_results_sorted = filtered_results.sort_values(by='MAS', ascending=False)
# top_genes = filtered_results_sorted[['Row.names', 'ensembl_gene_id', 'logFC', 'p-value', 'MAS']].head(10)

# # Save and print the top genes
# # top_genes.to_csv("top_filtered_genes_with_MAS4.csv", index=False)
# # print("Top genes with MAS:", top_genes)

# # Step 4: ROC Curve and AUC
# # Now compute true_labels for only the filtered results
# true_labels = (filtered_results['p-value'] <= 0.05).astype(int)  # True labels: 1 for significant genes, 0 for others

# # Generate the ROC curve
# fpr, tpr, thresholds = roc_curve(true_labels, filtered_results['MAS'])
# roc_auc = auc(fpr, tpr)

# # Plot the ROC curve
# plt.figure()
# plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
# plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver Operating Characteristic')
# plt.legend(loc="lower right")
# plt.show()

# # Print the AUC score
# print(f"AUC (Area Under Curve) for MAS: {roc_auc:.2f}")

# # Step 5: Calculate other metrics
# # Predict class labels based on a threshold for MAS (e.g., if MAS > median, classify as 1, else 0)
# threshold_mas = filtered_results['MAS'].median()
# predicted_labels = (filtered_results['MAS'] > threshold_mas).astype(int)

# # Calculate accuracy, precision, recall, F1-score
# accuracy = accuracy_score(true_labels, predicted_labels)
# precision = precision_score(true_labels, predicted_labels)
# recall = recall_score(true_labels, predicted_labels)
# f1 = f1_score(true_labels, predicted_labels)

# # Print confusion matrix and classification report
# conf_matrix = confusion_matrix(true_labels, predicted_labels)
# class_report = classification_report(true_labels, predicted_labels)

# # Output the evaluation metrics
# print(f"Accuracy: {accuracy:.2f}")
# print(f"Precision: {precision:.2f}")
# print(f"Recall: {recall:.2f}")
# print(f"F1-Score: {f1:.2f}")
# print("Confusion Matrix:")
# print(conf_matrix)
# print("Classification Report:")
# print(class_report)







# import pandas as pd
# import numpy as np
# from sklearn.metrics import confusion_matrix, accuracy_score, classification_report

# # Load the results file
# results_file_path = "gene_expression_p_values_with_ensembl_ids.csv"  # Update with your file path
# results = pd.read_csv(results_file_path)

# # Ensure numeric data for 'logFC' and 'p-value'
# results['logFC'] = pd.to_numeric(results['logFC'], errors='coerce')
# results['p-value'] = pd.to_numeric(results['p-value'], errors='coerce')

# # Drop rows with missing or invalid values in 'logFC' or 'p-value'
# results = results.dropna(subset=['logFC', 'p-value'])

# # Add ground truth (create a 'target' column if not present)
# # Replace this with the actual ground truth data if available
# # Here, assume genes with 'logFC' > 0 and 'p-value' < 0.05 are significant
# results['target'] = ((results['logFC'] > 0) & (results['p-value'] < 0.05)).astype(int)

# # Apply thresholds for filtering
# logFC_threshold = -0.2  # Filter for substantial expression differences
# pvalue_threshold = 0.3  # Significant p-value threshold

# filtered_results = results[
#     (results['logFC'].abs() >= logFC_threshold) &
#     (results['p-value'] <= pvalue_threshold)
# ]

# # Calculate SMAS
# M = 1  # Example value for M
# A = 1  # Example value for A

# filtered_results['log2FC_abs_M'] = np.abs(filtered_results['logFC']) ** M
# filtered_results['log10_p_abs_A'] = np.abs(np.log10(filtered_results['p-value'])) ** A
# filtered_results['MAS'] = filtered_results['log2FC_abs_M'] * filtered_results['log10_p_abs_A']

# # Define MAS threshold (choose based on analysis or experimentation)
# mas_threshold = 1.0  # Example threshold (adjust as needed)

# # Predict significant genes based on MAS
# filtered_results['predicted'] = (filtered_results['MAS'] >= mas_threshold).astype(int)

# # Compare predictions to the ground truth
# y_true = filtered_results['target']
# y_pred = filtered_results['predicted']

# # Confusion Matrix
# conf_matrix = confusion_matrix(y_true, y_pred)
# print("Confusion Matrix:")
# print(conf_matrix)

# # Accuracy Score
# accuracy = accuracy_score(y_true, y_pred)
# print(f"Accuracy: {accuracy:.2f}")

# # Classification Report (Precision, Recall, F1-Score)
# class_report = classification_report(y_true, y_pred)
# print("Classification Report:")
# print(class_report)








import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
import matplotlib.pyplot as plt

# Load the results file
results_file_path = "gene_expression_p_values_with_ensembl_ids.csv"  # Update with your file path
results = pd.read_csv(results_file_path)

# Ensure numeric data for 'logFC' and 'p-value'
results['logFC'] = pd.to_numeric(results['logFC'], errors='coerce')
results['p-value'] = pd.to_numeric(results['p-value'], errors='coerce')

# Drop rows with missing or invalid values in 'logFC' or 'p-value'
results = results.dropna(subset=['logFC', 'p-value'])

# Step 1: Apply thresholds for filtering
logFC_threshold = -0.2  # Filter for substantial expression differences
pvalue_threshold = 0.3  # Significant p-value threshold

filtered_results = results[
    (results['logFC'].abs() >= logFC_threshold) & 
    (results['p-value'] <= pvalue_threshold)
]

# Step 2: Calculate SMAS (Supervised Magnitude-Altitude Scoring)
M, A = 1, 1  # Set example values for M and A
filtered_results['log2FC_abs_M'] = np.abs(filtered_results['logFC']) ** M
filtered_results['log10_p_abs_A'] = np.abs(np.log10(filtered_results['p-value'])) ** A
filtered_results['SMAS'] = filtered_results['log2FC_abs_M'] * filtered_results['log10_p_abs_A']

# Debug: Check the top genes based on SMAS scores
filtered_results_sorted = filtered_results.sort_values(by='SMAS', ascending=False)

# Step 3: Perform clustering using logFC and p-value as features
features = filtered_results[['logFC', 'p-value']].values
kmeans = KMeans(n_clusters=2, random_state=42)
clusters = kmeans.fit_predict(features)

# Add cluster labels to the DataFrame
filtered_results['cluster_label'] = clusters

# Step 4: Assign pseudo-ground-truth labels
# Assume one cluster corresponds to "significant" (1) and the other to "not significant" (0)
# Use cluster centroids to determine which cluster is significant
cluster_centroids = kmeans.cluster_centers_
significant_cluster = np.argmin(cluster_centroids[:, 1])  # Cluster with lower p-value
filtered_results['pseudo_target'] = (filtered_results['cluster_label'] == significant_cluster).astype(int)

# Step 5: Evaluate SMAS rankings against pseudo-target labels
# Rank genes by SMAS and classify top N as "significant"
N = len(filtered_results[filtered_results['pseudo_target'] == 1])  # Match number of true significant genes
filtered_results['SMAS_rank'] = filtered_results['SMAS'].rank(ascending=False)
filtered_results['SMAS_predicted'] = (filtered_results['SMAS_rank'] <= N).astype(int)

# Calculate accuracy
accuracy = accuracy_score(filtered_results['pseudo_target'], filtered_results['SMAS_predicted'])
print(f"Accuracy of SMAS model against pseudo-ground truth: {accuracy:.2f}")

# Confusion Matrix and Classification Report
conf_matrix = confusion_matrix(filtered_results['pseudo_target'], filtered_results['SMAS_predicted'])
print("Confusion Matrix:")
print(conf_matrix)
print("Classification Report:")
print(classification_report(filtered_results['pseudo_target'], filtered_results['SMAS_predicted']))

# Step 6: Plot clusters and SMAS rankings
plt.figure(figsize=(10, 6))
plt.scatter(
    filtered_results['logFC'], 
    -np.log10(filtered_results['p-value']), 
    c=filtered_results['cluster_label'], 
    cmap='coolwarm', 
    label='Clusters'
)
plt.colorbar(label="Cluster Label")
plt.title('Clustering of Genes (LogFC vs -log10(P-value))')
plt.xlabel('LogFC')
plt.ylabel('-log10(P-value)')
plt.show()
