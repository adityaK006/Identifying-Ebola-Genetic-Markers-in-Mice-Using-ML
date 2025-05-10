


import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, classification_report, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

gene_expression_data = pd.read_csv("gene_expression_p_values_with_ensembl_ids.csv")

gene_expression_data = gene_expression_data[['ensembl_gene_id', 'logFC', 'p-value']]

imputer = SimpleImputer(strategy='mean')
gene_expression_data[['logFC', 'p-value']] = imputer.fit_transform(gene_expression_data[['logFC', 'p-value']])

gene_expression_data['target'] = gene_expression_data['logFC'].apply(lambda x: 1 if x > 0 else 0)

X = gene_expression_data[['logFC', 'p-value']]
y = gene_expression_data['target']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

model = LogisticRegression()
model.fit(X_train_scaled, y_train)

y_pred = model.predict(X_test_scaled)

conf_matrix = confusion_matrix(y_test, y_pred)

class_report = classification_report(y_test, y_pred)

auc_score = roc_auc_score(y_test, model.predict_proba(X_test_scaled)[:, 1])


print("Classification Report:")
print(class_report)
print(f"AUC Score: {auc_score:.2f}")

# Plot Confusion Matrix
plt.figure(figsize=(6, 6))
plt.matshow(conf_matrix, cmap='Blues', fignum=1)
plt.title("Confusion Matrix")
plt.colorbar()
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.xticks([0, 1], ['Class 0 (Suppressed)', 'Class 1 (Activated)'])
plt.yticks([0, 1], ['Class 0 (Suppressed)', 'Class 1 (Activated)'])
plt.show()

# Visualize the AUC score
plt.figure(figsize=(6, 6))
plt.text(0.5, 0.5, f"AUC = {auc_score:.2f}", fontsize=20, ha='center', va='center', bbox=dict(facecolor='lightblue', alpha=0.5))
plt.axis('off')
plt.title("Area Under the Curve (AUC)")
plt.show()

# Plot scatter plot to visualize relationship between logFC and p-value
plt.figure(figsize=(8, 6))
plt.scatter(gene_expression_data['logFC'], gene_expression_data['p-value'], c=gene_expression_data['target'], cmap='coolwarm', alpha=0.7)
plt.xlabel('LogFC')
plt.ylabel('P-value')
plt.title('Scatter plot of LogFC vs P-value with Classification')
plt.colorbar(label='Target (1=Activated, 0=Suppressed)')
plt.show()
