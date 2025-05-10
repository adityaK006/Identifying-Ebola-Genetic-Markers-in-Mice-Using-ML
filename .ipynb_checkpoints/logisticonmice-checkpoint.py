# import pandas as pd
# import matplotlib.pyplot as plt
# from sklearn.model_selection import train_test_split
# from sklearn.linear_model import LogisticRegression
# from sklearn.metrics import confusion_matrix, classification_report, roc_curve, auc
# from sklearn.preprocessing import StandardScaler
# from sklearn.impute import SimpleImputer

# # Load your gene expression data (make sure to adjust the file path)
# gene_expression_data = pd.read_csv("gene_expression_p_values_with_ensembl_ids.csv")

# # Selecting relevant columns: gene_id, logFC, p-value
# gene_expression_data = gene_expression_data[['ensembl_gene_id', 'logFC', 'p-value']]

# # Handle missing values: drop rows with missing values (Option 1)
# # gene_expression_data.dropna(inplace=True)

# # Or, replace missing values with the column mean (Option 2: Imputation)
# imputer = SimpleImputer(strategy='mean')  # Impute missing values with the mean
# gene_expression_data[['logFC', 'p-value']] = imputer.fit_transform(gene_expression_data[['logFC', 'p-value']])

# # Create a binary target variable based on logFC (1 for activated genes, 0 for suppressed genes)
# gene_expression_data['target'] = gene_expression_data['logFC'].apply(lambda x: 1 if x > 0 else 0)

# # Drop the gene_id column as it's not used in the model
# X = gene_expression_data[['logFC', 'p-value']]
# y = gene_expression_data['target']

# # Split data into training and testing sets
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# # Standardize the data (important for models like Logistic Regression)
# scaler = StandardScaler()
# X_train_scaled = scaler.fit_transform(X_train)
# X_test_scaled = scaler.transform(X_test)

# # Initialize and train the Logistic Regression model
# model = LogisticRegression()
# model.fit(X_train_scaled, y_train)

# # Make predictions
# y_pred = model.predict(X_test_scaled)

# # Confusion Matrix
# conf_matrix = confusion_matrix(y_test, y_pred)

# # Classification Report
# class_report = classification_report(y_test, y_pred)

# # ROC Curve
# fpr, tpr, _ = roc_curve(y_test, model.predict_proba(X_test_scaled)[:, 1])
# roc_auc = auc(fpr, tpr)

# # Plot Confusion Matrix
# plt.figure(figsize=(6, 6))
# plt.matshow(conf_matrix, cmap='Blues', fignum=1)
# plt.title("Confusion Matrix")
# plt.colorbar()
# plt.xlabel("Predicted Label")
# plt.ylabel("True Label")
# plt.xticks([0, 1], ['Class 0 (Suppressed)', 'Class 1 (Activated)'])
# plt.yticks([0, 1], ['Class 0 (Suppressed)', 'Class 1 (Activated)'])
# plt.show()

# # Print classification report
# print("Classification Report:")
# print(class_report)

# # Plot ROC Curve
# plt.figure(figsize=(8, 6))
# plt.plot(fpr, tpr, color='darkorange', lw=2, label=f"ROC curve (area = {roc_auc:.2f})")
# plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver Operating Characteristic')
# plt.legend(loc='lower right')
# plt.show()

# # Plot scatter plot to visualize relationship between logFC and p-value
# plt.figure(figsize=(8, 6))
# plt.scatter(gene_expression_data['logFC'], gene_expression_data['p-value'], c=gene_expression_data['target'], cmap='coolwarm', alpha=0.7)
# plt.xlabel('LogFC')
# plt.ylabel('P-value')
# plt.title('Scatter plot of LogFC vs P-value with Classification')
# plt.colorbar(label='Target (1=Activated, 0=Suppressed)')
# plt.show()






# import pandas as pd
# import matplotlib.pyplot as plt
# from sklearn.model_selection import train_test_split
# from sklearn.linear_model import LogisticRegression
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.svm import SVC
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn.metrics import confusion_matrix, classification_report, roc_curve, auc
# from sklearn.preprocessing import StandardScaler
# from sklearn.impute import SimpleImputer
# from sklearn.model_selection import learning_curve

# # Load your gene expression data (make sure to adjust the file path)
# gene_expression_data = pd.read_csv("gene_expression_p_values_with_ensembl_ids.csv")

# # Selecting relevant columns: gene_id, logFC, p-value
# gene_expression_data = gene_expression_data[['ensembl_gene_id', 'logFC', 'p-value']]

# # Handle missing values: Imputation (replacing missing values with the mean)
# imputer = SimpleImputer(strategy='mean')
# gene_expression_data[['logFC', 'p-value']] = imputer.fit_transform(gene_expression_data[['logFC', 'p-value']])

# # Create a binary target variable based on logFC (1 for activated genes, 0 for suppressed genes)
# gene_expression_data['target'] = gene_expression_data['logFC'].apply(lambda x: 1 if x > 0 else 0)

# # Drop the gene_id column as it's not used in the model
# X = gene_expression_data[['logFC', 'p-value']]
# y = gene_expression_data['target']

# # Split data into training and testing sets
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# # Standardize the data (important for models like Logistic Regression, SVM, etc.)
# scaler = StandardScaler()
# X_train_scaled = scaler.fit_transform(X_train)
# X_test_scaled = scaler.transform(X_test)

# # Models to try
# models = {
#     "Logistic Regression": LogisticRegression(),
#     "Random Forest": RandomForestClassifier(),
#     "Support Vector Machine": SVC(probability=True),
#     "K-Nearest Neighbors": KNeighborsClassifier()
# }

# # Function to evaluate models
# def evaluate_model(model, X_train, y_train, X_test, y_test):
#     # Train the model
#     model.fit(X_train, y_train)
    
#     # Make predictions
#     y_pred = model.predict(X_test)
    
#     # Calculate metrics
#     conf_matrix = confusion_matrix(y_test, y_pred)
#     class_report = classification_report(y_test, y_pred)
    
#     # ROC Curve
#     fpr, tpr, _ = roc_curve(y_test, model.predict_proba(X_test)[:, 1])
#     roc_auc = auc(fpr, tpr)
    
#     return conf_matrix, class_report, fpr, tpr, roc_auc

# # Plot results for each model
# for model_name, model in models.items():
#     print(f"Evaluating {model_name}...")
    
#     # Evaluate the model
#     conf_matrix, class_report, fpr, tpr, roc_auc = evaluate_model(model, X_train_scaled, y_train, X_test_scaled, y_test)
    
#     # Plot Confusion Matrix
#     plt.figure(figsize=(6, 6))
#     plt.matshow(conf_matrix, cmap='Blues', fignum=1)
#     plt.title(f"Confusion Matrix for {model_name}")
#     plt.colorbar()
#     plt.xlabel("Predicted Label")
#     plt.ylabel("True Label")
#     plt.xticks([0, 1], ['Class 0 (Suppressed)', 'Class 1 (Activated)'])
#     plt.yticks([0, 1], ['Class 0 (Suppressed)', 'Class 1 (Activated)'])
#     plt.show()

#     # Print Classification Report
#     print(f"Classification Report for {model_name}:")
#     print(class_report)
    
#     # Plot ROC Curve
#     plt.figure(figsize=(8, 6))
#     plt.plot(fpr, tpr, color='darkorange', lw=2, label=f"ROC curve (area = {roc_auc:.2f})")
#     plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
#     plt.xlim([0.0, 1.0])
#     plt.ylim([0.0, 1.05])
#     plt.xlabel('False Positive Rate')
#     plt.ylabel('True Positive Rate')
#     plt.title(f'ROC Curve for {model_name}')
#     plt.legend(loc='lower right')
#     plt.show()

#     # Plot Learning Curves for each model (train and test score as a function of the training size)
#     train_sizes, train_scores, test_scores = learning_curve(model, X_train_scaled, y_train, cv=5)
#     plt.figure(figsize=(8, 6))
#     plt.plot(train_sizes, train_scores.mean(axis=1), label=f"Training score ({model_name})")
#     plt.plot(train_sizes, test_scores.mean(axis=1), label=f"Test score ({model_name})")
#     plt.xlabel('Training Set Size')
#     plt.ylabel('Score')
#     plt.title(f"Learning Curves for {model_name}")
#     plt.legend()
#     plt.show()

# # For Random Forest, plot feature importance
# rf_model = RandomForestClassifier()
# rf_model.fit(X_train_scaled, y_train)
# importances = rf_model.feature_importances_

# # Plot Feature Importance
# plt.figure(figsize=(8, 6))
# plt.barh(X.columns, importances, color='steelblue')
# plt.xlabel('Feature Importance')
# plt.title('Feature Importance (Random Forest)')
# plt.show()


# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# from sklearn.model_selection import train_test_split
# from sklearn.linear_model import LogisticRegression
# from sklearn.preprocessing import StandardScaler

# # Load your gene expression data (adjust the file path as needed)
# gene_expression_data = pd.read_csv("gene_expression_p_values_with_ensembl_ids.csv")

# # Select relevant columns: gene_id, logFC, p-value
# gene_expression_data = gene_expression_data[['ensembl_gene_id', 'logFC', 'p-value']]

# # Drop rows with missing values
# gene_expression_data = gene_expression_data.dropna()

# # Create a binary target variable based on logFC (logFC > 0 = Activated, else Suppressed)
# gene_expression_data['target'] = gene_expression_data['logFC'].apply(lambda x: 1 if x > 0 else 0)

# # Select features and target variable
# X = gene_expression_data[['logFC', 'p-value']]
# y = gene_expression_data['target']

# # Split data into training and testing sets
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# # Standardize the features
# scaler = StandardScaler()
# X_train_scaled = scaler.fit_transform(X_train)
# X_test_scaled = scaler.transform(X_test)

# # Train a Logistic Regression model
# model = LogisticRegression()
# model.fit(X_train_scaled, y_train)

# # Create a meshgrid for plotting decision boundaries
# x_min, x_max = X_train_scaled[:, 0].min() - 1, X_train_scaled[:, 0].max() + 1
# y_min, y_max = X_train_scaled[:, 1].min() - 1, X_train_scaled[:, 1].max() + 1
# xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.01),
#                      np.arange(y_min, y_max, 0.01))

# # Use the model to predict on the grid points
# Z = model.predict(np.c_[xx.ravel(), yy.ravel()])
# Z = Z.reshape(xx.shape)

# # Plot the decision boundary
# plt.figure(figsize=(10, 6))
# plt.contourf(xx, yy, Z, alpha=0.8, cmap="coolwarm")
# sns.scatterplot(x=X_train_scaled[:, 0], y=X_train_scaled[:, 1], 
#                 hue=y_train, palette="coolwarm", edgecolor="k", s=40)

# # Add labels and title
# plt.title("Decision Boundary of Logistic Regression")
# plt.xlabel("Standardized logFC")
# plt.ylabel("Standardized p-value")
# plt.legend(title="Gene Class", labels=["Suppressed (0)", "Activated (1)"])
# plt.show()


# # Import necessary libraries
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns

# # Step 1: Load your gene expression data (adjust the file path as needed)
# gene_expression_data = pd.read_csv("gene_expression_p_values_with_ensembl_ids.csv")

# # Step 2: Select relevant columns: ensembl_gene_id, logFC, p-value
# gene_expression_data = gene_expression_data[['ensembl_gene_id', 'logFC', 'p-value']]

# # Step 3: Drop rows with missing values
# gene_expression_data = gene_expression_data.dropna()

# # Step 4: Create a binary target variable based on logFC
# # logFC > 0 = Activated (1), logFC <= 0 = Suppressed (0)
# gene_expression_data['target'] = gene_expression_data['logFC'].apply(lambda x: 1 if x > 0 else 0)

# # Step 5: Compute the correlation matrix for selected columns
# correlation_matrix = gene_expression_data[['logFC', 'p-value', 'target']].corr()

# # Step 6: Plot the heatmap
# plt.figure(figsize=(8, 6))
# sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", fmt=".2f", cbar=True)
# plt.title("Correlation Heatmap: logFC, p-value, and Target")
# plt.show()

# # Optional: Print the correlation matrix for verification
# print("Correlation Matrix:")
# print(correlation_matrix)


# # Plot pairwise relationships between logFC, p-value, and target (colored by target)
# sns.pairplot(gene_expression_data[['logFC', 'p-value', 'target']], hue='target', palette='coolwarm')
# plt.suptitle("Pairplot of logFC, p-value, and Target", y=1.02)
# plt.show()







import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, classification_report, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

# Load your gene expression data (adjust the file path as needed)
gene_expression_data = pd.read_csv("gene_expression_p_values_with_ensembl_ids.csv")

# Selecting relevant columns: gene_id, logFC, p-value
gene_expression_data = gene_expression_data[['ensembl_gene_id', 'logFC', 'p-value']]

# Handle missing values: Impute missing values with the column mean
imputer = SimpleImputer(strategy='mean')
gene_expression_data[['logFC', 'p-value']] = imputer.fit_transform(gene_expression_data[['logFC', 'p-value']])

# Create a binary target variable based on logFC (1 for activated genes, 0 for suppressed genes)
gene_expression_data['target'] = gene_expression_data['logFC'].apply(lambda x: 1 if x > 0 else 0)

# Drop the gene_id column as it's not used in the model
X = gene_expression_data[['logFC', 'p-value']]
y = gene_expression_data['target']

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Standardize the data (important for models like Logistic Regression)
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Initialize and train the Logistic Regression model
model = LogisticRegression()
model.fit(X_train_scaled, y_train)

# Make predictions
y_pred = model.predict(X_test_scaled)

# Confusion Matrix
conf_matrix = confusion_matrix(y_test, y_pred)

# Classification Report
class_report = classification_report(y_test, y_pred)

# Compute AUC Score
auc_score = roc_auc_score(y_test, model.predict_proba(X_test_scaled)[:, 1])

# Print results
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
