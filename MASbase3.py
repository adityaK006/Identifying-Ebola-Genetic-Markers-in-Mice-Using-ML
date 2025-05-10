
#random forest
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import GridSearchCV

# Step 1: Load all three datasets (adjust the paths as needed)
sheet1 = pd.read_excel('ML base code/aaq1016_Table S1.xlsx', sheet_name=None)  # Load all sheets in Table S1
sheet2 = pd.read_excel('ML base code/aaq1016_Table S2.xlsx')
sheet3 = pd.read_excel('ML base code/aaq1016_Table S3.xlsx')

# Step 2: Preprocess `sheet3` (assume it contains normalized expression counts)
gene_col_sheet3 = 'Gene'  # Adjust if the gene column name differs

# Step 3: Process `sheet1` and `sheet2` as before
fold_columns1 = []
pval_columns1 = []
processed_data1 = []

for name, df in sheet1.items():
    fold_columns1 = [col for col in df.columns if 'fold' in col.lower()]
    pval_columns1 = [col for col in df.columns if 'pval' in col.lower()]
    for fold_col, pval_col in zip(fold_columns1, pval_columns1):
        condition_data = {
            'Gene': df['Gene'],  # Adjust if 'Gene' column name differs
            'Fold_Change': df[fold_col],
            'P_Value': df[pval_col],
            'Condition': fold_col.split('_')[1]  # Adjust as needed
        }
        processed_data1.append(pd.DataFrame(condition_data))

fold_columns2 = [col for col in sheet2.columns if 'fold' in col.lower()]
pval_columns2 = [col for col in sheet2.columns if 'pval' in col.lower()]
processed_data2 = []

for fold_col, pval_col in zip(fold_columns2, pval_columns2):
    condition_data = {
        'Gene': sheet2['Gene'],  # Adjust if 'Gene' column name differs
        'Fold_Change': sheet2[fold_col],
        'P_Value': sheet2[pval_col],
        'Condition': fold_col.split('_')[1]  # Adjust as needed
    }
    processed_data2.append(pd.DataFrame(condition_data))

# Combine data from `sheet1` and `sheet2`
final_df1 = pd.concat(processed_data1, ignore_index=True)
final_df2 = pd.concat(processed_data2, ignore_index=True)
combined_df = pd.concat([final_df1, final_df2], ignore_index=True)

# Step 4: Integrate `sheet3` (expression data)
combined_with_expression = pd.merge(combined_df, sheet3, on='Gene', how='outer')

# Step 5: Filter genes based on selection criteria: p-value < 0.05 and logFC > 1
combined_with_expression['Log2_FC'] = np.log2(combined_with_expression['Fold_Change'].abs())
combined_with_expression['Log10_Pval'] = np.log10(combined_with_expression['P_Value'].abs())

# Step 6: Create a binary target variable based on p-value and Log2 fold change
combined_with_expression['Target'] = np.where(
    (combined_with_expression['Log2_FC'].abs() > 1) & (combined_with_expression['P_Value'] < 0.05),
    'Positive',  # Classifying as Positive if these conditions are met
    'Negative'   # Otherwise, Negative
)

# Step 7: Handle missing values
imputer = SimpleImputer(strategy='mean')  # You can change the strategy if needed
combined_with_expression[['Log2_FC', 'Log10_Pval']] = imputer.fit_transform(
    combined_with_expression[['Log2_FC', 'Log10_Pval']]
)

# Step 8: Feature scaling
scaler = StandardScaler()
X_scaled = scaler.fit_transform(combined_with_expression[['Log2_FC', 'Log10_Pval']])

# Step 9: Encode target labels (Positive/Negative) into 0/1 for classification
y = combined_with_expression['Target'].map({'Positive': 1, 'Negative': 0})

# Step 10: Split data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

# Step 11: Train a Random Forest Classifier
rf_classifier = RandomForestClassifier(random_state=42)

# Step 12: Hyperparameter tuning with Grid Search
param_grid = {
    'n_estimators': [100, 200],
    'max_depth': [None, 10, 20],
    'min_samples_split': [2, 5],
    'min_samples_leaf': [1, 2]
}

grid_search = GridSearchCV(rf_classifier, param_grid, cv=5, n_jobs=-1)
grid_search.fit(X_train, y_train)

# Step 13: Get the best model from Grid Search
best_rf_model = grid_search.best_estimator_

# Step 14: Make predictions and evaluate the model
y_pred = best_rf_model.predict(X_test)

# Step 15: Calculate accuracy
accuracy = accuracy_score(y_test, y_pred)
print(f"Random Forest Classifier Accuracy: {accuracy}")

# Optionally, use cross-validation to get a better estimate of the model's performance
cv_scores = cross_val_score(best_rf_model, X_scaled, y, cv=5)
print(f"Cross-validation accuracy: {cv_scores.mean()}")

print(combined_with_expression['Target'].value_counts())
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)
print("Confusion Matrix:")
print(cm)
