import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.impute import SimpleImputer

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
filtered_df = combined_with_expression[
    (combined_with_expression['P_Value'] < 0.05) &
    (np.log2(combined_with_expression['Fold_Change'].abs()) > 1)
]

# Step 6: Calculate MAS scores for filtered genes
M = 1  # Hyperparameter for scaling log2FC
A = 1  # Hyperparameter for scaling log10(p-value)

filtered_df['Log2_FC'] = np.log2(filtered_df['Fold_Change'].abs())
filtered_df['Log10_Pval'] = np.log10(filtered_df['P_Value'].abs())

filtered_df['MAS'] = (
    filtered_df['Log2_FC'] * M *
    filtered_df['Log10_Pval'].abs() * A
)

# Step 7: Sort and display the top genes by MAS score
top_genes_df = filtered_df.sort_values(by='MAS', ascending=False)

print("Top genes selected by MAS scoring for classification:")
print(top_genes_df.head(10)[['Gene', 'Fold_Change', 'P_Value', 'MAS', 'Condition']])

# Optionally, save the top genes to a file
top_genes_df.to_csv('top_genes_MAS.csv', index=False)

# Check for a specific gene (e.g., 'IFI27')
print(top_genes_df[top_genes_df['Gene'] == 'IFI27'])

# Step 8: Handle missing values (either impute or drop)
# Option 1: Impute missing values using SimpleImputer
imputer = SimpleImputer(strategy='mean')

# Impute missing values for features X
X_imputed = imputer.fit_transform(filtered_df[['Log2_FC', 'Log10_Pval']])

# Impute missing values for the target variable y (if necessary)
y_imputed = SimpleImputer(strategy='most_frequent').fit_transform(filtered_df['Condition'].values.reshape(-1, 1))

# Option 2: Alternatively, drop rows with missing values (if you prefer)
# Drop rows with NaN values in the feature set or target
# filtered_df.dropna(subset=['Log2_FC', 'Log10_Pval', 'Condition'], inplace=True)

# Now, split the data into features (X) and target (y)
X = X_imputed  # Features after imputation
y = y_imputed.ravel()  # Flatten y for classification (if imputed)

# Step 9: Split into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Step 10: Train the logistic regression model
classifier = LogisticRegression(max_iter=10000)  # Increase max_iter if convergence is an issue
classifier.fit(X_train, y_train)

# Step 11: Make predictions and calculate accuracy
y_pred = classifier.predict(X_test)

# Step 12: Calculate accuracy
accuracy = accuracy_score(y_test, y_pred)
print(f'Accuracy of the classifier: {accuracy}')
