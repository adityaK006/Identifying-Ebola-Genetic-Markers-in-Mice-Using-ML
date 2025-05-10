import pandas as pd

# Load Day 3 and Day 6 data from Excel files
day3_file = "ML base code/Miceday03.xlsx"  # Replace with the actual path to your Day 3 file
day6_file = "ML base code/Miceday06.xlsx"  # Replace with the actual path to your Day 6 file

# Assuming the data is in the first sheet of the Excel files
day3_data = pd.read_excel(day3_file)
day6_data = pd.read_excel(day6_file)

# Display the first few rows to identify the column names
print("Day 3 Data Columns:", day3_data.columns)
print("Day 6 Data Columns:", day6_data.columns)

# Merge based on a common column, e.g., 'ensembl_gene_id' or 'gene_symbol'
# Replace 'common_column_name' with the actual column used for merging
merged_data = pd.merge(day3_data, day6_data, on='ensembl_transcript_id', suffixes=('_day3', '_day6'))

# Display the shape and head of the merged dataframe
print("Merged Data Shape:", merged_data.shape)
print(merged_data.head())

# Save the merged data to a new Excel file
output_file = "merged_data.xlsx"
merged_data.to_excel(output_file, index=False)
print(f"Merged data saved to {output_file}")
# Save the merged dataset to a CSV file in the current directory
merged_data.to_csv("merged_dataset.csv", index=False)

print("Merged dataset saved as 'merged_dataset.csv' in the current directory.")