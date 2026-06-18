import sys
import pandas as pd

# Parameters
input_file = sys.argv[1] #'input.tsv'   # replace with your file name
target_individual = 'BIC_family3_F_female'
threshold = 0.05

# Read file using header
df = pd.read_csv(input_file, sep='\t')

# Make sure Mean_Scaled_DP is numeric
df['Mean_Scaled_DP'] = pd.to_numeric(df['Mean_Scaled_DP'], errors='coerce')

# Extract gene name from Feature_ID (before '.t')
df['Gene'] = df['Feature_ID'].str.split('.t').str[0]

# Process each gene
output_lines = []

for gene, gene_df in df.groupby('Gene'):
    df_exon = gene_df[gene_df['Feature_Type'] == 'exon']
    unique_absent_exons = []

    for feature_id in df_exon['Feature_ID'].unique():
        subset = df_exon[df_exon['Feature_ID'] == feature_id]
        target_row = subset[subset['Sample'] == target_individual]

        if target_row.empty:
            continue

        target_dp = target_row['Mean_Scaled_DP'].iloc[0]

        if target_dp < threshold:
            others = subset[subset['Sample'] != target_individual]
            if (others['Mean_Scaled_DP'] >= threshold).all():
                unique_absent_exons.append(feature_id)

    if unique_absent_exons:
        output_lines.append(f"{gene}\t1")
    else:
        output_lines.append(f"{gene}\t0")

# Print output as tab-separated lines
for line in output_lines:
    print(line)
