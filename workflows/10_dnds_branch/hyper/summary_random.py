import pandas as pd
import sys
# Parameters
file_list_path = sys.argv[1]#'input_file_list.txt'  # Your file list file
species = sys.argv[2]
output_file = sys.argv[3]#'output.tsv'              # Your final output file
target_individual = sys.argv[4]#'BIC_family3_F_female'
threshold = 0.1

# Open output file for writing
with open(output_file, 'w') as out_f:
    out_f.write("Gene\tUniqueAbsence\n")  # write header

    # Read file list
    with open(file_list_path, 'r') as list_f:
        for line in list_f:
            gene_id = line.strip("\n")
            input_file = f"/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/steps/hyper/random/vcfs/{species}/{gene_id}/{species}_{gene_id}_summary.tsv"
            if not input_file:
                continue  # skip empty lines

            # Read one input file
            df = pd.read_csv(input_file, sep='\t')
            df['Mean_Scaled_DP'] = pd.to_numeric(df['Mean_Scaled_DP'], errors='coerce')

            # Extract gene name from Feature_ID
            df['Gene'] = df['Feature_ID'].str.split('.t').str[0]

            # There should be only one gene in this file
            gene_name = df['Gene'].unique()[0]

            df_exon = df[df['Feature_Type'] == 'exon']
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
                out_f.write(f"{gene_name}\t1\n")
            else:
                out_f.write(f"{gene_name}\t0\n")

print("Batch processing completed!")

