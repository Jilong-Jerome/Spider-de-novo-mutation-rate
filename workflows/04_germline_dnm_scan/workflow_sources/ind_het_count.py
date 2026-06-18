import sys
import pysam
import pandas as pd

# Define input parameters
vcf_file = sys.argv[1]#"your_file.vcf"  # Replace with the path to your VCF file
species = sys.argv[2]#"your_species"  # Replace with the species name
chromosome_id = sys.argv[3]#"your_chromosome_id"  # Replace with the chromosome ID
minDP = sys.argv[4]#10  # Replace with the minimum depth value (integer)

# Open VCF file
vcf = pysam.VariantFile(vcf_file)

# Initialize a dictionary to store counts
genotype_counts = {}

# Loop through each variant and each sample to count genotypes
for record in vcf:
    for sample in record.samples:
        genotype = record.samples[sample]['GT']  # Get the genotype tuple (e.g., (0, 0), (0, 1), (1, 1))
        # Initialize counts if sample is encountered for the first time
        if sample not in genotype_counts:
            genotype_counts[sample] = {"hom_ref": 0, "het": 0, "hom_alt": 0}
        
        # Increment appropriate genotype count
        if genotype == (0, 0):
            genotype_counts[sample]["hom_ref"] += 1
        elif genotype == (0, 1) or genotype == (1, 0):
            genotype_counts[sample]["het"] += 1
        elif genotype == (1, 1):
            genotype_counts[sample]["hom_alt"] += 1

# Convert the counts to a DataFrame for easy output
df = pd.DataFrame([
    {
        "species": species,
        "chromosome_id": chromosome_id,
        "minDP": minDP,
        "Individual ID": sample,
        "hom_ref": counts["hom_ref"],
        "het": counts["het"],
        "hom_alt": counts["hom_alt"]
    }
    for sample, counts in genotype_counts.items()
])

# Save to a tab-separated file
output_file = sys.argv[5]#"genotype_summary_with_metadata.tsv"
df.to_csv(output_file, sep="\t", index=False)

print(f"Genotype summary with metadata saved to {output_file}")

