#!/bin/bash
# vcf_summary.sh
# Usage: ./vcf_summary.sh <vcf_file> <species_code> <individual_id> <chromosome_id> <output_file>

# Check for exactly 5 arguments.
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <vcf_file> <species_code> <individual_id> <chromosome_id> <output_file>"
    exit 1
fi

VCF_FILE="$1"
SPECIES="$2"
INDIVIDUAL="$3"
CHROM="$4"
OUTPUT_FILE="$5"

# Check if the VCF file exists.
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: VCF file '$VCF_FILE' not found."
    exit 1
fi

# Process the VCF file:
# - Only consider lines for the specified chromosome (column 1).
# - Count the reference alleles (column 4) that are A, T, C, or G.
awk -v species="$SPECIES" -v individual="$INDIVIDUAL" -v chrom_id="$CHROM" '
BEGIN {
    # Initialize counts for A, T, C, and G.
    A=0; T=0; C=0; G=0;
}
# Ignore header lines and process only lines with the specified chromosome.
!/^#/ && $1 == chrom_id {
    if ($4=="A") A++;
    else if ($4=="T") T++;
    else if ($4=="C") C++;
    else if ($4=="G") G++;
}
END {
    # Output header and counts.
    printf "species\tindividual\tchrom\tA\tT\tC\tG\n%s\t%s\t%s\t%d\t%d\t%d\t%d\n", species, individual, chrom_id, A, T, C, G;
}
' "$VCF_FILE" > "$OUTPUT_FILE"

echo "Summary written to $OUTPUT_FILE"

