#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 vcf_list_file individuals_list_file output_vcf_file"
    exit 1
fi

# Read arguments
vcf_list_file=$1
individuals_list_file=$2
output_vcf_file=$3

# Read VCF files into an array
readarray -t vcf_files < "$vcf_list_file"

# Check and re-compress with bgzip if needed, then index with tabix
for vcf in "${vcf_files[@]}"; do
    if file "$vcf" | grep -q 'gzip compressed data'; then
        echo "$vcf is gzipped, recompressing with bgzip..."
        zcat "$vcf" | bgzip > "${vcf}.bgz"
        mv "${vcf}.bgz" "$vcf"
    fi
    echo "Indexing $vcf..."
    tabix -p vcf "$vcf"
done

# Join array elements into a string separated by space
vcf_files_joined=$(IFS=" "; echo "${vcf_files[*]}")

# Read individuals into an array
readarray -t individuals < "$individuals_list_file"

# Join array elements into a string separated by comma
individuals_joined=$(IFS=","; echo "${individuals[*]}")

# Merge VCF files
bcftools merge $vcf_files_joined -Oz -o temp_merged.vcf.gz

# Filter for specified individuals
bcftools view -Oz -s "$individuals_joined" temp_merged.vcf.gz -o temp_filtered_individuals.vcf.gz

# Filter for SNP variable sites
bcftools view -Oz -g het -v snps temp_filtered_individuals.vcf.gz -o "$output_vcf_file"

# Clean up intermediate files
rm temp_merged.vcf.gz temp_filtered_individuals.vcf.gz

