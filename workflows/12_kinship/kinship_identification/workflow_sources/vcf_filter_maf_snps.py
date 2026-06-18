import sys
import pysam

def is_snp(ref, alt):
    """Check if the variant is a SNP (single nucleotide polymorphism)."""
    return len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT"

def extract_af(info):
    """Extract the allele frequency (AF) from the INFO field."""
    if 'AF' in info:
        try:
            return float(info['AF'][0])  # Assuming it's a list and we take the first allele frequency
        except (ValueError, IndexError):
            return None
    return None

def is_autosome(chrom):
    """Check if the chromosome is an autosome (second part of chromosome ID is numeric)."""
    parts = chrom.split('_')
    if len(parts) == 2 and parts[1].isdigit():
        return True
    return False

def calculate_maf(af):
    """Calculate Minor Allele Frequency (MAF) given allele frequency AF."""
    return min(af, 1.0 - af)

def filter_vcf(input_vcf, output_vcf, maf_threshold=0.1):
    """Filter autosomal SNPs from input VCF and save to output VCF.

    Keeps biallelic autosomal SNPs. If maf_threshold >= 0, additionally requires
    MAF > maf_threshold (population MAF read from INFO/AF). If maf_threshold < 0
    the MAF gate is skipped entirely (no-MAF mode) and every biallelic autosomal
    SNP is kept regardless of allele frequency -- used by the trio Mendelian
    check, where the informative-site filter is the per-trio opposite-homozygote
    condition, not a population frequency cutoff.
    """
    # Open input VCF file using pysam
    vcf_in = pysam.VariantFile(input_vcf, 'r')

    # Open output VCF file
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

    skip_maf = maf_threshold < 0

    # Iterate through variants
    for record in vcf_in:
        ref = record.ref
        alts = record.alts
        chrom = record.chrom

        # Check if it's a SNP and if the chromosome is autosomal
        if alts and all(is_snp(ref, alt) for alt in alts) and is_autosome(chrom):
            if skip_maf:
                # No-MAF mode: keep every biallelic autosomal SNP.
                vcf_out.write(record)
                continue
            af = extract_af(record.info)
            if af is not None:
                maf = calculate_maf(af)
                if maf > maf_threshold:
                    # Write autosomal SNPs with MAF > maf_threshold to the new VCF file
                    vcf_out.write(record)

    # Close files
    vcf_in.close()
    vcf_out.close()

# Example usage
input_vcf = sys.argv[1]#"input.vcf"
output_vcf = sys.argv[2]#"filtered_autosomal_snps.vcf"
# Optional 3rd arg overrides the MAF threshold; a negative value disables the
# MAF gate (no-MAF mode). Default preserves the original MAF > 0.1 behaviour.
maf_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
filter_vcf(input_vcf, output_vcf, maf_threshold=maf_threshold)

