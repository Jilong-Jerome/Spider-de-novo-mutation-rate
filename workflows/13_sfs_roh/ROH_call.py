import vcf
import sys
from collections import defaultdict

def parse_fai_file(fai_file):
    """
    Parse the FAI file to get chromosome lengths.
    Returns a dictionary {chromosome: length}
    """
    chrom_lengths = {}
    with open(fai_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            chrom_lengths[parts[0]] = int(parts[1])
    return chrom_lengths

def calculate_window_heterozygosity(vcf_file, fai_file, detail_file, window_size=100000, step_size=10000):
    # Parse the FAI file
    chrom_lengths = parse_fai_file(fai_file)

    # Open the VCF file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    
    # Initialize variables for sliding window
    current_chrom = None
    current_window_start = 1
    current_window_end = window_size
    window_heterozygosity_count = defaultdict(int)  # {sample: heterozygosity_count}

    # Store detailed results
    detailed_results = []

    # Process each chromosome
    for record in vcf_reader:
        if record.CHROM != current_chrom:
            # Reset for new chromosome
            current_chrom = record.CHROM
            current_window_start = 1
            current_window_end = window_size
            window_heterozygosity_count.clear()

        while record.POS > current_window_end:
            # Record results for the current window and reset for new window
            for sample in vcf_reader.samples:
                count = window_heterozygosity_count[sample]
                detailed_results.append((sample, current_chrom, current_window_start, current_window_end, count))
            current_window_start += step_size
            current_window_end = current_window_start + window_size
            window_heterozygosity_count.clear()

        for sample in record.samples:
            # Count heterozygous genotypes within the window
            if sample.gt_type is not None and sample.gt_type == 1:
                window_heterozygosity_count[sample.sample] += 1

    # Process the last window of the last chromosome
    for sample in vcf_reader.samples:
        count = window_heterozygosity_count.get(sample, 0)
        detailed_results.append((sample, current_chrom, current_window_start, current_window_end, count))

    # Write detailed results to file
    with open(detail_file, "w") as det_file:
        det_file.write("Individual\tChromosome\tWindowStart\tWindowEnd\tHeterozygosityCount\n")
        for detail in detailed_results:
            det_file.write("\t".join(map(str, detail)) + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py <vcf_file> <fai_file> <detail_file> [window_size] [step_size]")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    fai_file = sys.argv[2]
    detail_file = sys.argv[3]
    window_size = int(sys.argv[4]) if len(sys.argv) > 4 else 100000
    step_size = int(sys.argv[5]) if len(sys.argv) > 5 else 10000

    calculate_window_heterozygosity(vcf_file, fai_file, detail_file, window_size, step_size)

