import pysam
import matplotlib.pyplot as plt
import sys
from collections import defaultdict


# Input VCF file
vcf_file = sys.argv[1]# "gene_variants.vcf"

depth_file = sys.argv[2]#"ind_cov.tsv"

gff_file = sys.argv[3] # "gene.gff"

# Output files
per_base_plot = sys.argv[4]#"normalized_coverage_plot.png"
per_feature_summary = sys.argv[5]#"feature_coverage.tsv"

############################################
# STEP 1: Read mean depth file
############################################
mean_depth = {}
with open(depth_file) as f:
    header = f.readline()  # skip header
    for line in f:
        fields = line.strip().split()
        sample = fields[0]
        mean_cov = float(fields[2])
        mean_depth[sample] = mean_cov

############################################
# STEP 2: Parse GFF and collect features
############################################
features = []
with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom = fields[0]
        feature_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        attributes = fields[8]
        feature_id = "NA"
        for attr in attributes.split(";"):
            if attr.startswith("ID="):
                feature_id = attr.split("=")[1]
                break
        features.append({'chrom': chrom, 'start': start, 'end': end,
                         'feature_type': feature_type, 'feature_id': feature_id})

############################################
# STEP 3: Parse VCF and collect data
############################################
vcf = pysam.VariantFile(vcf_file)
samples = list(vcf.header.samples)

# Storage for per-site plotting
sample_data = {sample: {'pos': [], 'scaled_dp': []} for sample in samples}
all_scaled_dp = []

# Storage for per-feature aggregation
agg_data = defaultdict(lambda: defaultdict(list))

for record in vcf:
    chrom = record.contig
    pos = record.pos

    # Per-sample per-site coverage
    for sample in samples:
        call = record.samples[sample]
        dp = call.get("DP", None)
        if dp is not None and sample in mean_depth:
            scaled_dp = dp / mean_depth[sample]
            sample_data[sample]['pos'].append(pos)
            sample_data[sample]['scaled_dp'].append(scaled_dp)
            all_scaled_dp.append(scaled_dp)

            # Check whether this site falls into any feature
            for feat in features:
                if chrom == feat['chrom'] and feat['start'] <= pos <= feat['end']:
                    feature_id = feat['feature_id']
                    feature_type = feat['feature_type']
                    agg_data[(feature_id, feature_type)][sample].append(scaled_dp)

############################################
# STEP 4: Per-base normalized coverage plot
############################################
n_samples = len(samples)
fig, axs = plt.subplots(n_samples, 1, figsize=(10, 3 * n_samples), sharex=True)

if n_samples == 1:
    axs = [axs]

y_min = 0
y_max = max(all_scaled_dp) + 0.1 if all_scaled_dp else 1.0

for i, sample in enumerate(samples):
    axs[i].scatter(sample_data[sample]['pos'], sample_data[sample]['scaled_dp'], s=10)
    axs[i].set_ylabel(f"{sample}\nScaled DP")
    axs[i].set_ylim(y_min, y_max)
    axs[i].grid(True)

axs[-1].set_xlabel("Genomic Position")
plt.tight_layout()
plt.savefig(per_base_plot, dpi=300)
plt.close()

print(f"Per-base plot saved to {per_base_plot}")

############################################
# STEP 5: Per-feature normalized coverage summary
############################################
with open(per_feature_summary, "w") as out:
    out.write("Sample\tFeature_ID\tFeature_Type\tMean_Scaled_DP\n")
    for (feature_id, feature_type), sample_dict in agg_data.items():
        for sample, dps in sample_dict.items():
            mean_scaled_dp = sum(dps) / len(dps)
            out.write(f"{sample}\t{feature_id}\t{feature_type}\t{mean_scaled_dp:.4f}\n")

print(f"Per-feature summary saved to {per_feature_summary}")

