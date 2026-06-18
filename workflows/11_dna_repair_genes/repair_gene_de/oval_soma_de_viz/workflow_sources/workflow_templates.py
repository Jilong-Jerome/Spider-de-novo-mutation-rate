import json
import shlex

from gwf import AnonymousTarget


def pathway_bias_bar_template(
    input_csv,
    out_dir,
    log_path,
    scripts_path,
    comparison,
    positive_label,
    negative_label,
    positive_color,
    negative_color,
    unbiased_color,
    pathway_order,
    account,
    log2fc_column="log2fc",
    log2fc_threshold=1.0,
    count_mode="log2fc_threshold",
    is_sig_column="is_sig",
    consistency_column=None,
    output_suffix="",
):
    done_file = f"{log_path}/{comparison}_pathway_bar{output_suffix}.DONE"
    pathway_order_str = ",".join(pathway_order)
    consistency_arg = f"    --consistency-column {consistency_column} \\\n" if consistency_column else ""

    inputs = {"csv": input_csv}
    outputs = {
        "tsv": f"{out_dir}/pathway_bias_counts{output_suffix}.tsv",
        "pdf": f"{out_dir}/pathway_bias_bar{output_suffix}.pdf",
        "png": f"{out_dir}/pathway_bias_bar{output_suffix}.png",
        "done": done_file,
    }
    options = {
        "cores": 1,
        "memory": "2g",
        "walltime": "00:30:00",
        "account": account,
    }

    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts_path}/plot_pathway_bias_bar.py \\
    --in {input_csv} \\
    --out-dir {out_dir} \\
    --comparison {comparison} \\
    --positive-label {positive_label} \\
    --negative-label {negative_label} \\
    --log2fc-column {log2fc_column} \\
    --log2fc-threshold {log2fc_threshold} \\
    --count-mode {count_mode} \\
    --is-sig-column {is_sig_column} \\
    --output-suffix "{output_suffix}" \\
    --positive-color "{positive_color}" \\
    --negative-color "{negative_color}" \\
    --unbiased-color "{unbiased_color}" \\
{consistency_arg}    --pathway-order "{pathway_order_str}"

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def combined_pathway_bias_bar_template(
    panels,
    out_dir,
    log_path,
    scripts_path,
    pathway_order,
    account,
    width=12.0,
    height=4.0,
    inconsistent_alpha=0.2,
    output_suffix="",
):
    done_file = f"{log_path}/combined_pathway_bar{output_suffix}.DONE"
    pathway_order_str = ",".join(pathway_order)
    panels_json = shlex.quote(json.dumps(panels))

    inputs = {f"tsv_{i}": p["tsv"] for i, p in enumerate(panels)}
    outputs = {
        "pdf": f"{out_dir}/combined_pathway_bias_bar{output_suffix}.pdf",
        "png": f"{out_dir}/combined_pathway_bias_bar{output_suffix}.png",
        "done": done_file,
    }
    options = {
        "cores": 1,
        "memory": "2g",
        "walltime": "00:30:00",
        "account": account,
    }

    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts_path}/combine_pathway_bars.py \\
    --pathway-order "{pathway_order_str}" \\
    --out-dir {out_dir} \\
    --width {width} --height {height} \\
    --inconsistent-alpha {inconsistent_alpha} \\
    --output-suffix "{output_suffix}" \\
    --panels {panels_json}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
