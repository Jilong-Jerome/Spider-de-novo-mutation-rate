import os
import yaml
from workflow_templates import (
    pathway_bias_bar_template,
    combined_pathway_bias_bar_template,
)


def social_subsocial_count_modes(config):
    if "count_modes" in config:
        return config["count_modes"]
    if config["comparison"].startswith("social_vs_subsocial"):
        return ["log2fc_threshold", "is_sig"]
    return ["log2fc_threshold"]


def mode_suffix(count_mode):
    return "" if count_mode == "log2fc_threshold" else f"_{count_mode}"


def pathway_bias_bar_workflow(config_file, gwf):
    with open(config_file) as f:
        config = yaml.safe_load(f)

    comparison = config["comparison"]
    out_dir = os.path.join(config["output_directory_path"], "pathway_bias_bar", comparison)
    log_path = config["log_directory_path"]

    for count_mode in social_subsocial_count_modes(config):
        output_suffix = mode_suffix(count_mode)
        gwf.target_from_template(
            name=f"{comparison}_pathway_bar{output_suffix}",
            template=pathway_bias_bar_template(
                input_csv=config["input_csv"],
                out_dir=out_dir,
                log_path=log_path,
                scripts_path=config["scripts_path"],
                comparison=comparison,
                positive_label=config["positive_label"],
                negative_label=config["negative_label"],
                positive_color=config["positive_color"],
                negative_color=config["negative_color"],
                unbiased_color=config["unbiased_color"],
                pathway_order=config["pathway_order"],
                account=config["account"],
                log2fc_column=config.get("log2fc_column", "log2fc"),
                log2fc_threshold=config.get("log2fc_threshold", 1.0),
                count_mode=count_mode,
                is_sig_column=config.get("is_sig_column", "is_sig"),
                consistency_column=config.get("consistency_column"),
                output_suffix=output_suffix,
            ),
        )
    return gwf


def combined_pathway_bias_bar_workflow(config_files, gwf):
    configs = []
    for cf in config_files:
        with open(cf) as f:
            configs.append(yaml.safe_load(f))

    configs.sort(key=lambda c: c.get("panel_index", 999))
    base = configs[0]
    out_dir = os.path.join(base["output_directory_path"], "pathway_bias_bar", "combined")

    def panels_for_mode(count_mode):
        panels = []
        for c in configs:
            comparison = c["comparison"]
            panel_mode = count_mode if comparison.startswith("social_vs_subsocial") else "log2fc_threshold"
            output_suffix = mode_suffix(panel_mode)
            tsv = os.path.join(c["output_directory_path"], "pathway_bias_bar",
                               comparison, f"pathway_bias_counts{output_suffix}.tsv")
            panels.append({
                "tsv": tsv,
                "title": c.get("panel_title", comparison),
                "positive_label": c["positive_label"],
                "negative_label": c["negative_label"],
                "positive_color": c["positive_color"],
                "negative_color": c["negative_color"],
                "unbiased_color": c["unbiased_color"],
            })
        return panels

    for count_mode in ["log2fc_threshold", "is_sig"]:
        output_suffix = mode_suffix(count_mode)
        gwf.target_from_template(
            name=f"combined_pathway_bar{output_suffix}",
            template=combined_pathway_bias_bar_template(
                panels=panels_for_mode(count_mode),
                out_dir=out_dir,
                log_path=base["log_directory_path"],
                scripts_path=base["scripts_path"],
                pathway_order=base["pathway_order"],
                account=base["account"],
                output_suffix=output_suffix,
            ),
        )
    return gwf
