import glob
import sys

sys.path.insert(0, './workflow_sources')

from gwf import Workflow
from workflow_sources import final_individual_variation_workflow


gwf = Workflow()

species_config_files = sorted(glob.glob('../dnm_variance/configurations/[A-Z]*.config.yaml'))
cv_config_file = '../dnm_variance/configurations/cv_comparison.config.yaml'

gwf = final_individual_variation_workflow(
    config_files=species_config_files,
    cv_config_file=cv_config_file,
    gwf=gwf,
)
