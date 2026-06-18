import sys
import glob

sys.path.insert(0, './workflow_sources')

from gwf import Workflow
from workflow_sources import dnm_variance_workflow, dnm_cv_comparison_workflow, dnm_individual_cv_comparison_workflow

gwf = Workflow()

all_config_files = glob.glob('./configurations/*.config.yaml')

species_config_files = [f for f in all_config_files if 'cv_comparison' not in f]
cv_config_file = './configurations/cv_comparison.config.yaml'

for config in species_config_files:
    gwf = dnm_variance_workflow(config_file=config, gwf=gwf)

gwf = dnm_cv_comparison_workflow(
    config_files=species_config_files,
    cv_config_file=cv_config_file,
    gwf=gwf,
)

gwf = dnm_individual_cv_comparison_workflow(
    config_files=species_config_files,
    cv_config_file=cv_config_file,
    gwf=gwf,
)
