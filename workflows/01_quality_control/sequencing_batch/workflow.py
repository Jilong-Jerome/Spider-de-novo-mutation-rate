"""
GWF workflow for the sequencing batch summary.
"""

import sys
import os
from glob import glob

from gwf import Workflow

# Make workflow_sources importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'workflow_sources'))

from workflow_sources import sequencing_batch_workflow  # noqa: E402

gwf = Workflow(defaults={'account': 'spider2'})

for config_file in glob('configurations/*.config.yaml'):
    sequencing_batch_workflow(config_file=config_file, gwf=gwf)
