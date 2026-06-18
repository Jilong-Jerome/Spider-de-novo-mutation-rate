from workflow_targets import *
from gwf import *

gwf = Workflow()

alphabet_list = generate_alphabet_list("aa", "cr")
for part in alphabet_list:
    merge_and_filter(gwf,part)
