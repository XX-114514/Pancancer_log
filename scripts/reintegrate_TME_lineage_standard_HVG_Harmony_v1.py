"""Run the existing standard recipe on an accepted full-parent lineage view.

Only the explicit, inherited clone-variable feature exclusion is added to the
parent pipeline's technical feature filter. Scanpy HVG/PCA/Harmony stay intact.
"""
import argparse,hashlib,json,re
from pathlib import Path
import reintegrate_TME_standard_HVG_Harmony_v1 as standard
def main(R):
 cfg=standard.load(R/'config/spec_v1.json');A=Path(cfg['lineage_scope_auditor']);assert standard.load(A/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;accepted=standard.load(A/'outputs/accepted_lineage_model_input_scope_v1.json');assert accepted['bridge_path']==cfg['working_acceptance'] and accepted['bridge_SHA256']==cfg['working_acceptance_SHA256'];assert cfg['additional_excluded_clone_variable_gene_pattern']==r'^(IG[HKL][VDJC]|TR[ABGD][VJ])'
 for name,digest in cfg['standard_adapter_code_SHA256'].items():assert standard.sha(R/'scripts'/name)==digest
 original=standard.technical;pattern=re.compile(cfg['additional_excluded_clone_variable_gene_pattern']);standard.technical=lambda gene:original(gene) or bool(pattern.match(gene))
 standard.save(R/'provenance/lineage_standard_recipe_adapter_v1.json',dict(lineage_scope_acceptance_SHA256=standard.sha(A/'outputs/accepted_lineage_model_input_scope_v1.json'),HVG_PCA_Harmony_recipe_unchanged=True,additional_excluded_clone_variable_gene_pattern=pattern.pattern,feature_exclusion_inherited_from_existing_lineage_protocol=True,transductive_joint_parent_refinement_not_ABC_winner=True,no_new_identity_or_QC_assigned=True))
 standard.main(R)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  standard.save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
