"""Freeze existing HNSCC native axis before any programme/clinical tests."""
import json,hashlib,os
from datetime import datetime,timezone
from pathlib import Path
import h5py,numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem


def main():
    project=Path(os.environ['PANCANCER_PROJECT_ROOT']);root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    contract=json.loads((atlas/'config/source_contract_v2_private.json').read_text())
    rec=next(r for r in contract['records'] if r['dataset']=='GSE286827')
    out=root/'HNSCC_frozen_MCP_gene_axis_20261014_v1';assert not out.exists();out.mkdir()
    g=h5py.File(rec['path'],'r') if str(rec['path']).endswith('.h5ad') else zarr.open_group(rec['path'],mode='r')
    var=read_elem(g[rec['var_key']]);columns=np.arange(len(var))
    if rec.get('new_parent_dataset'):
        source_order=list(read_elem(g['uns'])['source_order'])
        measured=np.asarray(g['varm/measured_in_source'])[:,source_order.index(rec['new_parent_dataset'])]
        columns=np.flatnonzero(measured);var=var.iloc[columns].copy()
    aliases={}
    for i,(original,row) in enumerate(var.iterrows()):
        if str(original).startswith('AMBIGUOUS::'):continue
        names={str(original)}
        for col in ['gene_symbol','feature_name','gene_name_original']:
            if col in var and str(row[col]) not in ['','nan','None','__SOURCE_VALUE_MISSING__']:names.add(str(row[col]))
        for name in names:aliases.setdefault(name,set()).add(int(columns[i]))
    rows=[]
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    for endpoint in ['Tnode9','Macnode18']:
        path=crc/'private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_genes.tsv')
        w=pd.read_csv(path,sep='\t');measured=np.array([len(aliases.get(g,set()))==1 for g in w.gene])
        fraction=float(np.abs(w.weight[measured]).sum()/np.abs(w.weight).sum())
        rows.append(dict(endpoint=endpoint,frozen_genes=len(w),measured_genes=int(measured.sum()),absolute_weight_retained=fraction,
                         complete_score_measurable=bool(measured.all()),partial95_measurable=bool(fraction>=.95)))
        mapped=w.copy();mapped['native_column']=[next(iter(aliases[g])) if len(aliases.get(g,set()))==1 else -1 for g in w.gene]
        mapped.to_csv(out/(endpoint+'_private_frozen_gene_columns.tsv'),sep='\t',index=False)
    pd.DataFrame(rows).to_csv(out/'aggregate_frozen_program_gene_coverage.tsv',sep='\t',index=False)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source='GSE286827',
      endpoint_definition='within-new-atlas multigene CD8_Effector and Macro_ candidates; not accepted original172 node mapping',
      sample_selection='first primary pre and post per exact GEO source-patient; exclude post2; both endpoints>=20',
      scoring='original frozen CRC genes/weights/mean/SD/quality coefficients; no refitting or gene selection',
      primary_scope='partial95 measurement sensitivity if complete score missing, retain complete-score missingness',
      planned_tests=['paired shift T score','paired shift myeloid score','paired shift equal-mean CRC-training-z score',
                     'pre cross-type Spearman','primary post cross-type Spearman','paired change cross-type Spearman'],
      family='six hypotheses*two fixed original program slots=12, unsupported MCP2 internally p1',
      random_seed=20261014,permutations=20000,patient_bootstrap=2000,
      inference_min_paired_patients=10,identity_frozen_before_programme_scores=True,
      actual_response_labels_missing=True,response_prediction_not_tested=True,
      therapy_interpretation='cohort treatment-time association, no randomized untreated paired control; pre biopsy/post resection sampling differences remain',
      no_imputed_clinical_outcomes=True,no_formal_identity_or_state_promotion=True)
    (out/'prospective_paired_cohort_plan_before_scores.json').write_text(json.dumps(plan,indent=2)+'\n')
    (out/'private_input_manifest.json').write_text(json.dumps(dict(status='HNSCC_GENE_AXIS_AND_PAIRED_COHORT_PLAN_FROZEN_NO_SCORES_READ',
      finished_UTC=datetime.now(timezone.utc).isoformat(),source_record=rec,
      source_contract_sha256=hashlib.sha256((atlas/'config/source_contract_v2_private.json').read_bytes()).hexdigest(),
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    print('HNSCC_GENE_AXIS_AND_PAIRED_COHORT_PLAN_FROZEN_NO_SCORES_READ')


if __name__=='__main__':main()
