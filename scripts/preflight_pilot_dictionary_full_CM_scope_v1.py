"""All frozen source axes and CM capture denominators; no new biological scores."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=Path(os.environ['PANCANCER_PILOT_FULL_SCOPE_PREFLIGHT_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256']
I=R/'all_CM_state_native_input_preflight_20261014_v1';assert sha(I/'inventory.json')==plan['full_scope_inventory_sha256']
assert json.loads((R/'all_CM_state_row_preflight_audit_20261014_v1/independent_audit.json').read_text())['status'].startswith('PASS_ALL4532933')
D=R/'state_programme_dictionary_20261014_v2/frozen_programme_dictionary.json';assert sha(D)==plan['pilot_dictionary_sha256'];dictionary=json.loads(D.read_text());programmes=dictionary['programmes'];assert len(programmes)==98
G=Path(os.environ['PANCANCER_CM_CORE_GRAPH']);assert sha(G)==plan['CM_core_graph_sha256'];graph=pd.read_csv(G,sep='\t');assert graph[['cm','node']].drop_duplicates().shape[0]==len(graph) and graph.cm.nunique()==9
units=pd.read_csv(I/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False);assert len(units)==35567 and units.sample_uid.nunique()==1290
assert set(graph.node).issubset(set(units.subCluster))
source=pd.read_csv(I/'all_source_native_scope.tsv',sep='\t');U=AR/'20261010_uniform_expression_annotation_v1';coverage=[];capture=[];member=[]
for row in source.itertuples():
 sid=row.gse_id;path=U/'by_gse'/sid/'genes.tsv';assert sha(path)==row.native_gene_axis_sha256;genes=set(pd.read_csv(path,sep='\t').gene);sub=units[units.gse_id.eq(sid)];bymajor=sub.groupby('major',observed=True).agg(cells=('n_cells','sum'),libraries=('sample_uid','nunique'))
 for p in programmes:
  names=p['primary_genes'];present=sorted(set(names)&genes);missing=sorted(set(names)-genes);fraction=len(present)/len(names);gene_gate=fraction>=plan['min_gene_fraction'] and len(present)>=plan['min_genes']
  for major in p['allowed_majors']:
   u=sub[sub.major.eq(major)];coverage.append(dict(gse_id=sid,major=major,programme_id=p['programme_id'],stability_main_candidate=p['stability_main_candidate'],signature_genes=len(names),native_measurable_genes=len(present),native_gene_fraction=fraction,native_gene_axis_gate=gene_gate,signature_missing_genes=';'.join(missing),original_major_captured=bool(len(u)),cells=int(u.n_cells.sum()),libraries=u.sample_uid.nunique(),original_nodes=u.subCluster.nunique(),sample_node_units_at_least20=int(u.n_cells.ge(20).sum()),expression_score_not_computed=True,programme_activity_not_assumed=True))
 for cm,g in graph.groupby('cm',sort=True):
  core=g.set_index('node');u=sub[sub.subCluster.isin(core.index)];observed=set(u.subCluster);eligible=set(u.loc[u.n_cells.ge(plan['min_cells_per_member_library']),'subCluster']);weight=float(core.core_weight.sum());assert np.isclose(weight,1,atol=1e-9)
  both=u[u.n_cells.ge(plan['min_cells_per_member_library'])].groupby('sample_uid',observed=True).subCluster.nunique();twomember=int(both.ge(2).sum())
  allmajors=set(g.major);pilotmajors=set(major for p in programmes for major in p['allowed_majors']);potentialpilot=eligible&set(g.loc[g.major.isin(pilotmajors),'node'])
  capture.append(dict(gse_id=sid,cm=cm,core_members=len(g),core_weight_total=weight,members_captured_any_cells=len(observed),captured_core_weight_any_cells=float(core.loc[list(observed),'core_weight'].sum()) if observed else 0,members_at_least20_in_any_library=len(eligible),captured_core_weight_at_least20=float(core.loc[list(eligible),'core_weight'].sum()) if eligible else 0,libraries_with_at_least2_members20=twomember,pilot_T_Mac_members_at_least20=len(potentialpilot),pilot_T_Mac_core_weight_at_least20=float(core.loc[list(potentialpilot),'core_weight'].sum()) if potentialpilot else 0,all_major_state_models_pending=True,expression_state_activity_not_tested=True,sample_units_not_independent_patients=True))
  for node,r in core.iterrows():
   uu=u[u.subCluster.eq(node)];member.append(dict(gse_id=sid,cm=cm,node=node,major=r.major,core_weight=r.core_weight,cells=int(uu.n_cells.sum()),libraries=uu.sample_uid.nunique(),libraries_at_least20=int(uu.n_cells.ge(20).sum()),has_fixed_pilot_programmes=r.major in pilotmajors))
O=R/'pilot_dictionary_full_CM_scope_preflight_20261014_v1';assert not O.exists();O.mkdir()
for name,rows in [('all52_source_pilot_programme_gene_capture.tsv',coverage),('all9CM_by52source_native_capture.tsv',capture),('all_CM_core_member_source_capture.tsv',member)]:pd.DataFrame(rows).to_csv(O/name,sep='\t',index=False)
c=pd.DataFrame(coverage);m=pd.DataFrame(capture)
out=dict(status='COMPLETE_FROZEN_PILOT_GENES_AND_ALL9CM_NATIVE_CAPTURE_DENOMINATORS_NOT_STATE_ASSOCIATIONS',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),pilot_dictionary_sha256=sha(D),full_scope_inventory_sha256=sha(I/'inventory.json'),CM_core_graph_sha256=sha(G),sources=52,CMs=9,original_cells=4532933,original_nodes=172,original_libraries=1290,programme_ids=98,source_major_programme_records=len(c),gene_axis_gate_pass_records=int(c.native_gene_axis_gate.sum()),gene_gate_and_original_major_captured_records=int((c.native_gene_axis_gate&c.original_major_captured).sum()),CM_source_records=len(m),core_member_records=len(graph),source_member_records=len(member),all9CM_scope_included=True,no_expression_copied_or_scores_computed=True,all_gene_gate_pass_does_not_imply_state_activation=True,formal_states=False,no_LR_or_clinical_tests=True,output_sha256={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
