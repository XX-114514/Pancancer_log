"""Technical finisher test using cropped real CPDB outputs; never a scientific High/Low result."""
from pathlib import Path
import tempfile,shutil,json,subprocess,sys,pandas as pd,numpy as np
r=Path(__file__).resolve().parents[1];root=Path(tempfile.mkdtemp(prefix='summary_smoke_',dir=r/'audit'));(root/'scripts').mkdir();(root/'audit').mkdir();(root/'results').mkdir()
for name in ['summarize_cpdb.py','summarize_environment_drivers.py']:shutil.copy2(r/'scripts'/name,root/'scripts'/name)
s=pd.read_csv(next((r/'results/pilot').glob('statistical_analysis_significant_means_*.txt')),sep='\t');m=pd.read_csv(next((r/'results/pilot').glob('statistical_analysis_means_*.txt')),sep='\t');nodes=sorted({x for c in s if '|' in c for x in c.split('|')})[:2];cols=[x+'|'+y for x in nodes for y in nodes];metadata=[c for c in s if '|' not in c];members=[]
for i in range(1,10):
 cm=f'CM{i:02d}'
 for rank,node in enumerate(nodes,1):members.append(dict(cm=cm,subCluster=node,network_member=True,top10_rank=rank))
 for group in ['High','Low']:
  d=root/'results'/(cm+'_'+group);d.mkdir();keep=[c for c in cols if not (i==2 and group=='Low' and c==nodes[0]+'|'+nodes[1])]
  s[metadata+keep].to_csv(d/'statistical_analysis_significant_means_fixture.txt',sep='\t',index=False);m[[c for c in m if '|' not in c]+keep].to_csv(d/'statistical_analysis_means_fixture.txt',sep='\t',index=False);(d/'complete.json').write_text('{}')
pd.DataFrame(members).to_csv(root/'audit/member_comparison.tsv',sep='\t',index=False)
res=subprocess.run([sys.executable,str(root/'scripts/summarize_cpdb.py')],capture_output=True,text=True);(root/'console.txt').write_text(res.stdout+res.stderr);assert res.returncode==0,res.stderr
pairs=pd.read_csv(root/'results/summary/all_CM_CCI_pairs.tsv',sep='\t');missing=pairs[(pairs.cm=='CM02')&(pairs.source==nodes[0])&(pairs.target==nodes[1])];assert missing.Low_pooled.isna().all();other=pairs[pairs.status=='evaluated'];assert (other.delta_pooled==0).all()
lr=pd.read_csv(root/'results/summary/CM02_LR_detection.tsv.gz',sep='\t');x=lr[(lr.source==nodes[0])&(lr.target==nodes[1])];assert (x.detection=='not_evaluable_missing_cellpair').all();assert not x.cellpair_present_both.any();tests=pd.read_csv(root/'results/summary/all_CM_CCI_tests.tsv',sep='\t');assert (tests.p_cellpair_Wilcoxon==1).all();assert (tests.q_extension_all_CM==1).all();driver=pd.read_csv(root/'results/summary/HighLow_concentration.tsv',sep='\t');assert (driver.net_candidate_change==0).all();receipt={'status':'technical_finisher_smoke_pass','fixture':'cropped_real_pilot_duplicated_groups; not_scientific_results','tested':['missing_pair_NA_not_zero','identical_groups_zero_effect','all_CM_family_p_and_q','driver_zero_for_identical_measurable_pairs'],'fixture_directory':root.name};(r/'audit/finisher_smoke_acceptance.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
