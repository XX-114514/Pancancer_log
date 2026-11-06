"""Enumerate all prespecified patient contrasts as unestimable; no model or p-value is fitted."""
from pathlib import Path
import csv,gzip,itertools,collections,json,hashlib
r=Path(__file__).resolve().parents[1];mem=list(csv.DictReader((r/'audit/network_members.tsv').open(),delimiter='\t'));db=list(csv.DictReader((r/'audit/LR_database_background.tsv').open(),delimiter='\t'));strata=list(csv.DictReader((r/'audit/metadata_v2/patient_strata_eligibility.tsv').open(),delimiter='\t'));assert not any(x['status'].startswith('candidate') for x in strata);d=r/'results/patient_extension';d.mkdir(exist_ok=True);out=d/'all_CM_pair_LR_not_estimable.tsv.gz';summary=[];total=0
with gzip.open(out,'wt',compresslevel=5,newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['cm','source','target','id_cp_interaction','LR_resource_measurable','categorical_CMT_effect','continuous_H_effect','CI_low','CI_high','p','q','status'])
 for cm in sorted({x['cm'] for x in mem}):
  nodes=[x['subCluster'] for x in mem if x['cm']==cm];n=0
  for a,b in itertools.product(nodes,repeat=2):
   for lr in db:w.writerow([cm,a,b,lr['id_cp_interaction'],lr['all_subunits_measured'],'NA','NA','NA','NA','NA','NA','not_estimable_no_eligible_patient_stratum']);n+=1
  total+=n;summary.append({'cm':cm,'ordered_pairs_including_self':len(nodes)**2,'database_LRs':len(db),'predefined_pair_LR_slots':n,'models_fitted':0,'status':'no_eligible_patient_stratum'})
with (r/'audit/patient_family_status.tsv').open('w') as f:w=csv.DictWriter(f,summary[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
receipt={'status':'eligibility_enumeration_only_not_fitted_models','rows':total,'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'logical_path':'${RUN_ROOT}/results/patient_extension/all_CM_pair_LR_not_estimable.tsv.gz'};(r/'audit/patient_family_receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
