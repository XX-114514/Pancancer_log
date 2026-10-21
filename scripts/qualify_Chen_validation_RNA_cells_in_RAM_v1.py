"""Complete public small RNA-ISH inputs in RAM, metadata qualification only.

No expression/coordinate/type row is written or returned. Image groups are not
patients. Numeric assay semantics, author type rules and clinical keys remain
separate from payload integrity and field availability.
"""
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import csv
import gc
import hashlib
import io
import json
import math
import os
import re
import resource
import subprocess
import sys
import urllib.request
import numpy as np
import pandas as pd

GENES=['IDO1','CCR7','FOXP3','CD4','PDCD1','CCL22','EPCAM','CCL19','CD8A','CXCL10','CD3E','TCF7']
FILES=[('ACD_stemimmunity.csv',6409466,'1dba4e9c4b2b8f1e9193a921dadefe38'),
       ('ACD_tumor.csv',97389572,'52adc008769b41af30c2faf4c1c1995b')]


def reclaimable_memory():
    """Conservative reported free/inactive/speculative, not an OS guarantee."""
    s=subprocess.run(['vm_stat'],capture_output=True,text=True,check=True).stdout
    page=int(re.search(r'page size of (\d+) bytes',s).group(1))
    return page*sum(int(re.search(r'Pages '+k+r':\s+(\d+)',s).group(1)) for k in ['free','inactive','speculative'])


def profile(blob):
    frame=pd.read_csv(io.BytesIO(blob),dtype={'slide':str,'type':str,'comp':str},keep_default_na=False)
    assert set(GENES+['slide','x','y','type'])<=set(frame.columns)
    assert frame.slide.ne('').all() and frame.type.ne('').all()
    values=frame[GENES+['x','y']].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    types=Counter(frame.type.tolist());groups=set(frame.slide)
    summary=dict(cell_observations=len(frame),image_groups_not_reviewed_patients=len(groups),
        author_type_counts=dict(types),numeric_nonfinite_entries=int(np.sum(~np.isfinite(values))),
        negative_gene_entries=int(np.sum(values[:,:12]<0)),
        fractional_gene_entries=int(np.sum(np.isfinite(values[:,:12])&(values[:,:12]!=np.floor(values[:,:12])))),
        coordinate_duplicate_rows_within_image=int(frame.duplicated(subset=['slide','x','y']).sum()),
        author_compartment_column_present='comp' in frame,
        compartment_label_cardinality=int(frame.comp.nunique()) if 'comp' in frame else None)
    # Independent scalar CSV parse of every row and every target/coordinate.
    raw=csv.DictReader(io.StringIO(blob.decode('utf-8-sig')))
    n=0;nonfinite=0;negative=0;fractional=0;tt=Counter();gg=set();keys=set()
    for r in raw:
        n+=1;tt[r['type']]+=1;gg.add(r['slide']);keys.add((r['slide'],r['x'],r['y']))
        for j,name in enumerate(GENES+['x','y']):
            try:v=float(r[name])
            except ValueError:v=float('nan')
            nonfinite+=not math.isfinite(v)
            if j<12:
                negative+=v<0
                fractional+=math.isfinite(v) and v!=math.floor(v)
    assert (n,nonfinite,negative,fractional)==tuple(summary[k] for k in ['cell_observations','numeric_nonfinite_entries','negative_gene_entries','fractional_gene_entries'])
    assert tt==types and gg==groups
    assert len(frame)-len(keys)==summary['coordinate_duplicate_rows_within_image']
    return summary,groups,keys


def main():
    out=Path(os.environ['PANCANCER_CHEN_TASK']);frozen=json.loads((out/'private_frozen_complete_validation_RNA_input_v1.json').read_text())
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==frozen['reader_SHA256']
    production=out/'validation_RNA_input_v1';production.mkdir(exist_ok=False)
    save=lambda n,v:(production/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    started=datetime.now(timezone.utc).isoformat()
    save('private_checkpoint.json',dict(stage='RAM_BUDGET',started_UTC=started))
    available=reclaimable_memory();assert available>=768*2**20+3*2**30,'RAM_HEADROOM_NOT_ACCEPTED'
    proxy=os.environ['PANCANCER_SOURCE_HTTPS_PROXY']
    op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':proxy}))
    summaries=[];sets=[];coordinate_sets=[]
    for name,size,digest in FILES:
        save('private_checkpoint.json',dict(stage='COMPLETE_RAM_INPUT',file=name,started_UTC=started))
        req=urllib.request.Request('https://zenodo.org/api/records/11198494/files/'+name+'/content',
            headers={'User-Agent':'Bounded-research-input-qualification/1','Accept-Encoding':'identity'})
        buf=io.BytesIO();md5=hashlib.md5();sha=hashlib.sha256()
        with op.open(req,timeout=120) as response:
            assert response.status==200 and int(response.headers.get('Content-Length'))==size
            while True:
                b=response.read(1024**2)
                if not b:break
                buf.write(b);md5.update(b);sha.update(b)
                assert buf.tell()<=size
        assert buf.tell()==size and md5.hexdigest()==digest
        blob=buf.getvalue();summary,groups,keys=profile(blob)
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
        assert rss<=768*2**20,'OBSERVED_RAM_BUDGET_EXCEEDED'
        summary.update(file=name,bytes=size,MD5=md5.hexdigest(),SHA256=sha.hexdigest(),
            complete_official_length_MD5_accepted=True,all_numeric_and_type_counts_scalar_replayed=True,
            peak_RSS_bytes=int(rss),individual_arrays_RAM_only=True)
        summaries.append(summary);sets.append(groups);coordinate_sets.append(keys)
        del blob,buf;gc.collect()
    save('private_complete_validation_RNA_input_scope.json',dict(source='Chen_NatImmunol2024_Zenodo11198494',files=summaries,
        image_group_union_not_patients=len(sets[0]|sets[1]),image_groups_shared_between_files=len(sets[0]&sets[1]),
        exact_image_xy_keys_shared_between_files=len(coordinate_sets[0]&coordinate_sets[1]),
        patients_clinical_linkage_not_accepted=True,author_type_and_numeric_assay_semantics_pending=True,
        manually_selected_ROIs_not_population_random_samples=True,no_new_state_classifier_effect_or_clinical_model=True))
    save('private_terminal.json',dict(status='COMPLETE_FULL_RNA_INPUT_SOURCE_METADATA_ACCEPTED',started_UTC=started,
        completed_UTC=datetime.now(timezone.utc).isoformat(),reader_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        frozen_SHA256=hashlib.sha256((out/'private_frozen_complete_validation_RNA_input_v1.json').read_bytes()).hexdigest(),
        complete_file_length_and_official_MD5_and_scalar_numeric_replay=True,
        no_raw_payload_expression_coordinate_case_or_type_rows_written=True,no_source_individual_keys_returned=True))
    print(json.dumps(dict(status='COMPLETE_FULL_RNA_INPUT_SOURCE_METADATA_ACCEPTED')))


if __name__=='__main__':main()
