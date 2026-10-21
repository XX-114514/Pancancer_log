"""Same complete RNA input gates with chunked parsing and lossless compact keys.

Keeps the v1 failure and unchanged scientific scope/768MiB+3GiB budget. Raw
strings are never expanded into an entire decoded StringIO/whole DataFrame.
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
import resource
import sys
import urllib.request
import numpy as np
import pandas as pd
from qualify_Chen_validation_RNA_cells_in_RAM_v1 import GENES,FILES,reclaimable_memory


def packed_key(slide,x,y):
    assert all('\0' not in a for a in [slide,x,y])
    return ('\0'.join([slide,x,y])).encode('utf-8')


def profile(blob):
    count=0;nonfinite=0;negative=0;fractional=0;types=Counter();groups=set();keys=set();comps=set()
    has_comp=None
    dtype={k:str for k in ['slide','type','comp','x','y']}
    for frame in pd.read_csv(io.BytesIO(blob),dtype=dtype,keep_default_na=False,chunksize=10000):
        assert set(GENES+['slide','x','y','type'])<=set(frame.columns)
        assert frame.slide.ne('').all() and frame.type.ne('').all()
        values=frame[GENES+['x','y']].apply(pd.to_numeric,errors='coerce').to_numpy(float)
        count+=len(frame);nonfinite+=int(np.sum(~np.isfinite(values)))
        negative+=int(np.sum(values[:,:12]<0))
        fractional+=int(np.sum(np.isfinite(values[:,:12])&(values[:,:12]!=np.floor(values[:,:12]))))
        types.update(frame.type.tolist());groups.update(frame.slide)
        keys.update(packed_key(a,b,c) for a,b,c in zip(frame.slide,frame.x,frame.y))
        if has_comp is None:has_comp='comp' in frame
        assert has_comp==('comp' in frame)
        if has_comp:comps.update(frame.comp)
    summary=dict(cell_observations=count,image_groups_not_reviewed_patients=len(groups),author_type_counts=dict(types),
        numeric_nonfinite_entries=nonfinite,negative_gene_entries=negative,fractional_gene_entries=fractional,
        coordinate_duplicate_rows_within_image=count-len(keys),author_compartment_column_present=has_comp,
        compartment_label_cardinality=len(comps) if has_comp else None)
    n=0;nf=0;ng=0;frac=0;tt=Counter();gg=set();kk=set()
    # TextIOWrapper decodes incrementally; no complete Unicode copy/StringIO.
    with io.TextIOWrapper(io.BytesIO(blob),encoding='utf-8-sig',newline='') as text:
        for r in csv.DictReader(text):
            n+=1;tt[r['type']]+=1;gg.add(r['slide']);kk.add(packed_key(r['slide'],r['x'],r['y']))
            for j,name in enumerate(GENES+['x','y']):
                try:v=float(r[name])
                except ValueError:v=float('nan')
                nf+=not math.isfinite(v)
                if j<12:
                    ng+=v<0;frac+=math.isfinite(v) and v!=math.floor(v)
    assert (n,nf,ng,frac)==(count,nonfinite,negative,fractional)
    assert tt==types and gg==groups and kk==keys
    return summary,groups,keys


def main():
    out=Path(os.environ['PANCANCER_CHEN_TASK']);fp=out/'private_frozen_complete_validation_RNA_input_v2.json'
    frozen=json.loads(fp.read_text());code=Path(__file__).parent
    for n,d in frozen['code_SHA256'].items():assert hashlib.sha256((code/n).read_bytes()).hexdigest()==d
    production=out/'validation_RNA_input_v2';production.mkdir(exist_ok=False)
    save=lambda n,v:(production/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    started=datetime.now(timezone.utc).isoformat()
    save('private_checkpoint.json',dict(stage='RAM_BUDGET',started_UTC=started))
    assert reclaimable_memory()>=768*2**20+3*2**30,'RAM_HEADROOM_NOT_ACCEPTED'
    op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':os.environ['PANCANCER_SOURCE_HTTPS_PROXY']}))
    summaries=[];sets=[];coordinate_sets=[]
    for name,size,digest in FILES:
        save('private_checkpoint.json',dict(stage='COMPLETE_RAM_INPUT',file=name,started_UTC=started))
        req=urllib.request.Request('https://zenodo.org/api/records/11198494/files/'+name+'/content',
            headers={'User-Agent':'Bounded-research-input-qualification/2','Accept-Encoding':'identity'})
        buf=io.BytesIO();md5=hashlib.md5();sha=hashlib.sha256()
        with op.open(req,timeout=120) as response:
            assert response.status==200 and int(response.headers.get('Content-Length'))==size
            while True:
                b=response.read(1024**2)
                if not b:break
                buf.write(b);md5.update(b);sha.update(b);assert buf.tell()<=size
        assert buf.tell()==size and md5.hexdigest()==digest
        blob=buf.getvalue()
        save('private_checkpoint.json',dict(stage='COMPLETE_SIGNATURE_ACCEPTED_CHUNKED_SCALAR_QUALIFICATION',file=name,started_UTC=started))
        summary,groups,keys=profile(blob)
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
        completed_UTC=datetime.now(timezone.utc).isoformat(),code_SHA256=frozen['code_SHA256'],
        frozen_SHA256=hashlib.sha256(fp.read_bytes()).hexdigest(),
        complete_file_length_and_official_MD5_and_scalar_numeric_replay=True,
        unchanged_768MiB_RSS_budget_and_3GiB_headroom=True,prior_v1_memory_failure_preserved=True,
        no_raw_payload_expression_coordinate_case_or_type_rows_written=True,no_source_individual_keys_returned=True))
    print(json.dumps(dict(status='COMPLETE_FULL_RNA_INPUT_SOURCE_METADATA_ACCEPTED')))


if __name__=='__main__':main()
