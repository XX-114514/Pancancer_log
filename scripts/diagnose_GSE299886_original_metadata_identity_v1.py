"""Original archived CSV identity/coordinate join audit; no spatial effects.

Only aggregate match diagnostics are emitted; original files remain read-only.
"""
import gzip
import io
import json
import os
import re
import tarfile
from datetime import datetime,timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def strings(values):return np.array([str(x) for x in np.asarray(values)],dtype=str)


def normalized(values):
    # Exact decimal integer serialization only; never cast IDs through float.
    return np.array([re.sub(r'\.0+$','',v) if re.fullmatch(r'[+-]?\d+\.0+',v) else v for v in values],dtype=str)


def readcsv(tar,name,**kwargs):
    f=tar.extractfile(name);assert f is not None
    stream=gzip.GzipFile(fileobj=f) if name.endswith('.gz') else f
    return pd.read_csv(stream,**kwargs)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    source=root/'CM02_existing_targeted_spatial_schema_20261018_v2/private_all_eligible_targeted_input_schema.json'
    terminal=json.loads((source.parent/'private_terminal.json').read_text());assert sha(source)==terminal['files'][source.name]['sha256']
    inputs=json.loads(source.read_text());out=root/'GSE299886_original_metadata_identity_diagnostic_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,schema_sha256=sha(source),
        scope='All three coverage-eligible archived sources. Inspect original member/header maps, exact raw-count ID to metadata-ID bijections and finite source coordinates. No cell IDs or rows exported; all summaries private.',
        normalization='Exact strings first, optional exact decimal integer .0 suffix removal; never binary float conversion or fuzzy matching.',
        restriction='No data repair, array overwrite, cell-state assignment or spatial test. Full original RNA-to-existing-matrix equality remains required before any sidecar reuse.')
    (out/'plan_before_original_join_diagnostic.json').write_text(json.dumps(plan,indent=2)+'\n')
    allrows=[];details=[]
    for item in inputs:
        p=Path(item['analysis_file']);assert sha(p)==item['input_sha256']
        with h5py.File(p,'r') as h:
            archive=Path(str(read_elem(h['uns']['archive_name'])))
            mapping=read_elem(h['uns']['archive_member_map']);mapping=json.loads(mapping) if isinstance(mapping,str) else mapping
            details.append(dict(sample_id=item['sample_id'],archive=str(archive),member_map=mapping))
            obs=h['obs'];obsids={}
            for name in ['_index','barcode','cell_id','EntityID']:
                if name in obs:obsids[name]=strings(read_elem(obs[name]))
        with tarfile.open(archive,'r:') as tar:
            # Associate by the exact archived GSM, preserving distinct TMA members.
            members=[m.name for m in tar.getmembers() if m.isfile() and item['sample_id'].split('|')[-1] in m.name]
            # sample_id need not equal GSM; use original member map when available.
            flatten=[]
            def walk(value):
                if isinstance(value,str):flatten.append(value)
                elif isinstance(value,dict):
                    for v in value.values():walk(v)
                elif isinstance(value,list):
                    for v in value:walk(v)
            walk(mapping)
            matched=[v for v in flatten if v in tar.getnames()]
            prefix=set(re.findall(r'GSM\d+', '\n'.join(matched)))
            if prefix:members=[m.name for m in tar.getmembers() if m.isfile() and any(g in m.name for g in prefix)]
            counts=[n for n in members if 'cell_by_gene' in n.lower() and n.endswith(('.csv','.csv.gz'))]
            metas=[n for n in members if 'cell_metadata' in n.lower() and n.endswith(('.csv','.csv.gz'))]
            assert len(counts)==1 and len(metas)==1, 'Original archive members ambiguous: '+str((len(counts),len(metas)))
            ch=readcsv(tar,counts[0],nrows=0).columns.tolist();mh=readcsv(tar,metas[0],nrows=0).columns.tolist()
            count_id=ch[0];cid=readcsv(tar,counts[0],usecols=[count_id],dtype=str,keep_default_na=False)[count_id].to_numpy(str)
            keys=[v for v in mh if re.search(r'entity|cell.?id|barcode|^id$|^unnamed',v,re.I)]
            if mh[0] not in keys:keys=[mh[0]]+keys
            geo=[v for v in ['center_x','center_y','fov','global_x','global_y','x','y'] if v in mh]
            meta=readcsv(tar,metas[0],usecols=sorted(set(keys+geo)),dtype=str,keep_default_na=False)
            assert len(cid)==item['shape'][0]
            for basis,convert in [('EXACT',lambda x:x),('EXACT_DECIMAL_INTEGER_SERIALIZATION',normalized)]:
                cc=convert(cid);validc=~np.isin(cc,['','nan','NaN','None'])
                for key in keys:
                    mm=convert(meta[key].to_numpy(str));validm=~np.isin(mm,['','nan','NaN','None'])
                    count_unique=len(set(cc[validc]))==len(cc) and validc.all();meta_unique=len(set(mm[validm]))==len(mm) and validm.all()
                    lookup={v:i for i,v in enumerate(mm) if validm[i]};idx=np.array([lookup.get(v,-1) for v in cc])
                    hit=idx>=0;xyok=0
                    if hit.any() and all(v in meta for v in ['center_x','center_y']):
                        xx=pd.to_numeric(meta.center_x,errors='coerce').to_numpy();yy=pd.to_numeric(meta.center_y,errors='coerce').to_numpy();xyok=int((np.isfinite(xx[idx[hit]])&np.isfinite(yy[idx[hit]])).sum())
                    allrows.append(dict(sample_id=item['sample_id'],basis=basis,count_id_column=count_id,metadata_id_column=key,count_rows=len(cid),metadata_rows=len(meta),count_ids_unique=bool(count_unique),metadata_ids_unique=bool(meta_unique),matched_count_rows=int(hit.sum()),matched_finite_coordinate_rows=xyok,complete_unique_bijection=bool(count_unique and meta_unique and hit.all() and len(set(idx))==len(idx)),source_header_coordinate_fields='|'.join(geo),old_obs_exact_matches={n:int(np.array_equal(values,cid)) for n,values in obsids.items()}))
            details[-1].update(count_member=counts[0],metadata_member=metas[0],count_header=ch,metadata_header=mh)
    pd.DataFrame(allrows).to_csv(out/'private_all_original_ID_join_diagnostics.tsv',sep='\t',index=False)
    (out/'private_original_source_member_headers.json').write_text(json.dumps(details,indent=2)+'\n')
    files=[p for p in out.iterdir() if p.is_file()]
    t=dict(status='COMPLETE_ORIGINAL_ID_JOIN_DIAGNOSTIC_NO_REPAIR_OR_SPATIAL_TEST',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},existing_inputs_modified=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],completed_UTC=t['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
