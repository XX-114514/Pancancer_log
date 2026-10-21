"""Public author annotation in RAM; exact cell-key opportunities, not adoption."""
from datetime import datetime,timezone
from pathlib import Path
from collections import Counter
import gzip,hashlib,io,json,os,subprocess,urllib.request
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists()
    url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_metadata_immunecells_broad.rds.gz'
    expected='7835c18873d51f8bd571f6bf21ade5d0f18c4eac78d942d26a752a03f5d419a6'
    out.mkdir(mode=0o700);save=lambda n,v:(out/n).write_text(json.dumps(v,indent=2)+'\n')
    save('plan_before_exact_cell_key_comparison.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),url=url,expected_metadata_SHA256=expected,metadata_bytes_and_individual_rows_RAM_only=True,R_values_whitelist=['row_names','Patient','Time','manual_celltype_broad'],no_expression_sequence_or_outcome_values_selected=True,key_policy='Exact literal current index/original/parent cell keys; report all opportunities, no fuzzy barcode or patient suffix mapping. No adoption until version, library/phase and count audits.',original_atlas_read_only=True))
    with urllib.request.urlopen(url,timeout=20) as r:b=r.read(4*1024**2+1)
    assert len(b)<=4*1024**2 and hashlib.sha256(b).hexdigest()==expected
    raw=b
    for _ in range(3):
        if raw[:2]!=b'\x1f\x8b':break
        raw=gzip.GzipFile(fileobj=io.BytesIO(raw)).read(256*1024**2+1);assert len(raw)<=256*1024**2
    command='d<-readRDS(file("stdin","rb")); stopifnot(is.data.frame(d)); w<-c("Patient","Time","manual_celltype_broad"); stopifnot(all(w %in% names(d))); v<-list(cell_key=rownames(d),patient=d$Patient,phase=d$Time,label=d$manual_celltype_broad); cat(jsonlite::toJSON(v,auto_unbox=FALSE))'
    r=subprocess.run(['Rscript','-e',command],input=raw,capture_output=True,timeout=30);assert r.returncode==0 and len(r.stdout)<128*1024**2
    data=json.loads(r.stdout);author=pd.DataFrame(data);assert author.cell_key.is_unique and not author.isna().any().any()
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    dataset=np.asarray(read_elem(z['obs/dataset'])).astype(str);mask=dataset=='GSE286827'
    a=author.set_index('cell_key');keys=set(a.index);rows=[];cross=[]
    major=np.asarray(read_elem(z['obs/major'])).astype(str)[mask]
    node=np.asarray(read_elem(z['obs/node_id'])).astype(str)[mask]
    phase=np.asarray(read_elem(z['obs/timepoint'])).astype(str)[mask]
    for field in ['_index','original_cell_id','parent_cell_id']:
        if field not in z['obs']:continue
        v=np.asarray(read_elem(z['obs'][field])).astype(str)[mask]
        unique=pd.Index(v).is_unique;hit=np.fromiter((s in keys for s in v),bool,len(v))
        rows.append(dict(source='GSE286827',current_key_field=field,current_source_cells=len(v),current_keys_unique=bool(unique),exact_matching_current_cells=int(hit.sum()),unique_author_keys_matched=len(set(v[hit])),metadata_author_cells=len(a),author_keys_globally_unique=True,formal_mapping_adopted=False))
        if unique and hit.any():
            labels=a.loc[v[hit],'label'].to_numpy();counter=Counter(zip(phase[hit],major[hit],node[hit],labels))
            cross.extend(dict(source='GSE286827',current_key_field=field,phase=str(p),current_major=str(m),current_node=str(n),author_broad_label=str(l),matched_cells=int(k)) for (p,m,n,l),k in sorted(counter.items()))
    label_counts=Counter(author.label);save('private_source_author_label_coverage.json',[dict(source='GSE286827',author_broad_label=str(k),author_cells=v) for k,v in sorted(label_counts.items())])
    save('private_all_exact_key_opportunities.json',rows);save('private_all_exact_key_node_author_crosscounts.json',cross)
    files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()};assert sum(v['bytes'] for v in files.values())<2*1024**2
    terminal=dict(status='COMPLETE_OFFICIAL_AUTHOR_BROAD_METADATA_EXACT_KEY_OPPORTUNITIES_NO_LABEL_ADOPTION',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),metadata_compressed_SHA256=expected,no_cell_rows_clinical_values_sequences_or_raw_metadata_persisted=True,original_atlas_unchanged=True,files=files)
    save('private_terminal.json',terminal);print(json.dumps(dict(status=terminal['status'],terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
