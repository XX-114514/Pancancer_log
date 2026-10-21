"""RAM-only author state bridge using reviewed library aliases and cell barcodes.

Only identity/state fields are selected from public annotation data frames.
TCR and clinical values are never selected, exported or used. This is a join
opportunity audit, not an authority-label update or clinical validation.
"""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import gzip,hashlib,io,json,os,re,subprocess,urllib.request
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

SOURCES=[('GSE286827_metadata_immunecells_broad.rds.gz','manual_celltype_broad','7835c18873d51f8bd571f6bf21ade5d0f18c4eac78d942d26a752a03f5d419a6'),('GSE286827_metadata_CD4T_TCR.rds.gz','cd4t_subcluster','b33adaaba779a271df11aa6834266133aa41a1952b745770124221909db28210'),('GSE286827_metadata_CD8T_TCR.rds.gz','cd8t_subcluster','dc78fc347c8b400e6ebf1add2492b1c7198ab68ccf833abae543db4dfbcf1343')]
PATIENT_FIELDS=['patient.tag','patient_title_token','patient_registry']
TIME_FIELDS=['time','title_timepoint','timepoint_registry']
BARCODE=re.compile(r'(?<![ACGT])([ACGT]{16})(?:-([0-9]+))?(?![ACGT])')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def norm(x):return str(x).strip().casefold()


def barcode(value,mode):
    m=list(BARCODE.finditer(str(value)))
    if len(m)!=1:return None
    return m[0][1] if mode=='within_reviewed_library_16nt' else m[0][1]+('-'+m[0][2] if m[0][2] else '')


def read_author(name,field,expected):
    url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/'+name
    with urllib.request.urlopen(url,timeout=20) as r:b=r.read(8*1024**2+1)
    assert len(b)<=8*1024**2 and hashlib.sha256(b).hexdigest()==expected
    raw=b
    for _ in range(3):
        if raw[:2]!=b'\x1f\x8b':break
        raw=gzip.GzipFile(fileobj=io.BytesIO(raw)).read(256*1024**2+1);assert len(raw)<=256*1024**2
    R='d<-readRDS(file("stdin","rb")); stopifnot(is.data.frame(d)); w<-c("Patient","Time",'+json.dumps(field)+'); stopifnot(all(w %in% names(d))); cat(jsonlite::toJSON(list(cell_key=rownames(d),patient=d$Patient,phase=d$Time,label=as.character(d[['+json.dumps(field)+']])),auto_unbox=FALSE))'
    p=subprocess.run(['Rscript','-e',R],input=raw,capture_output=True,timeout=30);assert p.returncode==0 and len(p.stdout)<128*1024**2
    f=pd.DataFrame(json.loads(p.stdout));assert f.cell_key.is_unique and not f.isna().any().any()
    return f


def fixture():
    assert barcode('sample_AAACCCAAGAAACCAT-1_12','full_barcode_suffix')=='AAACCCAAGAAACCAT-1'
    assert barcode('AAACCCAAGAAACCAT-2','within_reviewed_library_16nt')=='AAACCCAAGAAACCAT'
    assert barcode('unknown','within_reviewed_library_16nt') is None
    assert ('LIB1','AAACCCAAGAAACCAT')!=('LIB2','AAACCCAAGAAACCAT')
    return dict(barcode_is_never_joined_without_reviewed_library=True,suffix_and_no_suffix_modes_explicit=True,ambiguous_or_missing_barcode_not_guessed=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700);save=lambda n,v:(out/n).write_text(json.dumps(v,indent=2)+'\n')
    table=root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    save('plan_before_library_barcode_comparison.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),author_table_SHA256=sha(table),sources=[dict(file=a,label_field=b,SHA256=c) for a,b,c in SOURCES],patient_alias_fields=PATIENT_FIELDS,time_alias_fields=TIME_FIELDS,alias_normalization='strip whitespace and casefold only; no patient-number/suffix guesses',binding_selection='maximum uniquely mapped author groups; equal-coverage mappings must give identical assignments or fail; report all nine alias opportunities',cell_key_fields=['_index','original_cell_id','parent_cell_id'],barcode_modes=['full_barcode_suffix','within_reviewed_library_16nt'],barcode_join='Both author and current (reviewed exact library, synthetic barcode) must be unique. Preserve missing/colliding rows.',identity_checks='matched current patient and phase must agree with prior reviewed library registry; report conflicts without relabeling',metadata_and_cells_RAM_only=True,R_selected_fields=['rownames','Patient','Time','requested_subcluster_label'],no_TCR_clinical_or_expression_values_selected=True,no_formal_mapping_adoption=True))
    save('private_synthetic_fixture.json',fixture())
    cols=['library_id','patient_registry','timepoint_registry']+PATIENT_FIELDS+TIME_FIELDS;cols=list(dict.fromkeys(cols))
    table_frame=pd.read_csv(table,sep='\t',usecols=cols,dtype=str,keep_default_na=False);assert table_frame.library_id.is_unique
    lib_registry=table_frame.set_index('library_id')
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    dataset=np.asarray(read_elem(z['obs/dataset'])).astype(str);mask=dataset=='GSE286827'
    current={f:np.asarray(read_elem(z['obs'][f])).astype(str)[mask] for f in ['library_id','patient_id','timepoint','major','node_id','_index','original_cell_id','parent_cell_id']}
    assert set(current['library_id'])==set(lib_registry.index)
    expected_patient=pd.Series(current['library_id']).map(lib_registry.patient_registry).map(norm).to_numpy()
    expected_phase=pd.Series(current['library_id']).map(lib_registry.timepoint_registry).map(norm).to_numpy()
    patient_discord=np.asarray([norm(v) for v in current['patient_id']])!=expected_patient
    phase_discord=np.asarray([norm(v) for v in current['timepoint']])!=expected_phase
    opportunities=[];alias_rows=[];cross=[];labels=[]
    for name,field,expected in SOURCES:
        a=read_author(name,field,expected);pairs=[(norm(p),norm(t)) for p,t in zip(a.patient,a.phase)];unique_pairs=set(pairs)
        candidates=[]
        for pf in PATIENT_FIELDS:
            for tf in TIME_FIELDS:
                lookup=defaultdict(set)
                for row in table_frame.to_dict('records'):lookup[norm(row[pf]),norm(row[tf])].add(row['library_id'])
                mapping={k:next(iter(v)) for k,v in lookup.items() if len(v)==1}
                assignment=[mapping.get(k) for k in pairs];good=sum(k in mapping for k in unique_pairs)
                alias_rows.append(dict(annotation_source=name,patient_field=pf,time_field=tf,author_groups=len(unique_pairs),uniquely_bound_groups=good,ambiguously_bound_groups=sum(len(lookup[k])>1 for k in unique_pairs),bound_author_cells=sum(v is not None for v in assignment)))
                candidates.append((good,pf,tf,assignment))
        best=max(v[0] for v in candidates);ties=[v for v in candidates if v[0]==best]
        if best:assert all(v[3]==ties[0][3] for v in ties)
        author_lib=ties[0][3];a['reviewed_library']=author_lib
        labels.extend(dict(annotation_source=name,author_label=str(l),author_cells=int(n)) for l,n in Counter(a.label).items())
        for mode in ['full_barcode_suffix','within_reviewed_library_16nt']:
            ab=[barcode(v,mode) for v in a.cell_key];author_keys=[(l,b) if l is not None and b is not None else None for l,b in zip(author_lib,ab)]
            ac=Counter(k for k in author_keys if k is not None);lookup={k:i for i,k in enumerate(author_keys) if k is not None and ac[k]==1}
            for key_field in ['_index','original_cell_id','parent_cell_id']:
                cb=[barcode(v,mode) for v in current[key_field]];ck=[(l,b) if b is not None else None for l,b in zip(current['library_id'],cb)];cc=Counter(k for k in ck if k is not None)
                matched=[(i,lookup[k]) for i,k in enumerate(ck) if k is not None and cc[k]==1 and k in lookup]
                current_pos=np.asarray([i for i,j in matched],int);author_pos=np.asarray([j for i,j in matched],int)
                phase_bad=int(phase_discord[current_pos].sum());patient_bad=int(patient_discord[current_pos].sum())
                opportunities.append(dict(annotation_source=name,label_field=field,patient_alias_field=ties[0][1],time_alias_field=ties[0][2],barcode_mode=mode,current_cell_key_field=key_field,current_cells=len(ck),author_cells=len(a),uniquely_bound_author_groups=best,current_unparsable_cells=sum(v is None for v in cb),author_unparsable_cells=sum(v is None for v in ab),current_collision_keys=sum(v>1 for v in cc.values()),author_collision_keys=sum(v>1 for v in ac.values()),matched_cells=len(matched),matched_current_phase_registry_disagreements=int(phase_bad),matched_current_patient_registry_disagreements=int(patient_bad),mapping_adopted=False))
                counter=Counter(zip(current['timepoint'][current_pos],current['major'][current_pos],current['node_id'][current_pos],a.label.to_numpy()[author_pos]))
                cross.extend(dict(annotation_source=name,barcode_mode=mode,current_cell_key_field=key_field,phase=str(p),current_major=str(m),current_node=str(n),author_label=str(l),matched_cells=int(k)) for (p,m,n,l),k in sorted(counter.items()))
    assert len(alias_rows)==27 and len(opportunities)==18
    save('private_all27_reviewed_library_alias_opportunities.json',alias_rows);save('private_all18_library_barcode_opportunities.json',opportunities);save('private_all_author_label_coverage.json',labels);save('private_all_source_node_state_crosscounts.json',cross)
    files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()};assert sum(v['bytes'] for v in files.values())<4*1024**2
    terminal=dict(status='COMPLETE_REVIEWED_LIBRARY_AND_BARCODE_AUTHOR_STATE_BRIDGE_OPPORTUNITIES_NO_ADOPTION',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),all_metadata_cell_rows_and_keys_RAM_only=True,TCR_clinical_and_expression_values_not_selected=True,original_atlas_unchanged=True,files=files)
    save('private_terminal.json',terminal);print(json.dumps(dict(status=terminal['status'],terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
