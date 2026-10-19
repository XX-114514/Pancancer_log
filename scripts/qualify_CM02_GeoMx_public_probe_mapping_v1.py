"""Recover public DCC probe mapping and unfiltered workbook structures, no models."""
import argparse
import gzip
import hashlib
import json
import re
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()


def main(a):
    assert not a.output.exists();a.output.mkdir(mode=0o700)
    with gzip.open(a.source/'GSE292098_Hs_R_NGS_WTA_v1.0.pkc.gz','rt') as f:pkc=json.load(f)
    rows=[];probes={};ambiguous=set()
    for target in pkc['Targets']:
        for p in target['Probes']:
            r=dict(RTS_ID=p['RTS_ID'],target=target['DisplayName'],CodeClass=target['CodeClass'],
                   ProbeID=p['ProbeID'],GeneID=p['GeneID'])
            rows.append(r)
            if p['RTS_ID'] in probes and probes[p['RTS_ID']]['target']!=r['target']:ambiguous.add(p['RTS_ID'])
            probes[p['RTS_ID']]=r
    for key in ambiguous:probes.pop(key)
    genes=sorted({p['target'] for p in probes.values() if p['CodeClass'].startswith('Endogenous')})
    gi={g:i for i,g in enumerate(genes)}
    samples=[];vectors=[];unmapped=set();all_probes=set();records=0
    with tarfile.open(a.source/'GSE271689_RAW.tar') as tar:
        for member in tar:
            if not member.isfile():continue
            raw=tar.extractfile(member).read();raw=gzip.decompress(raw) if raw[:2]==b'\x1f\x8b' else raw
            lines=raw.decode().splitlines();inside=False;v=np.zeros(len(genes),np.int64);n=0;unmapped_n=0
            for line in lines:
                if line=='<Code_Summary>':inside=True;continue
                if line=='</Code_Summary>':inside=False
                if not inside:continue
                key,count=line.split(',')
                if key=='RTS_ID':continue
                value=int(count);assert value>=0;n+=1;records+=1;all_probes.add(key)
                if key not in probes:unmapped.add(key);unmapped_n+=1;continue
                target=probes[key]
                if target['CodeClass'].startswith('Endogenous'):v[gi[target['target']]]+=value
            samples.append(dict(source_file=member.name,GSM=re.match(r'GSM\d+',Path(member.name).name).group(),
                                probe_records=n,unmapped_probe_records=unmapped_n))
            vectors.append(v)
    assert len(samples)==586
    pd.DataFrame(rows).to_csv(a.output/'private_candidate_PKC_probe_target_map.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(samples).to_csv(a.output/'private_all586_DCC_gene_mapping_registry.tsv',sep='\t',index=False)
    np.savez_compressed(a.output/'private_candidate_mapped_DCC_gene_counts.npz',counts=np.array(vectors,np.int64),gene=np.array(genes,str))
    (a.output/'private_unmapped_or_ambiguous_probes.json').write_text(json.dumps(dict(unmapped=sorted(unmapped),ambiguous=sorted(ambiguous)),indent=2)+'\n')
    books=[]
    for p in sorted(a.source.glob('*_0_filter_Greek.xlsx')):
        w=openpyxl.load_workbook(p,read_only=True,data_only=True);tables={}
        for name in ['SegmentProperties','TargetProperties','TargetCountMatrix']:
            s=w[name];s.reset_dimensions();it=s.values;header=list(next(it));data=list(it)
            assert all(len(r)==len(header) for r in data)
            tables[name]=pd.DataFrame(data,columns=header)
        segments=tables['SegmentProperties'];targets=tables['TargetProperties'];counts=tables['TargetCountMatrix']
        assert segments.SegmentDisplayName.is_unique and targets.TargetName.is_unique and counts.TargetName.is_unique
        assert set(counts.columns[1:])<=set(segments.SegmentDisplayName)
        matrix=counts.iloc[:,1:].apply(pd.to_numeric).to_numpy(float)
        assert np.isfinite(matrix).all() and (matrix>=0).all()
        prefix=p.stem
        segments.to_csv(a.output/(prefix+'_private_segment_metadata.tsv.gz'),sep='\t',index=False,compression='gzip')
        targets.to_csv(a.output/(prefix+'_private_target_metadata.tsv.gz'),sep='\t',index=False,compression='gzip')
        np.savez_compressed(a.output/(prefix+'_private_matrix.npz'),values=matrix,target=counts.TargetName.astype(str).to_numpy(),
                            segment=np.array(counts.columns[1:],str))
        books.append(dict(source=p.name,sha256=sha(p),segments=len(segments),targets=len(counts),
                          measurement_nonnegative=True,integer_counts=bool(np.equal(matrix,np.rint(matrix)).all()),
                          raw_or_normalized_scale_not_inferred_from_filename=True,
                          clinical_field_names=[str(x) for x in segments.columns if re.search('patient|case|response|recist|pfs|surviv|follow|status|treatment|therapy',str(x),re.I)],
                          unique_target_to_segment_axis=True))
    result=dict(status='COMPLETE_PUBLIC_DCC_PKC_AND_GREEK_WORKBOOK_QUALIFICATION_NO_EFFECTS',
                finished_utc=datetime.now(timezone.utc).isoformat(),DCC_files=len(samples),DCC_probe_records=records,
                observed_unique_probes=len(all_probes),unmapped_observed_probes=len(unmapped),ambiguous_PKC_probes=len(ambiguous),
                candidate_PKC_version=pkc['Version'],candidate_PKC_name=pkc['Name'],
                same_author_Greek_PKC_not_proof_of_Yale_exact_panel_version=True,
                patient_ROI_linkage_and_assay_scale_review_required=True,Greek_tables=books,
                no_single_cell_or_new_clinical_effect=True,producer_sha256=sha(Path(__file__)),
                source_archive_sha256=sha(a.source/'GSE271689_RAW.tar'),PKC_sha256=sha(a.source/'GSE292098_Hs_R_NGS_WTA_v1.0.pkc.gz'))
    (a.output/'probe_workbook_technical_public.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
