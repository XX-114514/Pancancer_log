"""All-row source capture qualification, not state assignment or association.

Source payloads, keys, coordinates and clinical rows stay in RAM. Only cohort
aggregates and technical progress are emitted. No outcome-guided thresholds.
"""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path,PurePosixPath
import argparse,csv,hashlib,io,json,math,os,re,resource,time,zipfile
import numpy as np
import pandas as pd
from Danenberg_bounded_member_transport_v1 import RangeFile,IntegrityError
from qualify_Danenberg_clinical_metadata_v1 import binary,clean,positive_time

URL='https://zenodo.org/api/records/7324285/files/MBTMEStrPublicProcessedDataCode.zip/content'
SIZE=910666172
DIRECTORY_SHA='17911e718220292268a8cddfac32fb4bfa6c021a615ead0e011464c4a3658814'
HEADER_SHA='8b5fe9782f47b2f7a29e2e75c514fab6924353a39683f4782dcc5fbcbbb78375'
MARKERS=['CD3','CD4','CD8','FOXP3','Ki-67','PD-1','ICOS','OX40','GITR','Histone H3','DNA1','DNA2']
FLAGS=['is_tumour','is_hotAggregate','is_epithelial','is_normal','is_dcis']


def require(ok,reason):
    if not ok:raise IntegrityError(reason)


def directory(z):
    return hashlib.sha256(json.dumps([(x.filename,x.file_size,x.compress_size,x.CRC)
        for x in z.infolist()],separators=(',',':')).encode()).hexdigest()


def member(z,name,size=None,crc=None):
    found=[x for x in z.infolist() if PurePosixPath(x.filename).name==name]
    require(len(found)==1,'member_not_unique');m=found[0]
    require(not m.is_dir() and not m.flag_bits&1 and m.compress_type==zipfile.ZIP_DEFLATED,
        'member_type_or_encryption')
    require(not PurePosixPath(m.filename).is_absolute() and '..' not in PurePosixPath(m.filename).parts,
        'member_path_rejected')
    require((m.external_attr>>16)&0o170000!=0o120000,'member_symlink')
    if size is not None:require(m.file_size==size,'member_size_changed')
    if crc is not None:require(m.CRC==crc,'member_CRC_changed')
    return m


def author_alias(label):
    s=re.sub('MHC-','MHC ',label);s=s.replace('-med','med')
    s=re.sub('-high|-hi|high','hi',s);s=re.sub('-low|-lo|low','lo',s)
    s=re.sub(r'(\+|lo$|med$|hi$|lo |med |hi |-)',lambda m:'^{'+m[0]+'}',s)
    s=s.replace(' }','}').replace('reg','Reg').replace('ex','Ex')
    s=re.sub('(Reg|Ex)',lambda m:'_{'+m[0]+'}',s)
    return s.replace('CK8^{-}','CK8-')


def qc_mask(frame):
    flags={f:frame[f].map(binary) for f in FLAGS}
    main=flags['is_tumour'].eq(True);artifact_clear=flags['is_hotAggregate'].eq(False)
    vals=frame[MARKERS+['Location_Center_X','Location_Center_Y','AreaShape_Area']].apply(pd.to_numeric,errors='coerce')
    finite=np.isfinite(vals.to_numpy()).all(axis=1)
    area=vals.AreaShape_Area.gt(0).to_numpy()
    good=main.to_numpy()&artifact_clear.to_numpy()&finite&area
    reasons=Counter()
    for cond,name in [(~main.to_numpy(),'NOT_EXPLICIT_MAIN_INVASIVE_SAMPLE'),
                      (main.to_numpy()&~artifact_clear.to_numpy(),'ARTIFACT_OR_UNKNOWN_ARTIFACT'),
                      (main.to_numpy()&artifact_clear.to_numpy()&~finite,'NONFINITE_REQUIRED_MEASUREMENT'),
                      (main.to_numpy()&artifact_clear.to_numpy()&finite&~area,'NONPOSITIVE_AREA')]:
        reasons[name]=int(cond.sum())
    return good,vals,reasons


class HashReader(io.RawIOBase):
    def __init__(self,stream):self.stream=stream;self.sha=hashlib.sha256();self.bytes=0
    def readable(self):return True
    def readinto(self,b):
        raw=self.stream.read(len(b));self.sha.update(raw);self.bytes+=len(raw);b[:len(raw)]=raw
        return len(raw)


def run(output):
    require(not output.exists(),'refuse_overwrite')
    progress=output.with_name(output.stem+'_progress.json')
    started=datetime.now(timezone.utc).isoformat();last=[0]
    def checkpoint(values):
        if time.monotonic()-last[0]>10:
            progress.write_text(json.dumps(dict(status='RUNNING_NATIVE_QUALIFICATION',started_UTC=started,
                observation_UTC=datetime.now(timezone.utc).isoformat(),**values),indent=2)+'\n');last[0]=time.monotonic()
    rf=RangeFile(URL,SIZE,int(os.environ['PANCANCER_EXISTING_CLASH_PORT']),checkpoint)
    try:
        with zipfile.ZipFile(rf) as z:
            require(directory(z)==DIRECTORY_SHA,'directory_changed')
            clinical_raw=z.read(member(z,'IMCClinical.csv',51316,1165356369))
            clinical_rows=list(csv.DictReader(io.StringIO(clinical_raw.decode('utf-8-sig'))))
            clinical={clean(r['metabric_id']):r for r in clinical_rows}
            require(len(clinical)==len(clinical_rows) and '' not in clinical,'clinical_key_conflict')
            raw=z.read(member(z,'SI_metaclusters.csv',1113,1536880894))
            require(hashlib.sha256(raw).hexdigest()=='e6602949a3806698489048d01e20e241f512ce6da003cf339afeda66ef09699d','dictionary_changed')
            published={r['meta_description'].strip() for r in csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))}
            aliases={}
            for label in published:
                for key in (label,author_alias(label)):
                    require(key not in aliases or aliases[key]==label,'ambiguous_label_alias');aliases[key]=label
            m=member(z,'SingleCells.csv',849239430,3748430769)
            require(m.compress_size==381180527,'compressed_size_changed')
            all_counts=Counter();reason_counts=Counter();flag_unknown=Counter();label_counts=Counter()
            objects=defaultdict(set);image_patient={};image_scope={};patient_counts=Counter()
            parent_counts=defaultdict(Counter);parent_ranges={};marker_negative=Counter()
            kept_images=set();all_patients=set();matched_patients=set();unmapped=Counter()
            with z.open(m) as source:
                h=HashReader(source);buffer=io.BufferedReader(h,buffer_size=1024**2)
                header=buffer.readline();require(hashlib.sha256(header).hexdigest()==HEADER_SHA,'header_changed')
                names=next(csv.reader([header.decode('utf-8-sig').rstrip('\r\n')]))
                for frame in pd.read_csv(buffer,names=names,header=None,dtype=str,keep_default_na=False,
                                         chunksize=20000,on_bad_lines='error'):
                    require(frame.shape[1]==len(names) and frame.notna().all().all(),'malformed_cells_row')
                    all_counts['source_cells']+=len(frame)
                    for (image,patient),part in frame.groupby(['ImageNumber','metabric_id'],sort=False):
                        require(bool(image.strip()) and bool(patient.strip()),'missing_image_or_patient_key')
                        require(image not in image_patient or image_patient[image]==patient,'image_patient_conflict')
                        image_patient[image]=patient;all_patients.add(patient)
                        scope=set(zip(part.is_tumour.map(binary),part.is_normal.map(binary),part.is_dcis.map(binary)))
                        # Main/normal sample flags must agree; DCIS can vary across fields of view, not within an image.
                        require(len(scope)==1,'within_image_sample_scope_conflict')
                        scope=next(iter(scope));require(image not in image_scope or image_scope[image]==scope,'image_scope_conflict')
                        image_scope[image]=scope
                        obj=part.ObjectNumber.tolist();require(all(x.strip() for x in obj),'missing_object_key')
                        unique=set(obj);require(len(unique)==len(obj) and not objects[image].intersection(unique),'duplicate_compound_cell_key')
                        objects[image].update(unique)
                    for f in FLAGS:flag_unknown[f]+=int(frame[f].map(binary).isna().sum())
                    good,values,reasons=qc_mask(frame);reason_counts.update(reasons)
                    kept=frame.loc[good];measure=values.loc[good];all_counts['QC_cells']+=len(kept)
                    kept_images.update(kept.ImageNumber);patient_counts.update(kept.metabric_id)
                    labels=kept.cellPhenotype.str.strip().map(aliases)
                    unmapped.update(kept.loc[labels.isna(),'cellPhenotype'].map(lambda x:'UNMAPPED_ENCODING'))
                    label_counts.update(labels.dropna())
                    for marker in MARKERS:marker_negative[marker]+=int(measure[marker].lt(0).sum())
                    for label in ('Treg & Tex','CD4+ T cells','CD8+ T cells','CD4+ T cells & APCs'):
                        mask=labels.eq(label)
                        for patient,idx in kept.loc[mask].groupby('metabric_id',sort=False).groups.items():
                            parent_counts[patient][label]+=len(idx);a=measure.loc[idx,MARKERS].to_numpy()
                            key=(patient,label);lo=a.min(axis=0);hi=a.max(axis=0)
                            if key in parent_ranges:lo=np.minimum(lo,parent_ranges[key][0]);hi=np.maximum(hi,parent_ranges[key][1])
                            parent_ranges[key]=(lo,hi)
                    require(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<=8*1024**3,'RAM_budget_exceeded')
                require(buffer.read()==b'' and h.bytes==m.file_size,'incomplete_member_read')
            require(sum(objects_i.__len__() for objects_i in objects.values())==all_counts['source_cells'],'independent_key_count_mismatch')
            require(sum(patient_counts.values())==all_counts['QC_cells'],'independent_QC_count_mismatch')
            require(all_counts['QC_cells']+sum(reason_counts.values())==all_counts['source_cells'],'QC_partition_mismatch')
            coverage=defaultdict(Counter)
            for patient,n in patient_counts.items():
                if patient not in clinical:continue
                matched_patients.add(patient);r=clinical[patient];flag=binary(r['isValidation'])
                bucket=coverage['AUTHOR_VALIDATION_TRUE' if flag is True else 'AUTHOR_VALIDATION_FALSE' if flag is False else 'UNKNOWN_FLAG']
                bucket['QC_patient_keys']+=1
                follow=positive_time(r['yearsToStatus']);event=binary(r['DeathBreast'])
                opportunity=follow is not None and event is not None
                if opportunity:bucket['reported_followup_patients']+=1;bucket['reported_death_events']+=int(event)
                counts=parent_counts[patient]
                for label in ('Treg & Tex','CD4+ T cells','CD8+ T cells','CD4+ T cells & APCs'):
                    if counts[label]>=20:
                        bucket[label+'_AT_LEAST_20']+=1
                        lo,hi=parent_ranges[(patient,label)]
                        bucket[label+'_ALL_FIXED_MARKERS_VARY']+=int(np.all(hi>lo))
                pair=counts['Treg & Tex']>=20 and counts['CD8+ T cells']>=20
                if pair:
                    bucket['MIXED_REG_TEX_AND_CD8_AT_LEAST_20_EACH']+=1
                    if opportunity:bucket['PAIR_REPORTED_FOLLOWUP']+=1;bucket['PAIR_REPORTED_DEATH_EVENTS']+=int(event)
                    complete=all(clean(r[f]) for f in ['ERStatus','LymphNodesOrdinal','sizeOrdinal','Grade','ERBB2_pos']) and binary(r['ERBB2_pos']) is not None
                    if opportunity and complete:bucket['PAIR_FIVE_COVARIATE_OPPORTUNITY']+=1;bucket['PAIR_FIVE_COVARIATE_DEATH_EVENTS']+=int(event)
            with zipfile.ZipFile(rf) as replay:require(directory(replay)==DIRECTORY_SHA,'fresh_directory_changed')
            result=dict(status='COMPLETE_NATIVE_CAPTURE_QUALIFICATION_NOT_STATE_OR_ASSOCIATION',
                started_UTC=started,completed_UTC=datetime.now(timezone.utc).isoformat(),
                **all_counts,QC_exclusion_partition=dict(reason_counts),unknown_flags=dict(flag_unknown),
                source_images=len(objects),source_author_patient_keys=len(all_patients),QC_images=len(kept_images),
                QC_author_patient_keys=len(patient_counts),clinical_matched_QC_keys=len(matched_patients),
                QC_keys_without_clinical=len(patient_counts)-len(matched_patients),
                clinical_keys_without_QC=len(set(clinical)-set(patient_counts)),
                published_label_QC_counts=dict(label_counts),unmapped_QC_label_count=sum(unmapped.values()),
                phenotype_encoding_qualified=not unmapped,negative_compensated_MFI_counts_preserved=dict(marker_negative),
                by_author_validation_flag=[dict(author_flag=k,**v) for k,v in sorted(coverage.items())],
                full_selected_member_CRC_and_bytes=True,selected_member_SHA256=h.sha.hexdigest(),
                fresh_directory_SHA256=DIRECTORY_SHA,all_compound_keys_unique_and_images_single_patient=True,
                independent_count_reconciliation=True,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                requested_bytes_including_retries=rf.requested,successful_range_bytes=rf.transferred,requests=rf.requests,
                failure_classes=rf.failures,no_individual_payload_written=True,no_full_archive_MD5_claim=True,
                response_identity_headers_available=rf.identity!=(None,None),
                no_CM02_projection_or_pure_Treg_identity_claim=True,no_state_cutoffs_or_clinical_effects_computed=True,
                CD45RO_isoform_unqualified_excluded=True,
                producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
            output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
            progress.write_text(json.dumps(dict(status=result['status'],completed_UTC=result['completed_UTC']),indent=2)+'\n')
    except BaseException as error:
        progress.write_text(json.dumps(dict(status='FAILED_NATIVE_QUALIFICATION',error_class=type(error).__name__,
            qualification_reason=str(error) if isinstance(error,IntegrityError) else 'EXTERNAL_ERROR_VALUES_WITHHELD',
            observation_UTC=datetime.now(timezone.utc).isoformat(),requested_bytes=rf.requested),indent=2)+'\n')
        raise
    finally:rf.close()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    try:run(p.parse_args().output)
    except BaseException as error:
        # Never leak source row values through library exception rendering.
        print(json.dumps(dict(status='FAILED',error_class=type(error).__name__)))
        raise SystemExit(1)
