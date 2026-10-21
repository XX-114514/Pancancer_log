"""Pinned public single-member RAM input and complete marked-space family.

Only source/cohort aggregates leave RAM. No patient, ROI, cell, coordinate,
score, graph or permutation arrays are written. Original atlas is untouched.
"""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import argparse,csv,hashlib,io,json,os,resource,time,zipfile
import numpy as np
import pandas as pd
from Danenberg_bounded_member_transport_v1 import RangeFile,IntegrityError
from qualify_Danenberg_native_cells_v1 import (URL,SIZE,DIRECTORY_SHA,HEADER_SHA,MARKERS,
    HashReader,member,directory,require,qc_mask,author_alias)
from qualify_Danenberg_clinical_metadata_v1 import binary,clean
from CM02_marked_spatial_core_v1 import (AXES_X,AXES_Y,rankz,marker_residuals,graph_edges,
    weighted_corr,permutation_groups,permutations,permuted_corr,paired_roi_contrast)
from CM02_marked_spatial_case_inference_v2 import cohort_summary,bh_fixed_family

MEMBER_SHA='9f693514acc02bbb6b2dabd14c15c806348c5b882659aea071a83ed16e2d5df0'
SEED=20261022;B=2000
ROLE_LABELS={'RegMixed':'Treg & Tex','CD4control':'CD4+ T cells','CD8':'CD8+ T cells'}
CONFIGS=[('primary25_grid2',25,2,False),('sensitivity50_grid2',50,2,False),
         ('sensitivity25_grid4',25,4,False),('sensitivity25_grid2_interfaces',25,2,True)]


def now():return datetime.now(timezone.utc).isoformat()


def guard_memory():require(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<8*1024**3,'RAM_budget_exceeded')


def slot_name(parent,k,l):return parent+'__'+AXES_X[k]+'__CD8_'+AXES_Y[l]


def main(output,qualification,frozen,calibration):
    require(not output.exists(),'refuse_overwrite')
    plan=json.loads(frozen.read_text());require(plan['status']=='FROZEN_IMPLEMENTATION_AND_CHECKS_BEFORE_SCIENCE','implementation_not_frozen')
    root=Path(__file__).parent
    for name,sha in plan['code_SHA256'].items():require(hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,'frozen_code_changed')
    calibration_result=json.loads(calibration.read_text())
    require(calibration_result['status']=='COMPLETE_EXTENDED_SYNTHETIC_CASE_NULL_CALIBRATION' and
        all(x['patient_signflip_calibration_flag']=='QUALIFIED_EXTENDED_SYNTHETIC_SCOPE' for x in calibration_result['results']),
        'case_null_calibration_not_qualified')
    previous=json.loads(qualification.read_text());require(previous['selected_member_SHA256']==MEMBER_SHA and previous['full_selected_member_CRC_and_bytes'],'native_qualification_not_accepted')
    progress=output.with_name(output.stem+'_progress.json');started=now();phase=['INPUT'];last=[0]
    def checkpoint(values=None,force=False):
        if force or time.monotonic()-last[0]>10:
            progress.write_text(json.dumps(dict(status='RUNNING_MARKED_SPATIAL',phase=phase[0],started_UTC=started,
                observation_UTC=now(),**(values or {})),indent=2)+'\n');last[0]=time.monotonic()
    rf=RangeFile(URL,SIZE,int(os.environ['PANCANCER_EXISTING_CLASH_PORT']),checkpoint)
    checkpoint(force=True)
    try:
        with zipfile.ZipFile(rf) as z:
            require(directory(z)==DIRECTORY_SHA,'directory_changed')
            clinical_raw=z.read(member(z,'IMCClinical.csv',51316,1165356369))
            rows=list(csv.DictReader(io.StringIO(clinical_raw.decode('utf-8-sig'))))
            clinical={clean(r['metabric_id']):r for r in rows}
            require(len(rows)==len(clinical) and '' not in clinical,'clinical_key_conflict')
            aliases={key:role for role,label in ROLE_LABELS.items() for key in (label,author_alias(label))}
            m=member(z,'SingleCells.csv',849239430,3748430769);require(m.compress_size==381180527,'compressed_size_changed')
            pieces=[];bounds={};QC=0;total=0;epithelial_conflict=0;unknown_interfaces=Counter()
            numeric=MARKERS+['Location_Center_X','Location_Center_Y','AreaShape_Area']
            with z.open(m) as stream:
                hashed=HashReader(stream);buffer=io.BufferedReader(hashed,buffer_size=1024**2)
                line=buffer.readline();require(hashlib.sha256(line).hexdigest()==HEADER_SHA,'header_changed')
                names=next(csv.reader([line.decode('utf-8-sig').rstrip('\r\n')]))
                for frame in pd.read_csv(buffer,names=names,header=None,dtype=str,keep_default_na=False,chunksize=20000,on_bad_lines='error'):
                    total+=len(frame);good,values,reasons=qc_mask(frame);QC+=int(good.sum())
                    kept=frame.loc[good]
                    for image,index in kept.groupby('ImageNumber',sort=False).groups.items():
                        xy=values.loc[index,['Location_Center_X','Location_Center_Y']].to_numpy();lo=xy.min(axis=0);hi=xy.max(axis=0)
                        if image in bounds:lo=np.minimum(lo,bounds[image][0]);hi=np.maximum(hi,bounds[image][1])
                        bounds[image]=(lo,hi)
                    roles=kept.cellPhenotype.str.strip().map(aliases);target=roles.notna()
                    part=kept.loc[target,['ImageNumber','ObjectNumber','metabric_id','is_epithelial','is_interface','is_perivascular']].copy()
                    part['role']=roles.loc[target]
                    epithelial_conflict+=int(part.is_epithelial.map(binary).ne(False).sum())
                    for f in ['is_interface','is_perivascular']:
                        part[f]=part[f].map(binary).map({True:1,False:0}).fillna(-1).astype(int)
                        unknown_interfaces[f]+=int(part[f].eq(-1).sum())
                    part=part.drop(columns='is_epithelial')
                    for c in numeric:part[c]=values.loc[part.index,c]
                    pieces.append(part);guard_memory()
                require(buffer.read()==b'' and hashed.bytes==m.file_size,'incomplete_member')
            require(hashed.sha.hexdigest()==MEMBER_SHA,'full_member_SHA_changed')
            require(total==previous['source_cells'] and QC==previous['QC_cells'],'native_count_reconciliation_failed')
            require(epithelial_conflict==0,'author_target_label_epithelial_conflict')
            with zipfile.ZipFile(rf) as replay:require(directory(replay)==DIRECTORY_SHA,'fresh_directory_changed')
        rf.close();data=pd.concat(pieces,ignore_index=True);del pieces
        require(not data.duplicated(['ImageNumber','ObjectNumber']).any(),'target_compound_key_duplicate')
        case_counts=data.groupby(['metabric_id','role']).size().unstack(fill_value=0)
        eligible=set(case_counts.index[(case_counts.RegMixed>=20)&(case_counts.CD8>=20)])&set(clinical)
        expected=sum(r['MIXED_REG_TEX_AND_CD8_AT_LEAST_20_EACH'] for r in previous['by_author_validation_flag'])
        require(len(eligible)==expected,'eligible_case_reconciliation_failed')
        data=data.loc[data.metabric_id.isin(eligible)].reset_index(drop=True)
        scope=dict(clinical_linked_target_eligible_cases=len(eligible),
            control_parent_also_at_least20_cases=int(sum(case_counts.loc[k,'CD4control']>=20 for k in eligible)),
            unknown_target_interface_flags=dict(unknown_interfaces),
            no_author_T_target_epithelial_conflicts=True)
        case_metadata={k:dict(author_flag=binary(clinical[k]['isValidation']),
             ER=clean(clinical[k]['ERStatus']).lower(),event=binary(clinical[k]['DeathBreast'])) for k in eligible}
        phase[0]='FIXED_MARKER_AND_GRAPH_QUALIFICATION';checkpoint(force=True)
        # No per-case outputs are written. Resident cells are used for all fixed configurations.
        roi_inputs=[];component_values=defaultdict(list);reasons_global=Counter();numeric_audit=Counter()
        for image,frame in data.groupby('ImageNumber',sort=True):
            case=frame.metabric_id.unique();require(len(case)==1,'image_case_conflict');case=case[0]
            parents={role:part.reset_index(drop=True) for role,part in frame.groupby('role',sort=True)}
            residuals={};status={}
            for role,part in parents.items():
                residuals[role],status[role]=marker_residuals(part,role)
                for index,reason in enumerate(status[role]):
                    reasons_global[role+'__'+reason]+=1
                    if reason=='RESIDUAL_DF_BELOW_FIXED_MINIMUM':
                        is_cycle=(role=='CD8' and index==0) or (role!='CD8' and index==2)
                        require(len(part)-(5 if is_cycle else 6)<8,'real_DF_gate_replay_failed')
                        numeric_audit['REAL_UNESTIMABLE_DF_COUNT_REPLAY']+=1
                if role!='CD8':
                    for a,b in [('FOXP3','CD4'),('ICOS','OX40'),('ICOS','GITR'),('OX40','GITR')]:
                        ra=rankz(part[a]);rb=rankz(part[b]);v=float(np.mean(ra*rb))
                        if np.isfinite(v):component_values[(case,role,a+'__'+b)].append(v)
            roi_inputs.append((case,image,parents,residuals,status))
        results={};configuration_coverage={}
        rng=np.random.default_rng(SEED)
        for config,radius,grid,interfaces in CONFIGS:
            phase[0]=config;checkpoint(force=True)
            # Each slot maps case -> observed ROI list and complete null ROI list.
            slots=defaultdict(lambda:defaultdict(list));missing=Counter();graph_cov=Counter()
            for case,image,parents,residuals,status in roi_inputs:
                parent_permutations={};parent_groups={}
                for role,part in parents.items():
                    groups,exchangeable=permutation_groups(part,bounds[image],grid,interfaces)
                    parent_groups[role]=(groups,exchangeable)
                    if len(part)>=10 and exchangeable>=10:
                        parent_permutations[role]=permutations(groups,len(part),B,rng)
                for parent in ['RegMixed','CD4control']:
                    if parent=='CD4control' and case_counts.loc[case,'CD4control']<20:
                        missing['CONTROL_PARENT_CASE_CAPTURE_BELOW20']+=1;continue
                    if parent not in parents or 'CD8' not in parents:
                        missing['PARENT_NOT_CAPTURED_IN_ROI']+=1;continue
                    source=parents[parent];target=parents['CD8']
                    if min(len(source),len(target))<10:missing['ROI_PARENT_BELOW10']+=1;continue
                    if parent not in parent_permutations or 'CD8' not in parent_permutations:
                        missing['EXCHANGEABLE_PARENT_CELLS_BELOW10']+=1;continue
                    graph,reason=graph_edges(source[['Location_Center_X','Location_Center_Y']],target[['Location_Center_X','Location_Center_Y']],radius)
                    if graph is None:missing[reason]+=1;continue
                    graph_cov[parent+'_eligible_ROIs']+=1;graph_cov[parent+'_heterotypic_edges']+=len(graph[0])
                    observed=weighted_corr(residuals[parent],residuals['CD8'],graph)
                    null=np.full((B,3,2),np.nan)
                    for start in range(0,B,32):
                        null[start:start+32]=permuted_corr(residuals[parent],residuals['CD8'],graph,
                            parent_permutations[parent][start:start+32],parent_permutations['CD8'][start:start+32])
                    for sign,name in [(1,'POSITIVE_REAL_EDGE_REPLAY'),(-1,'NEGATIVE_REAL_EDGE_REPLAY')]:
                        candidates=np.argwhere(np.isfinite(observed)&(observed*sign>0))
                        if len(candidates) and not numeric_audit[name]:
                            k0,l0=candidates[0];i,j,w=graph
                            a=residuals[parent][i,k0];b=residuals['CD8'][j,l0]
                            am=sum(float(weight)*float(value) for weight,value in zip(w,a))
                            bm=sum(float(weight)*float(value) for weight,value in zip(w,b))
                            cov=sum(float(weight)*(float(x)-am)*(float(y)-bm) for weight,x,y in zip(w,a,b))
                            va=sum(float(weight)*(float(x)-am)**2 for weight,x in zip(w,a))
                            vb=sum(float(weight)*(float(y)-bm)**2 for weight,y in zip(w,b))
                            require(np.isclose(cov/np.sqrt(va*vb),observed[k0,l0],atol=1e-10),'independent_real_edge_replay_failed')
                            for draw in range(3):
                                replay=weighted_corr(residuals[parent][parent_permutations[parent][draw]],
                                    residuals['CD8'][parent_permutations['CD8'][draw]],graph)
                                require(np.allclose(replay,null[draw],equal_nan=True,atol=1e-10),'independent_null_batch_replay_failed')
                            numeric_audit[name]+=1
                    for k in range(3):
                        for l in range(2):
                            key=slot_name(parent,k,l)
                            if not np.isfinite(observed[k,l]) or not np.isfinite(null[:,k,l]).all():
                                missing[key+'__'+(status[parent][k] if status[parent][k]!='ESTIMABLE' else status['CD8'][l] if status['CD8'][l]!='ESTIMABLE' else 'DEGENERATE_EDGE_OR_PERMUTED_VARIANCE')]+=1
                                continue
                            slots[key][case].append((image,float(observed[k,l]),null[:,k,l].copy()))
                guard_memory();require(time.monotonic()-rf.started<=7200,'TOTAL_WALL_BUDGET_EXCEEDED')
                checkpoint(dict(completed_configuration=config))
            cases={key:{case:(float(np.mean([x[1] for x in values])),np.mean([x[2] for x in values],axis=0))
                        for case,values in by_case.items()} for key,by_case in slots.items()}
            summary={}
            for parent in ['RegMixed','CD4control']:
                for k in range(3):
                    for l in range(2):
                        key=slot_name(parent,k,l);values=cases.get(key,{})
                        summary[key]=summarize_cohorts(values,case_metadata,config=='primary25_grid2' and k<2)
            for k in range(2):
                for l in range(2):
                    a=slots.get(slot_name('RegMixed',k,l),{});b=slots.get(slot_name('CD4control',k,l),{})
                    contrast={}
                    for case in sorted(set(a)&set(b)):
                        paired=paired_roi_contrast(a[case],b[case])
                        if paired is not None:contrast[case]=paired
                    summary[slot_name('RegMinusCD4',k,l)]=summarize_cohorts(contrast,case_metadata,config=='primary25_grid2')
            if config=='primary25_grid2':
                main={key:v['ALL_ELIGIBLE'] for key,v in summary.items() if 'Ki67_reference' not in key}
                require(len(main)==12,'primary_family_slot_count');bh_fixed_family(main,12)
            results[config]=summary;configuration_coverage[config]=dict(graph_coverage=dict(graph_cov),ROI_exclusion_reasons=dict(missing))
        components=[]
        for role in ['RegMixed','CD4control']:
            for pair in ['FOXP3__CD4','ICOS__OX40','ICOS__GITR','OX40__GITR']:
                values=[float(np.mean(v)) for (case,r,p),v in component_values.items() if r==role and p==pair]
                components.append(dict(parent=role,marker_pair=pair,n_cases=len(values),
                    mean_within_ROI_rank_correlation=float(np.mean(values)) if values else None,
                    fraction_positive=float(np.mean(np.asarray(values)>0)) if values else None))
        phase[0]='FINAL_SOURCE_AGGREGATE_AUDIT';checkpoint(force=True)
        require(len(results['primary25_grid2'])==16,'full_family_and_reference_slots')
        for config,summary in results.items():
            require(len(summary)==16,'configuration_slots_incomplete')
            for key,cohorts in summary.items():
                for group,item in cohorts.items():
                    if config!='primary25_grid2' or group!='ALL_ELIGIBLE' or 'Ki67_reference' in key:
                        require('patient_signflip_P' not in item,'unplanned_test')
                    if 'BH_fixed_family' in item:require(0<=item['BH_fixed_family']<=1,'BH_bounds')
        result=dict(status='COMPLETE_FIXED_CONTINUOUS_MARKER_SPATIAL_CASE_FAMILY',started_UTC=started,
            completed_UTC=now(),source_scope=scope,selected_member_SHA256=hashed.sha.hexdigest(),full_member_CRC=True,
            fresh_directory_SHA256=DIRECTORY_SHA,complete_primary_BH_family_size=12,permutations=B,bootstrap=2000,seed=SEED,
            marker_measurement_design_reasons=dict(reasons_global),component_support=components,
            independent_real_positive_negative_and_null_replay=dict(numeric_audit),
            configuration_coverage=configuration_coverage,results=results,
            full_native_count_and_target_case_reconciliation=True,all_fixed_slots_and_sensitivities_reported=True,
            randomization='Cell residual-vector shuffle defines the spatial baseline only, with no inferential P. Whole-case sign-flip P requires independent cases and symmetric null case effects.',
            synthetic_case_calibration_SHA256=hashlib.sha256(calibration.read_bytes()).hexdigest(),
            case_signflips=2000,cell_residual_P_excluded_after_failed_null_calibration=True,
            no_pure_Treg_or_exhaustion_or_LR_or_exact_CM02_projection_claim=True,
            no_clinical_model_or_heldout_validation_claim=True,no_individual_or_source_payload_files=True,
            no_global_phenotype_mapping_completion_claim=True,response_identity_headers_available=rf.identity!=(None,None),
            range_requested_bytes=rf.requested,range_successful_bytes=rf.transferred,range_failures=rf.failures,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        progress.write_text(json.dumps(dict(status=result['status'],completed_UTC=result['completed_UTC']),indent=2)+'\n')
    except BaseException as error:
        progress.write_text(json.dumps(dict(status='FAILED_MARKED_SPATIAL',phase=phase[0],error_class=type(error).__name__,
            qualification_reason=str(error) if isinstance(error,IntegrityError) else 'EXTERNAL_ERROR_VALUES_WITHHELD',
            observation_UTC=now()),indent=2)+'\n');raise
    finally:rf.close()


def summarize_cohorts(values,metadata,inference):
    groups={'ALL_ELIGIBLE':set(metadata),
        'AUTHOR_VALIDATION_FALSE':{k for k,v in metadata.items() if v['author_flag'] is False},
        'AUTHOR_VALIDATION_TRUE':{k for k,v in metadata.items() if v['author_flag'] is True},
        'ER_POSITIVE':{k for k,v in metadata.items() if v['ER'] in ('positive','pos','+','1','true')},
        'ER_NEGATIVE':{k for k,v in metadata.items() if v['ER'] in ('negative','neg','-','0','false')}}
    result={}
    for group,members in groups.items():
        keys=sorted(set(values)&members);observed=[values[k][0] for k in keys]
        null=np.asarray([values[k][1] for k in keys])
        result[group]=cohort_summary(observed,null,keys,inference=inference and group=='ALL_ELIGIBLE')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--calibration',type=Path,required=True);p.add_argument('--qualification',type=Path,required=True);p.add_argument('--frozen',type=Path,required=True);a=p.parse_args()
    try:main(a.output,a.qualification,a.frozen,a.calibration)
    except BaseException as error:
        print(json.dumps(dict(status='FAILED',error_class=type(error).__name__)));raise SystemExit(1)
