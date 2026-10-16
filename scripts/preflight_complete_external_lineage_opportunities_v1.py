"""All inherited major lineages and frozen candidate pools, no expression/effects."""
import hashlib,itertools,json,os,subprocess
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from external_lineage_qualification_inputs_v1 import load,STRATA


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def summarize(g):
    reviewed=g.reviewed_patient
    return dict(libraries=len(g),reviewed_patient_blocks=int(g.loc[reviewed,'patient_block'].nunique()),unreviewed_library_fallback=int((~reviewed).sum()))


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);data=load(r);lib=data['libraries'];counts=data['counts'];pools=data['pool_nodes'];out=r/'complete_external_lineage_opportunities_20261016_v1';assert not out.exists()
    df=subprocess.check_output(['df','-P',str(r.resolve())],text=True,timeout=20);fs=df.splitlines()[-1].split()[0];quota=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True,timeout=20);qr=[x.split() for x in quota.splitlines() if x.split() and x.split()[0]==fs];assert len(qr)==1;remaining=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert remaining>500_000_000_000+1_073_741_824+268_435_456
    out.mkdir(mode=0o700);plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_helper_sha256=sha(Path(__file__).with_name('external_lineage_qualification_inputs_v1.py')),input_sha256=data['input_sha256'],
        scope='All52 frozen-original annotation source entries/1290 libraries/172nodes and all9 archived immune source entries/400 libraries/270nodes. Source-arm entries and cell records are not automatically unique sources/patients/cells.',
        pools='Every inherited18 raw majors in each arm; explicit disjoint multigene CD8/Macrophage and Treg/Macrophage candidates in immune atlas; B+Plasma broad third-lineage extension. No labels merged or promoted, no exact NSCLC endpoint bridge.',
        gate='Both lineages>=20 cells in the same physical library. Report all zero source/pair positions and all inherited tissue/timepoint/treatment/capture/exvivo/stage strata. Counts are opportunities, not model identifiability or clinical acceptance.',
        patients='Latest verified patient blocks counted separately from library fallbacks; newer HNSCC continuous clinical definition and official TNBC sampling sidecar preserve original labels/conflicts. Never sum source-arm entries as independent cohorts.',
        no_expression_scores_P_or_new_clinical_acceptance=True,quota_free_bytes=remaining,output_cap_bytes=268_435_456)
    (out/'plan_before_complete_capture_counts.json').write_text(json.dumps(plan,indent=2)+'\n')
    lib.to_csv(out/'private_all_library_latest_metadata.tsv',sep='\t',index=False);counts.to_csv(out/'private_all_original_node_library_counts.tsv.gz',sep='\t',index=False,compression='gzip')
    definitions=[dict(arm=arm,pool=pool,node_id=node) for arm,pp in pools.items() for pool,nodes in pp.items() for node in sorted(nodes)];pd.DataFrame(definitions).to_csv(out/'private_all_pool_node_definitions.tsv',sep='\t',index=False)
    source_rows=[];strata_rows=[];candidate_rows=[];pool_rows=[];paired=[]
    for arm,ll in lib.groupby('arm',sort=True):
        cc=counts[counts.arm.eq(arm)];wide=cc.pivot_table(index=['dataset','library_id'],columns='node_id',values='n_cells',aggfunc='sum',fill_value=0).reindex(pd.MultiIndex.from_frame(ll[['dataset','library_id']]),fill_value=0)
        totals={pool:wide.reindex(columns=sorted(nodes),fill_value=0).sum(axis=1).to_numpy(int) for pool,nodes in pools[arm].items()};ll=ll.copy().reset_index(drop=True)
        for pool,values in totals.items():
            for i,row in ll.iterrows():pool_rows.append(dict(arm=arm,dataset=row.dataset,library_id=row.library_id,pool=pool,n_cells=int(values[i])))
        questions=[('major::'+a,'major::'+b,'inherited_raw_major') for a,b in itertools.combinations(data['raw_majors'][arm],2)]
        if arm=='archived_immune_atlas':questions += [('candidate::CD8','candidate::Macrophage','frozen_candidate_extension'),('candidate::Treg','candidate::Macrophage','frozen_candidate_extension')]
        for a,b,kind in questions:
            assert not pools[arm][a]&pools[arm][b];ok=(totals[a]>=20)&(totals[b]>=20)
            for source,indices in ll.groupby('dataset',sort=True).indices.items():
                gg=ll.iloc[indices];eligible=gg.loc[ok[indices]];record=dict(arm=arm,dataset=source,lineage_A=a,lineage_B=b,scope=kind,all_source_libraries=len(gg),**summarize(eligible),new_model_or_independence_accepted=False)
                if kind=='inherited_raw_major':source_rows.append(record)
                else:
                    flags=ok&(totals['extension::B_plus_Plasma']>=20);triple=gg.loc[flags[indices]];record.update(B_plus_Plasma_libraries=len(triple),B_plus_Plasma_reviewed_patient_blocks=summarize(triple)['reviewed_patient_blocks']);candidate_rows.append(record)
                    if source=='GSE286827':
                        for label,subset in [('two_end',eligible),('two_end_plus_B_Plasma',triple)]:
                            valid=subset[subset.reviewed_patient & subset.timepoint.isin(['pre','post'])];phases=valid.groupby('patient_block').timepoint.agg(set);paired.append(dict(dataset=source,lineage_A=a,lineage_B=b,scope=label,reviewed_pre_post_patient_pairs=int(phases.map(lambda v:{'pre','post'}<=v).sum()),pre_libraries=int(valid.timepoint.eq('pre').sum()),post_libraries=int(valid.timepoint.eq('post').sum()),other_phase_libraries=int((~subset.timepoint.isin(['pre','post'])).sum())))
                for levels,sg in gg.groupby(STRATA,sort=True,dropna=False):
                    eg=sg.loc[ok[sg.index.to_numpy()]];strata_rows.append(dict(arm=arm,dataset=source,lineage_A=a,lineage_B=b,scope=kind,**dict(zip(STRATA,levels)),all_stratum_libraries=len(sg),**summarize(eg),new_model_or_independence_accepted=False))
    for name,rows in [('private_all_library_pool_counts.tsv.gz',pool_rows),('private_all_source_raw_major_pairs.tsv',source_rows),('private_all_source_stratum_pairs.tsv.gz',strata_rows),('private_all_frozen_candidate_pair_opportunities.tsv',candidate_rows),('private_all_HNSCC_reviewed_paired_opportunities.tsv',paired)]:pd.DataFrame(rows).to_csv(out/name,sep='\t',index=False,compression='gzip' if name.endswith('.gz') else None)
    assert sum(p.stat().st_size for p in out.iterdir())<plan['output_cap_bytes']
    terminal=dict(status='COMPLETE_ALL_FROZEN_AND_IMMUNE_MAJOR_CANDIDATE_CAPTURE_COUNTS_LATEST_SIDECARS_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),source_arm_entries=int(lib.groupby(['arm','dataset']).ngroups),physical_library_arm_entries=len(lib),raw_major_pair_rows=len(source_rows),stratum_rows=len(strata_rows),candidate_pair_rows=len(candidate_rows),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_expression_or_scores_or_P=False,new_identity_or_clinical_acceptance=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
