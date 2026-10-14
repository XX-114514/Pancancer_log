#!/usr/bin/env python3
"""Bounded per-library CNV rerun, exact merge, then authority archive gate."""
import argparse,hashlib,json,os,subprocess,time,traceback
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem
import zarr


def now():return datetime.now(timezone.utc).isoformat()


def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)


def main(a):
    r=a.run;s=a.scratch;c=a.code_dir;m=r/'manifests';version='logspace_v1'
    save(m/'CNV_logspace_v1_launch_private.json',dict(started_utc=now(),pid=os.getpid(),code_git_sha=a.code_git_sha,
        script_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in c.iterdir() if p.suffix in ['.py','.R']},
        original_library_workers=a.workers,infercnv_subprocesses_per_library=2,old_raw_integer_CNVs_not_reused=True))
    env=os.environ.copy()
    for name in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:env[name]='2'
    env['NUMBA_NUM_THREADS']='4';env['MPLBACKEND']='Agg'
    def stage(name,command):
        save(m/(name+'_logspace_v1_started_private.json'),dict(started_utc=now(),command=command))
        with (r/'logs'/(name+'_logspace_v1.log')).open('xb') as f:
            result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env)
        save(m/(name+'_logspace_v1_terminal_public.json'),dict(finished_utc=now(),exit_code=result.returncode))
        if result.returncode:raise RuntimeError(name+' failed; new attempt preserved')
    try:
        stage('core_build',[a.python,str(c/'build_immune_cnv_logspace_core.py'),'--original',str(a.core_dir/'malignant_inference.py'),
            '--output',str(c/'corrected_malignant_inference.py'),'--report',str(m/'CNV_logspace_core_patch_public.json')])
        coords=pd.read_csv(a.gene_annotation,sep=r'\s+',header=None,names=['gene_symbol','chromosome','start','end'])
        coords=coords.loc[coords.chromosome.str.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])')].drop_duplicates('gene_symbol')
        path=r/'config/CNV_autosomal_coordinate_contract_logspace_v1.tsv'
        if path.exists():raise ValueError('Existing new coordinate attempt preserved')
        coords.to_csv(path,sep='\t',index=False,header=False)
        children=[];streams=[]
        for i in range(a.workers):
            command=[a.python,'-u',str(c/'immune_native_cnv.py'),'--contract',str(r/'config/source_contract_v2_private.json'),
                '--run',str(r),'--scratch',str(s),'--core-dir',str(a.core_dir),'--gene-annotation',str(a.gene_annotation),
                '--corrected-core',str(c/'corrected_malignant_inference.py'),'--attempt-version',version,
                '--shard-index',str(i),'--shard-count',str(a.workers)]
            stream=(r/'logs'/('CNV_logspace_v1_shard_'+str(i)+'.log')).open('xb');streams.append(stream)
            child=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,env=env);children.append(child)
            save(m/('CNV_logspace_v1_shard_'+str(i)+'_start_private.json'),dict(started_utc=now(),pid=child.pid,command=command))
        while any(x.poll() is None for x in children):time.sleep(20)
        codes=[x.returncode for x in children]
        for stream in streams:stream.close()
        save(m/'CNV_logspace_v1_workers_terminal_public.json',dict(finished_utc=now(),exit_codes=codes))
        if any(codes):raise RuntimeError('Corrected CNV shard failed; preserve completed outputs for recovery')
        obs=read_elem(zarr.open_group(str(s/'immune_refined.zarr'),mode='r')['obs'])
        parts=[];libraries=[];summaries=[]
        for i in range(a.workers):
            suffix=version+'_shard_'+str(i)
            summaries.append(json.loads((m/('CNV_'+suffix+'_public.json')).read_text()))
            parts.append(pd.read_csv(r/'tables'/('all_core_CNV_'+suffix+'_private.tsv.gz'),sep='\t',index_col=0,keep_default_na=False))
            libraries.append(pd.read_csv(r/'tables'/('CNV_library_method_index_'+suffix+'_private.tsv'),sep='\t'))
        joined=pd.concat(parts);library=pd.concat(libraries,ignore_index=True)
        if not joined.index.is_unique or set(joined.index)!=set(obs.index):raise ValueError('Corrected CNV cell coverage/uniqueness differs')
        joined=joined.reindex(obs.index)
        if not np.array_equal(joined.sampleID.astype(str),obs.sampleID.astype(str)):raise ValueError('Corrected CNV sample mapping differs')
        if len(library)!=400 or library.duplicated(['dataset','library_token']).any():raise ValueError('Corrected physical library coverage differs')
        if library.status.eq('failed').any():raise ValueError('Corrected CNV method failures remain')
        destination=r/'tables/all_core_CNV_logspace_v1_private.tsv.gz'
        if destination.exists():raise ValueError('Corrected combined outputs preserved')
        joined.to_csv(destination,sep='\t');library.to_csv(r/'tables/CNV_library_method_index_logspace_v1_private.tsv',sep='\t',index=False)
        mapping=obs[['dataset','sampleID']].copy();mapping['cell_id']=obs.index
        mapping['mapped']=joined.sampleID.to_numpy()==obs.sampleID.astype(str).to_numpy()
        mapping.groupby(['dataset','sampleID'],observed=True).agg(n_input_cells=('cell_id','size'),n_mapped=('mapped','sum')).to_csv(r/'tables/CNV_exact_library_mapping_logspace_v1_private.tsv',sep='\t')
        public=dict(observed_utc=now(),n_cells=len(joined),n_libraries=len(library),all_cell_IDs_exact=True,
            all_sample_IDs_exact=True,duplicate_or_unmapped_cells=0,library_method_status_counts=library.status.value_counts().to_dict(),
            cell_malignant_status_counts=joined.malignant_status.value_counts().to_dict(),raw_counts_preserved=True,
            inference_input='sample_native_raw_counts_to_log1p_CP10K_full_native_denominator',same_original_library_references=True,
            infercnv_process_workers=2,library_workers=a.workers,old_raw_integer_CNVs_scientifically_quarantined=True,
            n_method_failed=0,malignancy_gold=False)
        save(m/'CNV_logspace_v1_public.json',public)
        CM_index=m/'CM_recovery_v3_public.json'
        if not json.loads(CM_index.read_text())['all_ready_scopes_completed']:raise ValueError('CM prerequisite not accepted')
        stage('authority_pack',[a.python,'-u',str(c/'pack_immune_authority.py'),'--contract',str(r/'config/source_contract_v2_private.json'),
            '--run',str(r),'--scratch',str(s),'--quota-group',a.quota_group,'--CM-discovery-index',str(CM_index),
            '--CNV-results',str(destination)])
        proof=json.loads((m/'immune_authority_release_public.json').read_text())
        if not proof['technical_authority_released'] or not proof['all_checks_pass']:raise RuntimeError('Authority artifact gate failed')
        save(m/'CNV_logspace_v1_pipeline_terminal_public.json',dict(finished_utc=now(),status='completed',n_cells=proof['n_cells'],n_features=proof['n_features'],technical_authority_released=True))
    except Exception as exc:
        (r/'logs/CNV_logspace_v1_pipeline_traceback_private.txt').write_text(traceback.format_exc())
        save(m/'CNV_logspace_v1_pipeline_terminal_public.json',dict(finished_utc=now(),status='failed',failure=str(exc),history_preserved=True,authoritative_release=False))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['run','scratch','code-dir','core-dir','gene-annotation']:p.add_argument('--'+name,type=Path,required=True)
    for name in ['python','code-git-sha','quota-group']:p.add_argument('--'+name,required=True)
    p.add_argument('--workers',type=int,default=4)
    args=p.parse_args()
    if not 1<=args.workers<=4:raise ValueError('Expected bounded 1 to 4 original-library workers')
    main(args)
