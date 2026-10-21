"""Independent DataFrame merge of every author state bridge opportunity."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,importlib.util,json,os
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);parent=Path(os.environ['PANCANCER_PARENT_OUTPUT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists()
    t=json.loads((parent/'private_terminal.json').read_text())
    for name,v in t['files'].items():assert sha(parent/name)==v['SHA256'] and (parent/name).stat().st_size==v['bytes']
    source=parent.parent/'probe_CM02_author_library_barcode_bridge_v1.py';assert sha(source)==t['source_SHA256']
    spec=importlib.util.spec_from_file_location('accepted_metadata_reader',source);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    aliases=json.loads((parent/'private_all27_reviewed_library_alias_opportunities.json').read_text())
    opportunities=json.loads((parent/'private_all18_library_barcode_opportunities.json').read_text());cross=json.loads((parent/'private_all_source_node_state_crosscounts.json').read_text());assert len(aliases)==27 and len(opportunities)==18
    table=root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    fields=list(dict.fromkeys(['library_id','patient_registry','timepoint_registry']+module.PATIENT_FIELDS+module.TIME_FIELDS))
    tab=pd.read_csv(table,sep='\t',usecols=fields,dtype=str,keep_default_na=False)
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');ds=np.asarray(read_elem(z['obs/dataset'])).astype(str);mask=ds=='GSE286827'
    cur=pd.DataFrame({k:np.asarray(read_elem(z['obs'][k])).astype(str)[mask] for k in ['library_id','patient_id','timepoint','major','node_id','_index','original_cell_id','parent_cell_id']});cur['current_position']=np.arange(len(cur))
    registry=tab.set_index('library_id');cohort_scope=[];checked=0;sets={}
    clean=lambda s:s.astype(str).str.strip().str.casefold()
    def extract(values,mode):
        text=values.astype(str);ok=text.str.findall(module.BARCODE).str.len().eq(1)
        parts=text.str.extract(module.BARCODE)
        result=parts[0] if mode=='within_reviewed_library_16nt' else parts[0]+parts[1].fillna('').map(lambda x:'-'+x if x else '')
        return result.where(ok)
    for name,label,expected_sha in module.SOURCES:
        a=module.read_author(name,label,expected_sha);a['author_position']=np.arange(len(a));a['p']=clean(a.patient);a['t']=clean(a.phase)
        chosen=next(r for r in opportunities if r['annotation_source']==name)
        for pf in module.PATIENT_FIELDS:
            for tf in module.TIME_FIELDS:
                lookup=tab[['library_id']].copy();lookup['p']=clean(tab[pf]);lookup['t']=clean(tab[tf]);sizes=lookup.groupby(['p','t']).library_id.nunique();good=sizes[sizes==1].reset_index()[['p','t']]
                lookup=lookup.merge(good,on=['p','t']).drop_duplicates(['p','t']);joined=a.merge(lookup,on=['p','t'],how='left',validate='many_to_one')
                record=next(r for r in aliases if r['annotation_source']==name and r['patient_field']==pf and r['time_field']==tf)
                assert record['bound_author_cells']==int(joined.library_id.notna().sum())
                assert record['uniquely_bound_groups']==len(joined.loc[joined.library_id.notna(),['p','t']].drop_duplicates())
                if pf==chosen['patient_alias_field'] and tf==chosen['time_alias_field']:bound=joined
        pair_hashes={}
        for mode in ['full_barcode_suffix','within_reviewed_library_16nt']:
            ar=bound.copy();ar['barcode']=extract(ar.cell_key,mode);author_bad=int(ar.barcode.isna().sum());ar=ar.dropna(subset=['library_id','barcode'])
            ac=int((ar.groupby(['library_id','barcode']).size()>1).sum());ar=ar[~ar.duplicated(['library_id','barcode'],keep=False)]
            for key in ['_index','original_cell_id','parent_cell_id']:
                c=cur.copy();c['barcode']=extract(c[key],mode);current_bad=int(c.barcode.isna().sum());c=c.dropna(subset=['barcode']);cc=int((c.groupby(['library_id','barcode']).size()>1).sum());c=c[~c.duplicated(['library_id','barcode'],keep=False)]
                m=c.merge(ar[['library_id','barcode','label','author_position']],on=['library_id','barcode'],how='inner',validate='one_to_one')
                old=next(r for r in opportunities if r['annotation_source']==name and r['barcode_mode']==mode and r['current_cell_key_field']==key)
                assert old['matched_cells']==len(m) and old['author_unparsable_cells']==author_bad and old['current_unparsable_cells']==current_bad
                assert old['author_collision_keys']==ac and old['current_collision_keys']==cc
                patient_bad=int((clean(m.patient_id)!=clean(m.library_id.map(registry.patient_registry))).sum());phase_bad=int((clean(m.timepoint)!=clean(m.library_id.map(registry.timepoint_registry))).sum())
                assert old['matched_current_patient_registry_disagreements']==patient_bad and old['matched_current_phase_registry_disagreements']==phase_bad
                got={(str(p),str(g),str(n),str(l)):int(v) for (p,g,n,l),v in m.groupby(['timepoint','major','node_id','label']).size().items()}
                expected={(r['phase'],r['current_major'],r['current_node'],r['author_label']):r['matched_cells'] for r in cross if r['annotation_source']==name and r['barcode_mode']==mode and r['current_cell_key_field']==key}
                assert got==expected
                pairs=m.sort_values('current_position')[['current_position','author_position']].to_numpy(dtype='<i8');pair_hashes[mode,key]=hashlib.sha256(pairs.tobytes()).hexdigest()
                if mode=='within_reviewed_library_16nt' and key=='original_cell_id':
                    sets[name]=set(m.current_position)
                    cohort_scope.append(dict(annotation_source=name,author_annotated_cells=len(a),matched_current_cells=len(m),current_patient_registry_disagreements=patient_bad,current_phase_registry_disagreements=phase_bad,full_author_population_or_functional_truth_not_accepted=True))
                checked+=1
        assert all(pair_hashes[mode,'_index']==pair_hashes[mode,'original_cell_id']==pair_hashes[mode,'parent_cell_id'] for mode in ['full_barcode_suffix','within_reviewed_library_16nt'])
    assert checked==18
    result=dict(status='ACCEPTED_ALL_LIBRARY_BARCODE_BRIDGE_OPPORTUNITIES_AND_SOURCE_NODE_STATE_CROSSCOUNTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),parent_terminal_SHA256=sha(parent/'private_terminal.json'),all18_join_positions_and_source_crosscounts_replayed=True,independent_Pandas_merge_and_regex_extraction=True,shared_pinned_R_metadata_decoder_not_second_RDS_parser=True,all_three_current_key_fields_give_identical_cell_pairings=True,source_scope=cohort_scope,CD4_CD8_shared_current_cells=len(sets[module.SOURCES[1][0]]&sets[module.SOURCES[2][0]]),individual_rows_keys_and_metadata_bytes_not_persisted=True,clinical_expression_and_sequence_values_not_selected=True,no_formal_label_adoption=True)
    out.mkdir(mode=0o700);(out/'private_independent_acceptance.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
