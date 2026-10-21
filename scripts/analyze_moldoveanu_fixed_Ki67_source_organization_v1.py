"""Three frozen source-conditional Ki67 comparisons, not patient inference.

Repeated/unreviewed patient identity does not become independent sample size.
Inference conditions on every observed material's fixed coordinates and parent
identities. No patient CI/bootstrap, clinical endpoint, new annotation or causal
claim. Shared CD8 shuffles preserve dependence for the direct contrast.
"""
from pathlib import Path
import csv
import hashlib
import importlib.util
import io
import json
import numpy as np
import pandas as pd

SEED=2026102102
PERMUTATIONS=2000
MIN_PARENT=20
MIN_EDGES=20
MIN_MATERIALS=5


def load(path,expected,name):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==expected
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def covariance(x,y,a,b):
    value=float(np.mean(x[a]*y[b])-np.mean(x[a])*np.mean(y[b]))
    mx=float(x[a].mean());my=float(y[b].mean())
    independent=sum((float(u)-mx)*(float(v)-my) for u,v in zip(x[a],y[b]))/len(a)
    assert abs(value-independent)<=1e-9
    return value


def means_by_group(values,groups):
    out=np.empty_like(values)
    for group in np.unique(groups):out[groups==group]=values[groups==group].mean()
    return out


def shuffled(values,groups,rng):
    out=values.copy()
    for group in np.unique(groups):
        take=np.flatnonzero(groups==group);out[take]=values[rng.permutation(take)]
    return out


def bh(p):
    order=np.argsort(p,kind='stable');q=np.ones(len(p));bound=1.
    for j in range(len(p)-1,-1,-1):
        k=order[j];bound=min(bound,p[k]*len(p)/(j+1));q[k]=bound
    return q


def finish(records):
    rows=[];probabilities=[]
    for slot,label in enumerate(['regulatory_Ki67_CD8_Ki67','conventional_CD4_Ki67_CD8_Ki67','regulatory_minus_conventional_shared_materials']):
        selected=[r for r in records if r['valid'][slot]]
        row=dict(comparison=label,source_materials=len(selected),patient_identity_reviewed=False)
        if len(selected)<MIN_MATERIALS:
            row.update(status='INSUFFICIENT_SOURCE_MATERIALS',effect=None,P=None,BH=None)
            probabilities.append(1.);rows.append(row);continue
        effects=[r['observed'][slot]-r['expected'][slot] for r in selected]
        point=sum(effects)/len(effects)
        null=np.mean(np.stack([r['draws'][:,slot]-r['expected'][slot] for r in selected]),axis=0)
        # Full independent scalar replay of every source permutation mean.
        for j in range(PERMUTATIONS):
            ordinary=sum(float(r['draws'][j,slot]-r['expected'][slot]) for r in selected)/len(selected)
            assert abs(float(null[j])-ordinary)<=1e-9
        p=(1+int((np.abs(null)>=abs(point)-1e-12).sum()))/(PERMUTATIONS+1)
        row.update(status='ESTIMABLE_FIXED_MATERIAL_CONDITIONAL_ORGANIZATION',effect=float(point),P=p,BH=None,
            P_scope='Exploratory material-unit conditional randomization; duplicate/patient independence unreviewed',
            positive_materials=int(sum(x>0 for x in effects)),negative_materials=int(sum(x<0 for x in effects)),
            leave_one_material_effect_range=[float((sum(effects)-max(effects))/(len(effects)-1)),float((sum(effects)-min(effects))/(len(effects)-1))],
            patient_bootstrap_95_CI=None,patient_bootstrap_status='NOT_ESTIMABLE_UNREVIEWED_PATIENT_BRIDGE',
            clinical_inference_not_performed=True)
        probabilities.append(p);rows.append(row)
    for row,q in zip(rows,bh(probabilities)):
        if row['P'] is not None:row['BH']=float(q)
    return rows


def source_statistics(frame,knn,grid,rank):
    roles=[['Treg'],['Th.ae','Th.naive'],['Tc.ae','Tc.naive']]
    outputs={mode:[] for mode in ['ROI_parent_shuffle','ROI_parent_grid4_shuffle']}
    rngs={mode:np.random.default_rng(SEED+offset) for mode,offset in [('ROI_parent_shuffle',0),('ROI_parent_grid4_shuffle',100000)]}
    scalar_covariances=0
    for _,part in frame.groupby('sample.id',sort=True):
        part=part.reset_index(drop=True)
        positions=[np.flatnonzero(part['Cluster'].isin(labels).to_numpy()) for labels in roles]
        if len(positions[2])<MIN_PARENT or max(len(positions[0]),len(positions[1]))<MIN_PARENT:continue
        coords=part[['coord.x','coord.y']].to_numpy(dtype=float);edges=knn(coords,10)
        pairs=[];values=[];valid=[]
        for context in range(2):
            a,b=edges[:,0],edges[:,1]
            forward=np.isin(a,positions[context])&np.isin(b,positions[2])
            reverse=np.isin(b,positions[context])&np.isin(a,positions[2])
            left=np.r_[a[forward],b[reverse]];right=np.r_[b[forward],a[reverse]]
            inverse_left=np.full(len(part),-1,dtype=int);inverse_right=inverse_left.copy()
            inverse_left[positions[context]]=np.arange(len(positions[context]));inverse_right[positions[2]]=np.arange(len(positions[2]))
            pairs.append((inverse_left[left],inverse_right[right]))
            scores,ok=rank(part.iloc[positions[context]][['Ki67']].to_numpy()) if len(positions[context]) else (np.empty((0,1)),np.array([False]))
            values.append(scores[:,0]);valid.append(len(positions[context])>=MIN_PARENT and len(left)>=MIN_EDGES and bool(ok[0]))
        right,ok=rank(part.iloc[positions[2]][['Ki67']].to_numpy());right=right[:,0]
        valid=[v and bool(ok[0]) for v in valid]
        if not any(valid):continue
        extent=(coords.min(0),coords.max(0))
        for mode in outputs:
            groups=[grid(coords[p],extent) if 'grid4' in mode else np.zeros(len(p),dtype=int) for p in positions]
            observed=np.zeros(3);expected=np.zeros(3);draws=np.zeros((PERMUTATIONS,3))
            for context in range(2):
                if valid[context]:
                    a,b=pairs[context]
                    observed[context]=covariance(values[context],right,a,b)
                    expected[context]=covariance(means_by_group(values[context],groups[context]),means_by_group(right,groups[2]),a,b)
                    scalar_covariances+=2
            for j in range(PERMUTATIONS):
                # Same CD8 permutation used for both roles and their contrast.
                perm_right=shuffled(right,groups[2],rngs[mode])
                for context in range(2):
                    if valid[context]:
                        a,b=pairs[context]
                        draws[j,context]=covariance(shuffled(values[context],groups[context],rngs[mode]),perm_right,a,b)
                        scalar_covariances+=1
            if all(valid):
                observed[2]=observed[0]-observed[1];expected[2]=expected[0]-expected[1]
                draws[:,2]=draws[:,0]-draws[:,1]
            outputs[mode].append(dict(valid=valid+[all(valid)],observed=observed,expected=expected,draws=draws))
    return {mode:finish(records) for mode,records in outputs.items()},scalar_covariances


def analyze(bundle,parameters):
    assert parameters=={'scope':'fixed_author_materials_no_patient_inference','permutations':2000,'seed':SEED}
    archive=bundle['archive'];assert hasattr(archive,'_own_ICI_qualified_cells')
    frame=archive._own_ICI_qualified_cells
    folder=Path(__file__).parent
    qualifier=load(folder/'qualify_moldoveanu_ICI_cells_in_RAM_v1.py',
        'a498faeb15053036e72003585ba82d5c9b44b2f8f9d3c16d392c3fcf8e693125','accepted_cell_reader')
    # Independently compare the published classes in both complete source tables.
    intensity=qualifier.decoded(qualifier.table_bytes(archive,'ICI_meanIntensity.tsv.gz'))
    labels=pd.read_csv(io.BytesIO(intensity),sep='\t',usecols=['Cell_ID','Sample_ID','Cluster.v2'],
        dtype={'Cell_ID':str,'Sample_ID':str})
    checked=frame.merge(labels,left_on=['sample.id','obj.id'],right_on=['Sample_ID','Cell_ID'],
        how='left',validate='one_to_one',suffixes=('','_replay'))
    assert len(checked)==len(frame) and checked['Cluster'].eq(checked['Cluster.v2']).all()
    # Metadata categories only; response, survival and follow-up age are unused.
    clinical=qualifier.decoded(qualifier.table_bytes(archive,'ST1_sample_info.txt'))
    _,separator=qualifier.header(clinical)
    fields=['Sample_ID','Cohort','Drug','CTLA4','PD1','Treatment','Type','Tissue_Source_Simplified']
    metadata=pd.read_csv(io.BytesIO(clinical),sep=separator,usecols=fields,dtype=str)
    assert not metadata['Sample_ID'].duplicated().any()
    selected=metadata.loc[metadata['Sample_ID'].isin(frame['sample.id'].unique())]
    assert set(selected['Sample_ID'])==set(frame['sample.id'])
    source_context=[dict(field=field,categories=[dict(value=str(value),source_materials=int(count))
        for value,count in selected[field].fillna('SOURCE_MISSING').value_counts().sort_index().items()]) for field in fields[1:]]
    helpers=load(folder/'analyze_melanoma_IMC_continuous_state_neighbors_v2.py',
        '51b8a4323d806ae011ecb7d9419d56548eba63ca947a6756d3ab1dc5f9baf7bf','native_rank_grid')
    neighbors=load(folder/'melanoma_IMC_numpy_neighbors_v1.py',
        'acbde69dae387ee79d6d8397899681bedf00553dd3da08878abd16b392858471','exact_native_neighbors')
    # Reuse the already qualified exact blocked kNN implementation.
    def knn(coords,k):
        _,multiplicity=np.unique(coords,axis=0,return_counts=True)
        assert multiplicity.max()<=50
        tree=neighbors.ExactEuclideanNeighbors(coords)
        distances,indices=tree.query(coords,k=min(k+1,len(coords)))
        pairs=set()
        radii=np.nextafter(distances[:,-1],np.inf)
        for i,nearby in enumerate(tree.query_ball_point(coords,radii)):
            for j in nearby:
                if i!=j:pairs.add((min(i,int(j)),max(i,int(j))))
        return np.asarray(sorted(pairs),dtype=int).reshape(-1,2)
    result,calls=source_statistics(frame,knn,helpers.tile_groups,helpers.rank_z)
    return dict(status='COMPLETE_FIXED_MATERIAL_KI67_CONDITIONAL_ORGANIZATION_NOT_PATIENT_OR_CLINICAL_VALIDATION',
        family=result,family_size=3,permutations=PERMUTATIONS,seed=SEED,
        direct_contrast_uses_only_shared_eligible_materials_and_same_CD8_draws=True,
        source_material_identity_not_promoted_to_patient=True,
        unreviewed_repeat_materials_may_affect_pooled_randomization=True,
        source_context=source_context,
        source_classes_identical_in_complete_ST4_and_raw_intensity_tables=True,
        metadata_response_survival_and_followup_age_not_selected=True,
        patient_CI_and_clinical_models_not_estimable_or_run=True,
        source_phase='Author-reported pretreatment melanoma; regimen/tissue stratification still unreviewed',
        source_sample_qualifier_SHA256='a498faeb15053036e72003585ba82d5c9b44b2f8f9d3c16d392c3fcf8e693125',
        independent_scalar_covariance_checks=calls,
        source_permutation_means_independently_replayed=True,
        graph='Within source material exact native k10 union with boundary ties; no micron/contact claim',
        no_Histone_background_or_segmentation_artifact_exclusion_claim=True,
        no_direct_communication_or_new_identity_claim=True,
        producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
