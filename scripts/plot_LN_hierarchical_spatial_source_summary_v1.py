"""Private source-level scientific plots; never patient/cell/graph arrays."""
from pathlib import Path
import json,math,os

def main():
    task=Path(os.environ['PANCANCER_LN_STRUCTURE_TASK'])
    os.environ.setdefault('MPLCONFIGDIR','/tmp/pancancer_LN_spatial_mpl_20261022')
    os.environ.setdefault('XDG_CACHE_HOME','/tmp/pancancer_LN_spatial_fontcache_20261022')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    snapshot=json.loads((task/'private_LN_full_spatial_and_audit_actual_snapshot_latest.json').read_text())
    files=snapshot['files'];prod=files['private_spatial_production_complete.json'];audit=files['private_spatial_audit_complete.json']
    assert audit['status']=='PASS_FULL_FROZEN_FAMILY_NUMERIC_AND_GRAPH_AUDIT'
    assert prod['family_contract_SHA256']==audit['family_contract_SHA256']
    measurement=json.loads((task/'private_LN_measurement_units_phase002_acceptance_v1.json').read_text())['scientific_qualification']
    geometry=measurement['geometry_under_candidate_corrected_patient_node_image_mapping']
    source_images=geometry['eligible_distinct_patient_node_ROI_tuples'];source_nodes=geometry['eligible_distinct_patient_node_tuples'];source_patients=geometry['eligible_distinct_literal_patients'];paired_patients=measurement['corrected_patient_labels_with_both_stations']
    rows=prod['primary_rows'];assert len(rows)==27
    pairs=list(dict.fromkeys(x['pair'] for x in rows));estimands=['whole','conditional','N1_minus_N2'];assert len(pairs)==9
    records={(x['pair'],x['estimand']):x for x in rows};assert len(records)==27
    out=task/'source_figures_v1';out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    labels=[f'{p.split("_")[0]}  /  {p.split("__")[1].split("_")[0]}' for p in pairs]
    colors=['#2166ac']*3+['#7b3294']*3+['#666666']*3
    fig,axes=plt.subplots(1,3,figsize=(15,7.9),sharey=True)
    titles=['Whole-neighborhood O/E','Conditional state-pair odds ratio','Paired N1/N2 odds-ratio ratio']
    for ax,est,title in zip(axes,estimands,titles):
        for j,pair in enumerate(pairs):
            r=records[pair,est];point=r.get('median_effect');lo=r.get('bootstrap_lower');hi=r.get('bootstrap_upper')
            if point is None or lo is None or hi is None:
                ax.text(.4,j,'not estimable',va='center',color='#777777');continue
            x,l,h=map(math.exp,(point,lo,hi));assert l<=x<=h
            by=r.get('BY');fill=colors[j] if by is not None and by<.05 else 'white'
            ax.errorbar(x,j,xerr=[[x-l],[h-x]],fmt='o',markersize=6,color=colors[j],mfc=fill,mec=colors[j],capsize=3,lw=1.4)
            q=r.get('BH');qt='NA' if q is None else f'{q:.3g}'
            ax.text(.98,j,f'n={r["cases"]}; q={qt}',transform=ax.get_yaxis_transform(),va='center',ha='right',fontsize=7.5)
        ax.axvline(1,color='#888888',ls='--',lw=.9);ax.axhline(5.5,color='#bbbbbb',lw=.8)
        ax.set_xscale('log');ax.set_xlim(.30,12);ax.set_xticks([.5,1,2,4,8]);ax.set_xticklabels(['0.5','1','2','4','8'])
        ax.set_title(title,fontsize=10);ax.set_xlabel('Ratio on log scale');ax.grid(axis='x',alpha=.14)
    axes[0].set_yticks(range(9));axes[0].set_yticklabels(labels);axes[0].invert_yaxis()
    fig.suptitle('NSCLC regional LN: prespecified relative protein-state spatial family',fontsize=13,y=.98)
    key=('R0: CD4-like FOXP3/CTLA4 relative-high    R1: R0 with MKI67 relative-high\n'
         'R2: general CD4-like MKI67 relative-high reference\n'
         'E0: CD8a-like PD1/TIM3 relative-high    E1: CD8a-like MKI67 relative-high\n'
         'E2: CD8a-like CCR7/CD127 relative-high\n'
         'Points: median registered-patient effect. Intervals: 95% percentile patient bootstrap (2000 draws), conditional on fixed annotation.\n'
         'q: BH across 27 planned tests; filled points: BY < 0.05. Nested/overlapping states are not independent replications.\n'
         'Nearest-neighbor direction does not establish signalling direction.\n'
         'Candidate marker states, not absolute positivity, functional suppression or clinical prediction.')
    fig.text(.07,.035,key,ha='left',va='bottom',fontsize=8.2)
    fig.tight_layout(rect=(.03,.25,.995,.94))
    for ext in ('png','svg'):fig.savefig(out/f'FROZEN_FAMILY_FOREST.{ext}',dpi=200)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,5.9),sharey=True)
    for ax,field,title in zip(axes,['eligible_images','cases'],['Eligible images','Registered-patient units']):
        z=np.array([[records[p,e][field] for e in estimands] for p in pairs])
        im=ax.imshow(z,aspect='auto',cmap='Blues',vmin=0,vmax=source_images if field=='eligible_images' else source_patients)
        for j in range(9):
            for k in range(3):ax.text(k,j,str(z[j,k]),ha='center',va='center',color='white' if z[j,k]>(source_images/2 if field=='eligible_images' else source_patients/2) else 'black')
        ax.set_xticks(range(3));ax.set_xticklabels(['Whole','Conditional','N1/N2 paired']);ax.set_title(title)
        fig.colorbar(im,ax=ax,fraction=.035,pad=.04)
    axes[0].set_yticks(range(9));axes[0].set_yticklabels(labels)
    fig.suptitle('Source-level coverage for every planned state pair',fontsize=12)
    fig.text(.08,.025,f'Captured source: {source_images} images / {source_nodes} nodes / {source_patients} registered-patient units; {paired_patients} units have both stations.\nImages are aggregated within nodes and patients. Missing measurement or coverage is never filled as negative.',fontsize=8)
    fig.tight_layout(rect=(.02,.09,1,.94))
    for ext in ('png','svg'):fig.savefig(out/f'SOURCE_COVERAGE.{ext}',dpi=200)
    plt.close(fig)
    choices=[('primary',5,'k5 vs k10'),('primary',20,'k20 vs k10'),('area_adjusted',10,'Area-adjusted states vs primary')]
    fig,axes=plt.subplots(1,3,figsize=(15,6),sharey=True)
    for ax,(mode,k,title) in zip(axes,choices):
        for j,pair in enumerate(pairs):
            r=next(x for x in prod['sensitivities'] if x['pair']==pair and x['mode']==mode and x['k']==k and x['estimand']=='conditional')
            a=r['primary_median_on_common_cases'];b=r['sensitivity_median_on_common_cases']
            if a is None or b is None:ax.text(.5,j,'not estimable',va='center');continue
            a,b=math.exp(a),math.exp(b);ax.plot([a,b],[j,j],color=colors[j],lw=1.5)
            ax.plot(a,j,'o',mfc='white',mec='#444444',ms=5);ax.plot(b,j,'s',color=colors[j],ms=5)
            ax.text(.98,j,f'{r["common_eligible_images"]} images; n={r["common_cases"]}',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7.2)
        ax.set_xscale('log');ax.set_xlim(.7,6);ax.set_xticks([1,2,3,4,6]);ax.set_xticklabels(['1','2','3','4','6']);ax.axvline(1,color='#888888',ls='--',lw=.8);ax.axhline(5.5,color='#bbbbbb',lw=.8);ax.set_title(title,fontsize=10);ax.set_xlabel('Conditional odds ratio');ax.grid(axis='x',alpha=.15)
    axes[0].set_yticks(range(9));axes[0].set_yticklabels(labels);axes[0].invert_yaxis()
    fig.suptitle('Prespecified sensitivities: matched eligible images and patient units',fontsize=12)
    fig.text(.07,.025,'Open circle: primary estimate on the common image/case set. Square: sensitivity estimate on that same set.\nPoint/coverage checks only; no extra P tests. Area adjustment can also remove biological size-associated variation.',fontsize=8)
    fig.tight_layout(rect=(.02,.09,1,.93))
    for ext in ('png','svg'):fig.savefig(out/f'MATCHED_SENSITIVITY.{ext}',dpi=200)
    plt.close(fig)
    selected=[]
    for pair in [pairs[0],pairs[3],pairs[2]]:
        r=records[pair,'conditional'];selected.append({'pair':pair,'conditional_odds_ratio':math.exp(r['median_effect']),'lower':math.exp(r['bootstrap_lower']),'upper':math.exp(r['bootstrap_upper']),'cases':r['cases'],'positive_cases':r['positive_cases'],'BH':r['BH'],'BY':r['BY']})
    result={'status':'SOURCE_FIGURES_CREATED_PENDING_VISUAL_REVIEW','family_contract_SHA256':prod['family_contract_SHA256'],'primary_tests':27,'source_specific_candidates':selected,'interpretation_limits':['relative marker candidates, not formal functional identities','general proliferation reference is also associated; regulatory specificity has not been directly contrasted','no consistent N1/N2 difference in this frozen family','not exact CM02 projection, direct LR binding or clinical-response validation','bootstrap intervals condition on frozen annotation; source selection and spatial background remain limitations'],'figures':[p.name for p in out.iterdir() if p.suffix in ('.png','.svg')]}
    (task/'private_LN_source_figures_and_candidate_summary_v1.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
