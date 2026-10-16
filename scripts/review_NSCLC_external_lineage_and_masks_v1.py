"""Complete source/stratum opportunity join; no external expression or inference."""
import itertools
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_clinical_inputs_v1 import matched_quota, sha
from NSCLC_external_lineage_role_rules_v1 import VIEWS, CHANNELS


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    capture = root/'NSCLC_native_MFA_external_lineage_qualification_20261016_v2'
    masks = root/'NSCLC_native_MFA_external_gene_mask_calibration_20261016_v1'
    proofs = {
        capture/'private_independent_lineage_case_audit.json': 'b07bf0d9b682e8b6cc2ef72192a0ff110a722ba99bd7a931c5e7745c7640699a',
        masks/'private_independent_native_mask_audit.json': 'f82d4106e187eaf433adf0ffafd9f6e3a03051d47d23a460fcf1292838f5f708',
    }
    for p, digest in proofs.items():
        assert sha(p) == digest
        a = json.loads(p.read_text())
        assert a['producer_terminal_sha256'] == sha(p.parent/'private_terminal.json')
    out = root/'NSCLC_native_MFA_external_opportunity_review_20261016_v1'
    assert not out.exists()
    quota = matched_quota(root)
    out.mkdir(mode=0o700)
    paths = [capture/'private_all2013_source_role_pair_capture.tsv',
             capture/'private_all_source_stratum_role_pair_capture.tsv.gz',
             masks/'private_all2440_source_mask_calibration.tsv']
    source, strata, cal = [pd.read_csv(p, sep='\t', keep_default_na=False) for p in paths]
    assert len(source) == 2013 and len(strata) == 10428 and len(cal) == 2440
    gates = cal.set_index(['arm', 'dataset', 'view', 'factor'])
    assert gates.index.is_unique
    scope = {v: [v] for v in VIEWS}
    scope.update({'pair::'+a+'::'+b: [a,b] for a,b in itertools.combinations(VIEWS,2)})
    scope['all4'] = VIEWS
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(),
                source_sha256=sha(__file__), quota=quota,
                input_sha256={p.name: sha(p) for p in paths},
                audit_sha256={p.name+str(i): sha(p) for i,p in enumerate(proofs)},
                scope='All61 source-arm entries,3 role channels,11 capture scopes and10 fixed factors; every source and recorded stratum retained.',
                gates='Join pre-existing physical and native-mask gates for ALL endpoints. Capture/20 reviewed blocks flags describe opportunities, not accepted external models.',
                labels='Frozen candidate pools differ from full native author views; author lexical roles remain candidates. No channel selected by opportunity count or effect.',
                strata='Preserve all original10 metadata fields. Full-stratum N is one diagnostic; a future model may prespecify within-source pooling over a covariate but may not count repeated libraries as patients.',
                limits='No expression scores, P, clinical acceptance, external independence or annotation changes.')
    (out/'plan_before_opportunity_join.json').write_text(json.dumps(plan,indent=2)+'\n')
    for name, table in [('source',source),('stratum',strata)]:
        records=[]
        for r in table.to_dict('records'):
            for k in range(1,11):
                factor='NSCLC_MFA_20261016_F%02d'%k
                ends=scope[r['capture_scope_id']]
                g=[gates.loc[(r['arm'],r['dataset'],v,factor)] for v in ends]
                physical=all(bool(v.physical_gene_gate) for v in g)
                calibrated=all(bool(v.native_mask_calibration_pass) for v in g)
                records.append(dict(r,factor=factor,all_endpoints_physical_gate=physical,
                    all_endpoints_native_mask_gate=calibrated,
                    calibrated_capture_any=bool(calibrated and r['captured_libraries']>0),
                    calibrated_reviewed20_opportunity=bool(calibrated and r['reviewed_patient_blocks']>=20),
                    new_external_expression_or_model_accepted=False))
        result=pd.DataFrame(records)
        assert len(result)==len(table)*10
        result.to_csv(out/('private_all_'+name+'_calibrated_opportunities.tsv.gz'),sep='\t',index=False,compression='gzip')
    keys=sorted(set(zip(source.arm,source.dataset)))
    labels=[('I: ' if a=='archived_immune_atlas' else 'O: ')+d for a,d in keys]
    # All three channels shown together; no author-label recovery is presumed.
    matrix=np.zeros((61,12),int)
    for i,(arm,dataset) in enumerate(keys):
        for ci,ch in enumerate(CHANNELS):
            for vi,view in enumerate(VIEWS):
                row=source[(source.arm==arm)&(source.dataset==dataset)&(source.channel==ch)&(source.capture_scope_id==view)]
                assert len(row)==1
                matrix[i,ci*4+vi]=int(row.captured_libraries.iloc[0])
    fig,ax=plt.subplots(figsize=(13,19))
    im=ax.imshow(np.log1p(matrix),aspect='auto',cmap='viridis',vmin=0)
    ax.set_yticks(range(61),labels,fontsize=7)
    ax.set_xticks(range(12),[c+'\n'+v.removeprefix('author_') for c in ['Frozen candidate','Author explicit','Author concordant'] for v in VIEWS],rotation=50,ha='right',fontsize=8)
    for i in range(61):
        for j in range(12):
            ax.text(j,i,str(matrix[i,j]),ha='center',va='center',fontsize=6,color='white' if matrix[i,j]<8 else 'black')
    ax.set_title('Captured libraries (each endpoint >=20 cells); all61 sources / three candidate channels\nCounts are not independent patients; O=original frozen scope, I=archived immune atlas',fontsize=10)
    fig.colorbar(im,ax=ax,fraction=.025,label='log(1 + captured libraries)')
    fig.tight_layout();fig.savefig(out/'private_all61_three_channel_capture.png',dpi=160);plt.close(fig)
    scores=np.zeros((61,40),int)
    for i,(arm,dataset) in enumerate(keys):
        for vi,view in enumerate(VIEWS):
            for k in range(10):
                r=gates.loc[(arm,dataset,view,'NSCLC_MFA_20261016_F%02d'%(k+1))]
                scores[i,vi*10+k]=2 if r.native_mask_calibration_pass else 1 if r.physical_gene_gate else 0
    fig,ax=plt.subplots(figsize=(14,19))
    from matplotlib.colors import ListedColormap,BoundaryNorm
    im=ax.imshow(scores,aspect='auto',cmap=ListedColormap(['#eeeeee','#e5a340','#257e63']),norm=BoundaryNorm([-.5,.5,1.5,2.5],3))
    ax.set_yticks(range(61),labels,fontsize=7)
    ax.set_xticks(range(40),['F%02d'%(i%10+1) for i in range(40)],rotation=90,fontsize=7)
    for j,v in enumerate(VIEWS):ax.text(j*10+4.5,-1.5,v.removeprefix('author_'),ha='center',fontsize=9)
    for x in [9.5,19.5,29.5]:ax.axvline(x,color='white',lw=2)
    ax.set_title('All61 source masks x40 fixed factor/view slots\nNative-mask calibration is a measurement gate, not external biological validation',pad=35,fontsize=10)
    bar=fig.colorbar(im,ax=ax,fraction=.025,ticks=[0,1,2]);bar.ax.set_yticklabels(['Physical fail','Physical only','Both pass'])
    fig.tight_layout();fig.savefig(out/'private_all2440_measurement_gates.png',dpi=160);plt.close(fig)
    terminal=dict(status='COMPLETE_ALL_SOURCE_AND_STRATUM_OPPORTUNITIES_AND_TWO_COVERAGE_FIGURES_NO_EXTERNAL_EFFECT',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        external_model_or_independent_validation_accepted=False)
    assert sum(v['bytes'] for v in terminal['files'].values())<268435456
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
