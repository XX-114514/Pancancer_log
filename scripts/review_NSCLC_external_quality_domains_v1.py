"""Outcome-free quality-domain diagnostics; not a scientific source filter."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_clinical_inputs_v1 import sha,VIEWS

FIELDS=['logN','mean_logRNA','mean_loggenes','symbol_pct_mito','symbol_pct_ribo']
FULL=['mean_log_full_source_RNA','mean_log_full_source_detected_genes','symbol_pct_mito_RNA','symbol_pct_ribo_RNA']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    src=root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    proof=src/'private_independent_native_count_input_audit.json'
    assert sha(proof)=='7bd028cc1c89ac7324ac25c65042b5b184cb377162a38f9544caf7bb24048c04'
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(src/'private_terminal.json')
    t=json.loads((src/'private_terminal.json').read_text());name='private_all6760_native_count_units_and_quality.tsv.gz'
    assert sha(src/name)==t['files'][name]['sha256']
    units=pd.read_csv(src/name,sep='\t');assert len(units)==6760
    out=root/'NSCLC_native_MFA_external_quality_domain_review_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),count_audit_sha256=sha(proof),
        scope='All61 source-arm x4 frozen pools x5 observed RNA quality variables. Eligible libraries have>=20 cells and positive source/common RNA; zero/small units remain counted separately.',
        reference='Original145 discovery view-specific quality mean/SD/min/max, never refit on outcome. Full SOURCE RNA quality compared descriptively to original whole native gene axis; common-native alternative RNA/detection also recorded.',
        thresholds='Original discovery range and absZ>5 are descriptive extrapolation diagnostics ONLY. No automatic source/feature/sample exclusion or new scientific P.',
        parent_fraction='External diagnostic parent definitions differ from original lung full T/NK etc. No exact parent bridge or exact latent-factor transfer accepted.',
        next='Use quality diagnostics to fix a clearly labelled frozen gene-programme assay extension before any external expression correlation. Do not silently apply original quality coefficients across incompatible definitions.',
        no_clinical_or_factor_effects=True)
    (out/'plan_before_quality_domain_diagnostics.json').write_text(json.dumps(plan,indent=2)+'\n')
    rows=[]
    for vi,view in enumerate(VIEWS):
        path=root/'NSCLC_native_MFA_inputs_20261016_v1'/('private_'+view+'_prepared.npz')
        with np.load(path,allow_pickle=False) as z:q=z['quality'][:,:5];training=z['train'];mean=z['quality_mean'][:5];sd=z['quality_sd'][:5]
        assert (sd>0).all();low=q[training].min(axis=0);high=q[training].max(axis=0)
        for (arm,source),g in units[units['view'].eq(view)].groupby(['arm','dataset'],sort=True):
            eligible=g.n_cells.ge(20)&g.full_source_RNA_sum.gt(0)&g.common_native_RNA_sum.gt(0)
            gg=g[eligible]
            external=np.column_stack([np.log1p(gg.n_cells),gg[FULL].to_numpy(float)])
            alternate=external.copy()
            if len(gg):
                alternate[:,1]=gg.mean_log_common_native_RNA;alternate[:,2]=gg.mean_log_common_native_detected_genes
            for j,field in enumerate(FIELDS):
                x=external[:,j];z=(x-mean[j])/sd[j]
                a=(alternate[:,j]-mean[j])/sd[j]
                rows.append(dict(arm=arm,dataset=source,view=view,quality_variable=field,all_library_pools=len(g),eligible_libraries=len(gg),
                    libraries_below20=int(g.n_cells.lt(20).sum()),zero_source_or_common_RNA_libraries=int((g.n_cells.ge(20)&~eligible).sum()),
                    native_discovery_mean=float(mean[j]),native_discovery_SD=float(sd[j]),
                    external_median=float(np.median(x)) if len(x) else np.nan,
                    outside_native_discovery_range=int(((x<low[j])|(x>high[j])).sum()),absZ_gt5=int((abs(z)>5).sum()),
                    median_absZ=float(np.median(abs(z))) if len(z) else np.nan,
                    common_RNA_alternative_median_absZ=float(np.median(abs(a))) if len(a) else np.nan,
                    source_or_identity_accepted=False))
    assert len(rows)==61*4*5
    pd.DataFrame(rows).to_csv(out/'private_all1220_source_view_quality_domain_descriptions.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL61_SOURCE4_VIEW5_QUALITY_DOMAIN_DESCRIPTIONS_NO_AUTOMATIC_FILTER',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        exact_native_factor_or_clinical_transfer_accepted=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
