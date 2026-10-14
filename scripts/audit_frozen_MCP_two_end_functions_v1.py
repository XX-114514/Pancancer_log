"""Independent replay of all two-end functional descriptions and LOO."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from audit_HNSCC_subtype_function_followup_v1 import partial,same


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out=root/'frozen_MCP_two_end_function_review_20261014_v1'
    result=pd.read_csv(out/'all1500_private_function_diagnostics.tsv',sep='\t')
    md=pd.read_csv(out/'private_patient_programme_function_means.tsv',sep='\t')
    assert len(result)==1500 and not result.duplicated(['source','endpoint','reference','variant']).any()
    assert not md.duplicated(['source','endpoint','block']).any()
    for (source,endpoint),family in result.groupby(['source','endpoint']):
        m=md[md.source.eq(source)&md.endpoint.eq(endpoint)]
        if m.empty:
            assert family.status.eq('no_reviewed_scored_patient_pair').all() and family.rho.isna().all();continue
        x=m.programme.to_numpy()
        covs=dict(unadjusted=np.empty((len(m),0)),two_native_RNA_quality_means=m[['logRNA','loggenes']].to_numpy(),
                  RNA_quality_and_cell_count=np.column_stack((m[['logRNA','loggenes']].to_numpy(),np.log1p(m.native_cell_count))))
        for row in family.itertuples():
            assert row.n_patients==len(m)
            if row.reference not in m:
                assert row.status=='reference_measurement_insufficient' and pd.isna(row.rho);continue
            y=m[row.reference].to_numpy();cv=covs[row.variant];point=partial(x,y,cv);same(point,row.rho)
            values=[]
            for i in range(len(x)):
                keep=np.arange(len(x))!=i;v=partial(x[keep],y[keep],cv[keep])
                if np.isfinite(v):values.append(v)
            assert len(values)==row.LOO_valid
            same(min(values) if values else np.nan,row.LOO_min);same(max(values) if values else np.nan,row.LOO_max)
    report=dict(status='PASS_ALL1500_FUNCTION_DESCRIPTORS_AND_COMPLETE_LOO',finished_UTC=datetime.now(timezone.utc).isoformat(),
      result_sha256=hashlib.sha256((out/'all1500_private_function_diagnostics.tsv').read_bytes()).hexdigest(),no_new_P_BH=True)
    (out/'private_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
