"""Private 24-slot/family/missingness audit; no scientific payload exports."""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out=root/'CRC_fine_pair_DIALOGUE_20261014_v1/private_MMR_cohort_test_v1'
    results=pd.read_csv(out/'all24_private_results.tsv',sep='\t')
    assert len(results)==24
    assert not results.duplicated(['partition','program','endpoint','method']).any()
    for _,group in results.groupby(['partition','method']):
        assert len(group)==4
        valid=group.p.notna()
        q=multipletests(group.p.fillna(1).to_numpy(),method='fdr_bh')[1]
        assert np.allclose(group.loc[valid,'BH_four_slot_family'],q[valid],atol=1e-14)
        assert group.loc[~valid,'BH_four_slot_family'].isna().all()
        assert (group.loc[valid,'n_MMRd']>=5).all() and (group.loc[valid,'n_MMRp']>=5).all()
        assert np.isfinite(group.loc[valid,'effect']).all()
        assert group.loc[valid,'p'].between(0,1).all()
    unsupported=results.status.eq('no_gene_supported_score')
    assert results.loc[unsupported,['effect','lower','upper','p','BH_four_slot_family']].isna().all().all()
    estimated=results.p.notna()
    assert (results.loc[estimated,'n_complete']==results.loc[estimated,'n_MMRd']+results.loc[estimated,'n_MMRp']).all()
    destination=out/'private_structure_audit.json'
    assert not destination.exists()
    destination.write_text(json.dumps(dict(status='PASS_ALL24_SLOT_KEYS_SIX_BH_FAMILIES_MISSINGNESS_AND_PATIENT_DENOMINATORS',
                                         results_private=True),indent=2)+'\n')
    print('PASS_PRIVATE_MMR_RESULT_STRUCTURE_NO_RESULTS_EXPORTED')


if __name__=='__main__':
    main()
