"""Full fixed source randomization, shared right-end draws and missing slots."""
import importlib.util
from pathlib import Path
import unittest
import json
import io
import tarfile
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
m=load(ROOT/'scripts/analyze_moldoveanu_fixed_Ki67_source_organization_v1.py','fixed_source')
h=load(ROOT/'scripts/analyze_melanoma_IMC_continuous_state_neighbors_v2.py','rank_grid')


class SourceOrganization(unittest.TestCase):
    def archive_bundle(self,corrupt=False):
        rows=[];metadata=[]
        for material in range(5):
            sample='synthetic_'+str(material)
            metadata.append(dict(Sample_ID=sample,Cohort='synthetic_cohort',Drug='synthetic_drug',
                CTLA4='0',PD1='1',Treatment='synthetic_regimen',Type='synthetic_type',Tissue_Source_Simplified='synthetic_tissue'))
            for typ,offset in [('Treg',0),('Th.ae',20),('Tc.ae',40)]:
                for i in range(20):rows.append({'sample.id':sample,'obj.id':str(offset+i),
                    'Cluster':typ,'coord.x':float(offset+i),'coord.y':0.,'Ki67':float(i)})
        frame=pd.DataFrame(rows)
        labels=frame.rename(columns={'sample.id':'Sample_ID','obj.id':'Cell_ID','Cluster':'Cluster.v2'})[['Sample_ID','Cell_ID','Cluster.v2']]
        if corrupt:labels.loc[0,'Cluster.v2']='wrong'
        out=io.BytesIO()
        with tarfile.open(fileobj=out,mode='w') as tar:
            for name,body in [('Data/ICI_meanIntensity.tsv.gz',labels.to_csv(sep='\t',index=False).encode()),
                ('Data/ST1_sample_info.txt',pd.DataFrame(metadata).to_csv(sep='\t',index=False).encode())]:
                item=tarfile.TarInfo(name);item.size=len(body);tar.addfile(item,io.BytesIO(body))
        archive=tarfile.open(fileobj=io.BytesIO(out.getvalue()),mode='r:');archive._own_ICI_qualified_cells=frame
        return dict(archive=archive)
    def test_complete_wrapper_real_neighbors_metadata_and_label_gates(self):
        bundle=self.archive_bundle()
        try:
            result=m.analyze(bundle,dict(scope='fixed_author_materials_no_patient_inference',permutations=2000,seed=m.SEED))
        finally:bundle['archive'].close()
        self.assertTrue(result['source_classes_identical_in_complete_ST4_and_raw_intensity_tables'])
        self.assertTrue(result['metadata_response_survival_and_followup_age_not_selected'])
        self.assertEqual(result['family_size'],3)
        self.assertNotIn('synthetic_0',json.dumps(result))
    def test_changed_author_classification_is_rejected_before_effects(self):
        bundle=self.archive_bundle(corrupt=True)
        try:
            with self.assertRaises(AssertionError):m.analyze(bundle,dict(scope='fixed_author_materials_no_patient_inference',permutations=2000,seed=m.SEED))
        finally:bundle['archive'].close()
    def test_full_2000_three_slots_positive_negative_direct_contrast(self):
        rows=[]
        for material in range(5):
            for typ,offset in [('Treg',0),('Th.ae',20),('Tc.ae',40)]:
                for i in range(20):
                    rows.append({'sample.id':'synthetic_'+str(material),'Cluster':typ,
                        'coord.x':offset+i,'coord.y':0.,'Ki67':float(19-i if typ=='Th.ae' else i)})
        frame=pd.DataFrame(rows)
        def knn(coords,k):
            self.assertEqual((len(coords),k),(60,10))
            return np.array([(i,40+i) for i in range(20)]+[(20+i,40+i) for i in range(20)])
        ordinary=m.shuffled;calls=[]
        def counted(values,groups,rng):calls.append(1);return ordinary(values,groups,rng)
        m.shuffled=counted
        try:result,checked=m.source_statistics(frame,knn,h.tile_groups,h.rank_z)
        finally:m.shuffled=ordinary
        primary=result['ROI_parent_shuffle'];self.assertEqual(len(primary),3)
        self.assertAlmostEqual(primary[0]['effect'],1.)
        self.assertAlmostEqual(primary[1]['effect'],-1.)
        self.assertAlmostEqual(primary[2]['effect'],2.)
        self.assertTrue(all(r['P']<=.01 and r['source_materials']==5 for r in primary))
        self.assertTrue(all(r['patient_bootstrap_95_CI'] is None and not r['patient_identity_reviewed'] for r in primary))
        self.assertEqual(len(calls),5*2*2000*3)
        self.assertEqual(checked,5*2*2*2002)
        self.assertNotIn('synthetic_',json.dumps(result))
    def test_missing_comparisons_keep_family_and_unknown_patient_identity(self):
        result=m.finish([]);self.assertEqual(len(result),3)
        self.assertTrue(all(r['P'] is None and r['BH'] is None for r in result))
        np.testing.assert_allclose(m.bh([.01,.04,1.]),[.03,.06,1.])
    def test_independent_interface_affinity_is_not_pair_covariance(self):
        values=np.arange(20,dtype=float)
        a,b=np.array([(i,j) for i in range(15,20) for j in range(15,20)]).T
        self.assertAlmostEqual(m.covariance(values,values,a,b),0.)


if __name__=='__main__':unittest.main()
