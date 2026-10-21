"""Synthetic identity-unit and missingness gates; no real records/network."""
import importlib.util
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

try:
    import pandas as pd
except ImportError:
    pd=None
try:
    import numpy as np
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree=None


@unittest.skipIf(pd is None,'Optional pandas/numpy analytical runtime unavailable')
class PatientCoverageSemantics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=Path(__file__).resolve().parents[1]/'scripts'/'summarize_melanoma_IMC_patient_coverage_v1.py'
        spec=importlib.util.spec_from_file_location('coverage_under_test',path)
        cls.module=importlib.util.module_from_spec(spec);spec.loader.exec_module(cls.module)
        path=path.parent/'qualify_melanoma_IMC_clinical_covariate_keys_v4.py'
        spec=importlib.util.spec_from_file_location('clinical_format_under_test',path)
        cls.clinical=importlib.util.module_from_spec(spec);spec.loader.exec_module(cls.clinical)
        path=path.parent/'wait_for_melanoma_IMC_patient_coverage_v1.py'
        spec=importlib.util.spec_from_file_location('dependency_under_test',path)
        cls.dependency=importlib.util.module_from_spec(spec);spec.loader.exec_module(cls.dependency)

    def setUp(self):
        self.registry=[dict(sampleROI2=k,region=v) for k,v in [
            ('roiA','invasive'),('roiB','invasive'),('roiC','invasive'),('roiD','invasive'),('roiCore','tumor')]]
        self.bridge=[dict(ROI=k,patient=v) for k,v in [
            ('roiA','demoA'),('roiB','demoA'),('roiC','demoB'),('roiD','demoC')]]
        rows=[]
        for roi,typ,n,value in [('roiA','Treg',10,1),('roiB','Treg',10,3),('roiA','CD8T',20,2),
                ('roiC','Treg',19,1),('roiC','CD8T',20,4),('roiD','Treg',20,2),('roiCore','CD8T',3,9)]:
            rows.extend(dict(sample_id=roi,celltype20=typ,**{m:float(value) for m in self.module.MARKERS}) for _ in range(n))
        self.frame=pd.DataFrame(rows)

    def test_patient_not_ROI_denominator_and_absent_mean(self):
        m=self.module;lookup,keys=m.explicit_margin_map(self.registry,self.bridge)
        result,cells,means=m.summarize(self.frame,lookup,keys)
        self.assertEqual(result['source_patient_units'],3)
        self.assertEqual(result['source_ROIs_with_explicit_patient'],4)
        self.assertEqual(result['all_author_type_pair_coverage'][0]['both_types_ge20_patients'],1)
        self.assertEqual(cells.loc['demoA','Treg'],20)
        self.assertEqual(cells.loc['demoB','Treg'],19)
        self.assertEqual(cells.loc['demoC','CD8T'],0)
        self.assertTrue(pd.isna(means['Ki-67'].unstack().loc['demoC','CD8T']))
        self.assertEqual(means.loc[('demoA','Treg'),'Ki-67'],2)
        self.assertEqual(result['known_nonmargin_cells'],3)
        self.assertTrue(all(x['mean_IQR_quartiles'] is None for x in result['all_marker_mean_coverage']))
        replay=m.scalar_replay(self.frame.to_csv(index=False).encode(),lookup,cells,means)
        self.assertTrue(replay['all_CSV_dictionary_patient_type_counts_and_marker_means_replayed'])

    def test_no_duplicate_or_incomplete_patient_bridge(self):
        for bad in [self.bridge+[dict(ROI='roiA',patient='other')],self.bridge[:-1]]:
            with self.assertRaises(AssertionError):self.module.explicit_margin_map(self.registry,bad)

    def test_nonfinite_and_unknown_ROI_not_zero_filled(self):
        lookup,keys=self.module.explicit_margin_map(self.registry,self.bridge)
        for col,value in [('Ki-67',float('nan')),('sample_id','unmapped_unknown')]:
            bad=self.frame.copy();bad.loc[0,col]=value
            with self.assertRaises(AssertionError):self.module.summarize(bad,lookup,keys)

    def test_author_patient_format_is_unique_and_not_guessed(self):
        result=self.clinical.author_patient_keys({'demo','otherF'},'X')
        self.assertEqual(result,{'demo':'XdemoF','otherF':'XotherF'})
        self.assertEqual(self.clinical.author_patient_keys({'demo'},'X.'),{'demo':'X.demoF'})
        with self.assertRaises(AssertionError):self.clinical.author_patient_keys({'demo','demoF'},'X')
        self.assertEqual(self.clinical.author_patient_keys({'demo'},'prefix123'),{'demo':'prefix123demoF'})
        with self.assertRaises(AssertionError):self.clinical.author_patient_keys({'demo'},'X\n')
        self.assertEqual(self.clinical.field_role('Age (years)'),'age')
        self.assertIsNone(self.clinical.field_role('PFS'))

    def test_progress_or_hash_mismatch_cannot_start_another_reader(self):
        m=self.dependency
        with tempfile.TemporaryDirectory() as directory:
            accepted=Path(directory)/'accepted.json'
            self.assertFalse(m.ready(accepted))
            data=dict(status=m.COMPLETE,producer_SHA256=m.PRODUCER,complete_member_CRC_checked=True,
                complete_independent_CSV_numeric_ROI_and_label_counter_replay=True,
                published_metadata_and_directory_replayed=True,all_literal_ROI_keys_exact=True,
                author_label_counts=dict(celltype20=dict(Treg=10)))
            accepted.write_text(json.dumps(data))
            self.assertFalse(m.ready(accepted))
            terminal=dict(status=m.COMPLETE,producer_SHA256=m.PRODUCER,
                files={accepted.name:dict(SHA256=hashlib.sha256(accepted.read_bytes()).hexdigest())})
            (accepted.parent/'private_terminal.json').write_text(json.dumps(terminal))
            self.assertTrue(m.ready(accepted))
            data['all_literal_ROI_keys_exact']=False;accepted.write_text(json.dumps(data))
            with self.assertRaises(AssertionError):m.ready(accepted)


@unittest.skipIf(pd is None or cKDTree is None,'Optional complete analytical runtime unavailable')
class AliasAndNativeGeometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=Path(__file__).resolve().parents[1]/'scripts'/'summarize_melanoma_IMC_patient_coverage_v2.py'
        spec=importlib.util.spec_from_file_location('native_coverage_under_test',path)
        cls.module=importlib.util.module_from_spec(spec);spec.loader.exec_module(cls.module)

    def test_author_ROI_aliases_are_bijective_and_complete(self):
        text="""if (stringr::str_detect(each, 'F')){
        splitID = str_split(each,'_')[[1]]
        a = paste0(splitID[1], 'F_', splitID[2])
        """
        m=self.module;infix=m.author_ROI_infix(text)
        result=m.canonical_ROIs({'demo_1','otherF_2'},{'demoF_1','otherF_2'},infix)
        self.assertEqual(result,{'demo_1':'demoF_1','otherF_2':'otherF_2'})
        for raw,keys in [({'demo_1','demoF_1'},{'demoF_1'}),({'demo_1'},{'otherF_2'})]:
            with self.assertRaises(AssertionError):m.canonical_ROIs(raw,keys,infix)

    def test_native_neighbor_ties_and_scale_match_full_distance_reference(self):
        m=self.module;xy=m.coordinate_array(pd.Series(['(0,0)','(1,0)','(0,1)','(1,1)']))
        d=((xy[:,None]-xy[None,:])**2).sum(axis=2);np.fill_diagonal(d,np.inf)
        expected=set()
        for i in range(len(xy)):
            for j in np.flatnonzero(d[i]==d[i].min()):expected.add(tuple(sorted((i,int(j)))))
        observed={tuple(x) for x in m.native_knn_edges(xy,1)}
        self.assertEqual(observed,expected)
        self.assertEqual(observed,{tuple(x) for x in m.native_knn_edges(xy*7+13,1)})
        with self.assertRaises(AssertionError):m.native_knn_edges(np.zeros((60,2)),1)

    def test_different_ROI_frames_never_create_neighbors(self):
        m=self.module;xy=np.tile([[0,0],[1,0],[0,1],[1,1]],(2,1))
        frame=pd.DataFrame(dict(sample_id=['roiA']*4+['roiB']*4,celltype20=['Treg','CD8T','Treg','CD8T']*2))
        cells=pd.DataFrame([[4,4]],index=['demo'],columns=['CD8T','Treg'])
        result=m.graph_coverage(frame,xy,{'roiA':'demo','roiB':'demo'},cells,k=1)
        self.assertEqual(result['source_margin_undirected_edges'],8)

    def test_next_module_cannot_return_arrays_or_record_values(self):
        m=self.module;sensitive={'demo','roiA'}
        m.source_summary_only(dict(n_patients=20,CI=[-.1,.3],rows=[dict(author_type='Treg',effect=.2)]),sensitive)
        for value in [np.arange(24),pd.DataFrame(dict(x=[1])),dict(patient='demo'),dict(x='demo'),list(range(24))]:
            with self.assertRaises(AssertionError):m.source_summary_only(value,sensitive)


if __name__=='__main__':unittest.main()
