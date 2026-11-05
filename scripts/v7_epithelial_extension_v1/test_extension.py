"""Semantic counterexamples and deterministic resource replay, not accuracy tests."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
RESOURCE=ROOT/'resources/annotation/v7_epithelial_extension_v1'
sys.path.insert(0,str(ROOT/'scripts/v7_panel_repair_v1'))
from evidence_engine import evaluate, validate_registry
from build import build

class ExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry=json.loads((RESOURCE/'registry_merged.json').read_text())
        cls.byid=validate_registry(cls.registry)

    def report(self, genes, contexts=None):
        expr=pd.DataFrame([dict(cluster='synthetic',gene=g,n_cells=100,n_measured=100,
                                pct_expressing_measured=v) for g,v in genes.items()])
        return evaluate(self.registry,expr,contexts=contexts).set_index('panel')

    def test_single_pax8_does_not_supply_thyroid_machinery(self):
        r=self.report({'PAX8':.8})
        self.assertFalse(r.loc['Epi_thyroid_follicular_like','each_core_component_has_measured_gene'])
        self.assertFalse(r.final_label_assigned.any())

    def test_collecting_marker_does_not_supply_podocyte_evidence(self):
        r=self.report({'AQP2':.8,'AQP3':.5})
        self.assertTrue(r.loc['Epi_kidney_principal_like','each_core_component_has_measured_gene'])
        self.assertEqual(r.loc['Epi_kidney_podocyte_like','review_queue'],'measurement_unavailable')
        self.assertEqual(r.loc['Epi_renal_tubular_like','review_queue'],'measurement_unavailable')

    def test_islet_union_does_not_require_all_hormones(self):
        r=self.report({'INS':.9,'CHGB':.5})
        self.assertTrue(r.loc['Epi_pancreatic_islet_like','each_core_component_has_measured_gene'])
        self.assertFalse(r.loc['Epi_pancreatic_islet_like','complete_core_measurement'])
        self.assertFalse(r.final_label_assigned.any())

    def test_missing_context_keeps_candidate_for_review(self):
        r=self.report({'PRSS1':.7,'PNLIP':.6})
        self.assertEqual(len(r),137)
        self.assertEqual(r.loc['Epi_pancreatic_acinar_like','review_queue'],'organ_context_review')
        self.assertFalse(r.final_label_assigned.any())

    def test_smoc2_alone_is_not_isc_panel_evidence(self):
        r=self.report({'SMOC2':.9})
        self.assertEqual(r.loc['Epi_intestinal_stem_like','review_queue'],'measurement_unavailable')

    def test_measured_zero_is_not_positive_or_missing(self):
        r=self.report({'TG':0.,'TPO':0.}, {'synthetic':{'organ_contexts':['thyroid'],'reviewed_parent_ids':['Epithelial']}})
        self.assertTrue(r.loc['Epi_thyroid_follicular_like','complete_core_measurement'])
        self.assertEqual(r.loc['Epi_thyroid_follicular_like','min_component_max_detection'],0.)
        self.assertFalse(r.final_label_assigned.any())

    def test_mixed_source_targets_do_not_cross_supply_core(self):
        r=self.report({'PGC':.7,'LIPF':.7,'SCGB1A1':.5,'UMOD':.8,'SLC12A1':.8})
        for pid in ['Epi_gastric_parietal_like','Epi_prostate_basal_like','Epi_renal_tubular_like']:
            self.assertFalse(r.loc[pid,'each_core_component_has_measured_gene'])
        self.assertTrue(r.loc['Epi_gastric_chief_like','each_core_component_has_measured_gene'])

    def test_full_source_accounting_and_original_base_preserved(self):
        source=json.loads((RESOURCE/'source_tissue_rows.json').read_text())
        trace=json.loads((RESOURCE/'original_gene_disposition.json').read_text())
        self.assertEqual({(r['label'],g) for r in source for g in r['genes'].split(',')},
                         {(r['source_label'],r['original_gene']) for r in trace})
        base=json.loads((ROOT/'resources/annotation/v7_panel_repair_v1/registry.json').read_text())
        changed=[p['id'] for p in base['panels'] if p!=self.byid[p['id']]]
        self.assertEqual(changed,['Epi_renal_tubular_like'])
        self.assertTrue(self.registry['policy']['malignancy_separate'])

    def test_rebuild_matches_published_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/'rebuilt'
            manifest=build(ROOT,out)
            for name in list(manifest['output_sha256'])+['manifest.json']:
                self.assertEqual((RESOURCE/name).read_bytes(),(out/name).read_bytes(),name)

if __name__=='__main__': unittest.main(verbosity=2)
