"""Counterexamples that previously produced misleading labels."""
import json,unittest
from pathlib import Path
import pandas as pd
import numpy as np
from evidence_engine import evaluate,resolve_symbols,validate_registry
REG=Path(__file__).resolve().parents[2]/'resources/annotation/v7_panel_repair_v1/registry.json'

class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=json.loads(REG.read_text())
    def expr(self,genes,missing=()):
        return pd.DataFrame([dict(cluster='case',gene=g,n_cells=100,n_measured=0 if g in missing else 100,pct_expressing_measured=v) for g,v in genes.items()])
    def get(self,genes,missing=(),de=None):
        return evaluate(self.r,self.expr(genes,missing),de).set_index('panel')
    def test_registry_is_acyclic_all_refs_and_signs(self):validate_registry(self.r)
    def test_cytotoxic_t_does_not_become_NK(self):
        d=self.get({'CD3D':.95,'CD3E':.9,'TRAC':.8,'NKG7':.95,'GNLY':.5,'KLRD1':.6,'GZMK':.45,'GZMA':.6,'CCL5':.99})
        self.assertTrue(d.loc['T','each_core_component_has_measured_gene'])
        self.assertTrue(d.loc['GZMK_program','each_core_component_has_measured_gene'])
        self.assertEqual(d.loc['T_CD8','review_queue'],'measurement_unavailable')
        self.assertFalse(d.final_label_assigned.any())
    def test_tigit_does_not_supply_CD8_or_Treg_core(self):
        d=self.get({'TIGIT':.9,'CTLA4':.8})
        self.assertFalse(d.loc['T_Treg','each_core_component_has_measured_gene'])
        self.assertFalse(d.loc['Exhaustion_like','each_core_component_has_measured_gene'])
    def test_all_candidates_retained_Treg_not_shadowed_by_TIGIT(self):
        genes={'FOXP3':.5,'IL2RA':.4,'IKZF2':.4,'TIGIT':.99,'CD3D':.9,'TRAC':.8,'MKI67':.3,'TOP2A':.4,'UBE2C':.3}
        a=self.get(genes);genes['TIGIT']=.9999;b=self.get(genes)
        self.assertEqual(a.loc['T_Treg','min_component_max_detection'],b.loc['T_Treg','min_component_max_detection'])
        self.assertEqual(len(a),len(self.r['panels']))
        self.assertTrue(a.loc['Cycling','each_core_component_has_measured_gene']);self.assertFalse(a.final_label_assigned.any())
    def test_measurement_missing_not_zero(self):
        a=self.get({'FOXP3':0.,'IL2RA':.5},missing=['FOXP3'])
        b=self.get({'FOXP3':0.,'IL2RA':.5})
        self.assertFalse(a.loc['T_Treg','each_core_component_has_measured_gene'])
        self.assertTrue(b.loc['T_Treg','each_core_component_has_measured_gene'])
        self.assertTrue(np.isnan(a.loc['T_Treg','min_component_max_detection']))
        self.assertEqual(b.loc['T_Treg','min_component_max_detection'],0)
    def test_partial_core_not_complete_measurement(self):
        d=self.get({'FOXP3':.5,'IL2RA':.6})
        self.assertTrue(d.loc['T_Treg','each_core_component_has_measured_gene'])
        self.assertFalse(d.loc['T_Treg','complete_core_measurement'])
    def test_named_marker_not_replaced_by_shared_genes(self):
        d=self.get({'HDC':.9,'LTC4S':.8,'KIT':.8})
        self.assertEqual(d.loc['Mast_IL1RL1','naming_status'],'naming_marker_unmeasured')
    def test_adaptive_direction_and_TCR(self):
        p=validate_registry(self.r)['NK_adaptive_like']
        self.assertIn('FCER1G',p['expected_lower_genes']);self.assertNotIn('FCER1G',p['positive_genes']);self.assertNotIn('TRDC',p['positive_genes'])
    def test_alias_collision_refused(self):
        with self.assertRaises(ValueError):resolve_symbols(['JCHAIN','IGJ'],self.r['aliases'])
        self.assertEqual(resolve_symbols(['NDUFB8.1','NDUFB8-1'],self.r['aliases']),['NDUFB8.1','NDUFB8-1'])
    def test_legacy_log2fc_is_not_Seurat_DEG(self):
        with self.assertRaises(ValueError):self.get({'CD3D':.9},de=pd.DataFrame([dict(cluster='case',gene='CD3D',log2fc=.8)]))
    def test_DEG_alias_uses_same_versioned_map(self):
        de=pd.DataFrame([{'cluster':'case','gene':'IGJ','avg_log2FC':1.,'pct.1':.8,'pct.2':.1,'p_val_adj':.001}])
        d=self.get({'JCHAIN':.8},de=de)
        self.assertIn('JCHAIN',d.loc['ASC','positive_DEG_genes'])
    def test_invalid_denominator_rejected(self):
        df=self.expr({'FOXP3':.8});df.loc[0,'n_measured']=np.nan
        with self.assertRaises(ValueError):evaluate(self.r,df)
    def test_organ_required_and_no_parent_inheritance(self):
        d=self.get({'SFTPC':.8,'SFTPB':.7,'ABCA3':.2})
        self.assertEqual(d.loc['Epi_AT2','organ_context'],'organ_context_missing')
        self.assertEqual(d.loc['Epi_AT2','parent_review'],'parent_not_independently_reviewed')
    def test_component_marginals_not_same_cell_coexpression(self):
        d=self.get({'FOXP3':.5,'IL2RA':.5})
        self.assertFalse(d.loc['T_Treg','component_detection_is_joint_cell_fraction'])
    def test_single_CD1C_does_not_complete_DC3(self):
        d=self.get({'CD1C':.9})
        self.assertFalse(d.loc['DC_DC3','each_core_component_has_measured_gene'])
if __name__=='__main__':unittest.main()
