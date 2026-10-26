import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from native_specificity_shadow_v1 import FAMILY_COMPONENTS,ANCHORS,shadow_support

class SpecificityShadowTests(unittest.TestCase):
    def calls(self,*active):
        return {p:p in active for ps in FAMILY_COMPONENTS.values() for p in ps}
    def run_case(self,calls,detected=(),measured=None,scenario='both_anchors'):
        if measured is None:measured={g for gs in ANCHORS.values() for g in gs}
        return shadow_support(calls,detected,measured,scenario)
    def test_shared_cytotoxic_Fc_signal_is_diagnostic_only(self):
        r=self.run_case(self.calls('NK_related','Myeloid'),['LYZ'])
        self.assertEqual(r['original_families'],{'T_NK','Myeloid'})
        self.assertEqual(r['shadow_families'],{'Myeloid'})
        self.assertFalse(r['identity_assigned']);self.assertFalse(r['doublet_or_ambient_diagnosed'])
    def test_T_support_is_not_removed_with_NK_component(self):
        r=self.run_case(self.calls('T','NK_related','Myeloid'),['LYZ'])
        self.assertEqual(r['shadow_families'],{'T_NK','Myeloid'})
    def test_DC_neutrophil_components_preserve_family(self):
        for panel in ['DC_related','Neutrophil_related']:
            self.assertEqual(self.run_case(self.calls('Myeloid',panel))['shadow_families'],{'Myeloid'})
    def test_missing_anchor_not_a_negative_measurement(self):
        r=self.run_case(self.calls('NK_related'),measured={'GNLY'})
        self.assertEqual(r['shadow_families'],{'T_NK'})
        self.assertFalse(r['anchor_review']['NK_related']['all_anchors_measured'])
    def test_true_two_family_support_is_preserved(self):
        r=self.run_case(self.calls('NK_related','Myeloid'),['GNLY','LYZ'])
        self.assertEqual(r['shadow_families'],{'T_NK','Myeloid'})
    def test_other_families_and_baseline_unchanged(self):
        c=self.calls('B','Plasma','Epithelial','NK_related','Myeloid')
        r=self.run_case(c,scenario='baseline')
        self.assertEqual(r['shadow_families'],r['original_families'])
        self.assertEqual(self.run_case(c)['shadow_families'],{'B_Plasma','Epithelial'})
    def test_single_scenarios_do_not_apply_both(self):
        c=self.calls('NK_related','Myeloid')
        self.assertEqual(self.run_case(c,scenario='NK_anchor')['shadow_families'],{'Myeloid'})
        self.assertEqual(self.run_case(c,scenario='myeloid_anchor')['shadow_families'],{'T_NK'})
        self.assertEqual(self.run_case(c)['shadow_families'],set())
    def test_invalid_measurement_or_incomplete_inputs_rejected(self):
        with self.assertRaises(ValueError):self.run_case(self.calls('NK_related'),['GNLY'],set())
        with self.assertRaises(ValueError):self.run_case({'NK_related':True})
        with self.assertRaises(ValueError):self.run_case(self.calls(),scenario='optimized_winner')

if __name__=='__main__':unittest.main()
