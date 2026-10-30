import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from TME_state_program_scope_policy_v1 import applicable

class ScopeCases(unittest.TestCase):
    def test_identity_and_state_scope_are_not_lexical_substrings(self):
        mast=dict(compartment='Mast',reference_parent='Mast',canonical_identity='Mast')
        t=dict(compartment='T_NK',reference_parent='T',canonical_identity='CD8_T_reference_related')
        md=dict(exact_Cytopus_gene_match=True,reference_parents=['mast'])
        self.assertTrue(applicable('mast',md,mast)[0]);self.assertFalse(applicable('mast',md,t)[0])
    def test_default_scope_cannot_be_widened_by_TNK_parent(self):
        md=dict(exact_Cytopus_gene_match=True,reference_parents=['TNK'])
        nk=dict(compartment='T_NK',reference_parent='NK',canonical_identity='NK_reference_related')
        cd8=dict(compartment='T_NK',reference_parent='T',canonical_identity='CD8_T_reference_related')
        self.assertFalse(applicable('CD8_T',md,nk)[0]);self.assertTrue(applicable('CD8_T',md,cd8)[0])
    def test_leukocyte_global_is_not_all_tissue(self):
        md=dict(exact_Cytopus_gene_match=True,reference_parents=['leukocyte'])
        self.assertFalse(applicable('global',md,dict(compartment='Epithelial'))[0])
        self.assertTrue(applicable('global',md,dict(compartment='Myeloid'))[0])
    def test_memory_B_does_not_include_plasma_or_naive(self):
        md=dict(exact_Cytopus_gene_match=True,reference_parents=['B'])
        for p,c in [('Plasma','Plasma_reference_related'),('B','Naive_B_reference_related')]:
            self.assertFalse(applicable('B_memory',md,dict(compartment='B_Plasma',reference_parent=p,canonical_identity=c))[0])
        self.assertTrue(applicable('B_memory',md,dict(compartment='B_Plasma',reference_parent='B',canonical_identity='Memory_B_reference_related'))[0])
    def test_unresolved_reference_is_held_not_guessed(self):
        c=dict(compartment='Myeloid',reference_parent='Myeloid',canonical_identity='Monocyte_reference_related')
        for parents,exact in [([],True),(['M'],True),(['Mac'],False)]:
            self.assertFalse(applicable('MDC',dict(reference_parents=parents,exact_Cytopus_gene_match=exact),c)[0])
    def test_all_cell_process_can_cover_compartment_only_context(self):
        self.assertTrue(applicable('global',dict(reference_parents=['all-cells'],exact_Cytopus_gene_match=True),dict(compartment='Mesenchymal_stromal'))[0])
    def test_primary_uses_explicit_resource_TNK_scope(self):
        md=dict(exact_Cytopus_gene_match=True,reference_parents=['TNK'])
        nk=dict(compartment='T_NK',reference_parent='NK',canonical_identity='NK_reference_related')
        self.assertTrue(applicable('CD8_T',md,nk,mode='resource_hierarchy_primary')[0])
        self.assertFalse(applicable('CD8_T',md,nk,mode='legacy_default_intersection')[0])
        md['reference_parents']=['CD8-T']
        self.assertFalse(applicable('CD8_T',md,nk,mode='resource_hierarchy_primary')[0])
    def test_primary_B_resource_scope_excludes_plasma(self):
        md=dict(exact_Cytopus_gene_match=True,reference_parents=['B'])
        b=dict(compartment='B_Plasma',reference_parent='B',canonical_identity='Naive_B_reference_related')
        self.assertTrue(applicable('B_memory',md,b,mode='resource_hierarchy_primary')[0])
        self.assertFalse(applicable('B_memory',md,b,mode='legacy_default_intersection')[0])
        b.update(reference_parent='Plasma',canonical_identity='Plasma_reference_related')
        self.assertFalse(applicable('B_memory',md,b,mode='resource_hierarchy_primary')[0])

if __name__=='__main__':unittest.main()
