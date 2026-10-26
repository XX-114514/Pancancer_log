from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from atlas_identity_release_keys_v1 import membership_key,assign_display_numbers,select_display_marker,display_label,version_crosswalk


class ReleaseKeyChecks(unittest.TestCase):
    def test_membership_order_namespace_and_exact_boundaries(self):
        a=membership_key('source1',['lib1:a','lib1:b'])
        self.assertEqual(a,membership_key('source1',['lib1:b','lib1:a']))
        self.assertNotEqual(a,membership_key('source2',['lib1:a','lib1:b']))
        self.assertNotEqual(membership_key('s',['ab','c']),membership_key('s',['a','bc']))
        self.assertNotEqual(a,membership_key('source1',['lib1:a']))
        for bad in [[],['x','x'],['']]:
            with self.assertRaises(ValueError):membership_key('s',bad)

    def test_registry_preserves_inactive_numbers_and_appends_new(self):
        a,b,c=[membership_key('s',[s]) for s in ['a','b','c']]
        old={a:5,b:9};r=assign_display_numbers([c,a],old)
        self.assertEqual(r['current'][a],5);self.assertEqual(r['current'][c],10)
        self.assertEqual(r['registry'][b],9);self.assertEqual(r['inactive_historical_keys'],[b]);self.assertEqual(old,{a:5,b:9})
        self.assertEqual(assign_display_numbers([a,c],old)['registry'],r['registry'])
        with self.assertRaises(ValueError):assign_display_numbers([a,a])
        with self.assertRaises(ValueError):assign_display_numbers([a],{a:1,b:1})

    def test_marker_precedence_membership_and_explicit_parent(self):
        key=membership_key('s',['a']);parent=membership_key('s',['a','b'])
        expr={'gene':'KRT17','membership_key':key,'review_status':'accepted','marker_evidence_type':'own_subset_expression_specificity'}
        deg={'gene':'KRT14','membership_key':key,'review_status':'accepted','marker_evidence_type':'own_subset_valid_DEG'}
        self.assertEqual(select_display_marker([expr,deg],key)['marker'],'KRT14')
        self.assertEqual(select_display_marker([expr],key)['marker_evidence_type'],'own_subset_expression_specificity')
        with self.assertRaises(ValueError):select_display_marker([expr],parent)
        borrowed={'gene':'KRT5','review_status':'accepted','marker_evidence_type':'parent_group_borrowed_explicit','parent_membership_key':parent,'borrow_reason':'own marker unavailable'}
        self.assertIsNone(select_display_marker([borrowed],key)['marker'])
        self.assertTrue(select_display_marker([borrowed],key,allow_parent=True)['parent_marker_borrowed'])
        pending={**expr,'review_status':'pending'}
        self.assertIsNone(select_display_marker([pending],key)['marker'])

    def test_display_is_one_gene_and_number_without_hash(self):
        self.assertEqual(display_label('Basal_epithelial_cell','KRT17',12),'Basal_epithelial_cell_KRT17_12')
        self.assertEqual(display_label('T','HLA-DRA',1),'T_HLA-DRA_1')
        self.assertEqual(display_label('unresolved',None,2),'unresolved_marker_unresolved_2')
        for marker in ['KRT17_KRT14','KRT17 KRT14','KRT17/KRT14']:
            with self.assertRaises(ValueError):display_label('T',marker,1)

    def test_crosswalk_many_to_many_unresolved_and_archive_exclusion(self):
        r=version_crosswalk('s',['l:a','l:b','l:c','l:d'],['oldA','oldA','oldB','oldC'],['newX','newY','newX','unresolved'],['accepted','accepted','accepted','archive_excluded'])
        self.assertEqual(sum(x['n_archive_records'] for x in r),4)
        self.assertEqual(next(x for x in r if x['old_label']=='oldA' and x['new_label']=='newX')['relation'],'many_to_many')
        self.assertEqual(next(x for x in r if x['old_label']=='oldA' and x['new_label']=='newY')['relation'],'split')
        self.assertEqual(next(x for x in r if x['old_label']=='oldB')['relation'],'merge')
        self.assertEqual(next(x for x in r if x['old_label']=='oldC')['new_identity_status'],'archive_excluded')
        with self.assertRaises(ValueError):version_crosswalk('s',['a','a'],['A','B'],['C','D'],['accepted']*2)


if __name__=='__main__':unittest.main()
