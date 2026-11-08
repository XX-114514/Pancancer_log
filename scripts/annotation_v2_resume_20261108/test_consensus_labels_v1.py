import unittest
import numpy as np
from consensus_labels_v1 import consensus

class ConsensusTests(unittest.TestCase):
 def test_identity_preserved(self):
  y=np.array(['A']*4+['B']*4);z=np.tile(np.array(['0']*4+['1']*4)[:,None],(1,10))
  lab,st,ties,_,_,comp=consensus(y,z)
  np.testing.assert_array_equal(lab,y);self.assertTrue((st==1).all());self.assertFalse(ties.any());self.assertEqual(len(comp),2)
 def test_bidirectional_denominator_and_ties(self):
  y=np.array(['A']*4+['B']*4);z=np.full((8,10),'0')
  _,_,_,pairs,mapping,comp=consensus(y,z)
  self.assertEqual(pairs[0]['bidirectional_mean'],.5);self.assertEqual(set(comp),{'A|B'});self.assertTrue(all(x['tied_majority'] for x in mapping))
 def test_unrelated_identity_not_merged(self):
  y=np.array(['A']*4+['B']*4+['C']*4);z=np.tile(np.array(['0']*8+['1']*4)[:,None],(1,10))
  *_,comp=consensus(y,z);self.assertEqual(set(comp),{'A|B','C'})
 def test_exact_threshold_transitive_union(self):
  y=np.array(['A']*2+['B']*4+['C']*4)
  z=np.array([['0']*6+['1']*4]*5+[['0']*2+['1']*4+['0']*4]*5).T
  _,_,_,pairs,_,comp=consensus(y,z)
  self.assertEqual(set(comp),{'A|B|C'});self.assertEqual(sum(x['merge'] for x in pairs),2)
if __name__=='__main__':unittest.main()
