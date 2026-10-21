"""Source-case selection, outcome orientation and chunked real-kernel replay."""
import unittest
import numpy as np
import GSE288199_state_pair_math_v1 as math
from NSCLC_observed_T_state_statistics_v1 import select,value,controls,agent,association,stratified_indices,replay


def fixture():
    rng=np.random.default_rng(103);rows=[]
    for qc in ['native_QC500_250_MT20','all_retained_native']:
        for parent in ['dual_regulatory_evidence','costate_regulatory_sensitivity']:
            for i in range(243):
                response=i%2;state=max(0,min(60,int(12+30*response+rng.normal(0,4))))
                rows.append(dict(case_index=i,QC=qc,parent=parent,metrics=[300,60,100,0,140,5,state,7,12,10,80-state,0,0],
                    clinical_eligible=True,response=response,cancer_type='LUAD' if i<80 else 'LUSC',drug='DrugA(C4)',
                    clinical=rng.normal(size=8).tolist(),quality={k:rng.uniform(1,4,size=4).tolist() for k in ['R','E','O']},
                    parent_counts=dict(R=60,E=100,O=50),other_T_states=[1,2,3]))
    return dict(rows=rows)


class StateStatisticsChecks(unittest.TestCase):
    def test_no_selected_discovery_subset_and_observed_zero_preserved(self):
        p=fixture();p['rows'][0]['metrics'][6]=0;p['rows'][1]['metrics'][1]=19
        rr,cov=select(p)
        self.assertEqual(len(rr),242);self.assertEqual(rr[0]['case_index'],0)
        self.assertEqual(value(rr,'R_LAG3')[0],0)

    def test_agent_suffix_does_not_merge_regimens(self):
        self.assertEqual(agent('DrugA(C4-C5)'),'DrugA')
        self.assertNotEqual(agent('DrugA+DrugB(C4)'),agent('DrugA(C4)'))
        self.assertEqual(agent('Yes'),'Yes')
        p=fixture()
        for r in p['rows']:
            if r['case_index']<5:r['drug']='RareDrug'
        rr,coverage=select(p);self.assertEqual(len(rr),238)
        self.assertEqual(coverage['removed_for_agent_group_below10'],5)

    def test_case_and_binary_group_gates(self):
        rows,_=select(fixture());r=rows[:49]
        m,f,why=association(r,np.arange(49),np.arange(49),np.empty((49,0)))
        self.assertIsNone(f);self.assertEqual(why,'FEWER_THAN50_AUTHOR_CASE_UNITS')
        r=[dict(v,response=1 if i<9 else 0) for i,v in enumerate(rows[:60])]
        _,f,why=association(r,np.arange(60),np.arange(60),np.empty((60,0)))
        self.assertIsNone(f);self.assertEqual(why,'FEWER_THAN10_IN_ONE_PATHOLOGY_GROUP')

    def test_clinical_direction_own_quality_and_complete_random_replay(self):
        rows,_=select(fixture());rows=rows[:80]
        c=controls(rows,['R']);x=np.asarray([r['response'] for r in rows]);y=value(rows,'R_LAG3')
        np.testing.assert_allclose(c[0,8:12],rows[0]['quality']['R'])
        _,f,status=association(rows,x,y,c);self.assertEqual(status,'ESTIMABLE');self.assertGreater(f['rho'],.5)
        _,negative,status=association(rows,x,value(rows,'E_memory'),controls(rows,['E']))
        self.assertLess(negative['rho'],-.5)
        indices=stratified_indices(x,30,'synthetic')
        for take in indices:self.assertEqual(int(x[take].sum()),int(x.sum()))
        rng=np.random.default_rng(708);signs=np.where(rng.random((1001,len(rows)))<.5,-1,1).astype(np.int8)
        m,_,_=association(rows,x,y,c);wild=math.wild_statistics(m,f,signs);boot,_=math.bootstrap(x,y,c,indices)
        proof=replay('synthetic_clinical',x,y,c,f,signs,wild,indices,boot)
        self.assertEqual(proof['full_wild_statistics_checked'],1001)
        self.assertEqual(proof['full_bootstrap_statistics_checked'],30)
        self.assertLess(proof['max_numeric_error'],1e-8)


if __name__=='__main__':unittest.main()
