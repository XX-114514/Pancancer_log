"""Clinical key and censoring semantics, without real patient fixtures."""
import unittest
from qualify_Chen_coregistered_clinical_scope_v1 import key,numeric,summarize


class MetadataChecks(unittest.TestCase):
    def test_no_prefix_or_fractional_ID_guess(self):
        for x in [True,False,'Pt12','01',1.1,-1,'ID12',None]:self.assertIsNone(key(x))
        for x in [12,12.,'12',' 12 ']:self.assertEqual(key(x),'12')

    def test_finite_duration_and_known_event_flag(self):
        for x in [None,True,float('nan'),float('inf'),'unknown']:self.assertIsNone(numeric(x))
        records=[]
        for response,time,event in [('CR',50,1),('SD',60,0),('PD',-2,1),('mPR',40,7)]:
            records.append({'RECIST':response,'PFS (days)':time,'PFS censored (0=yes)':event,
                'OS (days)':time,'OS Censored (0=not dead)':event})
        d=summarize(records)
        self.assertEqual(d['RECIST_evaluable'],3)
        self.assertEqual(d['CR_PR'],1)
        self.assertEqual(d['SD_PD'],2)
        self.assertEqual(d['reported_survival_metadata']['OS']['positive_reported_duration_and_binary_flag'],2)
        self.assertEqual(d['reported_survival_metadata']['OS']['events_code1'],1)


if __name__=='__main__':unittest.main()
