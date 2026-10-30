"""Regression cases for WEHI sample metadata: GEO CRLF text, patient/run parsing, replicate grouping."""
from pathlib import Path
import csv, os, sys, tempfile, unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts' / 'wehi_xenium'))
from fetch_geo_samples import parse
import make_samples

SOFT = ('^SAMPLE = GSM1\r\n!Sample_title = Patient 9, run 2\r\n!Sample_source_name_ch1 = HNSCC tumor\r\n'
        '!Sample_characteristics_ch1 = p16 status: Positive\r\n'
        '^SAMPLE = GSM2\r\n!Sample_title = Ameloblastoma section\r\n!Sample_source_name_ch1 = Ameloblastoma\r\n')


class GeoParseTests(unittest.TestCase):
    def test_crlf_does_not_leak_into_values(self):
        rows = parse(SOFT)
        self.assertEqual([r['gsm'] for r in rows], ['GSM1', 'GSM2'])
        for r in rows:
            for v in r.values():
                self.assertNotIn('\r', v)
        self.assertEqual((rows[0]['patient'], rows[0]['run']), ('Patient 9', '2'))
        self.assertEqual(rows[0]['char_p16_status'], 'Positive')

    def test_title_without_run_leaves_unit_empty(self):
        rows = parse(SOFT)
        self.assertEqual((rows[1]['patient'], rows[1]['run']), ('', ''))


class SampleTableTests(unittest.TestCase):
    def test_replicate_runs_share_patient_and_unknown_stays_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            meta = tmp / 'meta' / 'geo'
            meta.mkdir(parents=True)
            with open(meta / 'samples_GSE1.tsv', 'w', newline='') as fh:
                w = csv.DictWriter(fh, ['gsm', 'patient', 'run'], delimiter='\t')
                w.writeheader()
                w.writerows([{'gsm': 'GSM1', 'patient': 'Patient 9', 'run': '1'},
                             {'gsm': 'GSM2', 'patient': 'Patient 9', 'run': '2'},
                             {'gsm': 'GSM3', 'patient': '', 'run': ''}])
            raw = tmp / 'raw' / 'S18_GSE1_HNSCC'
            for gsm in ['GSM1', 'GSM2', 'GSM3']:
                (raw / gsm).mkdir(parents=True)
                (raw / gsm / f'{gsm}_cell_feature_matrix.h5').write_bytes(b'')
            out = tmp / 'samples.tsv'
            old = os.environ.get('META_ROOT')
            os.environ['META_ROOT'] = str(tmp / 'meta')
            try:
                make_samples.main(tmp / 'raw', out)
            finally:
                if old is None:
                    os.environ.pop('META_ROOT')
                else:
                    os.environ['META_ROOT'] = old
            rows = {r['sample_id'].split('__')[-1]: r for r in csv.DictReader(open(out), delimiter='\t')}
            self.assertEqual(rows['GSM1']['patient'], 'GSE1_P9')
            self.assertEqual(rows['GSM2']['patient'], 'GSE1_P9')
            self.assertEqual(rows['GSM2']['run'], '2')
            self.assertEqual(rows['GSM3']['patient'], '')


if __name__ == '__main__':
    unittest.main()
