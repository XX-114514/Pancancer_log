import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_native_author_label_evidence_v1 import author_context


class AuthorEvidenceTests(unittest.TestCase):
    def test_explicit_context_can_recover_nonimmune_author(self):
        self.assertEqual(author_context(' Hepatocyte '), 'Hepatocyte')
        self.assertEqual(author_context('natural killer cell'), 'T_NK')

    def test_malignancy_does_not_infer_epithelial_lineage(self):
        self.assertEqual(author_context('malignant cells'), 'author_malignancy_label_lineage_unspecified')
        self.assertEqual(author_context('epithelial cells'), 'Epithelial')

    def test_mixed_and_marker_suffix_labels_are_not_auto_collapsed(self):
        for value in ['T cells/B cells', 'Tcell_CCL5_3', 'Cycling cells', 'hepatocyte-like tumor']:
            self.assertEqual(author_context(value), 'unreviewed_author_value')

    def test_missing_and_unreviewed_are_distinct(self):
        self.assertEqual(author_context('NA'), 'missing_or_explicit_unknown')
        self.assertEqual(author_context('cluster 12'), 'unreviewed_author_value')

    def test_merged_source_missing_sentinels_are_not_annotation_coverage(self):
        for value in ['__SOURCE_VALUE_MISSING__', '__MISSING__', 'AUTHOR_LABEL_MISSING', 'SOURCE_MISSING']:
            self.assertEqual(author_context(value), 'missing_or_explicit_unknown')


if __name__ == '__main__':
    unittest.main()
