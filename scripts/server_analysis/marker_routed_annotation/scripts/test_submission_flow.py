"""Scheduler orchestration tests: every scheduler command is mocked."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


class SubmissionTests(unittest.TestCase):
    def run_case(self, *, failed_parent=False, active=False, auth_failure=False):
        with tempfile.TemporaryDirectory(prefix='tme-submit-test-') as tmp:
            root = Path(tmp); task = root / 'task'; scripts = task / 'scripts'
            scripts.mkdir(parents=True)
            for script in Path(__file__).parent.iterdir():
                if script.suffix in ['.py', '.R', '.sbatch']:
                    shutil.copy2(script, scripts / script.name)
            pc = root / 'pc'; shards = pc / 'integration_v2/run_v2/shards'
            shards.mkdir(parents=True); malig = pc / 'malignancy_v1'; malig.mkdir()
            repo = root / 'repo'; (repo / 'scripts').mkdir(parents=True)
            (repo / 'scripts/CellTypist_official_unicode_axis_v1.py').write_text('# fixture\n')
            pending = ['GSE125449', 'GSE183904', 'GSE253013', 'GSE269936']
            names = pending + [f'FIXTURE{i:03d}' for i in range(68)]
            paths = []
            for name in names:
                path = shards / (name + '.h5ad'); path.write_text('fixture')
                paths.append(str(path))
                if name not in pending:
                    (malig / (name + '.tsv.gz')).write_text('fixture')
            (task / 'sources.txt').write_text('\n'.join(paths) + '\n')
            spec = importlib.util.spec_from_file_location('submit_fixture', scripts / 'submit_all.py')
            module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
            calls = []

            def execute(cmd, **kwargs):
                calls.append(cmd)
                if cmd[0] == 'munge':
                    return subprocess.CompletedProcess(cmd, 1 if auth_failure else 0, '', 'fixture authentication failure' if auth_failure else '')
                if cmd[:2] == ['scontrol', 'ping']:
                    return subprocess.CompletedProcess(cmd, 0, 'Slurmctld(primary) at fixture is UP\n', '')
                if cmd[0] == 'squeue':
                    return subprocess.CompletedProcess(cmd, 0, '999|marker-routed-tme-v1|RUNNING\n' if active else '', '')
                if cmd[:3] == ['scontrol', 'show', 'job']:
                    state = 'FAILED' if failed_parent and cmd[3] == '30645863_0' else 'RUNNING'
                    return subprocess.CompletedProcess(cmd, 0, 'JobState=' + state, '')
                if cmd[0] == 'sbatch':
                    return subprocess.CompletedProcess(cmd, 0, str(1000 + len([c for c in calls if c[0] == 'sbatch'])) + '\n', '')
                raise AssertionError('Unexpected subprocess: ' + repr(cmd))

            argv = ['submit_all.py', '--submit', '--run-dir', str(task / 'run'),
                    '--repo', str(repo), '--output-root', str(root / 'out')]
            with patch.object(module.pwd, 'getpwuid', return_value=object()), \
                 patch.object(module.subprocess, 'run', side_effect=execute), \
                 patch.object(sys, 'argv', argv), contextlib.redirect_stdout(io.StringIO()):
                if active or auth_failure:
                    with self.assertRaises(RuntimeError):
                        module.main()
                    self.assertFalse((task / 'run').exists())
                    self.assertFalse(any(c[0] == 'sbatch' for c in calls))
                    return
                self.assertEqual(module.main(), 0)
            receipt = json.loads((task / 'run/submission_receipt.json').read_text())
            jobs = [c for c in calls if c[0] == 'sbatch']
            self.assertEqual(len(jobs), 5 if failed_parent else 6)
            self.assertIn('--array=0-67%4', jobs[0])
            dependencies = [a for c in jobs for a in c if a.startswith('--dependency=afterok:')]
            self.assertEqual(len(dependencies), 3 if failed_parent else 4)
            self.assertEqual(len(receipt['blocked_sources']), int(failed_parent))
            self.assertTrue(any(a.startswith('--dependency=afterany:') for a in jobs[-1]))
            self.assertFalse(any(c[0] in ['scancel', 'srun'] for c in calls))

    def test_ready_and_four_live_dependencies(self):
        self.run_case()

    def test_failed_parent_is_recorded_without_rerun(self):
        self.run_case(failed_parent=True)

    def test_existing_jobs_prevent_duplicate_submission(self):
        self.run_case(active=True)

    def test_auth_failure_prevents_all_submission(self):
        self.run_case(auth_failure=True)


if __name__ == '__main__':
    unittest.main()
