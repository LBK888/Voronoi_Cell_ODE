"""The public launcher exercises the same headless entry point as the EXE."""

import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class AutomationCliTests(unittest.TestCase):
    def test_cli_opt_in_inspection_and_jobs(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            environment = {**os.environ, 'VORONOI_CELL_ODE_DATA_DIR': str(folder / 'settings')}
            report = folder / 'report.json'
            model = folder / 'model.json'
            config = {
                'grid_mode': 3, 'cell_positions': [[0, 0], [2, 0], [2, 2], [0, 2], [1, 1]],
                'T': '0.25', 'params': 'dict(dT=0.1, a=1, init_Y=[0.0])',
                'ode': 'def ode(Y, t, params, cell_positions):\n    return np.full_like(Y, params["a"])',
                'experiments': {
                    'x': {'target': 'a', 'minimum': 1, 'maximum': 2, 'count': 2},
                    'y': {'kind': 'initial', 'target': 0, 'minimum': 0, 'maximum': 1, 'count': 2},
                },
            }
            model.write_text(json.dumps(config), encoding='utf-8')

            def call(*args):
                process = subprocess.run(
                    [sys.executable, str(ROOT / 'windows_launcher.py'), '--cli', '--report', str(report), *map(str, args)],
                    cwd=folder, env=environment, capture_output=True, text=True, timeout=30,
                )
                return process.returncode, json.loads(report.read_text(encoding='utf-8'))

            code, response = call('simulate', model, '--output', folder / 'blocked')
            self.assertNotEqual(code, 0)
            self.assertIn('Enable Agent CLI', response['error'])
            self.assertFalse((folder / 'blocked').exists())
            (folder / 'settings' / 'automation.json').write_text('{"cli_enabled": true}', encoding='utf-8')

            inspection = folder / 'inspect.json'
            inspection.write_text(json.dumps({**config, 'ode': 'raise RuntimeError("must not execute")'}), encoding='utf-8')
            self.assertEqual(call('inspect', inspection)[0], 0)
            self.assertTrue(call('list-templates')[1]['templates'])
            simulation = folder / 'simulation'
            self.assertEqual(call('simulate', model, '--output', simulation)[0], 0)
            with (simulation / 'cells.csv').open(encoding='utf-8-sig', newline='') as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 15)
            self.assertAlmostEqual(float(rows[-1]['Y[0]']), 0.2)
            self.assertNotEqual(call('simulate', model, '--output', simulation)[0], 0)

            experiment = folder / 'experiment'
            code, response = call('experiment', model, '--output', experiment, '--seed', 8)
            self.assertEqual(code, 0, response)
            self.assertEqual(response['shape'], [2, 2])
            self.assertEqual(response['seed'], 8)
            self.assertEqual(response['failed_runs'], 0)
            self.assertTrue((experiment / 'heatmap.pdf').read_bytes().startswith(b'%PDF'))
