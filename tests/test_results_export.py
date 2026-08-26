import json
import os
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
from PyQt6.QtCore import QCoreApplication

from gui.results_export_worker import ExportRequest, ResultsExportWorker
from gui.simulation_session import SimulationResult
from voronoi_grid import VoronoiGrid


class ResultsExportInterfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QCoreApplication.instance() or QCoreApplication([])

    def test_run_writes_one_consistent_result_bundle(self):
        cells = np.array([
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, 1.0],
        ])
        grid = VoronoiGrid.from_cells(cells)
        simulation = SimulationResult(
            history=[np.full((5, 2), value) for value in range(4)],
            cell_positions=[cells.copy() for _ in range(4)],
            grid=grid,
            params={'dT': 0.5, 'labels': ['A', 'B']},
            cell_ids=[np.arange(10, 15) for _ in range(4)],
        )

        with tempfile.TemporaryDirectory() as folder:
            request = ExportRequest(
                folder=folder,
                prefix='test_',
                config={'name': 'fixture'},
                simulation=simulation,
                color_func=lambda values, mode='polygon': ['red'] * len(values),
                color_func_source='',
            )
            worker = ResultsExportWorker(request)
            worker._save_video = lambda: Path(
                folder, 'test_simulation.mp4'
            ).write_bytes(b'video-stand-in')
            worker._save_concentration_plot = lambda: Path(
                folder, 'test_concentration.pdf'
            ).write_bytes(b'plot-stand-in')

            completed = []
            errors = []
            worker.export_finished.connect(completed.append)
            worker.error_occurred.connect(errors.append)
            worker.run()

            self.assertEqual(errors, [])
            self.assertEqual(completed, [folder])
            expected = {
                'test_config.json',
                'test_simulation.mp4',
                'test_concentration.pdf',
                'test_cells.xlsx',
                'test_cells.h5',
            }
            self.assertTrue(expected.issubset(set(os.listdir(folder))))
            self.assertNotIn('test_cells.csv', os.listdir(folder))
            with open(Path(folder, 'test_config.json'), encoding='utf-8') as handle:
                self.assertEqual(json.load(handle), {'name': 'fixture'})
            workbook = pd.read_excel(Path(folder, 'test_cells.xlsx'))
            self.assertEqual(
                list(workbook.columns),
                ['cell_id', 'T', 'step', 't', 'x', 'y', 'Y[0]', 'Y[1]'],
            )
            self.assertEqual(workbook['cell_id'].tolist(), list(range(10, 15)) * 4)
            self.assertEqual(
                workbook['t'].tolist(),
                [0] * 5 + [1] * 5 + [2] * 5 + [3] * 5,
            )
            self.assertEqual(workbook['T'].tolist(), [0] * 10 + [1] * 10)
            self.assertEqual(
                workbook['step'].tolist(),
                [0] * 5 + [1] * 5 + [0] * 5 + [1] * 5,
            )
            with h5py.File(Path(folder, 'test_cells.h5'), 'r') as handle:
                self.assertEqual(handle['cells'].shape, (20, 8))
                columns = [
                    value.decode() if isinstance(value, bytes) else value
                    for value in handle['cells'].attrs['columns']
                ]
                self.assertEqual(columns, list(workbook.columns))


if __name__ == '__main__':
    unittest.main()
