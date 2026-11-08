from dataclasses import replace
import unittest

import numpy as np

from gui.experiments import (
    ExperimentCancelled, ExperimentRequest, ExperimentRunner, SweepAxis,
    measure_response,
)
from gui.simulation_session import SimulationRequest
from voronoi_grid import VoronoiGrid


CELLS = np.array([[0., 0.], [2., 0.], [2., 2.], [0., 2.], [1., 1.]])
ODE = '''
def helper(a, b):
    return a + b

def ode(Y, t, params, cell_positions):
    dY = np.zeros_like(Y)
    dY[:, 0] = helper(params['a'], params['b'])
    return dY
'''


class ExperimentTests(unittest.TestCase):
    def setUp(self):
        self.simulation = SimulationRequest(
            VoronoiGrid.from_cells(CELLS.copy()), .4,
            'dict(dT=.1, a=1., b=2., init_Y=[0., 7.], labels=["A", "B"])', ODE,
        )
        self.request = ExperimentRequest(
            self.simulation, SweepAxis('parameter', 'a', 1, 3, 3),
            SweepAxis('parameter', 'b', 2, 4, 2),
        )

    def test_parameter_sweep_uses_real_ode_and_y_rows_x_columns(self):
        progress = []
        result = ExperimentRunner(self.request).run(lambda percent, msg: progress.append(percent))
        np.testing.assert_allclose(result.values, [[.9, 1.2, 1.5], [1.5, 1.8, 2.1]])
        np.testing.assert_allclose(result.std, 0)
        self.assertEqual(progress[-1], 100)
        self.assertTrue(all(a <= b for a, b in zip(progress, progress[1:])))
        np.testing.assert_equal(self.simulation.grid.cells, CELLS)
        self.assertEqual(len(result.records), 6)

    def test_initial_concentration_scan_broadcasts_only_the_chosen_species(self):
        request = replace(self.request, x=SweepAxis('initial', 0, 10, 20, 2))
        result = ExperimentRunner(request).run()
        np.testing.assert_allclose(result.values, [[10.9, 20.9], [11.5, 21.5]])
        other = ExperimentRunner(replace(request, observable=1)).run()
        np.testing.assert_allclose(other.values, 7)

    def test_noise_replicates_are_reproducible_and_restore_rng(self):
        simulation = replace(self.simulation, params_source='dict(dT=.1, a=1., b=2., n_var=1, sigma=1., allow_negative=True)',
                             ode_source='def ode(Y, t, params, cell_positions):\n    return np.zeros_like(Y)')
        request = replace(self.request, simulation=simulation, seed=42, repeats=3)
        state = np.random.get_state()
        first = ExperimentRunner(request).run()
        state_after = np.random.get_state()
        np.testing.assert_equal(state[1], state_after[1])
        self.assertEqual(state[2:], state_after[2:])
        second = ExperimentRunner(request).run()
        np.testing.assert_equal(first.values, second.values)
        self.assertGreater(first.std[0, 0], 0)
        self.assertEqual(len(first.records), 18)
        # Shared replicate seeds keep random initial conditions equal across pixels.
        np.testing.assert_allclose(first.values, first.values[0, 0])

    def test_invalid_pixels_remain_nan_and_report_errors(self):
        simulation = replace(self.simulation, ode_source='''
def ode(Y, t, params, cell_positions):
    if params['a'] > 2:
        raise ValueError('unstable fixture')
    return np.zeros_like(Y)
''')
        result = ExperimentRunner(replace(self.request, simulation=simulation)).run()
        self.assertTrue(np.isnan(result.values[:, 2]).all())
        self.assertTrue(np.isfinite(result.values[:, :2]).all())
        self.assertEqual(sum(bool(record['error']) for record in result.records), 2)
        self.assertIn('unstable fixture', result.records[2]['error'])

    def test_cancel_discards_partial_sweep_and_restores_rng(self):
        calls = []
        with self.assertRaises(ExperimentCancelled):
            ExperimentRunner(self.request).run(lambda *args: calls.append(args),
                                               lambda: len(calls) >= 1)
        self.assertTrue(calls)

    def test_log_axes_and_validation(self):
        np.testing.assert_allclose(SweepAxis('parameter', 'a', 1, 100, 3, 'log').values(), [1, 10, 100])
        with self.assertRaisesRegex(ValueError, 'positive'):
            SweepAxis('parameter', 'a', 0, 1, 2, 'log').values()
        with self.assertRaisesRegex(ValueError, 'different'):
            ExperimentRunner(replace(self.request, y=self.request.x))
        with self.assertRaisesRegex(ValueError, 'scalar'):
            ExperimentRunner(replace(self.request, x=SweepAxis('parameter', 'labels', 0, 1, 2))).run()

    def test_statistics_are_defined_for_variable_cell_counts(self):
        history = [np.array([[0.], [2.]]), np.array([[2.], [4.], [6.]]), np.array([[4.], [8.]])]
        self.assertEqual(measure_response(history, 0, 'final_mean', 1), 6)
        self.assertEqual(measure_response(history, 0, 'final_std', 1), 2)
        self.assertAlmostEqual(measure_response(history, 0, 'final_cv', 1), 1 / 3)
        self.assertEqual(measure_response(history, 0, 'tail_mean', 2 / 3), 5)
        self.assertEqual(measure_response(history, 0, 'tail_amplitude', 1), 5)
        with self.assertRaisesRegex(ValueError, 'undefined'):
            measure_response([np.zeros((5, 1))], 0, 'final_cv', 1)
        with self.assertRaisesRegex(ValueError, 'non-finite'):
            measure_response([np.full((5, 1), np.inf)], 0, 'final_mean', 1)


if __name__ == '__main__':
    unittest.main()
