import unittest

import numpy as np

from gui.simulation_session import (
    SimulationConfigurationError,
    SimulationRequest,
    SimulationSession,
)
from voronoi_grid import VoronoiGrid


ODE_SOURCE = '''
def ode(Y, t, params, cell_positions):
    return np.zeros_like(Y)
'''


class SimulationSessionInterfaceTests(unittest.TestCase):
    def setUp(self):
        self.grid = VoronoiGrid.from_cells(np.array([
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, 1.0],
        ]))

    def test_run_prepares_model_and_returns_one_result(self):
        request = SimulationRequest(
            grid=self.grid,
            duration=0.3,
            params_source='params = dict(dT=0.1, n_var=2, sigma=0.0)',
            ode_source=ODE_SOURCE,
        )
        result = SimulationSession(request).run()
        self.assertEqual(len(result.history), 3)
        self.assertEqual(len(result.cell_positions), 3)
        self.assertEqual(len(result.cell_ids), 3)
        self.assertEqual(result.cell_ids[0].tolist(), list(range(5)))
        self.assertEqual(result.cell_ids[-1].tolist(), list(range(5)))
        self.assertEqual(result.history[-1].shape, (5, 2))
        self.assertIs(result.grid, self.grid)
        self.assertEqual(result.params['n_var'], 2)

    def test_initial_state_vector_is_broadcast_to_all_cells(self):
        request = SimulationRequest(
            grid=self.grid,
            duration=0.2,
            params_source='dict(dT=0.1, init_Y=[1.0, 2.0])',
            ode_source=ODE_SOURCE,
        )
        result = SimulationSession(request).run()
        expected = np.tile([1.0, 2.0], (5, 1))
        self.assertTrue(np.allclose(result.history[0], expected))

    def test_new_cell_gets_a_new_persistent_id_after_division(self):
        division_grid = VoronoiGrid(grid_shape=(7, 7), cell_dist=1.0)
        division_target = division_grid.get_inner_cell_indices()[0]
        request = SimulationRequest(
            grid=division_grid,
            duration=0.2,
            params_source='dict(dT=0.1, n_var=1, sigma=0.0)',
            ode_source=f'''
def ode(Y, t, params, cell_positions):
    return np.zeros_like(Y)

def mechanics(Y, t, params, cell_positions):
    return {{'divide': [{division_target}], 'die': [], 'forces': None}}
''',
        )

        result = SimulationSession(request).run()

        initial_count = len(result.cell_ids[0])
        self.assertEqual(result.cell_ids[0].tolist(), list(range(initial_count)))
        self.assertEqual(
            result.cell_ids[1].tolist(), list(range(initial_count + 1))
        )

    def test_surviving_cells_keep_their_ids_after_apoptosis(self):
        request = SimulationRequest(
            grid=self.grid,
            duration=0.2,
            params_source='dict(dT=0.1, n_var=1, sigma=0.0)',
            ode_source='''
def ode(Y, t, params, cell_positions):
    return np.zeros_like(Y)

def mechanics(Y, t, params, cell_positions):
    return {'divide': [], 'die': [1], 'forces': None}
''',
        )

        result = SimulationSession(request).run()

        self.assertEqual(result.cell_ids[0].tolist(), [0, 1, 2, 3, 4])
        self.assertEqual(result.cell_ids[1].tolist(), [0, 2, 3, 4])

    def test_invalid_params_fail_at_the_session_interface(self):
        request = SimulationRequest(
            grid=self.grid,
            duration=0.2,
            params_source='dict(n_var=2)',
            ode_source=ODE_SOURCE,
        )
        with self.assertRaises(SimulationConfigurationError):
            SimulationSession(request)


if __name__ == '__main__':
    unittest.main()
