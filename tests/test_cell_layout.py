import unittest

import numpy as np

from voronoi_grid import MovementPlan, VoronoiGrid


class VoronoiGridInterfaceTests(unittest.TestCase):
    def setUp(self):
        self.cells = np.array([
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, 1.0],
        ])

    def test_mutations_keep_voronoi_state_synchronized(self):
        grid = VoronoiGrid.from_cells(self.cells)
        added_index = grid.add_cell([3.0, 1.0])
        self.assertEqual(added_index, 5)
        grid.move_cell(added_index, [3.0, 1.5])
        self.assertTrue(np.allclose(grid.cells[added_index], [3.0, 1.5]))
        grid.remove_cell(added_index)
        self.assertEqual(grid.cells.shape, (5, 2))
        self.assertEqual(len(grid.vor.points), len(grid.cells))

    def test_neighbor_query_is_symmetric(self):
        grid = VoronoiGrid.from_cells(self.cells)
        neighbors = grid.finite_ridge_neighbors(max_ridge_length=None)
        for source, targets in neighbors.items():
            for target, length in targets.items():
                self.assertAlmostEqual(length, neighbors[target][source])

    def test_movement_plan_updates_the_layout_through_one_interface(self):
        grid = VoronoiGrid.from_cells(self.cells)
        before = grid.copy_cells()
        grid.apply_movement(MovementPlan(away_strength=0.1))
        self.assertFalse(np.allclose(before, grid.cells))
        self.assertEqual(len(grid.vor.points), len(grid.cells))


if __name__ == '__main__':
    unittest.main()
