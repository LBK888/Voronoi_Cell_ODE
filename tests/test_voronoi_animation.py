import unittest
from contextlib import contextmanager
from unittest.mock import patch

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

from voronoi_animation import VoronoiAnimator
from voronoi_grid import VoronoiGrid


class VoronoiAnimatorColorRangeTests(unittest.TestCase):
    def setUp(self):
        self.cells = np.array([
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, 1.0],
        ])
        self.history = [
            np.array([[1.0], [2.0], [3.0], [2.5], [1.5]]),
            np.array([[100.0], [200.0], [300.0], [250.0], [150.0]]),
        ]
        self.received_bounds = []

        def color_func(values, mode='polygon', global_bounds=None):
            self.received_bounds.append(global_bounds)
            if mode != 'polygon':
                return None
            vmin, vmax = global_bounds[0]
            norm = (values[:, 0] - vmin) / (vmax - vmin)
            return plt.cm.Reds(norm)

        self.animator = VoronoiAnimator(
            VoronoiGrid.from_cells(self.cells),
            self.history,
            color_func=color_func,
            cell_positions_history=[self.cells.copy(), self.cells.copy()],
            params={'labels': ['signal'], 'n_var': 1},
            color_func_str='return plt.cm.Reds(values)',
        )

    def test_draw_frame_passes_only_that_frames_bounds_to_color_function(self):
        figure = Figure()
        FigureCanvasAgg(figure)
        axis = figure.add_subplot(111)

        self.animator._draw_frame(axis, 0)
        first_frame_bounds = list(self.received_bounds)
        self.received_bounds.clear()
        axis.clear()
        self.animator._draw_frame(axis, 1)

        self.assertTrue(first_frame_bounds)
        self.assertTrue(
            all(bounds == [(1.0, 3.0)] for bounds in first_frame_bounds)
        )
        self.assertTrue(self.received_bounds)
        self.assertTrue(
            all(bounds == [(100.0, 300.0)] for bounds in self.received_bounds)
        )

    def test_gui_colorbar_tracks_the_current_frame_bounds(self):
        figure = Figure()
        FigureCanvasAgg(figure)
        figure.add_subplot(111)
        animation = self.animator.animateGUI(figure)
        animation._draw_was_started = True
        colorbar = figure.axes[1]._colorbar

        animation._func(0)
        self.assertEqual((colorbar.norm.vmin, colorbar.norm.vmax), (1.0, 3.0))

        animation._func(1)
        self.assertEqual(
            (colorbar.norm.vmin, colorbar.norm.vmax),
            (100.0, 300.0),
        )

    def test_exported_video_has_labeled_colorbar_in_every_frame(self):
        captured = []

        class RecordingWriter:
            def __init__(self, **kwargs):
                pass

            @contextmanager
            def saving(self, figure, path, dpi):
                self.figure = figure
                yield

            def grab_frame(self):
                self.figure.canvas.draw()
                captured.append([
                    (axis.get_ylabel(), axis._colorbar.norm.vmin,
                     axis._colorbar.norm.vmax)
                    for axis in self.figure.axes if hasattr(axis, '_colorbar')
                ])

        with patch('matplotlib.animation.FFMpegWriter', RecordingWriter):
            self.animator.export_video('unused.mp4')

        self.assertEqual(captured, [
            [('signal', 1.0, 3.0)],
            [('signal', 100.0, 300.0)],
        ])

    def test_explicit_legends_support_reordered_columns_and_reused_colormaps(self):
        animator = VoronoiAnimator(
            self.animator.vor_grid,
            [np.column_stack((frame[:, 0], frame[:, 0] * 10)) for frame in self.history],
            color_func=lambda values: plt.cm.Reds(values[:, 0] / values[:, 0].max()),
            params={'labels': ['A', 'B'], 'colorbars': [
                dict(variable=1, cmap='Reds'), dict(variable=0, cmap='Reds'),
            ]},
        )
        figure = Figure()
        FigureCanvasAgg(figure)
        figure.add_subplot(111)
        animation = animator.animateGUI(figure)
        animation._draw_was_started = True
        animation._func(1)
        self.assertEqual([axis.get_ylabel() for axis in figure.axes[1:]], ['B', 'A'])
        self.assertEqual([(axis._colorbar.norm.vmin, axis._colorbar.norm.vmax)
                          for axis in figure.axes[1:]], [(1000., 3000.), (100., 300.)])


if __name__ == '__main__':
    unittest.main()
