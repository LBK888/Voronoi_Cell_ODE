"""Exercise the real Qt controls, worker, clipboard and vector export offscreen."""

import json
import logging
import os
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

import numpy as np
import pandas as pd
from PyQt6.QtCore import QPoint, QPointF, QCoreApplication, QEvent, QMimeData, QUrl, Qt
from PyQt6.QtGui import QDragEnterEvent, QDropEvent
from PyQt6.QtWidgets import QApplication

from gui.figure_export import save_figure_dialog
from gui.main_window import MainWindow
from gui.simulation_session import SimulationRequest, SimulationSession
from voronoi_grid import VoronoiGrid


class ExperimentUiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.data_patch = patch.dict(os.environ, {'VORONOI_CELL_ODE_DATA_DIR': self.temp.name})
        self.data_patch.start()
        self.config_patch = patch('gui.main_window.CONFIG_PATH', str(Path(self.temp.name, 'config.json')))
        self.example_patch = patch('gui.main_window.CONFIG_EXAMPLE_PATH', str(Path(self.temp.name, 'no-example.json')))
        self.config_patch.start()
        self.example_patch.start()
        # Each test owns a window; avoid carrying GUI log handlers referring
        # to deleted windows into the next test.
        self.logger_patch = patch('gui.main_window.setup_logger', return_value=logging.Logger('test-window'))
        self.logger_patch.start()
        self.window = MainWindow()
        self.window.apply_config({
            'T': '.4',
            'params': 'dict(dT=.1, a=1., b=2., init_Y=[0.], labels=["Signal"])',
            'ode': 'def ode(Y, t, params, cell_positions):\n    return np.full_like(Y, params["a"] + params["b"])',
            'color_func': 'def color_func(Y, mode="polygon"):\n    return plt.cm.Reds((Y[:,0]-Y[:,0].min())/(np.ptp(Y[:,0])+1e-8)) if mode=="polygon" else None',
        })
        self.cells = np.array([[0., 0.], [2., 0.], [2., 2.], [0., 2.], [1., 1.]])
        self.window.preview_canvas.set_grid(VoronoiGrid.from_cells(self.cells))

    def tearDown(self):
        if self.window.experiments_panel.worker is not None:
            self.window.experiments_panel.cancel()
            self.window.experiments_panel.worker.wait(5000)
            self.app.processEvents()
        self.window.close()
        # Release native widgets before the shared QApplication is torn down.
        # Keeping closed Matplotlib/Qt windows alive until interpreter shutdown
        # can otherwise crash the offscreen Windows test process after "OK".
        self.window.deleteLater()
        self.app.clipboard().clear()
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        self.app.processEvents()
        self.config_patch.stop()
        self.example_patch.stop()
        self.logger_patch.stop()
        self.data_patch.stop()
        self.temp.cleanup()

    def wait_until(self, predicate):
        deadline = time.monotonic() + 10
        while not predicate() and time.monotonic() < deadline:
            self.app.processEvents()
            time.sleep(.01)
        self.assertTrue(predicate(), self.window.status_bar.currentMessage())

    def test_save_menu_and_drop_model_with_imported_grid(self):
        self.window.grid_mode_combo.setCurrentIndex(3)
        path = Path(self.temp.name, 'templates', '我的模板.json')
        with patch('gui.main_window.QFileDialog.getSaveFileName', return_value=(str(path.with_suffix('')), '')):
            self.window.save_template_dialog()
        config = json.loads(path.read_text(encoding='utf-8'))
        self.assertEqual(config['cell_positions'], self.cells.tolist())
        self.window.refresh_templates_menu()
        self.assertIn(path.stem, [action.text() for action in self.window.template_actions])

        self.window.T_edit.setText('9')
        self.window.simulation_result = object()
        self.window.download_btn.setEnabled(True)
        mime = QMimeData()
        mime.setUrls([QUrl.fromLocalFile(str(path))])
        enter = QDragEnterEvent(QPoint(1, 1), Qt.DropAction.CopyAction, mime,
                                Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier)
        drop = QDropEvent(QPointF(1, 1), Qt.DropAction.CopyAction, mime,
                          Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier)
        self.app.sendEvent(self.window.params_edit, enter)
        self.app.sendEvent(self.window.params_edit, drop)
        self.assertTrue(drop.isAccepted())
        self.assertEqual(self.window.T_edit.text(), '.4')
        np.testing.assert_allclose(self.window.preview_canvas.get_grid().copy_cells(), self.cells)
        self.assertIsNone(self.window.simulation_result)
        self.assertFalse(self.window.download_btn.isEnabled())

        broken = path.with_name('invalid.json')
        broken.write_text('{"grid_shape_x": 2}', encoding='utf-8')
        self.assertFalse(self.window.open_model(broken))
        self.assertEqual(self.window.T_edit.text(), '.4')
        self.window.run_btn.click()
        self.wait_until(lambda: self.window.simulation_result is not None)
        np.testing.assert_allclose(self.window.simulation_result.history[-1], 0.9)

    def test_agent_menu_persists_cli_preference(self):
        from automation_cli import cli_enabled
        self.assertFalse(cli_enabled())
        self.window.enable_cli_action.setChecked(True)
        self.assertTrue(cli_enabled())
        self.window.enable_cli_action.setChecked(False)
        self.assertFalse(cli_enabled())

    def test_real_tab_worker_heatmap_exports_and_settings_roundtrip(self):
        tabs = self.window.settings_tabs
        self.assertEqual([tabs.tabText(i) for i in range(tabs.count())][-2:], ['Run', 'Experiments'])
        tabs.setCurrentIndex(self.window.experiments_tab_index)
        panel = self.window.experiments_panel
        self.assertEqual(self.window.result_stack.currentWidget(), panel.result_widget)
        panel.x_axis.target.setCurrentIndex(panel.x_axis.target.findData('a'))
        panel.y_axis.target.setCurrentIndex(panel.y_axis.target.findData('b'))
        panel.x_axis.minimum.setText('1')
        panel.x_axis.maximum.setText('2')
        panel.y_axis.minimum.setText('3')
        panel.y_axis.maximum.setText('4')
        panel.x_axis.count.setValue(2)
        panel.y_axis.count.setValue(2)
        panel.run_btn.click()
        self.assertIsNotNone(panel.worker)
        self.assertFalse(self.window.run_btn.isEnabled())
        self.wait_until(lambda: panel.worker is None)
        self.assertIsNotNone(panel.result)
        np.testing.assert_allclose(panel.result.values, [[1.2, 1.5], [1.5, 1.8]])
        self.assertTrue(panel.pdf_btn.isEnabled())
        self.assertEqual(len(panel.figure.axes), 2)
        csv = Path(self.temp.name, 'sweep.csv')
        panel.export_data(csv)
        table = pd.read_csv(csv)
        self.assertEqual(len(table), 4)
        np.testing.assert_allclose(table.pixel_mean, [1.2, 1.5, 1.5, 1.8])
        metadata = json.loads(csv.with_name('sweep_settings.json').read_text(encoding='utf-8'))
        self.assertEqual(metadata['x']['target'], 'a')
        self.assertEqual(metadata['initial_cell_positions'], self.cells.tolist())
        saved = self.window.save_config()
        panel.x_axis.minimum.setText('999')
        self.window.apply_config(saved)
        panel.refresh_targets()
        self.assertEqual(panel.x_axis.minimum.text(), '1')
        self.assertEqual(panel.x_axis.target.currentData(), 'a')

    def test_saved_targets_survive_loading_before_opening_experiments(self):
        panel = self.window.experiments_panel
        panel.refresh_targets()
        saved = panel.settings()
        saved['x']['target'] = 'b'
        saved['y']['target'] = 'a'
        panel.restore_settings(saved)
        # An ordinary save before the first refresh must retain the loaded targets.
        self.assertEqual(panel.settings()['x']['target'], 'b')
        panel.refresh_targets()
        self.assertEqual(panel.x_axis.target.currentData(), 'b')
        self.assertEqual(panel.y_axis.target.currentData(), 'a')

    def test_result_right_click_copies_image_without_opening_cell_trace(self):
        request = SimulationRequest(VoronoiGrid.from_cells(self.cells), .4,
                                    self.window.params_edit.toPlainText(), self.window.ode_edit.toPlainText())
        result = SimulationSession(request).run()
        self.window.on_simulation_finished(result)
        self.window.update_anim_frame(3)
        with patch.object(self.window, 'show_cell_plot') as trace:
            self.window.on_anim_click(SimpleNamespace(button=3))
            trace.assert_not_called()
        with patch('gui.figure_export.QMenu') as menu:
            copy_action, save_action = object(), object()
            menu.return_value.addAction.side_effect = [copy_action, save_action]
            menu.return_value.exec.return_value = copy_action
            self.window.anim_canvas.customContextMenuRequested.emit(QPoint(10, 10))
        image = self.app.clipboard().image()
        self.assertFalse(image.isNull())
        self.assertGreater(image.width(), 100)
        self.assertIn('copied', self.window.status_bar.currentMessage())
        self.assertEqual(self.window.figure.axes[0].get_title(), 'Cell state animation, frame 3')
        self.assertAlmostEqual(self.window.figure.axes[1]._colorbar.norm.vmin, .9)

    def test_save_frame_defaults_to_editable_pdf_and_adds_missing_extension(self):
        self.window.figure.add_subplot(111).plot([0, 1], [0, 1])
        target = str(Path(self.temp.name, 'current'))
        with patch('gui.figure_export.QFileDialog.getSaveFileName', return_value=(target, 'PDF Files (*.pdf)')) as dialog:
            self.window.on_save_frame()
        self.assertEqual(dialog.call_args.args[2], 'frame.pdf')
        self.assertTrue(dialog.call_args.args[3].startswith('PDF Files'))
        self.assertTrue(Path(target + '.pdf').read_bytes().startswith(b'%PDF-'))
        with patch('gui.figure_export.QFileDialog.getSaveFileName', return_value=(target + '_vector', 'SVG Files (*.svg)')):
            path = save_figure_dialog(self.window, self.window.figure)
        self.assertTrue(Path(path).read_text(encoding='utf-8').find('<svg') >= 0)


if __name__ == '__main__':
    unittest.main()
