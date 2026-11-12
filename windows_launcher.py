"""Windowed executable entry point with a real frozen-application smoke test."""

import argparse
import json
import multiprocessing
import os
from pathlib import Path
import sys
import traceback


def smoke_test(app, window, output, model):
    """Exercise bundled libraries, external templates, UI, simulation and exports."""
    import numpy as np
    import pandas as pd
    import h5py
    from matplotlib.animation import FFMpegWriter
    from gui.model_files import read_model, template_files, user_templates
    from gui.experiments import ExperimentRequest, ExperimentRunner, SweepAxis
    from gui.node_editor import NodeEditorDialog, NodeEditorDocument
    from gui.code_editor_window import CodeEditorDialog
    from gui.simulation_session import SimulationRequest, SimulationSession
    from voronoi_animation import VoronoiAnimator
    from voronoi_grid import VoronoiGrid

    config = read_model(model)
    assert window.open_model(model)
    window.grid_shape_x.setValue(5)
    window.grid_shape_y.setValue(5)
    window.T_edit.setText('0.2')
    window.on_preview()
    # Exercise the Numba-compiled spatial helper, not just its import.
    from gui.sim_utils import diffusion_weighted_mean
    cells = window.preview_canvas.get_grid().copy_cells()
    assert np.isfinite(diffusion_weighted_mean(np.ones(len(cells)), cells, 0, 1.))
    request = SimulationRequest(VoronoiGrid.from_cells(cells), .2, config['params'], config['ode'])
    result = SimulationSession(request).run()
    assert all(np.isfinite(frame).all() for frame in result.history)
    window.on_simulation_finished(result)
    window.update_anim_frame(len(result.history) - 1)
    window.figure.savefig(output / 'frame.pdf')
    from gui.figure_export import copy_figure_to_clipboard
    copy_figure_to_clipboard(window.figure)
    assert not app.clipboard().image().isNull()
    source = config['color_func']
    namespace = {'np': np, 'plt': __import__('matplotlib.pyplot', fromlist=['plt'])}
    exec(source, namespace)
    assert FFMpegWriter.isAvailable()
    animator = VoronoiAnimator(result.grid, result.history, namespace['color_func'],
                               result.cell_positions, params=result.params, color_func_str=source)
    animator.export_video(str(output / 'simulation.mp4'), fps=10, dpi=80)
    sweep = ExperimentRunner(ExperimentRequest(
        request, SweepAxis('parameter', 'alpha', 1., 2., 2),
        SweepAxis('initial', 1, .1, .5, 2), observable=0,
    )).run()
    assert np.isfinite(sweep.values).all()
    panel = window.experiments_panel
    panel._finished(sweep)
    panel.figure.savefig(output / 'heatmap.pdf')
    panel.export_data(output / 'heatmap.csv')
    pd.DataFrame(sweep.values).to_excel(output / 'data.xlsx', index=False)
    with h5py.File(output / 'data.h5', 'w') as handle:
        handle.create_dataset('heatmap', data=sweep.values)
    editor = CodeEditorDialog(window, target_edit=window.ode_edit, editor_type='ode')
    assert editor.editor is not None
    editor.deleteLater()
    node_editor = NodeEditorDialog(window, NodeEditorDocument(
        ode_source=config['ode'], params_source=config['params'], color_func_source=config['color_func'],
    ))
    node_editor.deleteLater()
    window.settings_tabs.setCurrentIndex(window.experiments_tab_index)
    app.processEvents()
    window.grab().save(str(output / 'window.png'))
    window.save_config()
    directory = user_templates()
    directory.mkdir(parents=True, exist_ok=True)
    saved = directory / 'Smoke Template.json'
    window.save_template(saved)
    assert window.open_model(saved)
    window.enable_cli_action.setChecked(True)
    report = {
        'status': 'passed', 'frozen': bool(getattr(sys, 'frozen', False)),
        'frames': len(result.history), 'heatmap_shape': list(sweep.values.shape),
        'templates': [str(path) for path in template_files()],
        'saved_template': str(saved),
        'ffmpeg': FFMpegWriter.bin_path(),
        'checks': ['Qt UI', 'Numba spatial helper', 'ODE', 'PDF', 'clipboard',
                   'MP4', 'Experiments', 'CSV', 'Excel', 'HDF5', 'code editor', 'node editor'],
    }
    (output / 'smoke-report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')


def main():
    multiprocessing.freeze_support()
    if sys.argv[1:2] == ['--cli']:
        from automation_cli import main as run_cli
        return run_cli(sys.argv[2:])
    parser = argparse.ArgumentParser(description='Voronoi Cell ODE')
    parser.add_argument('--smoke-test', type=Path, metavar='OUTPUT_FOLDER')
    parser.add_argument('--smoke-model', type=Path, metavar='MODEL_JSON')
    args = parser.parse_args()
    if args.smoke_test and not args.smoke_model:
        if getattr(sys, 'frozen', False):
            parser.error('--smoke-test requires --smoke-model in the EXE.')
        args.smoke_model = Path(__file__).parent / 'skills/paper-to-voronoi-ode/assets/hill-demo.json'
    if args.smoke_test:
        args.smoke_test = args.smoke_test.resolve()
        args.smoke_test.mkdir(parents=True, exist_ok=True)
        os.environ['VORONOI_CELL_ODE_DATA_DIR'] = str(args.smoke_test / 'settings')
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    try:
        from PyQt6.QtWidgets import QApplication
        from gui.main_window import MainWindow
        app = QApplication(sys.argv[:1])
        window = MainWindow()
        if args.smoke_test:
            smoke_test(app, window, args.smoke_test, args.smoke_model)
            window.close()
            return 0
        window.show()
        return app.exec()
    except Exception:
        details = traceback.format_exc()
        if args.smoke_test:
            (args.smoke_test / 'smoke-error.txt').write_text(details, encoding='utf-8')
        else:
            from gui.app_paths import data_path
            data_path('startup-error.txt').write_text(details, encoding='utf-8')
            from PyQt6.QtWidgets import QApplication, QMessageBox
            app = QApplication.instance() or QApplication(sys.argv[:1])
            QMessageBox.critical(None, 'Voronoi Cell ODE could not start', details)
        return 1


if __name__ == '__main__':
    sys.exit(main())
