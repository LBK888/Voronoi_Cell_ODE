import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QGroupBox, QFormLayout, QLineEdit, QTextEdit, QComboBox,
    QCheckBox, QFileDialog, QProgressBar, QStatusBar, QSpinBox, QSlider,
    QDoubleSpinBox, QScrollArea, QMenuBar, QMenu
)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QFont, QColor, QPalette
from gui.logger import setup_logger
from gui.preview_canvas import PreviewCanvas
import numpy as np
import ast
import os
import logging
import matplotlib
matplotlib.use('QtAgg')  # 使用QtAgg backend，與PyQt6兼容
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import json

from voronoi_grid import VoronoiGrid
from biophysics_model import BiophysicsModel
from voronoi_animation import VoronoiAnimator
from gui.simulation_worker import SimulationWorker

CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config.json')
CONFIG_EXAMPLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config_example.json')

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.logger = setup_logger()
        self.setWindowTitle("Voronoi Notch Simulation")
        self.resize(1280, 880)
        self.init_ui()
        self.load_config()
        self.sim_history = None
        self.cell_positions_history = None
        self.grid = None
        self.params = None
        self.worker = None
        self.logger.addHandler(GUIStatusHandler(self.status_bar))

    def init_ui(self):
        menubar = self.menuBar()
        template_menu = menubar.addMenu('Templates')
        
        template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'templates')
        if os.path.exists(template_dir) and os.path.isdir(template_dir):
            for filename in sorted(os.listdir(template_dir)):
                if filename.endswith('.json'):
                    name = filename[:-5]
                    action = template_menu.addAction(name)
                    action.triggered.connect(lambda checked, f=filename: self.load_template(f))

        main_widget = QWidget()
        main_layout = QHBoxLayout(main_widget)

        # Left: Settings
        settings_widget = QWidget()
        settings_layout = QVBoxLayout(settings_widget)

        # Cell grid generation mode
        grid_mode_box = QGroupBox("Cell Grid Generation Mode")
        grid_mode_layout = QVBoxLayout()
        self.grid_mode_combo = QComboBox()
        self.grid_mode_combo.addItems([
            "Honeycomb",
            "Random",
            "Regular Grid",
            "User define by import"
        ])
        grid_mode_layout.addWidget(self.grid_mode_combo)
        grid_mode_box.setLayout(grid_mode_layout)
        settings_layout.addWidget(grid_mode_box)

        # Cell grid parameters
        grid_param_box = QGroupBox("Cell Grid Parameters")
        grid_param_form = QFormLayout()
        grid_shape_row = QHBoxLayout()
        self.grid_shape_x = QSpinBox()
        self.grid_shape_x.setRange(1, 200)
        self.grid_shape_x.setValue(20)
        self.grid_shape_y = QSpinBox()
        self.grid_shape_y.setRange(1, 200)
        self.grid_shape_y.setValue(20)
        grid_shape_row.addWidget(QLabel("X:"))
        grid_shape_row.addWidget(self.grid_shape_x)
        grid_shape_row.addWidget(QLabel("Y:"))
        grid_shape_row.addWidget(self.grid_shape_y)
        grid_param_form.addRow("Grid Shape:", grid_shape_row)
        self.cell_dist_edit = QLineEdit("1.0")
        self.pos_rand_edit = QLineEdit("0.25")
        grid_param_form.addRow("Cell Distance:", self.cell_dist_edit)
        grid_param_form.addRow("Position Randomness:", self.pos_rand_edit)
        grid_param_box.setLayout(grid_param_form)
        settings_layout.addWidget(grid_param_box)

        # Cell movement modes
        move_mode_box = QGroupBox("Cell Movement Modes")
        move_mode_layout = QVBoxLayout()
        self.move_random = QCheckBox("Basic Random Movement str:")
        self.move_random_strength = QLineEdit("0.05")
        move_random_row = QHBoxLayout()
        move_random_row.addWidget(self.move_random)
        move_random_row.addWidget(self.move_random_strength)
        move_mode_layout.addLayout(move_random_row)
        self.move_away = QCheckBox("Move Away From Center str:")
        self.move_away_strength = QLineEdit("0.1")
        move_away_row = QHBoxLayout()
        move_away_row.addWidget(self.move_away)
        move_away_row.addWidget(self.move_away_strength)
        move_mode_layout.addLayout(move_away_row)
        self.move_ce = QCheckBox("Convergent Extension str:")
        self.move_ce_strength = QLineEdit("0.02")
        move_ce_row = QHBoxLayout()
        move_ce_row.addWidget(self.move_ce)
        move_ce_row.addWidget(self.move_ce_strength)
        move_mode_layout.addLayout(move_ce_row)
        self.move_repulsion = QCheckBox("Repulsion (Neighbors, No Outer) str:")
        self.move_repulsion_strength = QLineEdit("0.08")
        move_repulsion_row = QHBoxLayout()
        move_repulsion_row.addWidget(self.move_repulsion)
        move_repulsion_row.addWidget(self.move_repulsion_strength)
        move_mode_layout.addLayout(move_repulsion_row)
        self.move_vertex_model = QCheckBox("Vertex Model (Adhesion) str:")
        self.move_vertex_model_strength = QLineEdit("0.05")
        move_vertex_model_row = QHBoxLayout()
        move_vertex_model_row.addWidget(self.move_vertex_model)
        move_vertex_model_row.addWidget(self.move_vertex_model_strength)
        move_mode_layout.addLayout(move_vertex_model_row)
        self.move_division = QCheckBox("Cell Division n:")
        self.move_division_n = QSpinBox()
        self.move_division_n.setRange(1, 10)
        self.move_division_n.setValue(1)
        self.move_division_method = QComboBox()
        self.move_division_method.addItems(["area", "random"])
        move_division_row = QHBoxLayout()
        move_division_row.addWidget(self.move_division)
        move_division_row.addWidget(self.move_division_n)
        move_division_row.addWidget(self.move_division_method)
        move_mode_layout.addLayout(move_division_row)
        self.move_apoptosis = QCheckBox("Apoptosis n:")
        self.move_apoptosis_n = QSpinBox()
        self.move_apoptosis_n.setRange(1, 10)
        self.move_apoptosis_n.setValue(1)
        self.move_apoptosis_method = QComboBox()
        self.move_apoptosis_method.addItems(["area", "random"])
        move_apoptosis_row = QHBoxLayout()
        move_apoptosis_row.addWidget(self.move_apoptosis)
        move_apoptosis_row.addWidget(self.move_apoptosis_n)
        move_apoptosis_row.addWidget(self.move_apoptosis_method)
        move_mode_layout.addLayout(move_apoptosis_row)
        move_mode_box.setLayout(move_mode_layout)
        settings_layout.addWidget(move_mode_box)

        # ODE editor
        ode_box = QGroupBox("ODE Function")
        ode_layout = QVBoxLayout()
        
        self.open_node_editor_btn = QPushButton("Open Node Editor Prototype")
        ode_layout.addWidget(self.open_node_editor_btn)
        
        self.ode_edit = QTextEdit()
        self.ode_edit.setMinimumHeight(50)
        font = QFont("Consolas", 10)
        self.ode_edit.setFont(font)
        palette = self.ode_edit.palette()
        palette.setColor(QPalette.ColorRole.Base, QColor(30,30,30))
        palette.setColor(QPalette.ColorRole.Text, QColor(220,220,220))
        self.ode_edit.setPalette(palette)
        self.ode_edit.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.ode_edit.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        ode_layout.addWidget(self.ode_edit)
        
        self.open_ode_editor_btn = QPushButton("Open Code Editor")
        self.open_ode_editor_btn.clicked.connect(lambda: self.open_code_editor(self.ode_edit, 'ode'))
        ode_layout.addWidget(self.open_ode_editor_btn)

        ode_box.setLayout(ode_layout)
        settings_layout.addWidget(ode_box)

        # Color function editor
        color_func_box = QGroupBox("Color Function (color_func)")
        color_func_layout = QVBoxLayout()
        self.color_func_edit = QTextEdit()
        self.color_func_edit.setMinimumHeight(50)
        self.color_func_edit.setFont(font)
        self.color_func_edit.setPalette(palette)
        self.color_func_edit.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.color_func_edit.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        color_func_layout.addWidget(self.color_func_edit)
        
        self.open_color_func_editor_btn = QPushButton("Open Code Editor")
        self.open_color_func_editor_btn.clicked.connect(lambda: self.open_code_editor(self.color_func_edit, 'color_func'))
        color_func_layout.addWidget(self.open_color_func_editor_btn)

        color_func_box.setLayout(color_func_layout)
        settings_layout.addWidget(color_func_box)

        # Params editor
        params_box = QGroupBox("Params (dict)")
        params_layout = QVBoxLayout()
        self.params_edit = QTextEdit()
        self.params_edit.setMinimumHeight(50)
        self.params_edit.setFont(font)
        self.params_edit.setPalette(palette)
        self.params_edit.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.params_edit.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        params_layout.addWidget(self.params_edit)
        
        self.open_params_editor_btn = QPushButton("Open Code Editor")
        self.open_params_editor_btn.clicked.connect(lambda: self.open_code_editor(self.params_edit, 'params'))
        params_layout.addWidget(self.open_params_editor_btn)

        self.generate_params_btn = QPushButton("Generate Visual Controls")
        params_layout.addWidget(self.generate_params_btn)

        self.dynamic_params_scroll = QScrollArea()
        self.dynamic_params_scroll.setWidgetResizable(True)
        self.dynamic_params_widget = QWidget()
        self.dynamic_params_layout = QFormLayout()
        self.dynamic_params_widget.setLayout(self.dynamic_params_layout)
        self.dynamic_params_scroll.setWidget(self.dynamic_params_widget)
        self.dynamic_params_scroll.setMinimumHeight(150)
        params_layout.addWidget(self.dynamic_params_scroll)

        params_box.setLayout(params_layout)
        settings_layout.addWidget(params_box)

        # Simulation parameters
        sim_param_box = QGroupBox("Simulation Parameters")
        sim_param_form = QFormLayout()
        self.T_edit = QLineEdit("30.0")
        self.replicate_edit = QLineEdit("5")
        self.repeats_edit = QLineEdit("5")
        sim_param_form.addRow("T:", self.T_edit)
        sim_param_form.addRow("Replicate:", self.replicate_edit)
        sim_param_form.addRow("Repeats:", self.repeats_edit)
        sim_param_box.setLayout(sim_param_form)
        settings_layout.addWidget(sim_param_box)

        # Preview button
        self.preview_btn = QPushButton("Preview")
        settings_layout.addWidget(self.preview_btn)
        # Apply movement button
        self.apply_move_btn = QPushButton("Apply Movement")
        settings_layout.addWidget(self.apply_move_btn)

        settings_widget.setLayout(settings_layout)
        main_layout.addWidget(settings_widget, 2)

        # Right: Preview
        preview_widget = QWidget()
        preview_layout = QVBoxLayout(preview_widget)
        self.preview_canvas = PreviewCanvas()
        self.preview_canvas.connect_events()
        preview_layout.addWidget(self.preview_canvas)
        # 連接 cell_count_changed signal
        self.preview_canvas.cell_count_changed.connect(lambda count: self.status_bar.showMessage(f"New cell number: {count}"))
        # 新增動畫嵌入區
        self.figure = Figure()
        self.anim_canvas = FigureCanvas(self.figure)
        self.anim_canvas.setVisible(False)
        preview_layout.addWidget(self.anim_canvas)
        # --- 新增動畫控制元件 ---
        anim_ctrl_layout = QHBoxLayout()
        self.play_btn = QPushButton("Play")
        self.play_btn.setEnabled(False)
        self.frame_slider = QSlider(Qt.Orientation.Horizontal)
        self.frame_slider.setEnabled(False)
        self.save_frame_btn = QPushButton("Download Pic")
        self.save_frame_btn.setEnabled(False)
        anim_ctrl_layout.addWidget(self.play_btn)
        anim_ctrl_layout.addWidget(self.frame_slider)
        anim_ctrl_layout.addWidget(self.save_frame_btn)
        preview_layout.addLayout(anim_ctrl_layout)
        # ---
        self.download_btn = QPushButton("Download Results")
        preview_layout.addWidget(self.download_btn)
        self.run_btn = QPushButton("Run Simulation")
        preview_layout.addWidget(self.run_btn)
        preview_widget.setLayout(preview_layout)
        main_layout.addWidget(preview_widget, 5)

        # Status bar and progress bar
        self.status_bar = QStatusBar()
        self.progress_bar = QProgressBar()
        self.status_bar.addPermanentWidget(self.progress_bar)
        self.setStatusBar(self.status_bar)

        self.setCentralWidget(main_widget)

        self.logger.info("GUI initialized.")

        self.preview_btn.clicked.connect(self.on_preview)
        self.apply_move_btn.clicked.connect(self.on_apply_movement)
        self.run_btn.clicked.connect(self.on_run_simulation)
        self.download_btn.clicked.connect(self.on_download_results)
        self.generate_params_btn.clicked.connect(self.on_generate_params)
        self.open_node_editor_btn.clicked.connect(self.open_node_editor)
        # 新增動畫控制事件
        self.play_btn.clicked.connect(self.on_play_pause)
        self.frame_slider.valueChanged.connect(self.on_frame_slider_changed)
        self.save_frame_btn.clicked.connect(self.on_save_frame)
        # QTimer for animation playback
        self.anim_timer = QTimer(self)
        self.anim_timer.timeout.connect(self._on_anim_timer_tick)
        self.anim_timer.setInterval(50)  # 20 FPS, can be adjusted
        # 設定變動時自動儲存
        for widget in [self.grid_shape_x, self.grid_shape_y, self.cell_dist_edit, self.pos_rand_edit,
                       self.move_random, self.move_random_strength, self.move_away, self.move_away_strength,
                       self.move_ce, self.move_ce_strength, self.move_repulsion, self.move_repulsion_strength,
                       self.move_vertex_model, self.move_vertex_model_strength,
                       self.move_division, self.move_division_n, self.move_division_method,
                       self.move_apoptosis, self.move_apoptosis_n, self.move_apoptosis_method,
                       self.ode_edit, self.params_edit, self.color_func_edit, self.T_edit, self.replicate_edit, self.repeats_edit]:
            if hasattr(widget, 'editingFinished'):
                widget.editingFinished.connect(self.save_config)
            elif hasattr(widget, 'valueChanged'):
                widget.valueChanged.connect(self.save_config)
            elif hasattr(widget, 'currentIndexChanged'):
                widget.currentIndexChanged.connect(self.save_config)
            elif hasattr(widget, 'stateChanged'):
                widget.stateChanged.connect(self.save_config)

        # Real-time IDE-like Syntax Checking
        self.ode_edit.textChanged.connect(lambda: self.validate_code(self.ode_edit, 'ode'))
        self.color_func_edit.textChanged.connect(lambda: self.validate_code(self.color_func_edit, 'color_func'))

    def validate_code(self, edit_widget, name):
        code = edit_widget.toPlainText().strip()
        if not code: return
        try:
            ast.parse(code)
            edit_widget.setStyleSheet("")
            self.status_bar.showMessage(f"{name} syntax OK", 2000)
        except SyntaxError as e:
            edit_widget.setStyleSheet("background-color: #4a1c1c;")
            self.status_bar.showMessage(f"{name} Syntax Error (Line {e.lineno}): {e.msg}")

    def on_generate_params(self):
        import ast
        try:
            params_text = self.params_edit.toPlainText().strip()
            if params_text.startswith('params ='):
                params_text = params_text[len('params ='):].strip()
            if params_text.startswith('dict'):
                params_dict = eval(params_text, {}, {})
            else:
                params_dict = ast.literal_eval(params_text)
            self.params = params_dict
        except Exception as e:
            self.status_bar.showMessage(f"Parse error: {e}")
            return

        while self.dynamic_params_layout.count():
            item = self.dynamic_params_layout.takeAt(0)
            if item.widget() is not None:
                item.widget().deleteLater()
            elif item.layout() is not None:
                while item.layout().count():
                    child = item.layout().takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()
                item.layout().deleteLater()

        self.dynamic_widgets = {}

        def update_dict_text():
            import pprint
            new_dict = self.params.copy()
            for k, w in self.dynamic_widgets.items():
                new_dict[k] = w.value()
            self.params = new_dict
            formatted = pprint.pformat(new_dict, width=60)
            self.params_edit.blockSignals(True)
            self.params_edit.setPlainText("params = " + formatted)
            self.params_edit.blockSignals(False)
            self.save_config()

        for k, v in self.params.items():
            if isinstance(v, float):
                w = QDoubleSpinBox()
                w.setDecimals(4)
                w.setRange(-float('inf'), float('inf'))
                w.setValue(v)
                w.setSingleStep(max(abs(v)/10.0, 0.001))
                w.valueChanged.connect(lambda val, key=k: update_dict_text())
                self.dynamic_params_layout.addRow(QLabel(k), w)
                self.dynamic_widgets[k] = w
            elif isinstance(v, int) and not isinstance(v, bool):
                w = QSpinBox()
                w.setRange(-9999999, 9999999)
                w.setValue(v)
                w.valueChanged.connect(lambda val, key=k: update_dict_text())
                self.dynamic_params_layout.addRow(QLabel(k), w)
                self.dynamic_widgets[k] = w
        
        self.status_bar.showMessage(f"Generated controls for {len(self.dynamic_widgets)} parameters.")

    def get_grid_shape(self):
        return (self.grid_shape_x.value(), self.grid_shape_y.value())

    def load_template(self, filename):
        template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'templates')
        file_path = os.path.join(template_dir, filename)
        if not os.path.exists(file_path):
            return
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                config = json.load(f)
            if config:
                self.apply_config(config)
                self.on_generate_params()
                self.status_bar.showMessage(f"Loaded template: {filename[:-5]}")
        except Exception as e:
            self.logger.error(f"Failed to load template {filename}: {e}")
            self.status_bar.showMessage(f"Failed to load template {filename}: {e}")

    def open_node_editor(self):
        from gui.node_editor import NodeEditorDialog
        dlg = NodeEditorDialog(self)
        dlg.exec()

    def open_code_editor(self, target_edit, editor_type):
        from gui.code_editor_window import CodeEditorDialog
        dlg = CodeEditorDialog(self, target_edit=target_edit, editor_type=editor_type)
        dlg.exec()

    def get_default_ode(self):
        from .default_ode import get_default_ode
        return get_default_ode()

    def get_default_params(self):
        from .default_params import get_default_params
        return get_default_params()

    def get_default_color_func(self):
        from .default_color_func import get_default_color_func
        return get_default_color_func()

    def on_preview(self):
        # 取得設定
        grid_shape = self.get_grid_shape()
        cell_dist = float(self.cell_dist_edit.text())
        pos_rand = float(self.pos_rand_edit.text())
        mode_map = {
            'Honeycomb': 'honeycomb',
            'Random': 'random',
            'Regular Grid': 'regular',
            'User define by import': 'import',
        }
        mode = mode_map[self.grid_mode_combo.currentText()]
        from voronoi_grid import VoronoiGrid
        try:
            if mode == 'import':
                file_path, _ = QFileDialog.getOpenFileName(self, "Import Cell Coordinates", "", "Table Files (*.xlsx *.xls *.csv);;Text Files (*.txt);;All Files (*)")
                if not file_path:
                    self.status_bar.showMessage("Import cancelled.")
                    return
                grid = VoronoiGrid(grid_shape=grid_shape, cell_dist=cell_dist, pos_rand=pos_rand, mode=mode, import_path=file_path)
            else:
                grid = VoronoiGrid(grid_shape=grid_shape, cell_dist=cell_dist, pos_rand=pos_rand, mode=mode)
            cell_count = len(grid.cells)
            self.logger.info(f"Generated grid with mode={mode}, shape={grid_shape}, cell_count={cell_count}")
            self.plot_voronoi(grid)
            self.status_bar.showMessage(f"Preview generated: {mode}, cell number: {cell_count}")
        except Exception as e:
            self.logger.error(f"Grid generation failed: {e}")
            self.status_bar.showMessage(f"Error: {e}")
        self.save_config()

    def plot_voronoi(self, grid):
        def plot_func(ax):
            from scipy.spatial import voronoi_plot_2d
            voronoi_plot_2d(grid.vor, ax=ax, show_vertices=False, line_colors='black')
            ax.plot(grid.cells[:,0], grid.cells[:,1], 'o', color='red')
            ax.set_aspect('equal')
            ax.set_title('Voronoi Preview')
        self.preview_canvas.set_cells(grid.cells)
        self.preview_canvas.plot(plot_func)

    def on_apply_movement(self):
        # 取得目前細胞座標
        cells = self.preview_canvas.get_cells()
        # 根據勾選的運動模式依序套用
        from gui.sim_utils import move_away_from_center, covergent_extension, repulsion_move_neighbors_no_outer, vertex_model_forces
        if self.move_random.isChecked():
            random_strength = float(self.move_random_strength.text())
            cells = cells + (np.random.rand(*cells.shape)-0.5)*2*random_strength
        if self.move_away.isChecked():
            away_strength = float(self.move_away_strength.text())
            cells = move_away_from_center(cells, away_strength)
        if self.move_ce.isChecked():
            ce_strength = float(self.move_ce_strength.text())
            cells = covergent_extension(cells, ce_strength)
        if self.move_repulsion.isChecked():
            repulsion_strength = float(self.move_repulsion_strength.text())
            cells = repulsion_move_neighbors_no_outer(cells, repulsion_strength)
        if getattr(self, 'move_vertex_model', None) and self.move_vertex_model.isChecked():
            vertex_strength = float(self.move_vertex_model_strength.text())
            cells = vertex_model_forces(cells, adhesion_strength=vertex_strength)
        if self.move_division.isChecked() or self.move_apoptosis.isChecked():
            from voronoi_grid import VoronoiGrid
            grid = VoronoiGrid(grid_shape=(len(cells), 1), mode='custom', custom_cells=cells)
            if self.move_division.isChecked():
                n = int(self.move_division_n.text())
                mode = self.move_division_method.currentText()
                grid.cell_proliferation(n=n, mode=mode)
            if self.move_apoptosis.isChecked():
                n = int(self.move_apoptosis_n.text())
                mode = self.move_apoptosis_method.currentText()
                grid.cell_apoptosis(n=n, mode=mode)
            cells = grid.cells

        self.preview_canvas.set_cells(cells)
        # 顯示最新細胞數
        cell_count = len(cells) if hasattr(cells, '__len__') else 0
        self.status_bar.showMessage(f"Movement applied. Cell number: {cell_count}")

    def on_run_simulation(self):
        if self.worker is not None and self.worker.isRunning():
            self.worker.cancel()
            self.run_btn.setText("Run Simulation")
            self.status_bar.showMessage("Simulation cancelled.")
            return

        from voronoi_grid import VoronoiGrid
        from biophysics_model import BiophysicsModel
        import numpy as np
        import ast
        # 每次 simulation 前，將 self.anim 設為 None
        if hasattr(self, 'anim') and self.anim is not None:
            self.anim = None
        # 清除前一次的 figure，避免重疊
        if hasattr(self, 'anim_canvas') and hasattr(self.anim_canvas, 'figure'):
            self.anim_canvas.figure.clf()
        cells = self.preview_canvas.get_cells()
        if len(cells) < 3:
            self.on_preview()
            cells = self.preview_canvas.get_cells()
            if len(cells) < 3:
                self.status_bar.showMessage("Cell number not enough, please generate cell first")
                return
        grid = VoronoiGrid(grid_shape=(len(cells), 1), mode='custom', custom_cells=cells)
        try:
            params_text = self.params_edit.toPlainText().strip()
            if params_text.startswith('params ='):
                params_text = params_text[len('params ='):].strip()
            if params_text.startswith('dict'):
                self.params = eval(params_text, {}, {})
            else:
                self.params = ast.literal_eval(params_text)
        except Exception as e:
            self.status_bar.showMessage(f"Params error: {e}")
            print(f"Params error: {e}")
            return
        
        # 如果沒有定義sigma，init_Y等於全部放0
        sigma = self.params.get('sigma', 0.0) or 0.0

        # 自動偵測物質數 n_var, 建立 init_Y
        n_var = None
        init_Y = None
        # 1. params['n_var']
        if 'n_var' in self.params:
            try:
                n_var = int(self.params['n_var'])
            except:
                n_var = None
        # 2. params['init_Y']
        if n_var is None and 'init_Y' in self.params:
            try:
                init_Y = np.array(self.params['init_Y'])
                if init_Y.ndim == 2:
                    n_var = init_Y.shape[1]
                elif init_Y.ndim == 1:
                    n_var = 1
            except:
                n_var = None
        # 3. ODE 預設值（略，需進階分析）
        # 4. fallback 預設 4
        if n_var is None:
            n_var = 4
        if init_Y is None:
            init_Y = np.random.randn(len(grid.cells), n_var) * sigma

        from gui.sim_utils import move_away_from_center, covergent_extension, repulsion_move_neighbors_no_outer, vertex_model_forces
        move_rules = []
        if self.move_away.isChecked():
            away_strength = float(self.move_away_strength.text())
            move_rules.append(lambda cells: move_away_from_center(cells, away_strength))
        if self.move_ce.isChecked():
            ce_strength = float(self.move_ce_strength.text())
            move_rules.append(lambda cells: covergent_extension(cells, ce_strength))
        if self.move_repulsion.isChecked():
            repulsion_strength = float(self.move_repulsion_strength.text())
            move_rules.append(lambda cells: repulsion_move_neighbors_no_outer(cells, repulsion_strength))
        if getattr(self, 'move_vertex_model', None) and self.move_vertex_model.isChecked():
            vertex_strength = float(self.move_vertex_model_strength.text())
            move_rules.append(lambda cells: vertex_model_forces(cells, adhesion_strength=vertex_strength))
        def combined_move_rule(cells):
            for rule in move_rules:
                cells = rule(cells)
            return cells
        move_rule = combined_move_rule if move_rules else None
        random_strength = 0.05 if self.move_random.isChecked() else 0.0
        try:
            T = float(self.T_edit.text())
        except:
            T = 30.0
            self.T_edit.setText('30.0')

        division_n = int(self.move_division_n.text()) if self.move_division.isChecked() else 0
        division_mode = self.move_division_method.currentText() if self.move_division.isChecked() else 'area'
        apoptosis_n = int(self.move_apoptosis_n.text()) if self.move_apoptosis.isChecked() else 0
        apoptosis_mode = self.move_apoptosis_method.currentText() if self.move_apoptosis.isChecked() else 'area'

        self.run_btn.setText("Cancel Simulation")
        self.progress_bar.setValue(0)
        self.status_bar.showMessage("Running simulation...")
        
        ode_func, mech_func = self.get_ode_func()
        model = BiophysicsModel(cell_count=len(grid.cells), params=self.params, ode_func=ode_func, mechanics_func=mech_func, init_Y=init_Y, vor_grid=grid, move_rule=move_rule, random_strength=random_strength)
        
        t_steps = int(T/self.params['dT']) if 'dT' in self.params and self.params['dT'] else 1
        step_interval = int(1/self.params['dT']) if 'dT' in self.params and self.params['dT'] else 1
        step_indices = list(range(0, t_steps, step_interval))
        
        if division_n:
            proliferation_steps = step_indices
        else:
            proliferation_steps = []
        if apoptosis_n:
            apoptosis_steps = step_indices
        else:
            apoptosis_steps = []
            
        self.worker = SimulationWorker(
            model=model,
            T=T,
            proliferation_steps=proliferation_steps,
            apoptosis_steps=apoptosis_steps,
            proliferation_n=division_n,
            apoptosis_n=apoptosis_n,
            proliferation_mode=division_mode,
            apoptosis_mode=apoptosis_mode
        )
        self.worker.progress_updated.connect(self.progress_bar.setValue)
        self.worker.simulation_finished.connect(self.on_simulation_finished)
        self.worker.error_occurred.connect(self.on_simulation_error)
        self.worker.warning_occurred.connect(lambda msg: self.status_bar.showMessage(msg))
        
        self.worker.start()

    def on_simulation_finished(self, history_arr, cell_positions_arr, grid):
        self.run_btn.setText("Run Simulation")
        self.progress_bar.setValue(100)
        self.sim_history = list(history_arr)
        self.cell_positions_history = list(cell_positions_arr)
        self.grid = grid
        self.status_bar.showMessage("Simulation finished. Previewing animation...")
        from voronoi_animation import VoronoiAnimator
        animator = VoronoiAnimator(grid, self.sim_history, self.get_color_func(), self.cell_positions_history, show_ticks=False, dynamic_range=False, params=self.params, color_func_str=self.color_func_edit.toPlainText().strip())
        self.anim_canvas.figure.clf()
        self.anim_canvas.setVisible(True)
        ax = self.figure.add_subplot(111)
        ax.set_xticks([])
        ax.set_yticks([])
        self.anim_canvas.draw()
        
        if hasattr(self, 'anim_cid') and self.anim_cid is not None:
            self.anim_canvas.mpl_disconnect(self.anim_cid)
        self.anim_cid = self.anim_canvas.mpl_connect('button_press_event', self.on_anim_click)
        
        # 建立動畫
        self.anim = animator.animateGUI(self.figure, interval=20, line_colors='black', line_width=1, line_alpha=0.5, point_size=15) #self.anim_canvas.figure
        
        # 初始化動畫狀態
        self.current_frame = 0
        self.anim._draw_frame(self.current_frame)
        self.anim_canvas.draw()
        self.anim.resume()

        self.is_anim_playing = False
        
        if hasattr(self.anim, 'frame_seq'):
            self.total_frames = len(list(self.anim.frame_seq))
        else:
            self.total_frames = len(self.sim_history)
        self.frame_slider.setEnabled(True)
        self.frame_slider.setMinimum(0)
        self.frame_slider.setMaximum(self.total_frames-1)
        self.frame_slider.setValue(0)
        self.play_btn.setEnabled(True)
        self.save_frame_btn.setEnabled(True)
        self.play_btn.setText("Play")
        self.update_anim_frame(0)
        self.status_bar.showMessage("Animation previewed.")
        self.save_config()

    def on_simulation_error(self, err_msg):
        self.run_btn.setText("Run Simulation")
        self.status_bar.showMessage(f"Simulation error: {err_msg}")
        self.logger.error(err_msg)

    def on_download_results(self):
        import json
        import datetime
        import numpy as np
        import os
        from openpyxl import Workbook
        if self.sim_history is None or self.cell_positions_history is None or self.grid is None:
            self.status_bar.showMessage("Please run simulation first.")
            return
        folder = QFileDialog.getExistingDirectory(self, "Select Download Folder")
        if not folder:
            self.status_bar.showMessage("Download cancelled.")
            return
        self.progress_bar.setValue(0)
        self.status_bar.showMessage("Saving results...")

        # 產生日期時間序號
        dt_prefix = datetime.datetime.now().strftime('%y%m%d_%H%M_')

        # 儲存參數（JSON格式，與config.json相同）
        param_path = os.path.join(folder, f"{dt_prefix}config.json")
        config = self.save_config(save_path=param_path)
        self.progress_bar.setValue(5)

        # 儲存動畫
        video_path = os.path.join(folder, f"{dt_prefix}simulation.mp4")
        animator = VoronoiAnimator(self.grid, self.sim_history, self.get_color_func(), self.cell_positions_history, show_ticks=False, dynamic_range=False, params=self.params, color_func_str=self.color_func_edit.toPlainText().strip())
        animator.animate(interval=5, save_path=video_path, line_colors='black', line_width=1, line_alpha=0.5, point_size=15)
        self.progress_bar.setValue(50)

        # 儲存濃度圖
        pdf_path = os.path.join(folder, f"{dt_prefix}concentration.pdf")
        labels = self.params.get('labels', [f'Y[{i}]' for i in range(self.sim_history[0].shape[1])])
        VoronoiAnimator.plot_concentration_over_time(self.sim_history, save_path=pdf_path, labels=labels, cell_indices=[0])
        self.progress_bar.setValue(55)

        # 儲存 cells.csv 與 cells.h5
        csv_path = os.path.join(folder, f"{dt_prefix}cells.csv")
        h5_path = os.path.join(folder, f"{dt_prefix}cells.h5")
        
        nY = self.sim_history[0].shape[1] if len(self.sim_history) > 0 else 0
        headers = ['T', 'step', 'x', 'y'] + [f'Y[{i}]' for i in range(nY)]
        
        try:
            T_val = float(self.T_edit.text())
        except:
            self.status_bar.showMessage("Unable to save cell grid, T is not a number or not available.")
            return
        try:
            dT_step = int(1/float(self.params.get('dT', None)))
        except:
            self.status_bar.showMessage("Unable to save cell grid, dT is not a number or not available.")
            return
        
        try:
            import pandas as pd
            records = []
            for t_idx, (Y, pos) in enumerate(zip(self.sim_history, self.cell_positions_history)):
                T_write = t_idx // dT_step
                dT_write = t_idx - (T_write * dT_step)
                self.status_bar.showMessage(f"saving {t_idx},{T_write},{dT_write}")
                for i in range(Y.shape[0]):
                    row = [T_write, dT_write, pos[i,0], pos[i,1]] + list(Y[i, :])
                    records.append(row)
                self.progress_bar.setValue(55 + int(35*(t_idx/(T_val*dT_step))))
            df = pd.DataFrame(records, columns=headers)
            df.to_csv(csv_path, index=False)
            try:
                df.to_hdf(h5_path, key='cells', mode='w')
            except Exception as e:
                self.logger.warning(f"Failed to save HDF5: {e}")
        except Exception as e:
            self.status_bar.showMessage(f"Error saving to CSV/HDF5: {e}")
            
        self.progress_bar.setValue(100)
        self.status_bar.showMessage(f"Results saved to {folder}")

    def save_config(self, save_path=None):
        ode_text = self.ode_edit.toPlainText().strip()
        params_text = self.params_edit.toPlainText().strip()
        color_func_text = self.color_func_edit.toPlainText().strip()
        #只有儲存預設config.json才會自動補上預設值
        if not ode_text and save_path is None:
            ode_text = self.get_default_ode()
        if not params_text and save_path is None:
            params_text = self.get_default_params()
        if not color_func_text and save_path is None:
            color_func_text = self.get_default_color_func()
        config = {
            'grid_shape_x': self.grid_shape_x.value(),
            'grid_shape_y': self.grid_shape_y.value(),
            'cell_dist': self.cell_dist_edit.text(),
            'pos_rand': self.pos_rand_edit.text(),
            'move_random': self.move_random.isChecked(),
            'move_random_strength': self.move_random_strength.text(),
            'move_away': self.move_away.isChecked(),
            'move_away_strength': self.move_away_strength.text(),
            'move_ce': self.move_ce.isChecked(),
            'move_ce_strength': self.move_ce_strength.text(),
            'move_repulsion': self.move_repulsion.isChecked(),
            'move_repulsion_strength': self.move_repulsion_strength.text(),
            'move_vertex_model': getattr(self, 'move_vertex_model', None) and self.move_vertex_model.isChecked(),
            'move_vertex_model_strength': getattr(self, 'move_vertex_model_strength', None) and self.move_vertex_model_strength.text(),
            'move_division': self.move_division.isChecked(),
            'move_division_n': self.move_division_n.value(),
            'move_division_method': self.move_division_method.currentIndex(),
            'move_apoptosis': self.move_apoptosis.isChecked(),
            'move_apoptosis_n': self.move_apoptosis_n.value(),
            'move_apoptosis_method': self.move_apoptosis_method.currentIndex(),
            'ode': ode_text,
            'params': params_text,
            'color_func': color_func_text,
            'T': self.T_edit.text(),
            'replicate': self.replicate_edit.text(),
            'repeats': self.repeats_edit.text(),
            'grid_mode': self.grid_mode_combo.currentIndex()
        }
        if save_path is None:
            save_path = CONFIG_PATH
        try:
            with open(save_path, 'w', encoding='utf-8') as f:
                json.dump(config, f, ensure_ascii=False, indent=2)
        except Exception as e:
            print(f"[Config] 儲存失敗: {e}")
        return config

    def load_config(self):
        # 若 config.json 不存在則用 config_example
        config = None
        if os.path.exists(CONFIG_PATH):
            try:
                with open(CONFIG_PATH, 'r', encoding='utf-8') as f:
                    config = json.load(f)
            except Exception as e:
                print(f"[Config] 載入失敗: {e}")
        elif os.path.exists(CONFIG_EXAMPLE_PATH):
            try:
                with open(CONFIG_EXAMPLE_PATH, 'r', encoding='utf-8') as f:
                    config = json.load(f)
            except Exception as e:
                print(f"[Config] 載入預設失敗: {e}")
        if config:
            self.apply_config(config)

    def apply_config(self, config):
        self.grid_shape_x.setValue(config.get('grid_shape_x', 20))
        self.grid_shape_y.setValue(config.get('grid_shape_y', 20))
        self.cell_dist_edit.setText(config.get('cell_dist', '1.0'))
        self.pos_rand_edit.setText(config.get('pos_rand', '0.25'))
        self.move_random.setChecked(config.get('move_random', False))
        self.move_random_strength.setText(config.get('move_random_strength', '0.05'))
        self.move_away.setChecked(config.get('move_away', False))
        self.move_away_strength.setText(config.get('move_away_strength', '0.1'))
        self.move_ce.setChecked(config.get('move_ce', False))
        self.move_ce_strength.setText(config.get('move_ce_strength', '0.02'))
        self.move_repulsion.setChecked(config.get('move_repulsion', False))
        self.move_repulsion_strength.setText(config.get('move_repulsion_strength', '0.08'))
        if hasattr(self, 'move_vertex_model'):
            self.move_vertex_model.setChecked(config.get('move_vertex_model', False))
            self.move_vertex_model_strength.setText(config.get('move_vertex_model_strength', '0.05'))
        self.move_division.setChecked(config.get('move_division', False))
        self.move_division_n.setValue(config.get('move_division_n', 1))
        self.move_division_method.setCurrentIndex(config.get('move_division_method', 0))
        self.move_apoptosis.setChecked(config.get('move_apoptosis', False))
        self.move_apoptosis_n.setValue(config.get('move_apoptosis_n', 1))
        self.move_apoptosis_method.setCurrentIndex(config.get('move_apoptosis_method', 0))
        ode_val = config.get('ode', '').strip()
        params_val = config.get('params', '').strip()
        color_func_val = config.get('color_func', '').strip()
        if not ode_val:
            ode_val = self.get_default_ode()
        if not params_val:
            params_val = self.get_default_params()
        if not color_func_val:
            color_func_val = self.get_default_color_func()
        self.ode_edit.setPlainText(ode_val)
        self.params_edit.setPlainText(params_val)
        self.color_func_edit.setPlainText(color_func_val)
        self.T_edit.setText(config.get('T', '30.0'))
        self.replicate_edit.setText(config.get('replicate', '5'))
        self.repeats_edit.setText(config.get('repeats', '5'))
        self.grid_mode_combo.setCurrentIndex(config.get('grid_mode', 0))

    def get_color_func(self):
        import numpy as np
        import matplotlib.pyplot as plt

        code = self.color_func_edit.toPlainText().strip()
        if not code:
            code = self.get_default_color_func()
        global_vars = {'np': np, 'plt': plt}
        local_vars = {}
        try:
            exec(code, global_vars, local_vars)
            func = local_vars.get('color_func', None)
            if func is None:
                raise Exception('color_func 未定義')
            return func
        except Exception as e:
            self.status_bar.showMessage(f"color_func 錯誤: {e}")
            exec(self.get_default_color_func(), global_vars, local_vars)
            return local_vars['color_func']

    def get_ode_func(self):
        import numpy as np
        import matplotlib.pyplot as plt
        try:
            from gui.sim_utils import get_voronoi_neighbors_wo_outer, diffusion_weighted_mean
        except ImportError:
            get_voronoi_neighbors_wo_outer = None

        code = self.ode_edit.toPlainText().strip()
        if not code:
            code = self.get_default_ode()
        global_vars = {'np': np, 'plt': plt, 'get_voronoi_neighbors_wo_outer': get_voronoi_neighbors_wo_outer, 'diffusion_weighted_mean': diffusion_weighted_mean}
        local_vars = {}
        try:
            exec(code, global_vars, local_vars)
            func = local_vars.get('ode', None)
            mech_func = local_vars.get('mechanics', None)
            if func is None:
                raise Exception('ode not defined')
            try:
                from numba import jit
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    func = jit(forceobj=True)(func)
            except ImportError:
                pass
            return func, mech_func
        except Exception as e:
            self.status_bar.showMessage(f"ODE error: {e}")
            exec(self.get_default_ode(), global_vars, local_vars)
            return local_vars['ode'], local_vars.get('mechanics', None)

    def closeEvent(self, event):
        self.save_config()
        super().closeEvent(event)

    def on_play_pause(self):
        if not hasattr(self, 'anim') or self.anim is None:
            return
        if self.is_anim_playing:
            self.anim_timer.stop()
            self.is_anim_playing = False
            self.play_btn.setText("Play")
        else:
            # If at last frame, restart
            if self.current_frame >= self.total_frames - 1:
                self.current_frame = 0
                self.update_anim_frame(0)
            self.anim_timer.start()
            self.is_anim_playing = True
            self.play_btn.setText("Pause")

    def _on_anim_timer_tick(self):
        if self.current_frame < self.total_frames - 1:
            self.current_frame += 1
            self.update_anim_frame(self.current_frame)
        else:
            self.anim_timer.stop()
            self.is_anim_playing = False
            self.play_btn.setText("Play")

    def on_frame_slider_changed(self, value):
        if not hasattr(self, 'anim') or self.anim is None:
            return
        self.current_frame = value
        self.update_anim_frame(value)
        # If playing, pause
        if self.is_anim_playing:
            self.anim_timer.stop()
            self.is_anim_playing = False
            self.play_btn.setText("Play")

    def update_anim_frame(self, frame_idx):
        if hasattr(self, 'anim') and hasattr(self.anim, '_draw_frame'):
            self.anim._draw_frame(frame_idx)
            self.anim_canvas.draw()
        self.frame_slider.blockSignals(True)
        self.frame_slider.setValue(frame_idx)
        self.frame_slider.blockSignals(False)
        self.current_frame = frame_idx

    def on_save_frame(self):
        if not hasattr(self, 'anim_canvas'):
            return
        from PyQt6.QtWidgets import QFileDialog
        file_path, _ = QFileDialog.getSaveFileName(self, "Save current image", "frame.png", "PNG Files (*.png);;All Files (*)")
        if file_path:
            self.anim_canvas.figure.savefig(file_path)
            self.status_bar.showMessage(f"Image saved: {file_path}")

    def on_anim_click(self, event):
        if getattr(self, 'is_anim_playing', False):
            return
        if event.inaxes is None:
            return
            
        x, y = event.xdata, event.ydata
        
        # Get current frame cell positions
        try:
            cells = self.cell_positions_history[self.current_frame]
            # Find nearest cell
            dists = np.sum((cells - np.array([x, y]))**2, axis=1)
            nearest_idx = np.argmin(dists)
            self.show_cell_plot(nearest_idx, self.current_frame)
        except Exception as e:
            print(f"Error handling click: {e}")

    def show_cell_plot(self, cell_idx, target_frame):
        from PyQt6.QtWidgets import QDialog, QVBoxLayout, QPushButton, QFileDialog
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
        import numpy as np
        
        # Trace backwards by position to handle index shifts (apoptosis/proliferation)
        trace = []
        current_pos = self.cell_positions_history[target_frame][cell_idx]
        
        for t in range(target_frame, -1, -1):
            positions = self.cell_positions_history[t]
            dists = np.sum((positions - current_pos)**2, axis=1)
            best_idx = np.argmin(dists)
            trace.append(self.sim_history[t][best_idx])
            current_pos = positions[best_idx]
            
        trace.reverse()
        trace_arr = np.array(trace) # shape: (target_frame+1, n_var)
        
        print(f"--- Cell Trace Debug (Clicked Index: {cell_idx}) ---")
        print(f"t=0   Concentration: {trace_arr[0]}")
        print(f"t=end Concentration: {trace_arr[-1]}")
        print("---------------------------------------")
        
        dlg = QDialog(self)
        dlg.setWindowTitle(f"Cell Concentration History (Target Frame: {target_frame}, Index: {cell_idx})")
        dlg.resize(600, 450)
        layout = QVBoxLayout(dlg)
        
        fig = Figure()
        canvas = FigureCanvasQTAgg(fig)
        layout.addWidget(canvas)
        
        btn_download = QPushButton("Download PDF")
        layout.addWidget(btn_download)
        
        ax = fig.add_subplot(111)
        labels = self.params.get('labels', [f"Var{i+1}" for i in range(trace_arr.shape[1])])
        
        for i in range(trace_arr.shape[1]):
            label = labels[i] if i < len(labels) else f"Var{i+1}"
            ax.plot(trace_arr[:, i], label=label)
            
        ax.set_xlabel("Time step")
        ax.set_ylabel("Concentration")
        ax.set_title(f"Concentration trace backwards from frame {target_frame} (position tracked)")
        ax.legend()
        
        def save_pdf():
            path, _ = QFileDialog.getSaveFileName(dlg, "Save Plot", f"cell_trace_{target_frame}_{cell_idx}.pdf", "PDF Files (*.pdf);;All Files (*)")
            if path:
                fig.savefig(path)
                self.status_bar.showMessage(f"Saved cell trace plot to {path}")
                
        btn_download.clicked.connect(save_pdf)
        
        dlg.exec()

class GUIStatusHandler(logging.Handler):
    def __init__(self, status_bar):
        super().__init__()
        self.status_bar = status_bar
    def emit(self, record):
        msg = self.format(record)
        self.status_bar.showMessage(msg)
        print(msg)  # 顯示到 terminal

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec()) 