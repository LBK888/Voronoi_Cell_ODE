import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QGroupBox, QFormLayout, QLineEdit, QTextEdit, QComboBox,
    QCheckBox, QFileDialog, QProgressBar, QStatusBar, QSpinBox, QSlider,
    QDoubleSpinBox, QScrollArea, QTabWidget, QSplitter,
    QSizePolicy, QFrame, QDialog, QStackedWidget, QDialogButtonBox, QTextBrowser
)
from PyQt6.QtCore import Qt, QTimer, QEvent, QUrl
from PyQt6.QtGui import QFont, QColor, QPalette, QAction, QKeySequence, QDesktopServices
from pathlib import Path
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
from voronoi_animation import VoronoiAnimator
from gui.simulation_session import (
    SimulationSession,
)
from gui.simulation_worker import SimulationWorker
from gui.results_export_worker import ExportRequest, ResultsExportWorker
from gui.figure_export import install_figure_menu, save_figure_dialog
from gui.experiments_panel import ExperimentsPanel
from gui.input_wheel_guard import disable_wheel_value_changes
from gui.app_paths import data_path, resource_path, configure_ffmpeg
from gui.model_files import model_grid, model_request, normalize_model, read_model, template_files, user_templates, write_json
from automation_cli import cli_enabled, set_cli_enabled

CONFIG_PATH = str(data_path('config.json'))
CONFIG_EXAMPLE_PATH = str(resource_path('config_example.json'))

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        disable_wheel_value_changes()
        configure_ffmpeg()
        self.logger = setup_logger()
        self.setWindowTitle("Voronoi Notch Simulation")
        self.resize(1280, 880)
        self._applying_config = False
        self._imported_cells = None
        self.setAcceptDrops(True)
        self.init_ui()
        self.load_config()
        self.sim_history = None
        self.cell_positions_history = None
        self.simulation_result = None
        self.grid = None
        self.params = None
        self.worker = None
        self.export_worker = None
        self.logger.addHandler(GUIStatusHandler(self.status_bar))
        QApplication.instance().installEventFilter(self)

    def init_ui(self):
        menubar = self.menuBar()
        self.template_menu = menubar.addMenu('Templates')
        open_model = self.template_menu.addAction('Open Model JSON...')
        open_model.setShortcut(QKeySequence('Ctrl+O'))
        open_model.triggered.connect(self.open_model_dialog)
        self.template_menu.addAction('Save Current Model as Template...', self.save_template_dialog)
        self.template_menu.addAction('Open User Templates Folder', self.open_templates_folder)
        self.template_menu.addSeparator()
        self.template_actions = []
        self.template_menu.aboutToShow.connect(self.refresh_templates_menu)
        self.refresh_templates_menu()
        model_menu = menubar.addMenu('Model')
        self.open_node_editor_action = QAction('Open Visual Node Editor...', self)
        self.open_node_editor_action.setShortcut(QKeySequence('Ctrl+Shift+N'))
        self.open_node_editor_action.setStatusTip(
            'Edit the current ODE model as a visual node graph'
        )
        self.open_node_editor_action.triggered.connect(self.open_node_editor)
        model_menu.addAction(self.open_node_editor_action)
        agent_menu = menubar.addMenu('Agent')
        self.enable_cli_action = agent_menu.addAction('Enable Agent CLI')
        self.enable_cli_action.setCheckable(True)
        self.enable_cli_action.setChecked(cli_enabled())
        self.enable_cli_action.toggled.connect(self.toggle_agent_cli)
        agent_menu.addAction('Connection Instructions...', self.show_agent_instructions)

        main_widget = QWidget()
        main_layout = QVBoxLayout(main_widget)
        main_layout.setContentsMargins(10, 10, 10, 8)
        workspace_splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left: Settings
        settings_widget = QWidget()
        settings_widget.setMinimumWidth(320)
        settings_widget.setMaximumWidth(500)
        settings_layout = QVBoxLayout(settings_widget)
        settings_layout.setContentsMargins(0, 0, 6, 0)
        settings_title = QLabel("Simulation setup")
        settings_title.setObjectName("sectionTitle")
        settings_layout.addWidget(settings_title)
        self.settings_tabs = QTabWidget()
        grid_tab = QWidget()
        grid_tab_layout = QVBoxLayout(grid_tab)
        movement_tab = QWidget()
        movement_tab_layout = QVBoxLayout(movement_tab)
        model_tab = QWidget()
        model_tab_layout = QVBoxLayout(model_tab)
        simulation_tab = QWidget()
        simulation_tab_layout = QVBoxLayout(simulation_tab)

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
        grid_tab_layout.addWidget(grid_mode_box)

        # Cell grid parameters
        grid_param_box = QGroupBox("Cell Grid Parameters")
        grid_param_form = QFormLayout()
        grid_param_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        grid_param_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )
        grid_shape_row = QHBoxLayout()
        self.grid_shape_x = QSpinBox()
        self.grid_shape_x.setRange(1, 200)
        self.grid_shape_x.setValue(20)
        self.grid_shape_x.setMinimumWidth(60)
        self.grid_shape_y = QSpinBox()
        self.grid_shape_y.setRange(1, 200)
        self.grid_shape_y.setValue(20)
        self.grid_shape_y.setMinimumWidth(60)
        grid_shape_row.addWidget(QLabel("X:"))
        grid_shape_row.addWidget(self.grid_shape_x)
        grid_shape_row.addWidget(QLabel("Y:"))
        grid_shape_row.addWidget(self.grid_shape_y)
        grid_shape_row.setStretch(1, 1)
        grid_shape_row.setStretch(3, 1)
        grid_param_form.addRow("Grid Shape:", grid_shape_row)
        self.cell_dist_edit = QLineEdit("1.0")
        self.pos_rand_edit = QLineEdit("0.25")
        grid_param_form.addRow("Cell Distance:", self.cell_dist_edit)
        grid_param_form.addRow("Position Randomness:", self.pos_rand_edit)
        grid_param_box.setLayout(grid_param_form)
        grid_tab_layout.addWidget(grid_param_box)

        # Cell movement modes
        move_mode_box = QGroupBox("Cell Movement Modes")
        move_mode_layout = QVBoxLayout()
        self.move_random = QCheckBox("Random movement")
        self.move_random_strength = QLineEdit("0.05")
        move_random_row = QHBoxLayout()
        move_random_row.addWidget(self.move_random)
        move_random_row.addWidget(self.move_random_strength)
        move_mode_layout.addLayout(move_random_row)
        self.move_away = QCheckBox("Move away from center")
        self.move_away_strength = QLineEdit("0.1")
        move_away_row = QHBoxLayout()
        move_away_row.addWidget(self.move_away)
        move_away_row.addWidget(self.move_away_strength)
        move_mode_layout.addLayout(move_away_row)
        self.move_ce = QCheckBox("Convergent extension")
        self.move_ce_strength = QLineEdit("0.02")
        move_ce_row = QHBoxLayout()
        move_ce_row.addWidget(self.move_ce)
        move_ce_row.addWidget(self.move_ce_strength)
        move_mode_layout.addLayout(move_ce_row)
        self.move_repulsion = QCheckBox("Neighbor repulsion")
        self.move_repulsion_strength = QLineEdit("0.08")
        move_repulsion_row = QHBoxLayout()
        move_repulsion_row.addWidget(self.move_repulsion)
        move_repulsion_row.addWidget(self.move_repulsion_strength)
        move_mode_layout.addLayout(move_repulsion_row)
        self.move_vertex_model = QCheckBox("Vertex-model adhesion")
        self.move_vertex_model_strength = QLineEdit("0.05")
        move_vertex_model_row = QHBoxLayout()
        move_vertex_model_row.addWidget(self.move_vertex_model)
        move_vertex_model_row.addWidget(self.move_vertex_model_strength)
        move_mode_layout.addLayout(move_vertex_model_row)
        self.move_division = QCheckBox("Cell division")
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
        self.move_apoptosis = QCheckBox("Apoptosis")
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
        for strength_edit in (
            self.move_random_strength,
            self.move_away_strength,
            self.move_ce_strength,
            self.move_repulsion_strength,
            self.move_vertex_model_strength,
        ):
            strength_edit.setFixedWidth(88)
            strength_edit.setToolTip("Movement strength")
        move_mode_box.setLayout(move_mode_layout)
        movement_tab_layout.addWidget(move_mode_box)

        # ODE editor
        ode_box = QGroupBox("ODE Function")
        ode_layout = QVBoxLayout()
        ode_hint = QLabel(
            "Edit the model as Python source or open the same model as a visual graph."
        )
        ode_hint.setObjectName("fieldHint")
        ode_hint.setWordWrap(True)
        ode_layout.addWidget(ode_hint)
        
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
        self.open_ode_editor_btn.setToolTip("Edit the ODE in a larger Python code editor")
        self.open_node_editor_btn = QPushButton("Visual Node Editor...")
        self.open_node_editor_btn.setObjectName("accentButton")
        self.open_node_editor_btn.setToolTip(
            "Open the current ODE, parameters and color function as a node graph (Ctrl+Shift+N)"
        )
        ode_editor_actions = QHBoxLayout()
        ode_editor_actions.addWidget(self.open_ode_editor_btn)
        ode_editor_actions.addWidget(self.open_node_editor_btn)
        ode_layout.addLayout(ode_editor_actions)

        ode_box.setLayout(ode_layout)
        model_tab_layout.addWidget(ode_box)

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
        model_tab_layout.addWidget(color_func_box)

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

        self.generate_params_btn = QPushButton("Generate Parameter Controls")
        params_layout.addWidget(self.generate_params_btn)

        self.dynamic_params_scroll = QScrollArea()
        self.dynamic_params_scroll.setWidgetResizable(True)
        self.dynamic_params_widget = QWidget()
        self.dynamic_params_layout = QFormLayout()
        self.dynamic_params_widget.setLayout(self.dynamic_params_layout)
        self.dynamic_params_scroll.setWidget(self.dynamic_params_widget)
        self.dynamic_params_scroll.setMinimumHeight(120)
        params_layout.addWidget(self.dynamic_params_scroll)

        params_box.setLayout(params_layout)
        model_tab_layout.addWidget(params_box)

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
        simulation_tab_layout.addWidget(sim_param_box)

        # Preview button
        self.preview_btn = QPushButton("Generate Preview")
        # Apply movement button
        self.apply_move_btn = QPushButton("Apply Movement Once")
        grid_actions = QVBoxLayout()
        grid_actions.addWidget(self.preview_btn)
        grid_actions.addWidget(self.apply_move_btn)
        grid_tab_layout.addLayout(grid_actions)

        grid_tab_layout.addStretch(1)
        movement_tab_layout.addStretch(1)
        model_tab_layout.addStretch(1)
        simulation_tab_layout.addStretch(1)
        for tab, title in (
            (grid_tab, "Grid"),
            (movement_tab, "Dynamics"),
            (model_tab, "Model"),
            (simulation_tab, "Run"),
        ):
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QFrame.Shape.NoFrame)
            scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            scroll.setWidget(tab)
            self.settings_tabs.addTab(scroll, title)
        self.experiments_panel = ExperimentsPanel(self._experiment_simulation_request, self)
        experiment_scroll = QScrollArea()
        experiment_scroll.setWidgetResizable(True)
        experiment_scroll.setFrameShape(QFrame.Shape.NoFrame)
        experiment_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        experiment_scroll.setWidget(self.experiments_panel)
        self.experiments_tab_index = self.settings_tabs.addTab(experiment_scroll, 'Experiments')
        settings_layout.addWidget(self.settings_tabs)

        settings_widget.setLayout(settings_layout)
        workspace_splitter.addWidget(settings_widget)

        # Right: Preview
        preview_widget = QWidget()
        preview_layout = QVBoxLayout(preview_widget)
        preview_layout.setContentsMargins(6, 0, 0, 0)
        preview_title = QLabel("Preview & results")
        preview_title.setObjectName("sectionTitle")
        preview_layout.addWidget(preview_title)
        self.result_stack = QStackedWidget()
        preview_layout.addWidget(self.result_stack, 1)
        simulation_view = QWidget()
        preview_layout = QVBoxLayout(simulation_view)
        preview_layout.setContentsMargins(0, 0, 0, 0)
        self.preview_canvas = PreviewCanvas()
        self.preview_canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.preview_canvas.connect_events()
        preview_layout.addWidget(self.preview_canvas)
        # 連接 cell_count_changed signal
        self.preview_canvas.cell_count_changed.connect(lambda count: self.status_bar.showMessage(f"New cell number: {count}"))
        # 新增動畫嵌入區
        self.figure = Figure()
        self.anim_canvas = FigureCanvas(self.figure)
        self.anim_canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.anim_canvas.setVisible(False)
        preview_layout.addWidget(self.anim_canvas)
        install_figure_menu(self.anim_canvas, self, enabled=lambda: self.sim_history is not None)
        # --- 新增動畫控制元件 ---
        anim_ctrl_layout = QHBoxLayout()
        self.play_btn = QPushButton("Play")
        self.play_btn.setEnabled(False)
        self.frame_slider = QSlider(Qt.Orientation.Horizontal)
        self.frame_slider.setEnabled(False)
        self.save_frame_btn = QPushButton("Save Frame")
        self.save_frame_btn.setEnabled(False)
        anim_ctrl_layout.addWidget(self.play_btn)
        anim_ctrl_layout.addWidget(self.frame_slider)
        anim_ctrl_layout.addWidget(self.save_frame_btn)
        preview_layout.addLayout(anim_ctrl_layout)
        # ---
        self.run_btn = QPushButton("Run Simulation")
        self.run_btn.setObjectName("primaryButton")
        self.download_btn = QPushButton("Export Results...")
        self.download_btn.setEnabled(False)
        self.download_btn.setToolTip("Export MP4, PDF, CSV, HDF5 and configuration files in the background")
        primary_actions = QHBoxLayout()
        primary_actions.addStretch(1)
        primary_actions.addWidget(self.download_btn)
        primary_actions.addWidget(self.run_btn)
        preview_layout.addLayout(primary_actions)
        self.result_stack.addWidget(simulation_view)
        self.result_stack.addWidget(self.experiments_panel.result_widget)
        workspace_splitter.addWidget(preview_widget)
        workspace_splitter.setStretchFactor(0, 0)
        workspace_splitter.setStretchFactor(1, 1)
        workspace_splitter.setSizes([380, 900])
        main_layout.addWidget(workspace_splitter)

        # Status bar and progress bar
        self.status_bar = QStatusBar()
        self.progress_bar = QProgressBar()
        self.progress_bar.setFixedWidth(220)
        self.status_bar.addPermanentWidget(self.progress_bar)
        self.setStatusBar(self.status_bar)
        self.experiments_panel.status_changed.connect(self.status_bar.showMessage)
        self.experiments_panel.busy_changed.connect(self._on_experiment_busy)
        self.settings_tabs.currentChanged.connect(self._on_settings_tab_changed)

        self.setCentralWidget(main_widget)
        self.setStyleSheet("""
            QMainWindow { background: #f3f5f7; }
            QLabel#sectionTitle { font-size: 17px; font-weight: 600; color: #25313c; padding: 2px 0 6px 2px; }
            QTabWidget::pane { border: 1px solid #d5dbe1; border-radius: 6px; background: white; top: -1px; }
            QTabBar::tab { padding: 8px 11px; background: #e8edf2; color: #52606d; }
            QTabBar::tab:selected { background: white; color: #1261a0; font-weight: 600; }
            QGroupBox { font-weight: 600; border: 1px solid #d8dee5; border-radius: 6px; margin-top: 12px; padding: 10px 8px 8px 8px; }
            QGroupBox::title { subcontrol-origin: margin; left: 9px; padding: 0 4px; }
            QLabel#fieldHint { color: #52606d; font-weight: 400; padding: 0 1px 4px 1px; }
            QPushButton { min-height: 28px; padding: 2px 12px; }
            QPushButton#primaryButton { background: #1261a0; color: white; border: none; border-radius: 4px; min-width: 150px; font-weight: 600; }
            QPushButton#primaryButton:hover { background: #0d4f85; }
            QPushButton#accentButton { color: #1A6FB5; border: 1px solid #1A6FB5; border-radius: 4px; font-weight: 600; }
            QPushButton#accentButton:hover { background: #e8f2fb; }
            QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox { min-height: 25px; }
        """)

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
                widget.editingFinished.connect(lambda *_: self.save_config())
            elif hasattr(widget, 'valueChanged'):
                widget.valueChanged.connect(lambda *_: self.save_config())
            elif hasattr(widget, 'currentIndexChanged'):
                widget.currentIndexChanged.connect(lambda *_: self.save_config())
            elif hasattr(widget, 'stateChanged'):
                widget.stateChanged.connect(lambda *_: self.save_config())

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

        self.clear_parameter_controls()

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

    def clear_parameter_controls(self):
        while self.dynamic_params_layout.rowCount():
            self.dynamic_params_layout.removeRow(0)
        self.dynamic_widgets = {}

    def refresh_templates_menu(self):
        for action in self.template_actions:
            self.template_menu.removeAction(action)
            action.deleteLater()
        self.template_actions = []
        try:
            files = template_files()
        except OSError as error:
            action = self.template_menu.addAction(f'Cannot read templates: {error}')
            action.setEnabled(False)
            self.template_actions.append(action)
            return
        for path in files:
            action = self.template_menu.addAction(path.stem)
            action.setStatusTip(str(path))
            action.triggered.connect(lambda checked, file=path: self.open_model(file))
            self.template_actions.append(action)
        if not files:
            action = self.template_menu.addAction('No template JSON files found')
            action.setEnabled(False)
            self.template_actions.append(action)

    def open_model_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Open Model JSON', str(user_templates()), 'Model JSON (*.json)')
        if path:
            self.open_model(path)

    def save_template_dialog(self):
        try:
            directory = user_templates()
            directory.mkdir(parents=True, exist_ok=True)
            path, _ = QFileDialog.getSaveFileName(self, 'Save Current Model as Template',
                                                 str(directory / 'model.json'), 'Model JSON (*.json)')
            if not path:
                return
            path = Path(path)
            if path.suffix.lower() != '.json':
                path = path.with_name(path.name + '.json')
                if path.exists():
                    self.status_bar.showMessage('Choose the existing .json file explicitly to confirm replacement.')
                    return
            self.save_template(path)
            self.status_bar.showMessage(f'Saved template: {path}')
        except (OSError, ValueError, SyntaxError) as error:
            self.status_bar.showMessage(f'Cannot save template: {error}')

    def save_template(self, path):
        write_json(path, normalize_model(self.current_config()))

    def open_templates_folder(self):
        try:
            directory = user_templates()
            directory.mkdir(parents=True, exist_ok=True)
            if not QDesktopServices.openUrl(QUrl.fromLocalFile(str(directory.resolve()))):
                self.status_bar.showMessage(f'Templates folder: {directory}')
        except OSError as error:
            self.status_bar.showMessage(f'Cannot open templates folder: {error}')

    def open_model(self, path):
        if self.experiments_panel.worker is not None or any(
            worker is not None and worker.isRunning() for worker in (self.worker, self.export_worker)
        ):
            self.status_bar.showMessage('Wait for the current job to finish before opening a model.')
            return False
        try:
            config = read_model(path)
            grid = model_grid(config) if config['grid_mode'] == 3 and 'cell_positions' in config else None
        except Exception as error:
            self.status_bar.showMessage(f'Cannot open {Path(path).name}: {error}')
            return False
        self.settings_tabs.setCurrentIndex(0)
        self.apply_config(config)
        self.clear_parameter_controls()
        self.params = None
        self.grid = grid
        self.preview_canvas.set_grid(grid)
        self.preview_canvas.show()
        self.anim_canvas.hide()
        self.anim_timer.stop()
        self.anim = None
        self.is_anim_playing = False
        self.play_btn.setText('Play')
        self.simulation_result = self.sim_history = self.cell_positions_history = None
        self.anim_canvas.figure.clear()
        for widget in (self.play_btn, self.frame_slider, self.save_frame_btn, self.download_btn):
            widget.setEnabled(False)
        panel = self.experiments_panel
        panel.result = None
        panel.figure.clear()
        panel.canvas.draw_idle()
        panel.summary.setText('Configure the X and Y axes, then run an experiment.')
        panel.progress_label.clear()
        panel.pdf_btn.setEnabled(False)
        panel.data_btn.setEnabled(False)
        self.progress_bar.setValue(0)
        self.save_config()
        self.status_bar.showMessage(f'Loaded {Path(path).name}. Press Run Simulation when ready.')
        return True

    def eventFilter(self, watched, event):
        # Text editors also accept drops. Route model files before they can
        # insert a filesystem URL into a parameter or source-code field.
        if event.type() in (QEvent.Type.DragEnter, QEvent.Type.DragMove, QEvent.Type.Drop):
            if isinstance(watched, QWidget) and (watched is self or self.isAncestorOf(watched)):
                urls = event.mimeData().urls()
                if len(urls) == 1 and urls[0].isLocalFile():
                    path = Path(urls[0].toLocalFile())
                    if path.suffix.lower() == '.json':
                        event.acceptProposedAction()
                        if event.type() == QEvent.Type.Drop:
                            self.open_model(path)
                        return True
        return super().eventFilter(watched, event)

    def toggle_agent_cli(self, enabled):
        try:
            set_cli_enabled(enabled)
            self.status_bar.showMessage('Agent CLI enabled.' if enabled else 'Agent CLI disabled.')
        except OSError as error:
            blocked = self.enable_cli_action.blockSignals(True)
            self.enable_cli_action.setChecked(cli_enabled())
            self.enable_cli_action.blockSignals(blocked)
            self.status_bar.showMessage(f'Cannot change Agent CLI setting: {error}')

    def show_agent_instructions(self):
        dialog = QDialog(self)
        dialog.setWindowTitle('Agent CLI Connection Instructions')
        dialog.resize(820, 640)
        layout = QVBoxLayout(dialog)
        launcher = sys.executable if getattr(sys, 'frozen', False) else str(resource_path('windows_launcher.py'))
        location = QLabel(launcher)
        location.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        location.setWordWrap(True)
        layout.addWidget(location)
        guide = QTextBrowser()
        instructions = resource_path('docs/agent-cli.md').read_text(encoding='utf-8')
        guide.setMarkdown(instructions)
        layout.addWidget(guide)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        copy = buttons.addButton('Copy instructions', QDialogButtonBox.ButtonRole.ActionRole)
        copy.clicked.connect(lambda: QApplication.clipboard().setText(f'{launcher}\n\n{instructions}'))
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.exec()

    def open_node_editor(self):
        from gui.node_editor import NodeEditorDialog, NodeEditorDocument

        document = NodeEditorDocument(
            ode_source=self.ode_edit.toPlainText(),
            params_source=self.params_edit.toPlainText(),
            color_func_source=self.color_func_edit.toPlainText(),
        )
        dlg = NodeEditorDialog(parent=self, document=document)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            self.status_bar.showMessage("Node editor closed without applying changes", 3000)
            return

        result = dlg.result_document()
        if result is None:
            return

        self.ode_edit.setPlainText(result.ode_source)
        self.params_edit.setPlainText(result.params_source)
        self.color_func_edit.setPlainText(result.color_func_source)
        self.on_generate_params()
        self.save_config()
        self.status_bar.showMessage("Model updated from Visual Node Editor", 5000)

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

    def on_preview(self, *, use_saved_coordinates=False):
        self.anim_timer.stop()
        self.preview_canvas.setVisible(True)
        self.anim_canvas.setVisible(False)
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
            if mode == 'import' and use_saved_coordinates and self._imported_cells is not None:
                grid = VoronoiGrid.from_cells(self._imported_cells, cell_dist=cell_dist)
            elif mode == 'import':
                file_path, _ = QFileDialog.getOpenFileName(self, "Import Cell Coordinates", "", "Table Files (*.xlsx *.xls *.csv);;Text Files (*.txt);;All Files (*)")
                if not file_path:
                    self.status_bar.showMessage("Import cancelled.")
                    return
                grid = VoronoiGrid(grid_shape=grid_shape, cell_dist=cell_dist, pos_rand=pos_rand, mode=mode, import_path=file_path)
                self._imported_cells = grid.copy_cells().tolist()
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
        self.grid = grid
        self.preview_canvas.set_grid(grid)

    def on_apply_movement(self):
        # 取得目前細胞座標
        grid = self.preview_canvas.get_grid()
        if grid is None:
            self.on_preview(use_saved_coordinates=True)
            grid = self.preview_canvas.get_grid()
            if grid is None:
                return
        cells = grid.copy_cells()
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
        grid.replace_cells(cells)
        if self.move_division.isChecked():
            n = int(self.move_division_n.text())
            mode = self.move_division_method.currentText()
            grid.cell_proliferation(n=n, mode=mode)
        if self.move_apoptosis.isChecked():
            n = int(self.move_apoptosis_n.text())
            mode = self.move_apoptosis_method.currentText()
            grid.cell_apoptosis(n=n, mode=mode)

        self.grid = grid
        self.preview_canvas.set_grid(grid)
        # 顯示最新細胞數
        cell_count = len(grid.cells)
        self.status_bar.showMessage(f"Movement applied. Cell number: {cell_count}")

    def on_run_simulation(self):
        if self.experiments_panel.worker is not None:
            self.status_bar.showMessage('Wait for the experiment to finish or cancel it first.')
            return
        if self.worker is not None and self.worker.isRunning():
            self.worker.cancel()
            self.run_btn.setText("Run Simulation")
            self.status_bar.showMessage("Simulation cancelled.")
            return

        # 每次 simulation 前，將 self.anim 設為 None
        if hasattr(self, 'anim') and self.anim is not None:
            self.anim = None
        # 清除前一次的 figure，避免重疊
        if hasattr(self, 'anim_canvas') and hasattr(self.anim_canvas, 'figure'):
            self.anim_canvas.figure.clf()
        preview_grid = self.preview_canvas.get_grid()
        if preview_grid is None or len(preview_grid.cells) < 3:
            self.on_preview(use_saved_coordinates=True)
            preview_grid = self.preview_canvas.get_grid()
            if preview_grid is None or len(preview_grid.cells) < 3:
                self.status_bar.showMessage("Cell number not enough, please generate cell first")
                return
        grid = VoronoiGrid.from_cells(preview_grid.copy_cells(), cell_dist=preview_grid.cell_dist)
        
        # 如果沒有定義sigma，init_Y等於全部放0
        # 自動偵測物質數 n_var, 建立 init_Y
        # 3. ODE 預設值（略，需進階分析）
        # 4. fallback 預設 4
        try:
            float(self.T_edit.text())
        except ValueError:
            self.T_edit.setText('30.0')

        try:
            request = model_request(self.current_config(), grid)
            session = SimulationSession(request)
        except Exception as error:
            self.status_bar.showMessage(f"Simulation configuration error: {error}")
            self.logger.error(f"Simulation configuration failed: {error}")
            return

        self.run_btn.setText("Cancel Simulation")
        self.download_btn.setEnabled(False)
        self.progress_bar.setValue(0)
        self.status_bar.showMessage("Running simulation...")
        
        self.worker = SimulationWorker(session=session)
        self.worker.progress_updated.connect(self.progress_bar.setValue)
        self.worker.simulation_finished.connect(self.on_simulation_finished)
        self.worker.error_occurred.connect(self.on_simulation_error)
        self.worker.warning_occurred.connect(lambda msg: self.status_bar.showMessage(msg))
        
        self.worker.start()

    def on_simulation_finished(self, result):
        self.run_btn.setText("Run Simulation")
        self.progress_bar.setValue(100)
        self.simulation_result = result
        self.sim_history = list(result.history)
        self.cell_positions_history = list(result.cell_positions)
        self.grid = result.grid
        self.params = dict(result.params)
        self.status_bar.showMessage("Simulation finished. Previewing animation...")
        from voronoi_animation import VoronoiAnimator
        animator = VoronoiAnimator(self.grid, self.sim_history, self.get_color_func(), self.cell_positions_history, show_ticks=False, dynamic_range=False, params=self.params, color_func_str=self.color_func_edit.toPlainText().strip())
        self.anim_canvas.figure.clf()
        self.preview_canvas.setVisible(False)
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
        self.anim.pause()

        self.is_anim_playing = False
        
        self.total_frames = len(self.sim_history)
        self.frame_slider.setEnabled(True)
        self.frame_slider.setMinimum(0)
        self.frame_slider.setMaximum(self.total_frames-1)
        self.frame_slider.setValue(0)
        self.play_btn.setEnabled(True)
        self.save_frame_btn.setEnabled(True)
        self.download_btn.setEnabled(True)
        self.play_btn.setText("Play")
        self.update_anim_frame(0)
        self.status_bar.showMessage("Animation previewed.")
        self.save_config()

    def _experiment_simulation_request(self):
        if ((getattr(self, 'worker', None) is not None and self.worker.isRunning())
                or (getattr(self, 'export_worker', None) is not None and self.export_worker.isRunning())):
            raise ValueError('Wait for the current simulation or export to finish first')
        preview_grid = self.preview_canvas.get_grid()
        if preview_grid is None:
            self.on_preview(use_saved_coordinates=True)
            preview_grid = self.preview_canvas.get_grid()
        if preview_grid is None:
            raise ValueError('Generate a valid grid first')
        grid = VoronoiGrid.from_cells(preview_grid.copy_cells(), cell_dist=preview_grid.cell_dist)
        return model_request(self.current_config(), grid)

    def _on_settings_tab_changed(self, index):
        experiments = index == self.experiments_tab_index
        self.result_stack.setCurrentIndex(1 if experiments else 0)
        if experiments and self.experiments_panel.worker is None:
            self.experiments_panel.refresh_targets()

    def _on_experiment_busy(self, busy):
        # Model code uses NumPy's global RNG. Serialize simulations and sweeps
        # to preserve reproducibility and prevent conflicting result updates.
        self.run_btn.setEnabled(not busy)
        self.download_btn.setEnabled(not busy and self.simulation_result is not None)
        self.preview_btn.setEnabled(not busy)
        self.apply_move_btn.setEnabled(not busy)
        self.menuBar().setEnabled(not busy)
        for index in range(self.experiments_tab_index):
            self.settings_tabs.setTabEnabled(index, not busy)

    def on_simulation_error(self, err_msg):
        self.run_btn.setText("Run Simulation")
        self.download_btn.setEnabled(self.sim_history is not None)
        self.status_bar.showMessage(f"Simulation error: {err_msg}")
        self.logger.error(err_msg)

    def on_download_results(self):
        self._start_results_export()

    def _start_results_export(self):
        if self.experiments_panel.worker is not None:
            self.status_bar.showMessage('Wait for the experiment to finish or cancel it first.')
            return
        import datetime

        if getattr(self, 'export_worker', None) is not None and self.export_worker.isRunning():
            self.export_worker.cancel()
            self.download_btn.setEnabled(False)
            self.status_bar.showMessage("Cancelling export after the current frame...")
            return
        if self.simulation_result is None:
            self.status_bar.showMessage("Please run simulation first.")
            return
        folder = QFileDialog.getExistingDirectory(self, "Select Download Folder")
        if not folder:
            self.status_bar.showMessage("Download cancelled.")
            return

        self.progress_bar.setValue(0)
        self.status_bar.showMessage("Preparing background export...")
        prefix = datetime.datetime.now().strftime('%y%m%d_%H%M_')
        request = ExportRequest(
            folder=folder,
            prefix=prefix,
            config=self.save_config(),
            simulation=self.simulation_result,
            color_func=self.get_color_func(),
            color_func_source=self.color_func_edit.toPlainText().strip(),
        )
        self.export_worker = ResultsExportWorker(
            request=request,
            parent=self,
        )
        self.export_worker.progress_updated.connect(self.on_export_progress)
        self.export_worker.export_finished.connect(self.on_export_finished)
        self.export_worker.error_occurred.connect(self.on_export_error)
        self.export_worker.cancelled.connect(self.on_export_cancelled)
        self.anim_timer.stop()
        self.is_anim_playing = False
        self.play_btn.setText("Play")
        self.play_btn.setEnabled(False)
        self.frame_slider.setEnabled(False)
        self.save_frame_btn.setEnabled(False)
        self.download_btn.setText("Cancel Export")
        self.run_btn.setEnabled(False)
        self.export_worker.start()

    def on_export_progress(self, value, message):
        self.progress_bar.setValue(value)
        self.status_bar.showMessage(message)

    def _reset_export_controls(self):
        self.download_btn.setText("Export Results...")
        self.download_btn.setEnabled(True)
        self.run_btn.setEnabled(True)
        has_results = self.sim_history is not None
        self.play_btn.setEnabled(has_results)
        self.frame_slider.setEnabled(has_results)
        self.save_frame_btn.setEnabled(has_results)

    def on_export_finished(self, folder):
        self.progress_bar.setValue(100)
        self._reset_export_controls()
        self.status_bar.showMessage(f"Results saved to {folder}")

    def on_export_cancelled(self):
        self._reset_export_controls()
        self.status_bar.showMessage("Export cancelled. Partially written files may remain.")

    def on_export_error(self, error_message):
        self._reset_export_controls()
        self.status_bar.showMessage("Export failed. See app.log for details.")
        self.logger.error(error_message)

    def current_config(self):
        ode_text = self.ode_edit.toPlainText().strip()
        params_text = self.params_edit.toPlainText().strip()
        color_func_text = self.color_func_edit.toPlainText().strip()
        if not ode_text:
            ode_text = self.get_default_ode()
        if not params_text:
            params_text = self.get_default_params()
        if not color_func_text:
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
        config['experiments'] = self.experiments_panel.settings()
        if config['grid_mode'] == 3:
            grid = self.preview_canvas.get_grid()
            if grid is not None and grid.mode in ('custom', 'import'):
                config['cell_positions'] = grid.copy_cells().tolist()
            elif self._imported_cells is not None:
                config['cell_positions'] = self._imported_cells
        return config

    def save_config(self, save_path=None):
        if self._applying_config:
            return
        if save_path is not None and not isinstance(save_path, (str, os.PathLike)):
            save_path = None
        config = self.current_config()
        if save_path is None:
            save_path = CONFIG_PATH
        try:
            write_json(save_path, config)
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
        self._applying_config = True
        try:
            self._apply_config_fields(config)
        finally:
            self._applying_config = False

    def _apply_config_fields(self, config):
        self._imported_cells = config.get('cell_positions')
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
        self.experiments_panel.restore_settings(config.get('experiments', {}))

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
            return func, mech_func
        except Exception as e:
            self.status_bar.showMessage(f"ODE error: {e}")
            exec(self.get_default_ode(), global_vars, local_vars)
            return local_vars['ode'], local_vars.get('mechanics', None)

    def closeEvent(self, event):
        experiment_worker = self.experiments_panel.worker
        if experiment_worker is not None and experiment_worker.isRunning():
            experiment_worker.cancel()
            if not experiment_worker.wait(3000):
                self.status_bar.showMessage('Finishing experiment cancellation; close again shortly.')
                event.ignore()
                return
        if self.worker is not None and self.worker.isRunning():
            self.worker.cancel()
            if not self.worker.wait(3000):
                self.status_bar.showMessage('Finishing simulation cancellation; close again shortly.')
                event.ignore()
                return
        if getattr(self, 'export_worker', None) is not None and self.export_worker.isRunning():
            self.export_worker.cancel()
            if not self.export_worker.wait(3000):
                self.status_bar.showMessage("Finishing export cancellation; close again shortly.")
                event.ignore()
                return
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
        try:
            path = save_figure_dialog(self, self.anim_canvas.figure)
            if path:
                self.status_bar.showMessage(f'Image saved: {path}')
        except Exception as error:
            self.status_bar.showMessage(f'Image export failed: {error}')

    def on_anim_click(self, event):
        if event.button != 1:
            return
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
        install_figure_menu(canvas, dlg, 'cell_trace.pdf')
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
