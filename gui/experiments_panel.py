"""Experiment controls and heatmap view, embedded after the Run settings tab."""

from pathlib import Path

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QPushButton, QSpinBox, QVBoxLayout, QWidget,
)

from gui.experiment_worker import ExperimentWorker
from gui.experiment_output import draw_heatmap, export_experiment_data
from gui.experiments import ExperimentRequest, METRICS, SweepAxis
from gui.figure_export import install_figure_menu, save_figure_dialog
from gui.simulation_session import SimulationSession


class AxisControls(QGroupBox):
    def __init__(self, title, parent=None):
        super().__init__(title, parent)
        form = QFormLayout(self)
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.kind = QComboBox()
        self.kind.addItem('Model parameter', 'parameter')
        self.kind.addItem('Initial concentration (all cells)', 'initial')
        self.target = QComboBox()
        self.minimum, self.maximum = QLineEdit('0'), QLineEdit('1')
        self.count = QSpinBox()
        self.count.setRange(2, 200)
        self.count.setValue(10)
        self.scale = QComboBox()
        self.scale.addItem('Linear', 'linear')
        self.scale.addItem('Logarithmic', 'log')
        for label, widget in [('Sweep', self.kind), ('Target', self.target),
                              ('Minimum', self.minimum), ('Maximum', self.maximum),
                              ('Samples', self.count), ('Spacing', self.scale)]:
            form.addRow(label, widget)

    def axis(self):
        if self.target.currentIndex() < 0:
            raise ValueError('Refresh the model and choose both axis targets')
        return SweepAxis(self.kind.currentData(), self.target.currentData(),
                         float(self.minimum.text()), float(self.maximum.text()),
                         self.count.value(), self.scale.currentData(), self.target.currentText())


class ExperimentsPanel(QWidget):
    busy_changed = pyqtSignal(bool)
    status_changed = pyqtSignal(str)

    def __init__(self, request_factory, parent=None):
        super().__init__(parent)
        self.request_factory = request_factory
        self.worker = None
        self.result = None
        self.parameter_options, self.species_options = [], []
        self._saved_targets = {}
        layout = QVBoxLayout(self)
        hint = QLabel('Sweep two independent inputs and measure one response. Uses the current Grid, Dynamics, Model and Run duration.')
        hint.setObjectName('fieldHint')
        hint.setWordWrap(True)
        layout.addWidget(hint)
        self.refresh_btn = QPushButton('Refresh targets from current model')
        layout.addWidget(self.refresh_btn)
        self.x_axis, self.y_axis = AxisControls('X axis'), AxisControls('Y axis')
        layout.addWidget(self.x_axis)
        layout.addWidget(self.y_axis)
        response_box = QGroupBox('Response & sampling')
        form = QFormLayout(response_box)
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.observable = QComboBox()
        self.metric = QComboBox()
        for key, label in METRICS.items():
            self.metric.addItem(label, key)
        self.repeats = QSpinBox()
        self.repeats.setRange(1, 100)
        self.seed = QSpinBox()
        self.seed.setRange(0, 2**31 - 1)
        self.tail = QDoubleSpinBox()
        self.tail.setRange(1, 100)
        self.tail.setValue(25)
        self.tail.setSuffix(' %')
        for label, widget in [('Species', self.observable), ('Statistic', self.metric),
                              ('Replicates / pixel', self.repeats), ('Random seed', self.seed),
                              ('Tail window', self.tail)]:
            form.addRow(label, widget)
        layout.addWidget(response_box)
        self.estimate = QLabel()
        self.estimate.setWordWrap(True)
        layout.addWidget(self.estimate)
        self.run_btn = QPushButton('Run Experiment')
        self.run_btn.setObjectName('primaryButton')
        self.cancel_btn = QPushButton('Cancel Experiment')
        self.cancel_btn.setEnabled(False)
        layout.addWidget(self.run_btn)
        layout.addWidget(self.cancel_btn)
        exports = QHBoxLayout()
        self.pdf_btn, self.data_btn = QPushButton('Save Heatmap...'), QPushButton('Export Data...')
        self.pdf_btn.setEnabled(False)
        self.data_btn.setEnabled(False)
        exports.addWidget(self.pdf_btn)
        exports.addWidget(self.data_btn)
        layout.addLayout(exports)
        layout.addStretch(1)

        self.result_widget = QWidget()
        result_layout = QVBoxLayout(self.result_widget)
        self.summary = QLabel('Configure the X and Y axes, then run an experiment. Each pixel is an independent simulation.')
        self.summary.setWordWrap(True)
        result_layout.addWidget(self.summary)
        self.figure = Figure(layout='constrained')
        self.canvas = FigureCanvasQTAgg(self.figure)
        result_layout.addWidget(self.canvas, 1)
        self.progress_label = QLabel('Statistics include all living cells. Failed runs appear in gray; details are included in the CSV export.')
        self.progress_label.setWordWrap(True)
        result_layout.addWidget(self.progress_label)
        install_figure_menu(self.canvas, self, 'experiment.pdf', enabled=lambda: self.result is not None)
        self.refresh_btn.clicked.connect(self.refresh_targets)
        self.run_btn.clicked.connect(self.start)
        self.cancel_btn.clicked.connect(self.cancel)
        self.pdf_btn.clicked.connect(self.save_plot)
        self.data_btn.clicked.connect(self.save_data)
        self.x_axis.kind.currentIndexChanged.connect(lambda: self._populate_axis(self.x_axis))
        self.y_axis.kind.currentIndexChanged.connect(lambda: self._populate_axis(self.y_axis))
        for widget in (self.x_axis.count, self.y_axis.count, self.repeats):
            widget.valueChanged.connect(self._update_estimate)
        self.metric.currentIndexChanged.connect(self._update_tail_enabled)
        self._update_estimate()
        self._update_tail_enabled()

    def settings(self):
        settings = {}
        for key, axis in (('x', self.x_axis), ('y', self.y_axis)):
            settings[key] = dict(kind=axis.kind.currentData(),
                                 target=self._saved_targets.get(key, axis.target.currentData()),
                                 minimum=axis.minimum.text(), maximum=axis.maximum.text(),
                                 count=axis.count.value(), scale=axis.scale.currentData())
        settings.update(observable=self._saved_targets.get('observable', self.observable.currentData()),
                        metric=self.metric.currentData(), repeats=self.repeats.value(),
                        seed=self.seed.value(), tail_percent=self.tail.value())
        return settings

    def restore_settings(self, settings):
        if not isinstance(settings, dict):
            return
        try:
            for key, axis in (('x', self.x_axis), ('y', self.y_axis)):
                saved = settings.get(key, {})
                axis.kind.setCurrentIndex(max(0, axis.kind.findData(saved.get('kind', 'parameter'))))
                axis.minimum.setText(str(saved.get('minimum', '0')))
                axis.maximum.setText(str(saved.get('maximum', '1')))
                axis.count.setValue(int(saved.get('count', 10)))
                axis.scale.setCurrentIndex(max(0, axis.scale.findData(saved.get('scale', 'linear'))))
                self._saved_targets[key] = saved.get('target')
            self._saved_targets['observable'] = settings.get('observable')
            self.metric.setCurrentIndex(max(0, self.metric.findData(settings.get('metric', 'final_mean'))))
            self.repeats.setValue(int(settings.get('repeats', 1)))
            self.seed.setValue(int(settings.get('seed', 0)))
            self.tail.setValue(float(settings.get('tail_percent', 25)))
        except (ValueError, TypeError, AttributeError, OverflowError):
            self.status_changed.emit('Some saved experiment settings are invalid; review the controls.')

    def _update_estimate(self):
        pixels = self.x_axis.count.value() * self.y_axis.count.value()
        self.estimate.setText(f'{pixels:,} pixels · {pixels * self.repeats.value():,} simulations. Start with a coarse scan to estimate runtime.')

    def _update_tail_enabled(self):
        self.tail.setEnabled(self.worker is None and self.metric.currentData().startswith('tail_'))

    def _populate_axis(self, axis):
        previous = axis.target.currentData()
        axis.target.clear()
        options = self.parameter_options if axis.kind.currentData() == 'parameter' else self.species_options
        for label, target in options:
            axis.target.addItem(label, target)
        index = axis.target.findData(previous)
        if index >= 0:
            axis.target.setCurrentIndex(index)

    def refresh_targets(self):
        try:
            request = self.request_factory()
            state = np.random.get_state()
            try:
                session = SimulationSession(request)
            finally:
                np.random.set_state(state)
            self.parameter_options = [
                (str(key), key) for key, value in session.params.items()
                if isinstance(value, (int, float)) and not isinstance(value, bool)
                and key != 'n_var' and np.isfinite(value)
            ]
            self.parameter_options.sort(key=lambda option: option[1] in ('dT', 'sigma'))
            labels = session.params.get('labels', [])
            self.species_options = [
                (f'{labels[i] if i < len(labels) else "Y"} [Y[{i}]]', i)
                for i in range(session.init_y.shape[1])
            ]
            for axis in (self.x_axis, self.y_axis):
                if not self.parameter_options:
                    axis.kind.setCurrentIndex(1)
                self._populate_axis(axis)
            for key, axis in (('x', self.x_axis), ('y', self.y_axis)):
                if key in self._saved_targets:
                    index = axis.target.findData(self._saved_targets[key])
                    if index >= 0:
                        axis.target.setCurrentIndex(index)
            if (self.x_axis.kind.currentData(), self.x_axis.target.currentData()) == (self.y_axis.kind.currentData(), self.y_axis.target.currentData()):
                if self.y_axis.target.count() > 1:
                    self.y_axis.target.setCurrentIndex(1)
                elif self.parameter_options:
                    self.y_axis.kind.setCurrentIndex(1)
            previous = self._saved_targets.get('observable', self.observable.currentData())
            self.observable.clear()
            for label, index in self.species_options:
                self.observable.addItem(label, index)
            self.observable.setCurrentIndex(max(0, self.observable.findData(previous)))
            self._saved_targets.clear()
            self.status_changed.emit('Experiment targets refreshed from the current model.')
        except Exception as error:
            self.status_changed.emit(f'Unable to refresh experiment targets: {error}')

    def start(self):
        if self.worker is not None:
            return
        try:
            request = ExperimentRequest(
                simulation=self.request_factory(), x=self.x_axis.axis(), y=self.y_axis.axis(),
                observable=self.observable.currentData(), observable_label=self.observable.currentText(),
                metric=self.metric.currentData(), repeats=self.repeats.value(), seed=self.seed.value(),
                tail_fraction=self.tail.value() / 100,
            )
            self.worker = ExperimentWorker(request, self)
        except Exception as error:
            self.status_changed.emit(f'Experiment configuration error: {error}')
            return
        self.worker.progress_updated.connect(self._progress)
        self.worker.experiment_finished.connect(self._finished)
        self.worker.error_occurred.connect(self._error)
        self.worker.cancelled.connect(lambda: self.status_changed.emit('Experiment cancelled. Previous heatmap retained.'))
        self.worker.finished.connect(self._idle)
        self._set_busy(True)
        self.busy_changed.emit(True)
        self.worker.start()

    def _set_busy(self, busy):
        for widget in (self.refresh_btn, self.x_axis, self.y_axis, self.observable,
                       self.metric, self.repeats, self.seed, self.tail, self.run_btn):
            widget.setEnabled(not busy)
        self.cancel_btn.setEnabled(busy)
        self.pdf_btn.setEnabled(not busy and self.result is not None)
        self.data_btn.setEnabled(not busy and self.result is not None)

    def _progress(self, value, message):
        self.progress_label.setText(f'{value}% · {message}')

    def cancel(self):
        if self.worker is not None:
            self.worker.cancel()
            self.cancel_btn.setEnabled(False)
            self.progress_label.setText('Cancelling experiment...')

    def _idle(self):
        old_worker = self.worker
        self.worker = None
        self._set_busy(False)
        self._update_tail_enabled()
        self.busy_changed.emit(False)
        old_worker.deleteLater()

    def _error(self, message):
        self.status_changed.emit(f'Experiment failed: {message}')
        self.progress_label.setText(message)

    def _finished(self, result):
        self.result = result
        request = result.request
        draw_heatmap(self.figure, result)
        failed = sum(bool(record['error']) for record in result.records)
        self.summary.setText(f'{len(result.x_values)} × {len(result.y_values)} pixels · {request.repeats} replicate(s) / pixel · seed {request.seed} · {failed} failed run(s)')
        self.progress_label.setText('Complete. Right-click to copy the heatmap. Export Data includes replicate values, errors and model settings.'
                                    + (f' {len(result.warnings)} model warning(s); see exported settings.' if result.warnings else ''))
        self.canvas.draw_idle()
        self.status_changed.emit(f'Experiment complete: {failed} failed run(s).')

    def save_plot(self):
        if self.result is None:
            return
        try:
            save_figure_dialog(self, self.figure, 'experiment.pdf')
        except Exception as error:
            self.status_changed.emit(f'Heatmap export failed: {error}')

    def save_data(self):
        if self.result is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, 'Export experiment data', 'experiment.csv', 'CSV Files (*.csv)')
        if not path:
            return
        if not Path(path).suffix:
            path += '.csv'
        try:
            self.export_data(path)
            self.status_changed.emit(f'Experiment data and settings saved: {path}')
        except Exception as error:
            self.status_changed.emit(f'Experiment data export failed: {error}')

    def export_data(self, path):
        export_experiment_data(self.result, path)
