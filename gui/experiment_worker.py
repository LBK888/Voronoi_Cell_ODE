"""Keep parameter sweeps off the GUI thread."""

import traceback
from PyQt6.QtCore import QThread, pyqtSignal
from gui.experiments import ExperimentCancelled, ExperimentRunner


class ExperimentWorker(QThread):
    progress_updated = pyqtSignal(int, str)
    experiment_finished = pyqtSignal(object)
    error_occurred = pyqtSignal(str)
    cancelled = pyqtSignal()

    def __init__(self, request, parent=None):
        super().__init__(parent)
        self.runner = ExperimentRunner(request)
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def run(self):
        try:
            result = self.runner.run(
                progress_callback=self.progress_updated.emit,
                check_cancel=lambda: self._cancelled,
            )
            self.experiment_finished.emit(result)
        except ExperimentCancelled:
            self.cancelled.emit()
        except Exception:
            self.error_occurred.emit(traceback.format_exc())
