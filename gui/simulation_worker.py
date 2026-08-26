import traceback
from PyQt6.QtCore import QThread, pyqtSignal

class SimulationWorker(QThread):
    progress_updated = pyqtSignal(int)
    simulation_finished = pyqtSignal(object)
    error_occurred = pyqtSignal(str)
    warning_occurred = pyqtSignal(str)

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self._is_cancelled = False

    def check_cancel(self):
        return self._is_cancelled

    def cancel(self):
        self._is_cancelled = True

    def progress_callback(self, current_step, total_steps):
        # Throttle progress updates a bit to not oversaturate the GUI thread
        if current_step % max(1, total_steps // 100) == 0:
            percentage = int((current_step / total_steps) * 100)
            self.progress_updated.emit(percentage)

    def run(self):
        try:
            def warning_cb(msg):
                self.warning_occurred.emit(msg)

            result = self.session.run(
                progress_callback=self.progress_callback,
                check_cancel=self.check_cancel,
                warning_callback=warning_cb
            )
            if self._is_cancelled:
                self.error_occurred.emit("Simulation cancelled.")
            else:
                self.progress_updated.emit(100)
                self.simulation_finished.emit(result)
        except Exception:
            error_msg = traceback.format_exc()
            self.error_occurred.emit(error_msg)
