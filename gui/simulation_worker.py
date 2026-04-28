import traceback
from PyQt6.QtCore import QThread, pyqtSignal

class SimulationWorker(QThread):
    progress_updated = pyqtSignal(int)
    simulation_finished = pyqtSignal(object, object, object) # history, cell_positions, grid
    error_occurred = pyqtSignal(str)
    warning_occurred = pyqtSignal(str)

    def __init__(self, model, T, proliferation_steps, apoptosis_steps, proliferation_n, apoptosis_n, proliferation_mode, apoptosis_mode):
        super().__init__()
        self.model = model
        self.T = T
        self.proliferation_steps = proliferation_steps
        self.apoptosis_steps = apoptosis_steps
        self.proliferation_n = proliferation_n
        self.apoptosis_n = apoptosis_n
        self.proliferation_mode = proliferation_mode
        self.apoptosis_mode = apoptosis_mode
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

            history, cell_positions = self.model.simulate(
                T=self.T,
                proliferation_steps=self.proliferation_steps,
                apoptosis_steps=self.apoptosis_steps,
                proliferation_n=self.proliferation_n,
                apoptosis_n=self.apoptosis_n,
                proliferation_mode=self.proliferation_mode,
                apoptosis_mode=self.apoptosis_mode,
                progress_callback=self.progress_callback,
                check_cancel=self.check_cancel,
                warning_callback=warning_cb
            )
            if self._is_cancelled:
                self.error_occurred.emit("Simulation cancelled.")
            else:
                self.progress_updated.emit(100)
                self.simulation_finished.emit(history, cell_positions, self.model.vor_grid)
        except Exception as e:
            error_msg = traceback.format_exc()
            self.error_occurred.emit(error_msg)
