from dataclasses import dataclass
import json
import os
import traceback

import h5py
import numpy as np
from openpyxl import Workbook
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from PyQt6.QtCore import QThread, pyqtSignal

from voronoi_animation import ExportCancelled, VoronoiAnimator
from gui.simulation_session import SimulationResult


@dataclass(frozen=True)
class ExportRequest:
    folder: str
    prefix: str
    config: dict
    simulation: SimulationResult
    color_func: object
    color_func_source: str


class ResultsExportWorker(QThread):
    """Export simulation artefacts without blocking the Qt event loop."""

    progress_updated = pyqtSignal(int, str)
    export_finished = pyqtSignal(str)
    error_occurred = pyqtSignal(str)
    cancelled = pyqtSignal()

    def __init__(self, request, parent=None):
        super().__init__(parent)
        if not isinstance(request, ExportRequest):
            raise TypeError('request must be an ExportRequest')
        self.request = request
        self.folder = request.folder
        self.prefix = request.prefix
        self.config = request.config
        self.grid = request.simulation.grid
        self.sim_history = request.simulation.history
        self.cell_positions_history = request.simulation.cell_positions
        self.cell_ids_history = request.simulation.cell_ids
        self.params = request.simulation.params or {}
        self.color_func = request.color_func
        self.color_func_str = request.color_func_source
        self._is_cancelled = False

    def cancel(self):
        self._is_cancelled = True

    def _check_cancel(self):
        return self._is_cancelled

    def _emit(self, value, message):
        self.progress_updated.emit(int(value), message)

    def _raise_if_cancelled(self):
        if self._is_cancelled:
            raise ExportCancelled()

    def run(self):
        try:
            os.makedirs(self.folder, exist_ok=True)
            self._save_config()
            self._emit(5, "Rendering video in background...")
            self._save_video()
            self._raise_if_cancelled()
            self._emit(72, "Saving concentration plot...")
            self._save_concentration_plot()
            self._raise_if_cancelled()
            self._emit(78, "Writing Excel and HDF5 data...")
            self._save_tabular_data()
            self._raise_if_cancelled()
            self._emit(100, "Export complete.")
            self.export_finished.emit(self.folder)
        except ExportCancelled:
            self.cancelled.emit()
        except Exception:
            self.error_occurred.emit(traceback.format_exc())

    def _save_config(self):
        path = os.path.join(self.folder, f"{self.prefix}config.json")
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(self.config, handle, ensure_ascii=False, indent=2)

    def _save_video(self):
        animator = VoronoiAnimator(
            self.grid, self.sim_history, self.color_func,
            self.cell_positions_history, show_ticks=False,
            dynamic_range=False, params=self.params,
            color_func_str=self.color_func_str,
        )
        video_path = os.path.join(self.folder, f"{self.prefix}simulation.mp4")

        def video_progress(current, total):
            progress = 5 + (current / max(total, 1)) * 65
            self._emit(progress, f"Rendering video frame {current}/{total}")

        animator.export_video(
            video_path, fps=20, dpi=120,
            progress_callback=video_progress,
            check_cancel=self._check_cancel,
            line_colors="black", line_width=1,
            line_alpha=0.5, point_size=15,
        )

    def _save_concentration_plot(self):
        labels = self.params.get(
            "labels", [f"Y[{i}]" for i in range(self.sim_history[0].shape[1])]
        )
        trace = np.asarray([frame[0] for frame in self.sim_history if len(frame) > 0])
        figure = Figure(figsize=(10, 6), tight_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.add_subplot(111)
        for index in range(trace.shape[1]):
            label = labels[index] if index < len(labels) else f"Y[{index}]"
            axis.plot(trace[:, index], label=label)
        axis.set_xlabel("Time step")
        axis.set_ylabel("Concentration")
        axis.set_title("Concentration over time (Cell 0)")
        axis.legend()
        figure.savefig(os.path.join(
            self.folder, f"{self.prefix}concentration.pdf"
        ), format="pdf")
        figure.clear()

    def _save_tabular_data(self):
        n_variables = self.sim_history[0].shape[1] if self.sim_history else 0
        headers = ["cell_id", "T", "step", "t", "x", "y"] + [
            f"Y[{i}]" for i in range(n_variables)
        ]
        d_t = float(self.params.get("dT", 1.0) or 1.0)
        max_step = max(1, round(1.0 / d_t))
        total_frames = min(len(self.sim_history), len(self.cell_positions_history))
        has_cell_ids = len(self.cell_ids_history) >= total_frames

        def frame_row_count(index):
            counts = [
                len(self.sim_history[index]),
                len(self.cell_positions_history[index]),
            ]
            if has_cell_ids:
                counts.append(len(self.cell_ids_history[index]))
            return min(counts)

        total_rows = sum(frame_row_count(index) for index in range(total_frames))
        excel_path = os.path.join(self.folder, f"{self.prefix}cells.xlsx")
        h5_path = os.path.join(self.folder, f"{self.prefix}cells.h5")
        chunk_rows = max(1, min(4096, total_rows))
        workbook = Workbook(write_only=True)
        worksheet = workbook.create_sheet(title="cells")
        worksheet.append(headers)

        with h5py.File(h5_path, "w") as h5_handle:
            dataset_options = {}
            if total_rows:
                dataset_options = {
                    "chunks": (chunk_rows, len(headers)),
                    "compression": "gzip",
                    "compression_opts": 1,
                }
            dataset = h5_handle.create_dataset(
                "cells", shape=(total_rows, len(headers)), dtype="f8",
                **dataset_options,
            )
            dataset.attrs["columns"] = np.asarray(
                headers, dtype=h5py.string_dtype("utf-8")
            )

            offset = 0
            for frame_index in range(total_frames):
                self._raise_if_cancelled()
                values = np.asarray(self.sim_history[frame_index])
                positions = np.asarray(self.cell_positions_history[frame_index])
                count = frame_row_count(frame_index)
                if count == 0:
                    continue
                if has_cell_ids:
                    cell_ids = np.asarray(
                        self.cell_ids_history[frame_index], dtype=int
                    )[:count]
                else:
                    cell_ids = np.arange(count, dtype=int)
                time_value, step_value = divmod(frame_index, max_step)
                t_value = time_value * max_step + step_value
                block = np.column_stack((
                    cell_ids,
                    np.full(count, time_value), np.full(count, step_value),
                    np.full(count, t_value),
                    positions[:count, :2], values[:count, :n_variables],
                ))
                for row in block.tolist():
                    worksheet.append(row)
                dataset[offset:offset + count] = block
                offset += count
                if frame_index % max(1, total_frames // 20) == 0:
                    progress = 78 + (frame_index / max(total_frames, 1)) * 21
                    self._emit(progress, "Writing simulation data...")
        workbook.save(excel_path)
