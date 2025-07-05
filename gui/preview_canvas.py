from PyQt6.QtWidgets import QWidget, QVBoxLayout, QApplication
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from PyQt6.QtCore import pyqtSignal

class PreviewCanvas(QWidget):
    cell_count_changed = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        layout = QVBoxLayout(self)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self.cells = np.array([])
        self.selected_idx = None
        self.dragging = False

    def plot(self, plot_func, *args, **kwargs):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        plot_func(ax, *args, **kwargs)
        self.canvas.draw()

    def set_cells(self, cells):
        self.cells = np.array(cells)
        self.selected_idx = None
        self.dragging = False
        self.update_plot()

    def get_cells(self):
        return self.cells.copy()

    def update_plot(self):
        def plot_func(ax):
            from scipy.spatial import Voronoi, voronoi_plot_2d
            if len(self.cells) < 3:
                ax.clear()
                ax.plot(self.cells[:,0], self.cells[:,1], 'o', color='red')
                ax.set_aspect('equal')
                return
            vor = Voronoi(self.cells)
            voronoi_plot_2d(vor, ax=ax, show_vertices=False, line_colors='black')
            ax.plot(self.cells[:,0], self.cells[:,1], 'o', color='red')
            ax.set_aspect('equal')
        self.plot(plot_func)

    # 只用 mpl 事件
    def connect_events(self):
        self.canvas.mpl_connect('button_press_event', self.on_mpl_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_mpl_motion)
        self.canvas.mpl_connect('button_release_event', self.on_mpl_release)

    def on_mpl_click(self, event):
        if event.inaxes is None:
            return
        pos = np.array([event.xdata, event.ydata])
        idx = self._find_nearest_cell(pos)
        if event.button == 1:
            if idx is not None and np.linalg.norm(self.cells[idx] - pos) < 0.5:
                self.selected_idx = idx
                self.dragging = True
            else:
                self.cells = np.vstack([self.cells, pos])
                self.update_plot()
                self.cell_count_changed.emit(len(self.cells))
        elif event.button == 3:
            if idx is not None and len(self.cells) > 3:
                self.cells = np.delete(self.cells, idx, axis=0)
                self.update_plot()
                self.cell_count_changed.emit(len(self.cells))

    def on_mpl_motion(self, event):
        if self.dragging and self.selected_idx is not None and event.inaxes is not None:
            pos = np.array([event.xdata, event.ydata])
            self.cells[self.selected_idx] = pos
            self.update_plot()

    def on_mpl_release(self, event):
        self.dragging = False
        self.selected_idx = None

    def _find_nearest_cell(self, pos):
        if len(self.cells) == 0:
            return None
        dists = np.linalg.norm(self.cells - pos, axis=1)
        idx = np.argmin(dists)
        return idx

    def canvas_mouseEventCoords(self, event):
        # 將 Qt event 轉換為資料座標
        x = event.position().x()
        y = event.position().y()
        inv = self.canvas.figure.axes[0].transData.inverted()
        return np.array(inv.transform((x, y)))