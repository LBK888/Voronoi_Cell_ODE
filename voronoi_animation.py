import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
import seaborn as sns
import pandas as pd
from scipy.spatial import voronoi_plot_2d, Voronoi
import numpy as np
import os


class ExportCancelled(Exception):
    """Raised when a background video export is cancelled."""

    pass

class VoronoiAnimator:
    def __init__(self, vor_grid, sim_history, color_func=None, cell_positions_history=None, show_ticks=False, dynamic_range=True, params=None, color_func_str=None):
        """
        vor_grid: VoronoiGrid 物件
        sim_history: 模擬結果 (time, cell, var)
        color_func: 根據狀態決定顏色的函數 (Y) -> color list
        cell_positions_history: 細胞位置歷史記錄 (time, cell, 2)
        show_ticks: 是否顯示x,y軸數值
        dynamic_range: 是否動態調整x,y軸範圍
        """
        self.vor_grid = vor_grid
        self.sim_history = sim_history
        self.color_func = color_func
        self.cell_positions_history = cell_positions_history
        self.show_ticks = show_ticks
        self.dynamic_range = dynamic_range
        self.params = params
        self.color_func_str = color_func_str
        self.xlim = None
        self.ylim = None
        self._static_cells = None
        self._static_voronoi = None
        if (self.cell_positions_history is not None
                and len(self.cell_positions_history) > 0):
            try:
                first_cells = np.asarray(self.cell_positions_history[0])
                if all(
                    np.array_equal(first_cells, np.asarray(frame_cells))
                    for frame_cells in self.cell_positions_history[1:]
                ):
                    self._static_cells = first_cells
                    if np.array_equal(first_cells, np.asarray(self.vor_grid.cells)):
                        self._static_voronoi = self.vor_grid.vor
                    else:
                        self._static_voronoi = Voronoi(first_cells)
            except (TypeError, ValueError):
                pass
        self._color_has_mode = False
        self._color_accepts_bounds = False
        if self.color_func:
            try:
                import inspect
                signature = inspect.signature(self.color_func)
                self._color_has_mode = 'mode' in signature.parameters
                self._color_accepts_bounds = (
                    'global_bounds' in signature.parameters
                    or any(p.kind == p.VAR_KEYWORD
                           for p in signature.parameters.values())
                )
            except (TypeError, ValueError):
                pass
        
        self.frame_bounds = [
            self._calculate_bounds(frame) for frame in self.sim_history
        ]
        # Retain the attribute for callers that inspect it, while rendering now
        # uses frame_bounds so the scale follows the currently displayed step.
        self.global_bounds = self._calculate_bounds_across_history()

    @staticmethod
    def _calculate_bounds(frame):
        values = np.asarray(frame)
        if values.ndim != 2:
            return []
        bounds = []
        for variable_index in range(values.shape[1]):
            column = np.asarray(values[:, variable_index], dtype=float)
            finite_values = column[np.isfinite(column)]
            if finite_values.size:
                bounds.append((float(finite_values.min()), float(finite_values.max())))
            else:
                bounds.append((0.0, 1.0))
        return bounds

    def _calculate_bounds_across_history(self):
        if not self.frame_bounds:
            return []
        variable_count = max((len(bounds) for bounds in self.frame_bounds), default=0)
        global_bounds = []
        for variable_index in range(variable_count):
            available = [
                bounds[variable_index]
                for bounds in self.frame_bounds
                if variable_index < len(bounds)
            ]
            if available:
                global_bounds.append((
                    min(lower for lower, _ in available),
                    max(upper for _, upper in available),
                ))
            else:
                global_bounds.append((0.0, 1.0))
        return global_bounds

    def _draw_frame(self, ax, frame, save_path=None, line_colors='black', line_width=2, line_alpha=0.7, point_size=10):
        sim_history = self.sim_history
        color_func = self.color_func
        cell_positions_history = self.cell_positions_history
        if self._static_voronoi is not None:
            cells = self._static_cells
            vor = self._static_voronoi
        elif cell_positions_history is not None:
            cells = cell_positions_history[frame]
            vor = Voronoi(cells)
        else:
            cells = self.vor_grid.cells
            vor = self.vor_grid.vor
        # 計算多邊形顏色
        poly_colors = None
        center_colors = None
        off_center_colors = None
        off_center2_colors = None
        membrane_colors = None
        if color_func:
            kwargs = {}
            if self._color_accepts_bounds:
                kwargs['global_bounds'] = self.frame_bounds[frame]
                
            if self._color_has_mode:
                try:
                    poly_colors = color_func(sim_history[frame], mode='polygon', **kwargs)
                except Exception as e:
                    poly_colors = None
                try:
                    center_colors = color_func(sim_history[frame], mode='center', **kwargs)
                except Exception:
                    center_colors = None
                try:
                    off_center_colors = color_func(sim_history[frame], mode='off_center', **kwargs)
                except Exception:
                    off_center_colors = None
                try:
                    off_center2_colors = color_func(sim_history[frame], mode='off_center2', **kwargs)
                except Exception:
                    off_center2_colors = None
                try:
                    membrane_colors = color_func(sim_history[frame], mode='membrane', **kwargs)
                except Exception:
                    membrane_colors = None
            else:
                try:
                    poly_colors = color_func(sim_history[frame], **kwargs)
                except Exception as e:
                    poly_colors = None
        # 畫每個細胞的 Voronoi polygon
        if poly_colors is not None:
            for i, region_index in enumerate(vor.point_region):
                region = vor.regions[region_index]
                if not -1 in region and len(region) > 0:
                    polygon = [vor.vertices[v] for v in region]  
                    edge_rgba = mcolors.to_rgba(line_colors, alpha=line_alpha)
                    poly_patch = Polygon(polygon, closed=True, facecolor=poly_colors[i], edgecolor=edge_rgba, linewidth=line_width)
                    ax.add_patch(poly_patch)
        
        voronoi_plot_2d(vor, ax=ax, show_vertices=False, line_colors=line_colors, line_width=line_width, line_alpha=line_alpha, point_size=0)
        
        # 畫中心點旁邊的蛋白質
        if off_center_colors is not None:
            ax.scatter(cells[:,0]+0.15, cells[:,1]-0.15, c=off_center_colors, s=point_size)
        if off_center2_colors is not None:
            ax.scatter(cells[:,0]-0.15, cells[:,1]+0.15, c=off_center2_colors, s=point_size)
        # 畫細胞中心點(細胞核)
        if center_colors is not None:
            ax.scatter(cells[:,0], cells[:,1], c=center_colors, s=point_size)

        # 畫細胞膜
        if membrane_colors is not None:
            line_rgba = mcolors.to_rgba(line_colors)
            for i, region_index in enumerate(vor.point_region):
                if np.allclose(membrane_colors[i], line_rgba):
                    continue
                region = vor.regions[region_index]
                if not -1 in region and len(region) > 0:
                    polygon = [vor.vertices[v] for v in region]
                    poly_patch = Polygon(
                        polygon,
                        closed=True,
                        edgecolor=membrane_colors[i],
                        linewidth=line_width+1,
                        alpha=1,
                        facecolor='none'  # 確保無填色
                    )
                    ax.add_patch(poly_patch)
        
        ax.set_title(f"Cell state animation, frame {frame}")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_aspect('equal')
        if not self.show_ticks:
            ax.set_xticks([])
            ax.set_yticks([])
        # Keep a fixed export viewport after the first frame. Recomputing outer
        # cell sets and ridge neighbours for every frame added avoidable work.
        if not self.dynamic_range and self.xlim and self.ylim:
            ax.set_xlim(self.xlim[0], self.xlim[1])
            ax.set_ylim(self.ylim[0], self.ylim[1])
        else:
            outer_set = set()
            for i, region_index in enumerate(vor.point_region):
                region = vor.regions[region_index]
                if not region or any(v < 0 for v in region):
                    outer_set.add(i)
            neighbor_set = set()
            for p1, p2 in vor.ridge_points:
                if p1 in outer_set:
                    neighbor_set.add(p2)
                if p2 in outer_set:
                    neighbor_set.add(p1)
            all_indices = set(range(len(vor.points)))
            inner_indices = list(all_indices - outer_set - neighbor_set)
            if not inner_indices:
                inner_indices = list(range(len(cells)))
            xlim = [cells[inner_indices, 0].min(), cells[inner_indices, 0].max()]
            ylim = [cells[inner_indices, 1].min(), cells[inner_indices, 1].max()]
            if self.dynamic_range:
                ax.set_xlim(xlim[0], xlim[1])
                ax.set_ylim(ylim[0], ylim[1])
            else:
                self.xlim = xlim
                self.ylim = ylim
                ax.set_xlim(self.xlim[0], self.xlim[1])
                ax.set_ylim(self.ylim[0], self.ylim[1])
        # 每50幀儲存一張png
        # Video export deliberately does not emit high-resolution key-frame
        # PNGs. The old behaviour rendered a second 480-DPI image every 50
        # frames and was a major source of latency and disk usage.

    def animate(self, save_path=None, interval=200, line_colors='black', line_width=2, line_alpha=0.7, point_size=10):
        if save_path:
            fps = max(1, round(1000 / max(interval, 1)))
            return self.export_video(
                save_path,
                fps=min(fps, 30),
                dpi=120,
                line_colors=line_colors,
                line_width=line_width,
                line_alpha=line_alpha,
                point_size=point_size,
            )
        fig, ax = plt.subplots()
        colorbars = self._create_colorbars(fig, ax)
        def animate_func(frame):
            ax.clear()
            self._draw_frame(ax, frame, save_path, line_colors, line_width, line_alpha, point_size)
            self._update_colorbars(colorbars, frame)
        from matplotlib.animation import FuncAnimation
        ani = FuncAnimation(fig, animate_func, frames=len(self.sim_history), interval=interval)
        plt.show()

    def export_video(self, save_path, fps=20, dpi=120,
                     progress_callback=None, check_cancel=None,
                     **draw_options):
        """Render an MP4 at screen resolution with progress and cancellation."""
        from matplotlib.animation import FFMpegWriter
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure

        figure = Figure(figsize=(8, 6), dpi=dpi, tight_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.add_subplot(111)
        colorbars = self._create_colorbars(figure, axis)
        writer = FFMpegWriter(
            fps=fps,
            bitrate=2400,
            metadata={"title": "Voronoi simulation"},
        )
        total = len(self.sim_history)
        try:
            with writer.saving(figure, save_path, dpi=dpi):
                for frame in range(total):
                    if check_cancel and check_cancel():
                        raise ExportCancelled()
                    axis.clear()
                    self._draw_frame(axis, frame, None, **draw_options)
                    self._update_colorbars(colorbars, frame)
                    writer.grab_frame()
                    if progress_callback and (
                        frame == total - 1
                        or frame % max(1, total // 100) == 0
                    ):
                        progress_callback(frame + 1, total)
        finally:
            figure.clear()

    def animateGUI(self, fig, interval=200, line_colors='black', line_width=2, line_alpha=0.7, point_size=10):
        ax = fig.gca()
        colorbar_entries = self._create_colorbars(fig, ax)

        def animate_func(frame):
            ax.clear()
            self._draw_frame(ax, frame, None, line_colors, line_width, line_alpha, point_size)
            self._update_colorbars(colorbar_entries, frame)
        from matplotlib.animation import FuncAnimation
        anim = FuncAnimation(fig, animate_func, frames=len(self.sim_history), interval=interval, blit=False, repeat=False)
        
        return anim

    def _create_colorbars(self, figure, axis):
        """Build the same legends for GUI, standalone animation and MP4."""
        import inspect
        import re
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        if not self.color_func or not self.frame_bounds:
            return []
        params = self.params or {}
        labels = params.get('labels', [])
        specs = params.get('colorbars')
        if specs is None:
            try:
                source = self.color_func_str or inspect.getsource(self.color_func)
            except (OSError, TypeError):
                source = ''
            cmaps = list(dict.fromkeys(re.findall(r'plt\.cm\.([a-zA-Z0-9_]+)', source)))
            specs = [dict(variable=i, cmap=cmap) for i, cmap in enumerate(cmaps)
                     if i < len(self.frame_bounds[0])]
        entries = []
        divider = make_axes_locatable(axis)
        for spec in specs:
            index = int(spec['variable'])
            if index < 0 or index >= len(self.frame_bounds[0]):
                raise ValueError('colorbars variable is outside the state vector')
            label = spec.get('label', labels[index] if index < len(labels) else f'Y[{index}]')
            mappable = plt.cm.ScalarMappable(cmap=spec['cmap'])
            mappable.set_array([])
            cax = divider.append_axes('right', size='5%', pad=0.8)
            colorbar = figure.colorbar(mappable, cax=cax)
            colorbar.set_label(label)
            entries.append((index, mappable, colorbar))
        self._update_colorbars(entries, 0)
        return entries

    def _update_colorbars(self, entries, frame):
        for index, mappable, colorbar in entries:
            vmin, vmax = self.frame_bounds[frame][index]
            if vmin == vmax:
                vmax = vmin + max(1e-8, abs(vmin) * 1e-8)
            # Replace the norm atomically, avoiding callbacks with transient
            # invalid limits when consecutive frames differ by orders of magnitude.
            mappable.set_norm(plt.Normalize(vmin=vmin, vmax=vmax))
            colorbar.update_normal(mappable)

    @staticmethod
    def plot_heatmap(results, title, save_path=None):
        fig = plt.figure(figsize=(10, 10))
        ax1 = fig.add_subplot(111)
        df1 = pd.DataFrame(results)
        sns.heatmap(df1, annot=False, fmt=".3f", linewidths=.5, ax=ax1, cmap='RdBu')
        ax1.set(xlabel="", ylabel="")
        ax1.xaxis.tick_top()
        ax1.set_title(title)
        if save_path:
            plt.savefig(save_path)
        plt.show()

    @staticmethod
    def plot_concentration_over_time(history, labels=None, save_path=None, cell_indices=None):
        """
        支援 history 為 object array（細胞數可變），可指定要畫哪些 cell index，
        若未指定則自動選存活最久的細胞（出現在所有幀的 index=0）。
        history: (time, cell, var) 或 object array
        cell_indices: list[int]，要畫的細胞索引（以每幀的陣列 index 為準）
        """
        arr = history
        # 預設只畫每一幀的第0個細胞（存活最久）
        if cell_indices is None:
            cell_indices = [0]
        # 準備每個 cell 的濃度歷程
        cell_traces = []
        for idx in cell_indices:
            trace = [frame[idx] for frame in arr if len(frame) > idx]
            cell_traces.append(np.array(trace))
        # 決定變數數量
        n_var = cell_traces[0].shape[1] if len(cell_traces[0].shape) > 1 else 1
        if labels is None:
            labels = [f"Var{i+1}" for i in range(n_var)]
        # 畫圖
        n_cells = len(cell_indices)
        if n_cells == 1:
            plt.figure(figsize=(12,8))
            for i in range(n_var):
                plt.plot(cell_traces[0][:,i], label=labels[i])
            plt.xlabel("Time step")
            plt.ylabel("Concentration")
            plt.legend()
            plt.title(f"Concentration over time (Cell {cell_indices[0]})")
            if save_path:
                plt.savefig(save_path)
            plt.show()
        else:
            fig, axes = plt.subplots(n_cells, 1, figsize=(12, 4*n_cells), sharex=True)
            if n_cells == 1:
                axes = [axes]
            for c, trace in enumerate(cell_traces):
                ax = axes[c]
                for i in range(n_var):
                    ax.plot(trace[:,i], label=labels[i])
                ax.set_ylabel("Concentration")
                ax.set_title(f"Cell {cell_indices[c]}")
                ax.legend()
            axes[-1].set_xlabel("Time step")
            fig.suptitle("Concentration over time (selected cells)")
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            if save_path:
                plt.savefig(save_path)
            plt.show() 

