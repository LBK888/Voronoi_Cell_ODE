import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
import seaborn as sns
import pandas as pd
from scipy.spatial import voronoi_plot_2d, Voronoi
import numpy as np
import os

class VoronoiAnimator:
    def __init__(self, vor_grid, sim_history, color_func=None, cell_positions_history=None, show_ticks=False, dynamic_range=True):
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
        self.xlim = None
        self.ylim = None

    def _draw_frame(self, ax, frame, save_path=None, line_colors='black', line_width=2, line_alpha=0.7, point_size=10):
        sim_history = self.sim_history
        color_func = self.color_func
        cell_positions_history = self.cell_positions_history
        if cell_positions_history is not None:
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
            if 'mode' in color_func.__code__.co_varnames:
                try:
                    poly_colors = color_func(sim_history[frame], mode='polygon')
                except Exception as e:
                    poly_colors = None
                try:
                    center_colors = color_func(sim_history[frame], mode='center')
                except Exception:
                    center_colors = None
                try:
                    off_center_colors = color_func(sim_history[frame], mode='off_center')
                except Exception:
                    off_center_colors = None
                try:
                    off_center2_colors = color_func(sim_history[frame], mode='off_center2')
                except Exception:
                    off_center2_colors = None
                try:
                    membrane_colors = color_func(sim_history[frame], mode='membrane')
                except Exception:
                    membrane_colors = None
            else:
                try:
                    poly_colors = color_func(sim_history[frame])
                except Exception as e:
                    poly_colors = None
        # 畫每個細胞的 Voronoi polygon
        if poly_colors is not None:
            for i, region_index in enumerate(vor.point_region):
                region = vor.regions[region_index]
                if not -1 in region and len(region) > 0:
                    polygon = [vor.vertices[v] for v in region]
                    poly_patch = Polygon(polygon, closed=True, facecolor=poly_colors[i], edgecolor=line_colors, linewidth=line_width, alpha=line_alpha)
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
        # 動態範圍
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
        if self.dynamic_range:
            ax.set_xlim(cells[inner_indices,0].min(), cells[inner_indices,0].max())
            ax.set_ylim(cells[inner_indices,1].min(), cells[inner_indices,1].max())
        else:
            if not self.xlim and not self.ylim:
                self.xlim=[cells[inner_indices,0].min(), cells[inner_indices,0].max()]
                self.ylim=[cells[inner_indices,1].min(), cells[inner_indices,1].max()]
            ax.set_xlim(self.xlim[0], self.xlim[1])
            ax.set_ylim(self.ylim[0], self.ylim[1])
        # 每50幀儲存一張png
        if save_path and frame % 50 == 0:
            base, ext = os.path.splitext(save_path)
            png_path = f"{base}_frame{frame}.png"
            plt.savefig(png_path, dpi=480)

    def animate(self, save_path=None, interval=200, line_colors='black', line_width=2, line_alpha=0.7, point_size=10):
        fig, ax = plt.subplots()
        def animate_func(frame):
            ax.clear()
            self._draw_frame(ax, frame, save_path, line_colors, line_width, line_alpha, point_size)
        from matplotlib.animation import FuncAnimation
        ani = FuncAnimation(fig, animate_func, frames=len(self.sim_history), interval=interval)
        if save_path:
            ani.save(save_path, dpi=480)
        else:
            plt.show()

    def animateGUI(self, fig, interval=200, line_colors='black', line_width=2, line_alpha=0.7, point_size=10):
        ax = fig.gca()
        def animate_func(frame):
            ax.clear()
            self._draw_frame(ax, frame, None, line_colors, line_width, line_alpha, point_size)
        from matplotlib.animation import FuncAnimation
        anim = FuncAnimation(fig, animate_func, frames=len(self.sim_history), interval=interval, blit=False, repeat=False)
        return anim

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

