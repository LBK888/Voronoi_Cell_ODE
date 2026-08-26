from collections import defaultdict
from dataclasses import dataclass

import numpy as np
from scipy.spatial import Voronoi, QhullError


@dataclass(frozen=True)
class MovementPlan:
    """Cell movement settings consumed by VoronoiGrid."""

    random_strength: float = 0.0
    away_strength: float = 0.0
    convergent_strength: float = 0.0
    repulsion_strength: float = 0.0
    vertex_strength: float = 0.0

class VoronoiGrid:
    def __init__(self, grid_shape=(2,1), cell_dist=1.0, pos_rand=0.0, mode='honeycomb', custom_cells=None, import_path=None):
        """
        grid_shape: (x, y) 格狀排列
        cell_dist: 細胞間距
        pos_rand: 位置隨機擾動幅度
        mode: 'honeycomb', 'random', 'regular', 'custom', 'import'
        custom_cells: np.ndarray, 若 mode='custom' 則用此
        import_path: str, 若 mode='import' 則從檔案讀取
        """
        self.grid_shape = grid_shape
        self.cell_dist = cell_dist
        self.pos_rand = pos_rand
        self.mode = mode
        self.custom_cells = custom_cells
        self.import_path = import_path
        self._safe_update_voronoi(self._init_cells())

    @classmethod
    def from_cells(cls, cells, cell_dist=1.0):
        """Create a layout from existing coordinates through the normal safety seam."""
        cells = np.asarray(cells, dtype=float)
        return cls(
            grid_shape=(len(cells), 1),
            cell_dist=cell_dist,
            mode='custom',
            custom_cells=cells,
        )

    def copy_cells(self):
        return self.cells.copy()

    def replace_cells(self, cells):
        cells = np.asarray(cells, dtype=float)
        if cells.ndim != 2 or cells.shape[1] != 2:
            raise ValueError('cells must have shape (n, 2)')
        if len(cells) < 3:
            raise ValueError('at least three cells are required for a Voronoi layout')
        self._safe_update_voronoi(cells)

    def add_cell(self, position):
        position = np.asarray(position, dtype=float)
        if position.shape != (2,):
            raise ValueError('position must have shape (2,)')
        self.replace_cells(np.vstack([self.cells, position]))
        return len(self.cells) - 1

    def move_cell(self, index, position):
        cells = self.copy_cells()
        cells[index] = np.asarray(position, dtype=float)
        self.replace_cells(cells)

    def remove_cell(self, index):
        if len(self.cells) <= 3:
            raise ValueError('a Voronoi layout must keep at least three cells')
        self.replace_cells(np.delete(self.cells, index, axis=0))

    def _safe_update_voronoi(self, new_cells, max_retries=3):
        import logging
        current_cells = new_cells.copy()
        for i in range(max_retries):
            try:
                self.vor = Voronoi(current_cells)
                self.cells = current_cells
                return
            except QhullError as e:
                logging.warning(f"QhullError encountered (attempt {i+1}/{max_retries}): {e}. Applying jitter.")
                current_cells += (np.random.rand(*current_cells.shape) - 0.5) * 1e-4
        
        # Last attempt, if it fails it will raise the error normally
        self.vor = Voronoi(current_cells)
        self.cells = current_cells

    def _init_cells(self):
        import os
        import numpy as np
        if self.mode == 'honeycomb':
            cells = []
            for x in range(self.grid_shape[0]):
                for y in range(self.grid_shape[1]):
                    if y % 2 == 0:
                        px = x + (1 + (np.random.random_sample() - 0.5) * self.pos_rand) * self.cell_dist
                        py = y + (1 + (np.random.random_sample() - 0.5) * self.pos_rand) * self.cell_dist
                    else:
                        px = x + (1 + (np.random.random_sample() - 0.5) * self.pos_rand) * self.cell_dist + 0.5 * self.cell_dist
                        py = y + (1 + (np.random.random_sample() - 0.5) * self.pos_rand) * self.cell_dist
                    cells.append([px, py])
            return np.array(cells)
        elif self.mode == 'random':
            n = self.grid_shape[0] * self.grid_shape[1]
            cells = np.random.rand(n, 2) * np.array([self.grid_shape[0], self.grid_shape[1]]) * self.cell_dist
            return cells
        elif self.mode == 'regular':
            xs = np.arange(self.grid_shape[0]) * self.cell_dist
            ys = np.arange(self.grid_shape[1]) * self.cell_dist
            grid = np.array(np.meshgrid(xs, ys)).T.reshape(-1, 2)
            return grid
        elif self.mode == 'custom':
            if self.custom_cells is not None:
                return np.array(self.custom_cells)
            else:
                raise ValueError('custom_cells must be provided for custom mode')
        elif self.mode == 'import':
            if self.import_path is not None:
                ext = os.path.splitext(self.import_path)[-1].lower()
                if ext == '.xlsx':
                    try:
                        from openpyxl import load_workbook
                    except ImportError:
                        raise ImportError('openpyxl is required to import xlsx files')
                    wb = load_workbook(self.import_path, read_only=True, data_only=True)
                    ws = wb.active
                    headers = [str(cell.value).strip().lower() if cell.value is not None else '' for cell in next(ws.iter_rows(min_row=1, max_row=1))]
                    # 優先找 T, dT, x, y
                    x_idx = None
                    y_idx = None
                    t_idx = None
                    step_idx = None
                    for idx, h in enumerate(headers):
                        if h == 'x':
                            x_idx = idx
                        elif h == 'y':
                            y_idx = idx
                        elif h == 't':
                            t_idx = idx
                        elif h == 'step':
                            step_idx = idx
                    if x_idx is not None and y_idx is not None:
                        # 若有 T/dT 欄，僅取 T=0/dT=0 的資料
                        rows = list(ws.iter_rows(min_row=2, values_only=True))
                        filtered = []
                        if t_idx is not None and step_idx is not None:
                            for row in rows:
                                if row[t_idx] == 0 and row[step_idx] == 0:      #the t=0 cell grid
                                    filtered.append([row[x_idx], row[y_idx]])
                        else:
                            # 沒有 T/dT 欄，全部取 x, y
                            filtered = [[row[x_idx], row[y_idx]] for row in rows if row[x_idx] is not None and row[y_idx] is not None]
                        if not filtered:
                            raise ValueError('No valid cell positions found in xlsx (T=0/step=0 or all x/y)')
                        return np.array(filtered)
                    else:
                        raise ValueError('xlsx 檔案需包含 x, y 欄位 (不分大小寫)')
                else:
                    # 預設 csv/txt
                    try:
                        return np.loadtxt(self.import_path, delimiter=',')
                    except Exception as e:
                        raise ValueError(f'Failed to load import file: {e}')
            else:
                raise ValueError('import_path must be provided for import mode')
        else:
            raise ValueError(f'Unknown mode: {self.mode}')

    def move_cells(self, rule=None, random_strength=0.0):
        """
        rule: function(cells) -> new_cells
        random_strength: 0~1, 隨機移動強度（相對 cell_dist）
        """
        new_cells = self.cells.copy()
        if rule:
            new_cells = rule(new_cells)
        if random_strength > 0:
            new_cells += (np.random.rand(*new_cells.shape)-0.5)*2*self.cell_dist*random_strength
        self._safe_update_voronoi(new_cells) 

    def finite_ridge_neighbors(self, max_ridge_length=2.0):
        """Return finite Voronoi neighbors weighted by their shared ridge length."""
        neighbors = defaultdict(dict)
        for (p1, p2), ridge_vertices in zip(self.vor.ridge_points, self.vor.ridge_vertices):
            if not all(vertex >= 0 for vertex in ridge_vertices):
                continue
            points = self.vor.vertices[ridge_vertices]
            length = np.linalg.norm(points[0] - points[1])
            if max_ridge_length is not None and length > max_ridge_length:
                continue
            neighbors[p1][p2] = length
            neighbors[p2][p1] = length
        return neighbors

    def apply_movement(self, plan):
        """Apply a complete movement plan while keeping topology updates local."""
        if not isinstance(plan, MovementPlan):
            raise TypeError('plan must be a MovementPlan')
        if plan.random_strength:
            self.move_cells(random_strength=plan.random_strength)
        if plan.away_strength:
            self.move_cells(lambda cells: self._move_away(cells, plan.away_strength))
        if plan.convergent_strength:
            self._safe_update_voronoi(self._convergent_extension(plan.convergent_strength))
        if plan.repulsion_strength:
            self._safe_update_voronoi(self._repel_neighbors(plan.repulsion_strength))
        if plan.vertex_strength:
            self._safe_update_voronoi(self._vertex_movement(plan.vertex_strength))

    def _outer_cell_indices(self):
        outer = set()
        for index, region_index in enumerate(self.vor.point_region):
            region = self.vor.regions[region_index]
            if not region or any(vertex < 0 for vertex in region):
                outer.add(index)
        return outer

    @staticmethod
    def _move_away(cells, strength):
        center = np.mean(cells, axis=0)
        vectors = cells - center
        norms = np.linalg.norm(vectors, axis=1, keepdims=True) + 1e-8
        return cells + (vectors / norms) * strength

    def _ridge_neighbor_sets(self):
        neighbors = defaultdict(set)
        for p1, p2 in self.vor.ridge_points:
            neighbors[p1].add(p2)
            neighbors[p2].add(p1)
        return neighbors

    def _convergent_extension(self, strength, max_force=1.0):
        cells = self.copy_cells()
        outer = self._outer_cell_indices()
        center_x, center_y = np.mean(cells, axis=0)
        new_cells = cells.copy()
        for index in range(len(cells)):
            if index in outer:
                continue
            force = np.array([
                np.sign(cells[index, 0] - center_x) * strength / 4,
                np.sign(center_y - cells[index, 1]) * strength,
            ])
            norm = np.linalg.norm(force)
            if norm > max_force:
                force = force / norm * max_force
            new_cells[index] += force
        return new_cells

    def _repel_neighbors(self, strength, min_dist=1.5, max_force=1.0):
        cells = self.copy_cells()
        neighbors = self._ridge_neighbor_sets()
        outer = self._outer_cell_indices()
        new_cells = cells.copy()
        for index in range(len(cells)):
            if index in outer:
                continue
            force = np.zeros(2)
            for neighbor in neighbors[index]:
                vector = cells[index] - cells[neighbor]
                distance = np.linalg.norm(vector)
                if 1e-8 < distance < min_dist:
                    force += (vector / distance) * (strength * (min_dist - distance))
            norm = np.linalg.norm(force)
            if norm > max_force:
                force = force / norm * max_force
            new_cells[index] += force
        return new_cells

    def _vertex_movement(self, adhesion_strength, repulsion_strength=0.1,
                         min_dist=1.5, max_force=1.0):
        cells = self.copy_cells()
        neighbors = self.finite_ridge_neighbors(max_ridge_length=None)
        outer = self._outer_cell_indices()
        new_cells = cells.copy()
        for index in range(len(cells)):
            if index in outer:
                continue
            force = np.zeros(2)
            for neighbor, edge_length in neighbors[index].items():
                vector = cells[neighbor] - cells[index]
                distance = np.linalg.norm(vector)
                if distance < 1e-8:
                    continue
                direction = vector / distance
                force += adhesion_strength * edge_length * direction
                if distance < min_dist:
                    force -= repulsion_strength * ((min_dist - distance) / distance) * direction
            norm = np.linalg.norm(force)
            if norm > max_force:
                force = force / norm * max_force
            new_cells[index] += force
        return new_cells

    def get_inner_cell_indices(self):
        """
        回傳內圈細胞的 index。
        1. 先排除 region 中帶有 -1 的細胞（外圈）
        2. 再排除與外圈細胞相鄰的細胞（利用 ridge_points）
        """
        vor = self.vor
        # 1. 找出所有外圈細胞（region 含 -1）
        outer_set = self._outer_cell_indices()
        # 2. 找出所有與外圈細胞相鄰的細胞
        neighbor_set = set()
        for p1, p2 in vor.ridge_points:
            if p1 in outer_set:
                neighbor_set.add(p2)
            if p2 in outer_set:
                neighbor_set.add(p1)
        # 3. 內圈細胞 = 全部細胞 - 外圈 - 鄰居
        all_indices = set(range(len(vor.points)))
        inner_indices = list(all_indices - outer_set - neighbor_set)
        return inner_indices

    def cell_proliferation(self, n=1, mode='area', concentrations=None, target_indices=None):
        from scipy.spatial import ConvexHull
        new_cells = []
        new_conc = []
        n_iters = len(target_indices) if target_indices is not None else n
        for iter_idx in range(n_iters):
            vor = self.vor  # 每次都用最新的
            inner_indices = self.get_inner_cell_indices()
            # 計算面積
            areas = []
            for i in inner_indices:
                region = vor.regions[vor.point_region[i]]
                polygon = [vor.vertices[v] for v in region]
                try:
                    hull = ConvexHull(polygon)
                    areas.append(hull.volume)
                except:
                    areas.append(0)
            if not areas:
                break
            if target_indices is not None:
                if target_indices[iter_idx] >= len(self.cells):
                    continue
                idx = target_indices[iter_idx]
            else:
                # 選最大面積
                if mode == 'area':
                    idx = inner_indices[np.argmax(areas)]
                elif mode == 'random':
                    idx = np.random.choice(inner_indices)
                else:
                    raise ValueError('mode 必須為 area 或 random')
            # 分裂
            region = vor.regions[vor.point_region[idx]]
            polygon = np.array([vor.vertices[v] for v in region])
            center = self.cells[idx]
            cov = np.cov(polygon - center, rowvar=False)
            eigvals, eigvecs = np.linalg.eigh(cov)
            axis = eigvecs[:, np.argmax(eigvals)]
            offset = axis * self.cell_dist * 0.75    # 0.75 is a magic number
            self.cells[idx] = self.cells[idx] - offset
            new_cell = center + offset
            new_cells_arr = np.vstack([self.cells, new_cell])
            
            if concentrations is not None:
                new_conc.append(concentrations[idx])
                concentrations = np.vstack([concentrations, concentrations[idx]])
                
            self._safe_update_voronoi(new_cells_arr)
            new_cells.append(len(self.cells)-1)
        if concentrations is not None:
            return new_cells, concentrations
        return new_cells

    def cell_apoptosis(self, n=1, mode='area', concentrations=None, target_indices=None):
        """
        細胞死亡：移除面積最小或隨機選的 n 個細胞
        mode: 'area'（預設，最小面積）或 'random'（隨機）
        concentrations: (cell, var) array，若有提供則同步移除濃度
        回傳：被移除的細胞索引list與新濃度陣列（若有）
        """
        from scipy.spatial import ConvexHull
        cells = self.cells
        vor = self.vor
        # 過濾外圈細胞
        inner_indices = self.get_inner_cell_indices()
        # 計算面積
        areas = []
        for i in inner_indices:
            region = vor.regions[vor.point_region[i]]
            polygon = [vor.vertices[v] for v in region]
            try:
                hull = ConvexHull(polygon)
                areas.append(hull.volume)
            except:
                areas.append(1e9)
        # 選擇要移除的細胞
        if target_indices is not None:
            chosen_indices = [idx for idx in target_indices if idx < len(self.cells)]
        elif mode == 'area':
            chosen = np.argsort(areas)[:n]
            chosen_indices = [inner_indices[i] for i in chosen]
        elif mode == 'random':
            chosen_indices = list(np.random.choice(inner_indices, n, replace=False))
        else:
            raise ValueError('mode 必須為 area 或 random')
        keep_mask = np.ones(len(cells), dtype=bool)
        keep_mask[chosen_indices] = False
        new_cells_arr = self.cells[keep_mask]
        self._safe_update_voronoi(new_cells_arr)
        if concentrations is not None:
            new_conc_arr = concentrations[keep_mask]
            return chosen_indices, new_conc_arr
        return chosen_indices
