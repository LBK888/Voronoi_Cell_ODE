from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from numba import njit
from scipy.spatial import Voronoi

from voronoi_grid import VoronoiGrid

def get_voronoi_neighbors_wo_outer(cells):
    """Compatibility adapter for user-authored ODE code."""
    return VoronoiGrid.from_cells(cells).finite_ridge_neighbors()

@njit
def diffusion_weighted_mean(values, positions, i, sigma):
    n_cells = positions.shape[0]
    
    # 提取中心細胞 i 的座標
    pos_i_x = positions[i, 0]
    pos_i_y = positions[i, 1]
    
    # 預先計算分母的常數部分，節省迴圈內的運算資源
    two_sigma_sq = 2.0 * sigma * sigma
    
    sum_weights = 0.0
    sum_weighted_values = 0.0
    
    # 執行單次迴圈遍歷
    for j in range(n_cells):
        # 1. 計算 X 與 Y 的座標差
        dx = positions[j, 0] - pos_i_x
        dy = positions[j, 1] - pos_i_y
        
        # 2. 取得距離的平方 (直接使用，不需要 np.sqrt)
        dist_sq = dx**2 + dy**2
        
        # 3. 計算高斯權重
        weight = np.exp(-dist_sq / two_sigma_sq)
        
        # 4. 累加至總和
        sum_weights += weight
        sum_weighted_values += weight * values[j]
        
    # 回傳加權平均值
    return sum_weighted_values / sum_weights


# below is not suitable for numba
#    dists = np.linalg.norm(positions - positions[i], axis=1)
#    weights = np.exp(-dists**2 / (2 * sigma**2))
#    return np.sum(weights * values) / np.sum(weights)

def sD_ode(Y, t, params, cell_positions):
    p = params
    dY = np.zeros_like(Y)
    D, R, sD3, sD4 = Y[:,0], Y[:,1], Y[:,2], Y[:,3]
    neighbors = get_voronoi_neighbors_wo_outer(cell_positions)
    sigma_diff_sD3 = p.get('sigma_diff_sD3', 2.0)
    sigma_diff_sD4 = p.get('sigma_diff_sD4', 2.0)
    for i in range(Y.shape[0]):
        neibs = neighbors[i]
        if neibs:
            weights = np.array(list(neibs.values()))
            dvals = np.array([D[j] for j in neibs.keys()])
            sD3_weighted = diffusion_weighted_mean(sD3, cell_positions, i, sigma_diff_sD3)
            sD4_weighted = diffusion_weighted_mean(sD4, cell_positions, i, sigma_diff_sD4)
            avgD = (
                np.sum(dvals * weights) / np.sum(weights)
                - p['Ktv3_inhib'] * sD3_weighted
                - p['Ktv4_inhib'] * sD4_weighted
            )
        else:
            dY[i,:] = 0
            continue
        dD = p['nu']*(1-p['sDtv3_ratio']-p['sDtv4_ratio'])*( p['betaD']/(1 + R[i]**p['h'])) - (1 + p['Dgr_Noise']*(np.random.randn()-0.5))*D[i]
        dR = (p['betaR']*(avgD**p['m']))/(1 + avgD**p['m']) - R[i]
        dsD3 = p['nu']*p['sDtv3_ratio']*(p['betaD']/(1 + R[i]**p['h'])) - p['Ktv3_Dgr']*sD3[i]
        dsD4 = p['nu']*p['sDtv4_ratio']*(p['betaD']/(1 + R[i]**p['h'])) - p['Ktv4_Dgr']*sD4[i]
        dY[i,0] = dD
        dY[i,1] = dR 
        dY[i,2] = dsD3
        dY[i,3] = dsD4
    Y_new = Y + dY * p['dT'] if 'dT' in p else Y + dY
    Y_new = np.clip(Y_new, 0, None)
    dY = (Y_new - Y) / (p['dT'] if 'dT' in p else 1)
    return dY

def color_func(Y, mode='polygon'):
    if mode == 'polygon':
        norm = (Y[:,0] - Y[:,0].min()) / (np.ptp(Y[:,0]) + 1e-8)
        return plt.cm.Reds(norm)
    elif mode == 'center':
        norm = (Y[:,1] - Y[:,1].min()) / (np.ptp(Y[:,1]) + 1e-8)
        return plt.cm.Blues(norm)
    elif mode == 'off_center':
        norm = (Y[:,2] - Y[:,2].min() + Y[:,3] - Y[:,3].min()) / (np.ptp(Y[:,2]) + np.ptp(Y[:,3])  + 1e-8)
        return plt.cm.Greens(norm)
    else:
        return ['gray'] * Y.shape[0]

def move_away_from_center(cells, strength=0.1, center_point=None):
    if center_point is None:
        center_point = np.mean(cells, axis=0)
    vectors = cells - center_point
    norms = np.linalg.norm(vectors, axis=1, keepdims=True) + 1e-8
    directions = vectors / norms
    return cells + directions * strength

def covergent_extension(cells, strength=0.02, max_force=1):
    vor = Voronoi(cells)
    neighbors = defaultdict(set)
    for p1, p2 in vor.ridge_points:
        neighbors[p1].add(p2)
        neighbors[p2].add(p1)
    outer_set = set()
    for i, region_index in enumerate(vor.point_region):
        region = vor.regions[region_index]
        if not region or any(v < 0 for v in region):
            outer_set.add(i)
    center_y = np.mean(cells[:, 1])
    center_x = np.mean(cells[:, 0])
    new_cells = cells.copy()
    for i in range(len(cells)):
        if i in outer_set:
            continue
        force = np.zeros(2)
        dy = center_y - cells[i, 1]
        force[1] += strength * (dy/abs(dy))
        dx = cells[i, 0] - center_x
        force[0] += strength * (dx/abs(dx)) /4
        norm = np.linalg.norm(force)
        if norm > max_force:
            force = force / norm * max_force
        new_cells[i] += force
    return new_cells

def repulsion_move_neighbors_no_outer(cells, strength=0.08, min_dist=1.5, max_force=1):
    vor = Voronoi(cells)
    neighbors = defaultdict(set)
    for p1, p2 in vor.ridge_points:
        neighbors[p1].add(p2)
        neighbors[p2].add(p1)
    outer_set = set()
    for i, region_index in enumerate(vor.point_region):
        region = vor.regions[region_index]
        if not region or any(v < 0 for v in region):
            outer_set.add(i)
    new_cells = cells.copy()
    for i in range(len(cells)):
        if i in outer_set:
            continue
        force = np.zeros(2)
        for j in neighbors[i]:
            if i == j:
                continue
            vec = cells[i] - cells[j]
            dist = np.linalg.norm(vec)
            if dist < 1e-8:
                continue
            if dist < min_dist:
                force += (vec / dist) * (strength * (min_dist - dist))
        if np.linalg.norm(force) > max_force:
            force = force / np.linalg.norm(force) * max_force
        new_cells[i] += force
    return new_cells

def vertex_model_forces(cells, adhesion_strength=0.05, repulsion_strength=0.1, min_dist=1.5, max_force=1.0):
    try:
        vor = Voronoi(cells)
    except:
        return cells  # Fallback if Voronoi fails

    neighbors = defaultdict(dict)
    for (p1, p2), ridge_vertices in zip(vor.ridge_points, vor.ridge_vertices):
        if all(v >= 0 for v in ridge_vertices):
            pts = vor.vertices[ridge_vertices]
            length = np.linalg.norm(pts[0] - pts[1])
            neighbors[p1][p2] = length
            neighbors[p2][p1] = length

    outer_set = set()
    for i, region_index in enumerate(vor.point_region):
        region = vor.regions[region_index]
        if not region or any(v < 0 for v in region):
            outer_set.add(i)

    new_cells = cells.copy()
    for i in range(len(cells)):
        if i in outer_set:
            continue
        force = np.zeros(2)
        for j, edge_len in neighbors[i].items():
            vec = cells[j] - cells[i]
            dist = np.linalg.norm(vec)
            if dist < 1e-8:
                continue
            unit_vec = vec / dist
            
            f_adhesion = adhesion_strength * edge_len * unit_vec
            
            f_repulsion = np.zeros(2)
            if dist < min_dist:
                f_repulsion = -repulsion_strength * ((min_dist - dist) / dist) * unit_vec
                
            force += f_adhesion + f_repulsion
            
        norm = np.linalg.norm(force)
        if norm > max_force:
            force = force / norm * max_force
        new_cells[i] += force
        
    return new_cells
