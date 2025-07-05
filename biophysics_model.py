import numpy as np
from tqdm import tqdm
from scipy.integrate import solve_ivp

class BiophysicsModel:
    def __init__(self, cell_count, params, ode_func, init_Y, vor_grid=None, move_rule=None, random_strength=0.0, LaterInhib_switch=None):
        self.cell_count = cell_count
        self.params = params
        if ode_func is None:
            raise Exception("ode_func must be provided!")
        self.ode_func = ode_func  # 必須傳入
        if init_Y is None:
            raise ValueError('init_Y must be provided，and shape=(cell_count, n_var)')
        self.init_Y = init_Y
        self.vor_grid = vor_grid
        self.move_rule = move_rule
        self.random_strength = random_strength
        self.reset()
        self.LaterInhib_switch = LaterInhib_switch

    def reset(self):
        self.Y = self.init_Y.copy()
        self.history = [self.Y.copy()]
        self.cell_positions_history = [self.vor_grid.cells.copy()] if self.vor_grid else []

    def _rhs(self, t, y_flat):
        y = y_flat.reshape((self.cell_count, -1))
        if self.ode_func is None:
            raise Exception("ode_func must be provided!")
        dy = self.ode_func(y, t, self.params)
        return dy.flatten()

    def simulate(self, T, proliferation_steps=None, apoptosis_steps=None, proliferation_n=5, apoptosis_n=3, proliferation_mode='area', apoptosis_mode='area'):
        """
        T: 總模擬時間
        proliferation_steps: list, 在哪些步驟進行細胞分裂
        apoptosis_steps: list, 在哪些步驟進行細胞死亡
        proliferation_n, apoptosis_n: 每次分裂/死亡的細胞數
        proliferation_mode, apoptosis_mode: 'area' 或 'random'
        """
        from tqdm import trange
        y0 = self.init_Y.flatten()
        t_eval = np.arange(0, T, self.params['dT'])
        Y = y0.reshape((self.cell_count, -1))
        history = [Y.copy()]
        cell_positions_history = [self.vor_grid.cells.copy()] if self.vor_grid else []
        for i, t in zip(trange(1, len(t_eval)), t_eval[1:]):
            # ODE
            dY = self.ode_func(Y, t, self.params, self.vor_grid.cells)
            Y = Y + dY * self.params['dT']

            # 細胞分裂
            if proliferation_steps == "all" or (isinstance(proliferation_steps, (list, tuple, set)) and i in proliferation_steps):
                new_cells, new_Y = self.vor_grid.cell_proliferation(n=proliferation_n, mode=proliferation_mode, concentrations=Y)
                Y = new_Y
            # 細胞死亡
            if apoptosis_steps == "all" or (isinstance(apoptosis_steps, (list, tuple, set)) and i in apoptosis_steps):
                removed_cells, new_Y = self.vor_grid.cell_apoptosis(n=apoptosis_n, mode=apoptosis_mode, concentrations=Y)
                Y = new_Y
            history.append(Y.copy())
            # 細胞移動
            if self.move_rule is not None:
                self.vor_grid.move_cells(self.move_rule, self.random_strength)
            if self.vor_grid is not None:
                cell_positions_history.append(self.vor_grid.cells.copy())
        
        # 由於細胞數會變動，history 需用 object array
        history_arr = np.empty(len(history), dtype=object)
        for idx, arr in enumerate(history):
            history_arr[idx] = arr
        cell_positions_arr = np.empty(len(cell_positions_history), dtype=object)
        for idx, arr in enumerate(cell_positions_history):
            cell_positions_arr[idx] = arr
        
        return history_arr, cell_positions_arr

    @staticmethod
    def parameter_scan(cell_count, params, scan_shape, T, ode_func, init_Y=None):
        results = np.zeros(scan_shape)
        for rR in tqdm(range(scan_shape[0]), desc='Parameter Scan (rR)'):
            for rT in range(scan_shape[1]):
                p = params.copy()
                p['betaDa'] += (np.random.rand()-0.5)*p['betaDa']*(1.25*rT/scan_shape[1])
                p['betaDb'] += (np.random.rand()-0.5)*p['betaDb']*(1.25*rT/scan_shape[1])
                p['betaAb'] += (np.random.rand()-0.5)*p['betaAb']*(1.25*rT/scan_shape[1])
                p['betaAr'] += (np.random.rand()-0.5)*p['betaAr']*(1.25*rT/scan_shape[1])
                p['betaBr'] += (np.random.rand()-0.5)*p['betaBr']*(1.25*rT/scan_shape[1])
                model = BiophysicsModel(cell_count, p, ode_func, init_Y)
                hist, _ = model.simulate(T)
                # hist 可能是 object array（細胞數可變），只取最後一幀的第0個細胞
                last = hist[-1]
                if len(last) > 0:
                    cell0 = last[0]
                    if cell0[2] > cell0[3]:
                        results[rR, rT] = cell0[3]/(cell0[2]+1e-7)
                    else:
                        results[rR, rT] = -cell0[2]/(cell0[3]+1e-7)
                else:
                    results[rR, rT] = np.nan  # 若最後一幀沒細胞
        return results 