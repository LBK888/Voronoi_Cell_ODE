import numpy as np
from tqdm import tqdm
from scipy.integrate import solve_ivp

class BiophysicsModel:
    def __init__(self, cell_count, params, ode_func, mechanics_func=None,
                 init_Y=None, vor_grid=None, move_rule=None,
                 random_strength=0.0, LaterInhib_switch=None,
                 movement_plan=None):
        self.cell_count = cell_count
        self.params = params
        if ode_func is None:
            raise Exception("ode_func must be provided!")
        self.ode_func = ode_func  # 必須傳入
        if init_Y is None:
            raise ValueError('init_Y must be provided，and shape=(cell_count, n_var)')
        self.init_Y = init_Y
        self.mechanics_func = mechanics_func
        self.vor_grid = vor_grid
        self.move_rule = move_rule
        self.random_strength = random_strength
        self.movement_plan = movement_plan
        self.reset()
        self.LaterInhib_switch = LaterInhib_switch

    def reset(self):
        self.Y = self.init_Y.copy()
        self.history = [self.Y.copy()]
        self.cell_positions_history = [self.vor_grid.cells.copy()] if self.vor_grid else []
        self.cell_ids_history = [np.arange(len(self.Y), dtype=int)]

    def _rhs(self, t, y_flat):
        y = y_flat.reshape((self.cell_count, -1))
        if self.ode_func is None:
            raise Exception("ode_func must be provided!")
        dy = self.ode_func(y, t, self.params)
        return dy.flatten()

    def simulate(self, T, proliferation_steps=None, apoptosis_steps=None, proliferation_n=5, apoptosis_n=3, proliferation_mode='area', apoptosis_mode='area', progress_callback=None, check_cancel=None, warning_callback=None):
        """
        T: 總模擬時間
        proliferation_steps: list, 在哪些步驟進行細胞分裂
        apoptosis_steps: list, 在哪些步驟進行細胞死亡
        proliferation_n, apoptosis_n: 每次分裂/死亡的細胞數
        proliferation_mode, apoptosis_mode: 'area' 或 'random'
        progress_callback: function(current_step, total_steps)
        check_cancel: function() -> bool (returns True if should cancel)
        warning_callback: function(msg)
        """
        import numpy as np
        try:
            from tqdm import trange
        except ImportError:
            trange = range
            
        y0 = self.init_Y.flatten()
        t_eval = np.arange(0, T, self.params['dT'])
        Y = y0.reshape((self.cell_count, -1))
        history = [Y.copy()]
        cell_positions_history = [self.vor_grid.cells.copy()] if self.vor_grid else []
        cell_ids = np.arange(len(Y), dtype=int)
        next_cell_id = len(cell_ids)
        cell_ids_history = [cell_ids.copy()]
        
        total_steps = len(t_eval) - 1
        
        # Use range instead of trange to avoid terminal clutter when using UI, but keep trange if progress_callback is None maybe?
        iterator = range(1, len(t_eval))
        
        warnings_emitted = 0
        
        for i, t in zip(iterator, t_eval[1:]):
            if check_cancel is not None and check_cancel():
                break

            # ODE
            dY = self.ode_func(Y, t, self.params, self.vor_grid.cells)
            Y = Y + dY * self.params['dT']

            # Biological Linter Check (加入開關機制)
            # 預設為 True (保護濃度不小於0)，但如果設定檔有寫 allow_negative=True 則不限制
            if not self.params.get('allow_negative', False):
                if np.any(Y < -1e-5):
                    if warning_callback is not None and warnings_emitted < 5:
                        warning_callback(f"Bio-Linter Warning: Negative concentrations at t={t:.2f}. Auto-clipped to 0.")
                        warnings_emitted += 1
                Y = np.clip(Y, 0, None)
                #print("allow_negative is false, Y clipped at t=", t)

            # Cell Mechanics API
            division_targets = []
            apoptosis_targets = []
            forces = None
            if self.mechanics_func is not None:
                mech_results = self.mechanics_func(Y, t, self.params, self.vor_grid.cells)
                if isinstance(mech_results, dict):
                    division_targets = mech_results.get('divide', [])
                    apoptosis_targets = mech_results.get('die', [])
                    forces = mech_results.get('forces', None)

            # 細胞分裂
            if division_targets:
                new_cells, new_Y = self.vor_grid.cell_proliferation(n=len(division_targets), mode=proliferation_mode, target_indices=division_targets, concentrations=Y)
                Y = new_Y
                if new_cells:
                    new_ids = np.arange(
                        next_cell_id, next_cell_id + len(new_cells), dtype=int
                    )
                    cell_ids = np.concatenate((cell_ids, new_ids))
                    next_cell_id += len(new_cells)
            elif proliferation_steps == "all" or (isinstance(proliferation_steps, (list, tuple, set)) and i in proliferation_steps):
                new_cells, new_Y = self.vor_grid.cell_proliferation(n=proliferation_n, mode=proliferation_mode, concentrations=Y)
                Y = new_Y
                if new_cells:
                    new_ids = np.arange(
                        next_cell_id, next_cell_id + len(new_cells), dtype=int
                    )
                    cell_ids = np.concatenate((cell_ids, new_ids))
                    next_cell_id += len(new_cells)
                
            # 細胞死亡
            if apoptosis_targets:
                removed_cells, new_Y = self.vor_grid.cell_apoptosis(n=len(apoptosis_targets), mode=apoptosis_mode, target_indices=apoptosis_targets, concentrations=Y)
                Y = new_Y
                if removed_cells:
                    cell_ids = np.delete(cell_ids, removed_cells)
            elif apoptosis_steps == "all" or (isinstance(apoptosis_steps, (list, tuple, set)) and i in apoptosis_steps):
                removed_cells, new_Y = self.vor_grid.cell_apoptosis(n=apoptosis_n, mode=apoptosis_mode, concentrations=Y)
                Y = new_Y
                if removed_cells:
                    cell_ids = np.delete(cell_ids, removed_cells)
                
            history.append(Y.copy())
            cell_ids_history.append(cell_ids.copy())
            
            # 細胞移動
            if forces is not None:
                self.vor_grid.move_cells(lambda cells: cells + forces, self.random_strength)
            elif self.movement_plan is not None:
                self.vor_grid.apply_movement(self.movement_plan)
            elif self.move_rule is not None:
                self.vor_grid.move_cells(self.move_rule, self.random_strength)
            elif self.random_strength > 0:
                self.vor_grid.move_cells(None, self.random_strength)
                
            if self.vor_grid is not None:
                cell_positions_history.append(self.vor_grid.cells.copy())
                
            if progress_callback is not None:
                progress_callback(i, total_steps)
        
        history_arr = np.empty(len(history), dtype=object)
        for idx, arr in enumerate(history):
            history_arr[idx] = arr
        cell_positions_arr = np.empty(len(cell_positions_history), dtype=object)
        for idx, arr in enumerate(cell_positions_history):
            cell_positions_arr[idx] = arr
        self.cell_ids_history = cell_ids_history
        
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
