# Voronoi 細胞 ODE 模擬器 - 腳本撰寫指南

[🇺🇸 English Version (英文版)](README.md)

本指南提供 Voronoi 細胞 ODE 模擬器介面中三個核心腳本元件的詳細說明：`ode()`、`params` 與 `color_func()`。這些元件定義了您細胞模型的數學運算規則、模擬參數以及視覺化呈現方式。

---

## 1. `ode()` 函數
`ode` 函數用於定義常微分方程式 (Ordinary Differential Equations)，決定了生化反應以及細胞內外的空間擴散。

### 函數簽章 (Signature)
```python
def ode(Y, t, params, cell_positions):
    ...
    return dY
```

### 參數說明
- **`Y`**: `numpy.ndarray`，形狀為 `(N, n_var)`
  - 代表所有 `N` 個細胞在當前時間點的變數濃度。
  - `Y[:, 0]` 是第一個變數，`Y[:, 1]` 是第二個變數，依此類推。
- **`t`**: `float`
  - 當前的模擬時間。
- **`params`**: `dict`
  - 來自 `params` 編輯區的參數字典，包含您定義的常數與擴散率等。
- **`cell_positions`**: `numpy.ndarray`，形狀為 `(N, 2)`
  - 每個細胞的空間 (x, y) 座標，可用於計算相鄰細胞或是空間梯度。

### 回傳值
- **`dY`**: `numpy.ndarray`，形狀為 `(N, n_var)`
  - 每個細胞中各個變數的濃度變化率 (導數)。

### 範例：簡單擴散模型
```python
def ode(Y, t, params, cell_positions):
    import numpy as np
    dY = np.zeros_like(Y)
    
    # 選用：取得相鄰細胞的拓樸關係以進行擴散計算
    try:
        from gui.sim_utils import get_voronoi_neighbors_wo_outer
        neighbors = get_voronoi_neighbors_wo_outer(cell_positions)
    except:
        neighbors = []
        
    for i in range(Y.shape[0]):
        # 計算拉普拉斯算子 (Laplacian)
        lap = 0.0
        if i < len(neighbors) and neighbors[i]:
            sum_w, sum_Yw = 0.0, 0.0
            for j, w in neighbors[i].items():
                sum_Yw += Y[j, 0] * w
                sum_w += w
            if sum_w > 0:
                lap = (sum_Yw / sum_w) - Y[i, 0]
                
        # 微分方程式 (反應 + 擴散)
        dY[i, 0] = -params['decay'] * Y[i, 0] + params['D'] * lap
        
    return dY
```

---

## 2. `params` 參數字典
`params` 字典定義了 ODE 以及 GUI 所使用的全域設定與物理常數。

### 格式
必須是合法的 Python 字典賦值語法。

### 標準鍵值 (Standard Keys)
- `n_var` (int): 生物變數的數量（例如 2 代表 Activator 與 Inhibitor）。
- `dT` (float): 模擬的時間步長。
- `sigma` (float): 加入至初始值 `Y=0` 的隨機雜訊標準差。
- `labels` (list of strings): 您的變數名稱列表（例如 `['Activator', 'Inhibitor']`）。這將用於自動生成 Colorbar 與圖表圖例。

### 範例
```python
params = {
    'n_var': 2,
    'dT': 0.02,
    'sigma': 0.01,
    'labels': ['Activator', 'Inhibitor'],
    'D_A': 0.005,
    'D_I': 0.2,
    'decay': 0.1
}
```

---

## 3. `color_func()` 函數
此函數決定了如何將數值濃度 (`Y`) 映射到動畫中的視覺色彩。

### 可用的 Colormaps (色彩對應表)
Matplotlib 提供了非常多種色彩對應表，您可以在以下連結中挑選並預覽：
👉 **[Matplotlib Colormaps 官方文件](https://matplotlib.org/stable/users/explain/colors/colormaps.html)**

### 單步標準化 (Per-Frame)
因為濃度的數值範圍變化很大，在對應到 Colormap 之前必須先被標準化到 `[0, 1]` 的區間。Preview 與 Results 的 Colorbar 現在會跟隨目前 step 的最小值與最大值。可使用 `global_bounds` 參數直接取得目前 step 的範圍；此參數保留舊名稱是為了相容既有設定檔。

```python
import numpy as np
import matplotlib.pyplot as plt

def color_func(Y, mode='polygon', global_bounds=None):
    if mode == 'polygon':
        vmin, vmax = global_bounds[0] # 取得變數 0 在目前 step 的極值
        norm = (Y[:,0] - vmin) / (vmax - vmin + 1e-8)
        return plt.cm.plasma(norm)
            
    elif mode == 'center':
        vmin, vmax = global_bounds[1] # 取得變數 1 在目前 step 的極值
        norm = (Y[:,1] - vmin) / (vmax - vmin + 1e-8)
        return plt.cm.viridis(norm)
            
    return None
```

也可以直接從目前 step 的 `Y` 計算相同範圍：

```python
def color_func(Y, mode='polygon', global_bounds=None):
    if mode == 'polygon':
        # 僅使用當前影格的極值進行標準化
        norm = (Y[:,0] - Y[:,0].min()) / (np.ptp(Y[:,0]) + 1e-8)
        return plt.cm.plasma(norm)
        
    elif mode == 'center':
        # 僅使用當前影格的極值進行標準化
        norm = (Y[:,1] - Y[:,1].min()) / (np.ptp(Y[:,1]) + 1e-8)
        return plt.cm.viridis(norm)
        
    return None
```

## 4. `mechanics(Y, t, params, cells)` 函數 (新增 - 細胞行為控制)
此函數用於在模擬循環中動態地控制細胞的增殖、凋亡與移動，而不僅僅是數值上的擴散。

### 函數簽章
```python
def mechanics(Y, t, params, cells):
    ...  
    return {"divide": [], "die": [], "forces": []}
```

### 參數說明
- **`Y`**: `numpy.ndarray`，形狀為 `(N, n_var)`
  - 所有細胞的當前濃度值。
- **`t`**: `float`
  - 當前模擬時間。
- **`params`**: `dict`
  - 模擬參數字典。
- **`cells`**: `numpy.ndarray`，形狀為 `(N, 2)`
  - 所有細胞的空間 (x, y) 座標。

### 回傳值
必須回傳一個字典，包含三個鍵：
1.  **`'divide'`**: List of integers，包含準備分裂的細胞索引。
2.  **`'die'`**: List of integers，包含準備凋亡的細胞索引。
3.  **`'forces'`**: numpy.ndarray，形狀為 `(N, 2)`，定義每個細胞受到的淨移動向量（可為零）。

### 範例程式碼
```python
def mechanics(Y, t, params, cells):
    division_targets = []  # 準備分裂的細胞 index 列表
    apoptosis_targets = [] # 準備凋亡的細胞 index 列表
    forces = np.zeros_like(cells) # 控制各細胞移動的力/位移向量 (大小為 N x 2)

    # 範例邏輯：
    for i in range(len(Y)):
        # 假設 Y[i, 0] 是某種生長因子的濃度
        concentration = Y[i, 0]
        
        # 1. 濃度過高 -> 觸發細胞分裂 (Proliferation)
        if concentration > 1.5:
            division_targets.append(i)
            
        # 2. 濃度過低 -> 觸發細胞凋亡 (Apoptosis)
        elif concentration < 0.1:
            apoptosis_targets.append(i)
            
        # 3. 控制細胞移動 (Movement)
        # 例如讓細胞朝右上方移動
        forces[i] = [0.01, 0.01] 
        
    # 最後必須回傳一個包含這三個 key 的字典 (字典 key 的名稱是固定的)
    return {
        'divide': division_targets,
        'die': apoptosis_targets,
        'forces': forces
    }
```

### 支援的繪製模式 (Modes)
您可以透過 `mode` 參數分別為 Voronoi 細胞的不同部位上色：
- `'polygon'`: Voronoi 細胞內部的填滿區域。
- `'center'`: 細胞的核心 / 中心點。
- `'membrane'`: Voronoi 細胞的邊界框線。
- `'off_center'`, `'off_center2'`: 中心點旁邊的附屬裝飾點。
