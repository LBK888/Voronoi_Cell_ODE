# Voronoi Cell ODE Simulation - Scripting Guide

[🇹🇼 繁體中文版本 (Traditional Chinese Version)](README_tw.md)

This guide provides detailed instructions on how to write the three core scripting components in the Voronoi Cell ODE Simulation GUI: `ode()`, `params`, and `color_func()`. These components define the mathematical rules, simulation parameters, and visual representation of your cell models.

For implementing a supplied paper, see the project-local [paper-to-voronoi-ode skill](../skills/paper-to-voronoi-ode/SKILL.md) and its executable model validator. For two-axis response maps, see [Experiments](experiments.md).

Colorbars are shared by GUI and MP4 rendering. Simple legacy color functions infer colormaps from unique `plt.cm.NAME` occurrences in source order. For repeated/reordered colormaps, explicitly add metadata to Params:

```python
colorbars=[dict(variable=0, cmap='Reds', label='A'),
           dict(variable=1, cmap='Blues', label='B')]
```

Use one state column per entry, with the same colormap and current-frame normalization used in `color_func`. An explicit list replaces automatic inference; `colorbars=[]` deliberately omits legends. Composite signals need a separate display interpretation rather than a misleading single-species label. Save Frame defaults to editable PDF; result canvases have a right-click copy/save menu.

---

## 1. The `ode()` Function
The `ode` function defines the ordinary differential equations governing the biochemical reactions and spatial diffusion within and between cells.

### Signature
```python
def ode(Y, t, params, cell_positions):
    ...
    return dY
```

### Arguments
- **`Y`**: `numpy.ndarray` of shape `(N, n_var)`
  - Represents the current concentration of variables for all `N` cells.
  - `Y[:, 0]` is the 1st variable, `Y[:, 1]` is the 2nd variable, etc.
- **`t`**: `float`
  - The current simulation time.
- **`params`**: `dict`
  - The parameter dictionary defined in the `params` editor. Contains your constants, diffusion rates, etc.
- **`cell_positions`**: `numpy.ndarray` of shape `(N, 2)`
  - The spatial (x, y) coordinates of each cell. Useful for calculating neighbors or spatial gradients.

### Return Value
- **`dY`**: `numpy.ndarray` of shape `(N, n_var)`
  - The calculated rate of change for each variable in each cell.

### Example: Simple Diffusion
```python
def ode(Y, t, params, cell_positions):
    import numpy as np
    dY = np.zeros_like(Y)
    
    # Optional: Get neighbor topology for diffusion
    try:
        from gui.sim_utils import get_voronoi_neighbors_wo_outer
        neighbors = get_voronoi_neighbors_wo_outer(cell_positions)
    except:
        neighbors = []
        
    for i in range(Y.shape[0]):
        # Calculate Laplacian for diffusion
        lap = 0.0
        if i < len(neighbors) and neighbors[i]:
            sum_w, sum_Yw = 0.0, 0.0
            for j, w in neighbors[i].items():
                sum_Yw += Y[j, 0] * w
                sum_w += w
            if sum_w > 0:
                lap = (sum_Yw / sum_w) - Y[i, 0]
                
        # Differential equation (Reaction + Diffusion)
        dY[i, 0] = -params['decay'] * Y[i, 0] + params['D'] * lap
        
    return dY
```

---

## 2. The `params` Dictionary
The `params` dictionary defines the global settings and physical constants used in the ODE and the GUI.

### Format
Must be a valid Python dictionary assignment.

### Standard Keys
- `n_var` (int): Number of biological variables (e.g., 2 for Activator/Inhibitor).
- `dT` (float): The time step size for the simulation.
- `sigma` (float): The standard deviation of the initial random noise added to `Y=0`.
- `labels` (list of strings): Names of your variables (e.g., `['Activator', 'Inhibitor']`). These are used for generating dynamic colorbars and plot legends.

### Example
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

## 3. The `color_func()` Function
This function determines how the numerical concentrations (`Y`) are mapped to visual colors in the animation. 

### Available Colormaps
Matplotlib offers a wide variety of colormaps. You can browse them here:
👉 **[Matplotlib Colormaps Documentation](https://matplotlib.org/stable/users/explain/colors/colormaps.html)**

### Per-Frame Normalization
Because concentration values can range drastically, they must be normalized to a `[0, 1]` scale before being mapped to a colormap. The preview and results colorbars now follow the minimum and maximum of the currently displayed step. Use the `global_bounds` parameter (the legacy name is retained for configuration compatibility) to receive those current-step bounds directly.

```python
import numpy as np
import matplotlib.pyplot as plt

def color_func(Y, mode='polygon', global_bounds=None):
    if mode == 'polygon':
        vmin, vmax = global_bounds[0] # Current-step bounds for Variable 0
        norm = (Y[:,0] - vmin) / (vmax - vmin + 1e-8)
        return plt.cm.plasma(norm)
            
    elif mode == 'center':
        vmin, vmax = global_bounds[1] # Current-step bounds for Variable 1
        norm = (Y[:,1] - vmin) / (vmax - vmin + 1e-8)
        return plt.cm.viridis(norm)
            
    return None
```

You can equivalently calculate the current-step bounds directly from `Y`:

```python
import numpy as np
import matplotlib.pyplot as plt

def color_func(Y, mode='polygon', global_bounds=None):
    if mode == 'polygon':
        # Normalize relative to the current frame only
        norm = (Y[:,0] - Y[:,0].min()) / (np.ptp(Y[:,0]) + 1e-8)
        return plt.cm.plasma(norm)
        
    elif mode == 'center':
        # Normalize relative to the current frame only
        norm = (Y[:,1] - Y[:,1].min()) / (np.ptp(Y[:,1]) + 1e-8)
        return plt.cm.viridis(norm)
        
    return None
```

### Supported Modes
The `mode` argument allows you to color different parts of the Voronoi cell independently:
- `'polygon'`: The interior fill of the Voronoi cell.
- `'center'`: The nucleus/center point of the cell.
- `'membrane'`: The boundary lines of the Voronoi cell.
- `'off_center'`, `'off_center2'`: Accessory points near the center.
