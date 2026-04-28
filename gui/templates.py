from gui.default_ode import get_default_ode
from gui.default_params import get_default_params
from gui.default_color_func import get_fulldemo_color_func, get_default_color_func

TEMPLATES = {
    "Default (Lateral Inhibition via sD)": {
        "ode": get_default_ode(),
        "params": get_default_params(),
        "color_func": get_fulldemo_color_func()
    },
    "Turing Pattern (Activator-Inhibitor)": {
        "ode": '''def ode(Y, t, params, cell_positions):
    p = params
    dY = np.zeros_like(Y)
    u, v = Y[:,0], Y[:,1]
    neighbors = get_voronoi_neighbors_wo_outer(cell_positions)
    
    for i in range(Y.shape[0]):
        neibs = neighbors[i]
        if neibs:
            weights = np.array(list(neibs.values()))
            u_vals = np.array([u[j] for j in neibs.keys()])
            v_vals = np.array([v[j] for j in neibs.keys()])
            
            # Simple discrete Laplacian approximation
            laplacian_u = np.sum(weights * (u_vals - u[i])) / np.sum(weights)
            laplacian_v = np.sum(weights * (v_vals - v[i])) / np.sum(weights)
        else:
            laplacian_u = laplacian_v = 0.0
            
        # Reaction-diffusion equations
        du = p['Du'] * laplacian_u + u[i]*(1 - u[i]**2) - v[i]
        dv = p['Dv'] * laplacian_v + p['c']*(u[i] - p['b']*v[i] - p['a'])
        dY[i,0] = du
        dY[i,1] = dv
        
    Y_new = Y + dY * p['dT'] if 'dT' in p else Y + dY
    dY = (Y_new - Y) / (p['dT'] if 'dT' in p else 1)
    return dY
''',
        "params": '''params = dict(
    Du=0.1,
    Dv=1.0,
    a=0.1,
    b=0.5,
    c=2.0,
    dT=0.01,
    n_var=2,
    sigma=0.1,
    labels=['Activator', 'Inhibitor']
)''',
        "color_func": '''def color_func(Y, mode='polygon'):
    if mode == 'polygon':
        norm = (Y[:,0] - Y[:,0].min()) / (np.ptp(Y[:,0]) + 1e-8)
        return plt.cm.viridis(norm)
    elif mode == 'center':
        norm = (Y[:,1] - Y[:,1].min()) / (np.ptp(Y[:,1]) + 1e-8)
        return plt.cm.plasma(norm)
    return None
'''
    },
    "Vertex Model Geometry Demo": {
        "ode": '''def ode(Y, t, params, cell_positions):
    # Pure mechanical model, no ODE
    return np.zeros_like(Y)
''',
        "params": '''params = dict(
    dT=1.0,
    n_var=1,
    sigma=0.0,
    labels=['None']
)''',
        "color_func": get_default_color_func()
    }
}
