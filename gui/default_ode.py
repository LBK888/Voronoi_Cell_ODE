def get_default_ode():
    return '''def ode(Y, t, params, cell_positions):
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

# Optional Mechanics Function
# def mechanics(Y, t, params, cell_positions):
#     """
#     Return a dictionary evaluating mechanical updates based on ODE outcomes.
#     Expected keys:
#       'divide': list of cell indices to proliferate
#       'die': list of cell indices to undergo apoptosis
#       'forces': np.ndarray of shape (N, 2) specifying X,Y force vectors for each cell.
#     """
#     return {}
'''

def simple_LI_ode():
    pass

def simple_foxi3_LI_ode():
    pass

def simple_turing_ode():
    pass
