def get_default_color_func():
    return '''def color_func(Y, mode='polygon'):
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
        None
'''

def get_fulldemo_color_func():
    return '''def color_func(Y, mode='polygon'):
    if mode == 'polygon':
        norm = (Y[:,0] - Y[:,0].min()) / (np.ptp(Y[:,0]) + 1e-8)
        return plt.cm.Reds(norm)
    elif mode == 'center':
        norm = (Y[:,1] - Y[:,1].min()) / (np.ptp(Y[:,1]) + 1e-8)
        return plt.cm.Blues(norm)
    elif mode == 'off_center':
        norm = (Y[:,2] - Y[:,2].min() + Y[:,3] - Y[:,3].min()) / (np.ptp(Y[:,2]) + np.ptp(Y[:,3])  + 1e-8)
        return plt.cm.Greens(norm)
    elif mode == 'off_center2':
        norm = (Y[:,3] - Y[:,3].min()) / (np.ptp(Y[:,3]) + 1e-8)
        return plt.cm.Greens(norm)
    elif mode == 'membrane':
        color = np.zeros((Y.shape[0], 4))  # RGBA
        mask = Y[:,0] > 40
        color[mask] = [1, 1, 0, 1]   # 黃色 (R,G,B,A)
        color[~mask] = [0, 0, 0, 1]  # 黑色
        return color
    else:
        None
'''
