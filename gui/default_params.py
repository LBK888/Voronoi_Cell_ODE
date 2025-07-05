def get_default_params():
    return '''params = dict(
    nu=1.0,
    betaD=50.0,
    betaR=50.0,
    h=3,
    m=3,
    sDtv3_ratio=0.,
    sDtv4_ratio=0.,
    Ktv3_Dgr=0.4,
    Ktv4_Dgr=0.5,
    Ktv3_inhib=0.15,
    Ktv4_inhib=0.3,
    sigma_diff_sD3=2.5,
    sigma_diff_sD4=4.0,
    LI_off_Dgr=1.0,
    Dgr_Noise=0.01,
    sigma = None,
    dT=1/20.0,
    n_var=4,  
    labels= ['Delta', 'Ractor', 'sD_tv3', 'sD_tv4']
)'''
