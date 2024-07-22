from myPackages import *

c0 = 3e5
H_0 = 70
Omega_l = 0.7
Omega_m = 0.3

def setup(work_path='.'):
    '''
    Set up all of the necessary directories
    '''
    for subdir in ('inputs', 'outputs', 'bin',  
                   'outputs/plots', 'outputs/eigenvalues', 'outputs/signals'
                   ):
        path = os.path.join(work_path, subdir)
        if not os.path.exists(path):
            os.mkdir(path)
            print(f'Built directory: {path}')
    outputs_dir = os.path.join(work_path, 'outputs')
    plots_dir = os.path.join(work_path, 'outputs/plots')
    eigenvalues_dir = os.path.join(work_path, 'outputs/eigenvalues')
    signals_dir = os.path.join(work_path, 'outputs/signals')
    inputs_dir = os.path.join(work_path, 'inputs')
    return outputs_dir, plots_dir, eigenvalues_dir, signals_dir, inputs_dir

def M_lim(z):
    """Fitting function for the mass completeness limit Weaver et. al 2022"""
    return np.log10(-1.51e6 * (1+z) + 6.8e7 * (1+z)**2)
    
def M_lim_ks(z):
    """Fitting function for mass completeness limit (on K_s) Weaver et al 2022"""
    return np.log10(-3.55e8 * (1+z) + 2.7e8 * (1+z)**2)

def slice_width(z, _physical_width):
    """Calculate the width of each redshift slice of size _physical_width_ (Mpc h^-1)"""
    return _physical_width * 100 / c0 * np.sqrt(Omega_m * (1+z)**3 + Omega_l)

def redshift_bins(zmin, zmax, Phys_width):
    """returns the slice centers and widths, given a physical length in (Mpc h^-1) """
    centers = []
    centers.append(zmin + 0.5 * slice_width(zmin, Phys_width))

    i = 0
    while (centers[i] + slice_width(centers[i], Phys_width) < zmax ):
        centers.append(centers[i] + slice_width(centers[i], Phys_width))
        i += 1

    centers = np.array(centers)

    "redshift edges"
    edges = np.zeros((len(centers), 2))

    for i in range(0, len(centers)):
        edges[i, 0] = centers[i] - slice_width(centers[i], Phys_width) / 2
        edges[i, 1] = centers[i] + slice_width(centers[i], Phys_width) / 2

    return (centers, edges)

def cartesian_from_polar(phi, theta):
    """ 
    phi, theta : float or numpy.array
        azimuthal and polar angle in radians.
    Returns
    -------
    nhat : numpy.array
        unit vector(s) in direction (phi, theta).
    """
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])

def cos_dist(alpha, delta, alpha0, delta0):
    """ gets all angles in [deg]"""
    phi = alpha * np.pi / 180
    theta = np.pi / 2 - delta * np.pi / 180
    phi0 = alpha0 * np.pi / 180
    theta0 = np.pi / 2 - delta0 * np.pi / 180
    
    x = cartesian_from_polar(phi, theta)
    x0 = cartesian_from_polar(phi0, theta0)
    cosdist = np.tensordot(x, x0, axes=[[0], [0]])
    return np.clip(cosdist, 0, 1)

def logsinh(x):
    if np.any(x < 0):
        raise ValueError("logsinh only valid for positive arguments")
    return x + np.log(1-np.exp(-2*x)) - np.log(2)

def Log_K(alpha, delta, alpha0, delta0, kappa):
    norm = -np.log(4 * np.pi / kappa) - logsinh(kappa)
    return norm + cos_dist(alpha, delta, alpha0, delta0) * kappa

def σ_k(X0, b, points):
    kappa = 1 / (b * np.pi / 180)**2
    X0_x = points[X0, 0]
    X0_y = points[X0, 1]
    rem = np.delete(points, X0, axis = 0)
    arr = rem[:, 2] * np.exp(Log_K(rem[:, 0], rem[:, 1], X0_x, X0_y, kappa))
    return np.sum(arr)
    
def LCV(b, points):
    N = len(points)
    arr1 = [np.log(σ_k(i, b, points)) for i in range(0, len(points))]
    return (1 / N) * np.sum(arr1)

def σ(alpha, delta, b_i, points):
    kappa = 1 / (b_i * np.pi / 180)**2
    arr2 = points[:, 2] * np.exp(Log_K(points[:, 0], points[:, 1], alpha, delta, kappa))
    return np.sum(arr2)

def Adaptive_b(b, points):
    g_i = np.array([np.log(points[i, 4] * σ(points[i, 0], points[i, 1], b, points)) for i in range(0, len(points))])
    log_g = 1 / len(points) * np.sum(g_i)
    b_i = np.array([(b * (points[i, 4] * σ(points[i, 0], points[i, 1], b, points) / np.exp(log_g))** -0.5) for i in tqdm_notebook(range(0, len(points)))])
    return b_i

def divider_NUV(rj):
    """Ilbert et al. 2013"""
    return (3*rj+1)