"""
Included are several models to use for Problems 1 and 2 in Homework 6.

These are not doing anything that could not have been managed with Numpy's `quad`,
but they are taking advantage of repeated calculations to improve speeds by 10x
"""

import numpy as np

def calc_distance(params, zs, num_samples=5000):
    """Computes the distance in Mpc given parameters and redshift.

    Inputs:
        params - list[float]: 
            a list of the cosmological parameters in H_0, Om_m, Om_k, Om_l order
        zs - np.array[float]: 
            a numpy array of redshift values
        num_samples - int (default = 5000):    
            the number of slices to compute the Riemann integral over
    Outputs:
        np.array[float]:      
            An array of distances in Mpc for each redshift

    Example:
        ```python
        params = [70, 0.3, 0, 0.7]
        redshifts = np.arange(0.1, 3, 0.1)

        distances = calc_distance(params, redshifts)
        ```
    """
    H0, om_m, om_k, om_l = params
    H = lambda z: H0 * np.sqrt(om_m*(1+z)**3 + om_k*(1+z)**2 + om_l)
    out = np.zeros(np.array(zs).shape)
    edges = np.linspace(0, np.max(zs), num_samples)
    bin_sizes = np.diff(edges)
    z_range = edges[:-1] + bin_sizes / 2
    Hs = H(z_range)
    integral = 1 / Hs * bin_sizes
    cum_int = np.cumsum(integral)
    for i,z in enumerate(zs):
        closest_z_idx = np.searchsorted(z_range, z, side='right') - 1
        out[i] = (1+z)*3E5*cum_int[closest_z_idx]
    return out

def test_distances():
    params = [70, 0.3, 0, 0.7]
    redshifts = np.arange(0.1, 3, 0.1)

    distances = calc_distance(params, redshifts)

    plt.plot(redshifts, distances, '.')
    plt.plot(redshifts, redshifts*3E5/70)
    plt.show()



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    test_distances()
