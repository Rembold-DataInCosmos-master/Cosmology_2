"""
Included are several models to use for Problems 1 and 2 in Homework 6.

These are not doing anything that could not have been managed with Numpy's `quad`,
but they are taking advantage of repeated calculations to improve speeds by 10x
"""

import numpy as np

def calc_velocity(params, dist, density, num_samples=5000, spacing='log'):
    """ Computes the velocity at various points given a density function.

    Computes the circular velocity of an orbiting object given a density function at
    a provided number of points. The density function should take two arguments: a
    list of parameters and a list of distances to compute the density at.

    Inputs:
        params - list[float]: 
            A list of the parameters necessary to compute the density
        distances - np.array[float]: 
            Array of distances (in km) to compute the velocity at
        density - function:
            A density function which returns an array of densities at various distances
            Densities should be in units of kg / m^3
        num_samples - integer (default=5000):
            The number of samples the Riemann sum should be divided into
        spacing - string (default='log'):
            How the Riemann rectangles should be sampled. The default is
            logarithmically, to better sample the quickly varying small radius 
            region. Can also be set to 'linear' for a more even sampling.

    Outputs:
        np.array[float]:
            An array of velocities at each indicated distance in units of m / s

    Example:
        ```python
        def density(params, distances):
            return np.where(distances < params[0], params[1], 0)

        distances = np.linspace(100, 1E8, 1000)
        params = [2E3, 100]

        vel = calc_velocity(params, distances, density)
        ```
    """
    G = 6.67E-11
    out = np.zeros(np.array(dist).shape)
    if spacing == 'log':
        edges = np.logspace(-6, np.log10(np.max(dist)), num_samples)
    elif spacing == 'linear':
        edges = np.linspace(1E-6, np.max(dist), num_samples)
    else:
        raise ValueError('Unknown value in spacing parameter')
    bin_sizes = np.diff(edges)
    distances = edges[:-1] + bin_sizes / 2
    densities = np.array(density(params, distances))
    masses = 4*np.pi * densities * distances ** 2 * bin_sizes
    cum_masses = np.cumsum(masses)
    for i,d in enumerate(dist):
        closest_dist_idx = np.searchsorted(distances, d, side='right') - 1
        out[i] = np.sqrt(G * cum_masses[closest_dist_idx] / d)
    return out


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

def test_velocities():
    def density(params, distances):
        return np.where(distances < params[0], params[1], 0)

    distances = np.linspace(100, 1E8, 1000)
    params = [2E6, 1000]

    vel = calc_velocity(params, distances, density)
    vel2 = calc_velocity(params, distances, density, spacing='linear')

    plt.plot(distances, vel, '.-')
    plt.plot(distances, vel2, '.-')
    plt.ylim([0,1200])
    plt.show()

def test_distances():
    params = [70, 0.3, 0, 0.7]
    redshifts = np.arange(0.1, 3, 0.1)

    distances = calc_distance(params, redshifts)

    plt.plot(redshifts, distances, '.')
    plt.plot(redshifts, redshifts*3E5/70)
    plt.show()



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    test_velocities()
    test_distances()
