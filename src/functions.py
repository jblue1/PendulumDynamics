'''
This file contains functions used for numerically solving ODEs describing
the motion of a pendulum subject to various forcing.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os
import numba as nb

g = 9.81


@nb.jit
def f(t, theta, A, L, d, omega):
    """Helper function for the parametrically forced pendulum with external
    magnetic forcing. Calculates -0.5*r^2, where r is the distance between the
    pendulum bob and the magnet

    :param type t: Time
    :param array_like theta: Angle between pendulum and the vertical
    :param float A: Amplitude of pivot oscillations
    :param float L: length of pendulum
    :param float d: Distance between origin and magnet
    :param float omega: Angular frequency of pivot oscillations
    :rtype: Array if theta is an array, float if theta is a float

    """

    z = np.cos(omega * t)
    return -0.5 * ((L**2) + (d**2) + ((A*z) ** 2) + 2*A*L*z*np.cos(theta)
                   - 2*A*d*z - 2*d*L*np.cos(theta))


def potential_energy_function(t, theta, A, L, d, omega, b, m):
    """Short summary.

    :param type t: Time
    :param array_like theta: Angle between pendulum and the vertical
    :param float A: Amplitude of pivot oscillations
    :param float L: length of pendulum
    :param float d: Distance between origin and magnet
    :param float omega: Angular frequency of pivot oscillations
    :param float b: Controls strength of magnet
    :param float m: Mass of Pendulum
    :return: Potential energy
    :rtype: Array if theta is an array, float if theta is a float

    """
    return -m*g*L*np.cos(theta) + m*A*L*(omega**2)*np.cos(theta) * \
        np.cos(omega*t) + b*np.exp(f(t, theta, A, L, d, omega))


@nb.jit
def difeqs_param_force_pend(t, angles, k, A, L, d, omega, b, m):
    """Right hand side of equations of motion for a parametrically forced
    pendulum with a repulsive magnet directly under the stable down position.
    Used by solve_ivp.

    :param float t: Time
    :param array_like angles: Array contaning theta and theta_dot
    :param float k: Linear friction coefficient
    :param float A: Amplitude of pivot oscillations
    :param float L: length of pendulum
    :param float d: Distance between origin and magnet
    :param float omega: Angular frequency of pivot oscillations
    :param float b: Controls strength of magnet
    :param float m: Mass of Pendulum
    :return: Array containing RHS of differential equation
    :rtype: Array

    """
    z = np.cos(omega*t)
    return [angles[1],
            -k*angles[1] - (np.sin(angles[0])/L)*(g - A*(omega**2)*z +
            (b/m)*(A*z - d)*np.exp(f(t, angles[0], A, L, d, omega)))]


def difeqs_simple_pendulum(t, angles, L):
    return [
        angles[1],
        -(g/L)*np.sin(angles[0])
    ]


def plot(x, y, save_path, scatter=True, labels=['y vs. x', 'x', 'y'],
         high_res=False):
    """Quick function for saving plots

    :param array_like or list of arrays x: x values to plot
    :param type y: y values to plot
    :param string save_path: Path to where plot will be saved
    :param list labels: Labels for the graph
    :return: None
    :rtype: None

    """
    if high_res:
        fig = plt.figure(figsize=[8, 5], dpi=300)
    else:
        fig = plt.figure(figsize=[5, 3], dpi=72)

    ax = fig.add_subplot(111)
    ax.set_title(labels[0])
    ax.set_xlabel(labels[1])
    ax.set_ylabel(labels[2])
    if isinstance(x, list):
        for i in range(len(x)):
            if scatter:
                ax.scatter(x[i], y[i], s=0.5)
            else:
                ax.plot(x[i], y[i])
    else:
        if scatter:
            ax.scatter(x, y, s=0.5)
        else:
            ax.plot(x, y)
    plt.savefig(save_path)
    plt.close()
    return


@nb.jit
def p_section(x, y, points_per_sec, period):
    """Generates points for a Poincare section

    :param array_like x: Coordinate in phase space
    :param array_like y: Time derivative of x
    :param int points_per_sec: Points per second in solv_ivp output
    :param type period: Period of orbit
    :return: Points for poincare section plot
    :rtype: Tuple of arrays

    """
    points_per_period = int(period*points_per_sec)
    num_points = int(len(x)/points_per_period)
    x_points = np.zeros(num_points)
    y_points = np.zeros(num_points)

    for i in range(num_points):
        x_points[i] = x[points_per_period*i+2]
        y_points[i] = y[points_per_period*i+2]

    return x_points, y_points


def solve_system(dif_eqs, init_vals, args, t_fin, points_per_sec, t_init=0):
    """Short summary.

    :param callable dif_eqs: RHS of equations of motion. Must take (t, y, args)
    :param array_like init_vals: Initial y values
    :param tuple args: Additional arguments for dif_eqs
    :param float t_init: Initial time value
    :param float t_fin: Final time value (in seconds)
    :param int points_per_sec: Points per second in solution output
    :return: Times and y values of solutions
    :rtype: tuple

    """
    num_points = int(points_per_sec*t_fin)
    times = np.linspace(t_init, t_fin, num_points)
    sol = solve_ivp(dif_eqs, (t_init, t_fin), [init_vals[0], init_vals[1]],
                    t_eval=times, args=args, vectorized=True, rtol=1e-10,
                    atol=1e-10, method='DOP853')
    return times, sol.y[0], sol.y[1]


def write_param_file(save_dir, args):
    fname = os.path.join(save_dir, 'params.txt')
    info_list = ['k = {}\n'.format(args[0]),
                 'A = {}\n'.format(args[1]),
                 'L = {}\n'.format(args[2]),
                 'd = {}\n'.format(args[3]),
                 'omega = {}\n'.format(args[4]),
                 'b = {}\n'.format(args[5]),
                 'm = {}\n'.format(args[6])]
    with open(fname, 'w') as f:
        f.writelines(info_list)

    return


def main():
    pass


if __name__ == '__main__':
    main()
