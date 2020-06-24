"""
This script numerically solves the equations of motion for a given system over
a range of initial conditions, and saves plots of the phase orbits and poincare
sections.
"""
import functions as f
from math import pi
# from datetime import date
# import time
# import os
import click
import numpy as np
# import matplotlib.pyplot as plt

g = 9.81


@click.command()
@click.argument('img_path', type=click.Path(exists=True, readable=True))
@click.option('--points_per_sec', default=300, help='Number of points the\
              solution will be evaluated at per second')
@click.option('--k', default=0.2, help='Linear friction coefficient')
@click.option('--a', default=0.0, help='Amplitude of pivot oscillations')
@click.option('--l', default=g/(pi**2), help='Length of pendulum')
@click.option('--d', default=5.0, help='Distance from origin to magnet')
@click.option('--omega', default=2*pi, help='Angular frequency of pivot\
              oscillations')
@click.option('--b', default=1, help='Controls strength of magnet')
@click.option('--m', default=0.1, help='Mass of pendulum')
@click.option('--t_fin', default=500, help='Length of time of solution')
@click.option('--run_number', default=1, help='ith run of the day')
def main(img_path, points_per_sec, k, a, l, d, omega, b, m, t_fin, run_number): # noqa
    A = a
    L = l
    """
    today = str(date.today())
    run_number = '_' + str(run_number)
    save_dir = img_path + '/' + today + run_number

    if os.path.exists(save_dir):
        ans = input(
            'The directory this run will write to already exists, would you like \
            to overwrite it? ([y/n])')
        if ans == 'y':
            pass
        else:
            return

    else:
        os.makedirs(save_dir)


    f.write_param_file(save_dir, args)
    """

    args = (k, A, L, d, omega, b, m)
    """
    xs = []
    xps = []
    ys = []
    yps = []
    """
    # count = 0
    for theta_init in np.linspace(-pi, pi, 5, endpoint=True):
        print(theta_init)
        for theta_dot_init in np.linspace(-3, 3, 5):
            # count += 1
            times, x, y = f.solve_system(f.difeqs_param_force_pend,
                                         [theta_init, theta_dot_init],
                                         args,
                                         t_fin,
                                         points_per_sec)
            """
            x = (x + pi) % (2*pi) - pi  # keep all values in range -pi to pi
            x = x[-len(x)//4:]
            y = y[-len(y)//4:]
            xp, yp = f.p_section(x, y, points_per_sec, 1)
            xs.append(x[-len(x)//4:])
            xps.append(xp)
            ys.append(y[-len(y)//4:])
            yps.append(yp)
            fname = save_dir + '/plot_' + str(count)
            fig = plt.figure(figsize=[5, 3], dpi=72)
            ax = fig.add_subplot(111)
            ax.set_title('Phase Space A = {}, Init conds = {:6.2f}, {:6.2f}'
                         .format(A, theta_init, theta_dot_init))
            ax.set_xlabel('Theta')
            ax.set_ylabel('Theta_dot')
            ax.scatter(x[-len(x)//4:], y[-len(y)//4:], s=0.1)
            ax.scatter(xp, yp, s=10)
            plt.savefig(fname)
            plt.close()
            """


"""
    f.plot(xs,
           ys,
           save_dir+'/PhaseSpace',
           labels=['Phase Space A={}'.format(A), 'Theta', 'Theta_dot'],
           high_res=True)
"""


if __name__ == '__main__':
    main()
