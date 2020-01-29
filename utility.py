import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_ci(array_of_numbers, as_string=True, ci=95):
    """
    Calculates a confidence interval and a median from an array of numbers
    :param array_of_numbers: The array to calculate
    :param as_string: Whether or not to return the output as a fancy string
    :param ci: The confidence interval percentage
    :return: the median, lower bound and upper bound
    """
    median_value = np.median(array_of_numbers)
    ci_frac = (100 - ci) / 2
    percentiles = np.percentile(array_of_numbers, [ci_frac, 100 - ci_frac])
    if as_string:
        return f'{median_value:.1f} wells ({percentiles[0]:.1f}-{percentiles[1]:.1f}, 95% CI)'
    else:
        return median_value, percentiles[0], percentiles[1]


def generate_3d_plot(x_label, y_label, z_label, x, y, z, filename):
    """
    Generates the 3d plots
    :param y_label: The label of the y axis
    :param x: The x axis values to plot
    :param y: The y axis values to plot
    :param z: The z axis values to plot
    :param filename: The filename to save it to
    :return: None
    """
    print(f'Generating plot for {x_label} vs {y_label} effect on {z_label} (may take some time).....')
    fig = plt.figure(figsize=[6, 6], dpi=100)
    ax = Axes3D(fig)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_zlabel(f' % {z_label}', fontsize=12)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)

    ax = fig.gca(projection='3d')

    def init():
        # Plot the surface.
        ax.plot_trisurf(x, y, z,
                        cmap=plt.cm.jet, linewidth=0.0, antialiased=True)
        return fig,

    def animate(i):
        # azimuth angle : 0 deg to 360 deg
        ax.view_init(elev=10, azim=i + 1)
        return fig,

    # Animate
    ani = animation.FuncAnimation(fig, animate, init_func=init,
                                  frames=360, interval=20, blit=True)

    ani.save(f'./plots/{filename}', writer='ffmpeg', fps=30)
