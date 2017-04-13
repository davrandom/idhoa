# from functions import sph2cart_acoustics
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import auxiliary as aux

'''
* This is not original work, but more a mix of online resources.
'''


def gainsPlot(title, angle, gains, ticks=[-180, -120, -60, 0, 60, 120, 180]):
    fig = plt.figure()
    ax = plt.gca()
    n = int(angle.size / 2)

    angle = np.mod(angle.T + np.pi, 2 * np.pi) - np.pi
    x = 180. / np.pi * np.roll(angle, n + 1)
    y = np.roll(gains, n + 1, 0)
    x = np.append(x, -x[0])
    y = np.vstack((y, y[0, :]))

    plt.xticks(ticks)
    rect = fig.patch
    rect.set_facecolor('w')

    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Gains (linear)')
    ax.set_title(title)
    ax.grid(linestyle=':')
    plt.plot([-180, 180], [0, 0], color='k', linestyle='-')
    plt.plot([0, 0], [-0.4, 1.2], color='k', linestyle='-')
    #    ax2 = ax.twinx()
    #    ax2.set_ylabel('Gains (dB)')
    #    ax2.set_yticks([0.501,0.251,0.126,0.063,-0.063,-0.126])
    #    ax2.set_yticklabels([-6, -12, -18, -24, -24, -18])
    plt.xlim(-180., 180.)
    plt.ylim(-0.4, 1.2)
    #    ax.legend(ticks.T)
    plt.plot(x, y, linewidth=4)
    fig.savefig(title + ".eps")


def sph2cart_acoustics(phi, theta, rho):
    """Acoustics convention!"""
    x = np.cos(phi) * np.cos(theta) * rho
    y = np.sin(phi) * np.cos(theta) * rho
    z = np.sin(theta) * rho
    return x, y, z


def threed_polar_plot(phi, theta, rho, numbers=False):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x, y, z = sph2cart_acoustics(phi, theta, rho)

    if not numbers:  # fast hack
        sc = ax.scatter(x, y, z)
    else:
        for i in range(len(x)):  # plot each point + it's index as text above
            ax.scatter(x[i], y[i], z[i], color='b')
            ax.text(x[i], y[i], z[i], '%s' % (str(i + 1)), size=14, zorder=1, color='k')

    # plt.show()
    fig.savefig("speakers_distribution.eps")


def polar_plot(config, title, angle, *variables):
    # plt.ion # uncomment to plot
    # radar black, solid grid lines
    plt.rc('grid', color='k', linewidth=1, linestyle='-')
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    # force square figure and square axes looks better for polar, IMO
    width, height = plt.rcParams['figure.figsize']
    size = min(width, height)
    # make a square figure
    fig = plt.figure(figsize=(size, size))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='1.0')  # axisbg is the background colour
    # white canvas colour
    rect = fig.patch
    rect.set_facecolor('w')

    colorlist = ['k', 'r', 'g', 'c', 'm', 'y', 'b', 'w']
    stylelist = ['dashdot', 'dashed', 'solid', 'solid', 'dashed', 'dashdot']
    i = 0
    maximum = 1.2
    for var in variables:
        if type(var) != tuple: ax.plot(angle, var, color=colorlist[i], linestyle=stylelist[i],
                                       linewidth=4); maximum = max(max(var) + 0.15, maximum)
        if type(var) == tuple: leg = var
        i += 1

    ax.set_rmax(maximum)
    plt.grid(True)

    plt.legend(leg, loc='upper right', bbox_to_anchor=(1.125, 1.13))

    ax.set_title(title, fontsize=20)
    plt.show()
    ti = title.split(" ", 1)[0]
    if len(title.split(" ", 1)) > 1:
        tle = title.split(" ", 1)[1]
        tle = tle.split(",", 1)
        tle = tle[0] + tle[1]
        fig.savefig(str(config.DEG) + "-" + str(config.DEC) + "-" + ti + "_" + tle + ".eps")
    else:
        fig.savefig(str(config.DEG) + "-" + str(config.DEC) + "-" + ti + ".eps")


def points_over_sphere_plotting(NP, cfg):
    # it is different from points_over_sphere because has some redundancy at 0 and 2*pi
    # you can increase redundancy to ease the plot task...
    # (plots in vertical plane now are a bit messy)
    thetaPrev = [((float(i) / (2 * NP)) * 2.0 * np.pi) for i in range(2 * NP + 1)]
    theta = []
    phi = []

    for i in range(len(thetaPrev)):
        n = max(int(2 * NP * np.cos(thetaPrev[i])), 1)
        phi.append([(float(jj) / n * 2) * np.pi for jj in range(n + 1)])
        temp = [thetaPrev[i]] * (n + 1)
        theta.append(temp)

    phiok = [item for sublist in phi for item in sublist]
    thetaok = [item for sublist in theta for item in sublist]

    if len(phiok) != len(thetaok):
        raise ValueError("Died generating points on the sphere")

    return phiok, thetaok
