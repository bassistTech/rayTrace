'''
Drawing functions for ray trace tools
Francis Deck 10-4-2023
MIT License
'''

import numpy as np
import matplotlib.pyplot as plt

import opticsRayTrace.rayTraceTools as rtt


def new_plot(axes=["x", "y"], **kwargs):
    '''
    axes is an optional list of axes that you actually want to plot
       ... the options are "x", "y", and "3d"

    The function plt.subplots() returns either a single axis, or a list of axes, which makes
    it hard to treat them the same in subsequent program. So I'm just reconciling this issue
    by always returning a list of axes, even if the list contains just one element.
    '''

    fig, axs = plt.subplots(len(axes), 1, sharex="col", sharey="row", **kwargs)

    if len(axes) == 1:
        return {"fig": fig, "axs": [axs]} | {axes[i]: i for i in range(len(axes))}
    else:
        return {"fig": fig, "axs": axs} | {axes[i]: i for i in range(len(axes))}


def _axdPlotGeneral(ax, x, h_axis, v_axis, color, **kwargs):
    '''
    General plot of 3d data projected onto 2d plane

    x = list of points in space
    h_axis = horizontal axis vector
    v_axis = vertical axis vector
    '''

    ha = np.sum(x*rtt.normalize_vec(h_axis), axis=1)
    va = np.sum(x*rtt.normalize_vec(v_axis), axis=1)
    ax.plot(ha, va, c=color, **kwargs)

rotate_matrix_3d = np.matmul(np.matmul(rtt.rotate_matrix([0, 1, 0], 30), 
                                       rtt.rotate_matrix([1, 0, 0], 30)), 
                             rtt.rotate_matrix([0, 0, 1], -10)).T

def _axdPlot(x, axd, label, color, **kwargs):
    '''
    Based on axis list, generate one of 3 kinds of plots
    '''

    # first graph shows the z-x coordinate axis in global space
    if "x" in axd:
        _axdPlotGeneral(axd["axs"][axd["x"]], x, np.array([0, 0, 1]), np.array([1, 0, 0]),
                        color, **kwargs)
        axd["axs"][axd["x"]].set_ylabel('x (mm)')
        #axd['fig'].supxlabel('Z (mm)')

    # second graph shows the z-y coordinate axis in global space
    if "y" in axd:
        _axdPlotGeneral(axd["axs"][axd["y"]], x, np.array([0, 0, 1]), np.array([0, 1, 0]),
                        color, **kwargs)
        axd["axs"][axd["y"]].set_ylabel('y (mm)')
        #axd['fig'].supxlabel('Z (mm)')

    if "3d" in axd:
        '''
        Form coordinates, then rotate in 2 directions to form basis for graph. For now, this
        is going to be quite crude, with only one option for orientation. It's enough to show
        that I can in fact produce a 3d drawing
        '''
        v1a = np.array([0, 0, 1])
        v2a = np.array([1, 0, 0])

        v1d = np.dot(v1a, rotate_matrix_3d)
        v2d = np.dot(v2a, rotate_matrix_3d)

        _axdPlotGeneral(axd["axs"][axd["3d"]], x, v1d, v2d,
                        color, **kwargs)
        
    axd["axs"][-1].set_xlabel('Z (mm)')

    [ax.set_aspect("equal") for ax in axd["axs"]]
    axd['fig'].tight_layout()


def plot_faces(axd, surface_list, **kwargs):
    '''
    axd = axis dictionary returned by new_plot()
    surface_list is from rtt.build_geometry()
    **kwargs are passed along to plotting function
    '''

    # This is where we will put the 3-d points for a given plot
    npts = 50
    x = np.empty((npts, 3))

    for s in surface_list:
        if s["draw"]:
            if "label" in s:
                label = s['label']
            else:
                label = None

            # create circular aperture in local coordinates, no option for rectangular drawing yet
            theta_ring = np.linspace(0, 2*np.pi, npts)
            x_ring = s["draw_radius"]*np.sin(theta_ring)
            y_ring = s["draw_radius"]*np.cos(theta_ring)

            # create values for x or y bar in local coordinates
            bar = np.linspace(-s["draw_radius"], s["draw_radius"], npts)

            # generate drawing curves in local coordinates of object
            match s["surf"]:
                case ("dummy" | "plane grating" | "plane grating 2" | "rotate" | "shift" | "paraxial"):
                    sag_ring = np.zeros((npts))
                    sag_xbar = np.zeros((npts))
                    sag_ybar = np.zeros((npts))

                case "conic":
                    r_ring = np.ones_like(x_ring)*s["draw_radius"]
                    sag_ring = rtt.sag_conic(r_ring, s["c"], s["k"])
                    sag_xbar = rtt.sag_conic(bar, s["c"], s["k"])
                    sag_ybar = sag_xbar

                case "cylindrical":
                    sag_ring = rtt.sag_conic(y_ring, s["c"], s["k"])
                    sag_xbar = np.zeros((npts))
                    sag_ybar = rtt.sag_conic(bar, s["c"], s["k"])

            local_points = np.empty((npts, 3))

            local_points[:, 0] = x_ring
            local_points[:, 1] = y_ring
            local_points[:, 2] = sag_ring
            global_points = np.dot(local_points, s["basis"]) + s["origin"]

            _axdPlot(global_points, axd, label, s["color"], **kwargs)

            local_points[:, 0] = bar
            local_points[:, 1] = np.zeros((npts))
            local_points[:, 2] = sag_xbar
            global_points = np.dot(local_points, s["basis"]) + s["origin"]

            _axdPlot(global_points, axd, label, s["color"], **kwargs)

            local_points[:, 0] = np.zeros((npts))
            local_points[:, 1] = bar
            local_points[:, 2] = sag_ybar
            global_points = np.dot(local_points, s["basis"]) + s["origin"]

            _axdPlot(global_points, axd, label, s["color"], **kwargs)

    axd['fig'].tight_layout()


def plot_rays(axd, geometry, ray_table, color_by='field', **kwargs):
    '''
    Plot the rays in the ray table.
    axd is from new_plot() above
    geometry is from rtt.build_geometry()
    ray_table = ray table
    **kwargs are passed along to plotting function
    '''

    keep = [s['draw_radius'] != 0 for s in geometry]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i in range(ray_table.shape[1]):
        if color_by == 'field':
            icolor = ray_table[0, i, 3, 0].astype(int)
        elif color_by == 'pupil':
            icolor = ray_table[0, i, 3, 1].astype(int)
        elif color_by == 'wavelength':
            icolor = ray_table[0, i, 3, 2].astype(int)
        else:
            icolor = i
        x = ray_table[:, i, 0, :][keep]
        color = colors[icolor % len(colors)]
        _axdPlot(x, axd, 'x', color, **kwargs)
    return

print('drawing loaded')