'''
Optical ray tracing formulas
Francis Deck 10-4-2023
MIT License
'''

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


def norm_vec(v):
    """
    Norm of a vector
    """
    return np.sqrt(np.sum(v*v))


def norm_vec_array(va):
    """
    Norms of an array of vectors
    va: [vectors]
    """
    return np.sqrt(np.sum(va*va, axis=1))


def normalize_vec(v):
    """
    Force vector to have unit length
    """
    return v / norm_vec(v)


def normalize_vec_array(va):
    """
    Force array of factors to have unit lengths
    """
    return va / norm_vec_array(va)[:, np.newaxis]

"""
Use Sympy to create the sag functions for conic section (standard) surfaces,
then turn them into Python functions using numpy math.
"""

r, c, k = sp.symbols('r c k')

sag_conic_sp = c*r**2/(1 + sp.sqrt(1 - (1 + k)*c**2*r**2))
sag_conic_diff_sp = sp.simplify(sp.diff(sag_conic_sp, r)/r)

sag_conic = sp.lambdify([r, c, k], sag_conic_sp, 'numpy')
sag_conic_diff = sp.lambdify([r, c, k], sag_conic_diff_sp, 'numpy')


def distance_point_line(a, n, p):
    """
    distance from point to line, possibly not used right now
    a = point on the line
    n = unit vector in direction of line
    p = point in space
    
    https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    """
    return norm_vec(p - a - n*np.dot(p - a, n))
    
def distance_point_line_array(a, n, pa):
    """
    distance from point to line
    a = point on the line
    n = unit vector in direction of line
    p = [point in space]
    
    https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    """
    return norm_vec_array(pa - a - np.outer(np.sum((pa - a)*n, axis=1), n))
    

def intersect_line_plane(x_surf, n_surf, x_line, v_line):
    """
    intersection of line and plane, possibly not used right now
    x_surf = point on plane
    n_surf = normal to plane
    x_line = point on line
    v_line = direction vector of line
    
    https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    """
    return x_line + v_line * np.dot((x_surf - x_line), n_surf) / np.dot(v_line, n_surf)
    
def intersect_line_plane_array(x_surf, n_surf, x_linea, v_line):
    """
    intersection of line and plane
    x_surf = point on plane
    n_surf = normal to plane
    x_line = [point on line]
    v_line = direction vector of line
    
    https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    """
    nps = (np.sum(
        (x_surf-x_linea)
        *n_surf, axis=1) /
           np.sum(v_line*n_surf, axis=1))[:, np.newaxis]
    return x_linea + v_line * nps


def rotate_matrix(axis, alpha_deg):
    '''
    generic rotation matrix
    
    https://www.haroldserrano.com/blog/rotations-in-computer-graphics
    '''
    alpha = alpha_deg * np.pi / 180
    a, b, c = axis
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    K = 1 - ca
    return np.matrix(
        [
            [a**2 * K + ca,       a * b * K - c * sa,  a * c * K + b * sa],
            [a * b * K + c * sa,  b**2 * K + ca,       b * c * K - a * sa],
            [a * c * K - b * sa,  b * c * K + a * sa,  c**2 * K + ca],
        ]
    )


def coord_rotate(v, axis, alpha_deg):
    '''
    coordinate rotation of a vector around an axis and angle
    v: vector or point = The thing to be rotated
    axis: vector = The rotation axis
    alpha:deg = Amount of rotation in degrees
    '''
    M = rotate_matrix(axis, alpha_deg)
    return np.array(np.matmul(M, v))[0, :]


def coord_shift(x, dx):
    '''
    coordinate shift of a vector by a distance vector
    Will probably get rid of, self explanatory
    '''
    return x + dx


def apply_geometry_defaults(d, defaults):
    """
    Helper function, this only overrides the values in d, if they are not
    already defined.
    """
    return d | {k: defaults[k] for k in defaults if k not in d}


def build_geometry(surface_list):
    geometry_defaults = {"dist": 0,
                         "draw_radius": 10, "draw_dx": 0, "draw_dy": 0}

    # This is the coordinate basis that will be cumulatively updated
    # at each surface
    origin = np.array((0, 0, 0))
    x_axis = np.array((1, 0, 0))
    y_axis = np.array((0, 1, 0))
    z_axis = np.array((0, 0, 1))

    # Run through the surfaces, applying the coordinate basis to e
    geometry = []
    for s in surface_list:
        g = apply_geometry_defaults(s, geometry_defaults)
        if g["surf"] == "rotate":
            x_axis = coord_rotate(x_axis, g["axis"], g["degrees"])
            y_axis = coord_rotate(y_axis, g["axis"], g["degrees"])
            z_axis = coord_rotate(z_axis, g["axis"], g["degrees"])
        elif g["surf"] == "shift":
            origin = coord_shift(
                origin,
                g["delta"][0] * x_axis
                + g["delta"][1] * y_axis
                + g["delta"][2] * z_axis,
            )
        # Apply the new coordinates to the current surface
        g["origin"] = np.array(origin)
        g["x_axis"] = np.array(x_axis)
        g["y_axis"] = np.array(y_axis)
        g["z_axis"] = np.array(z_axis)
        # Set the origin of the next surface
        origin = origin + g["dist"] * z_axis
        geometry.append(g)

    return geometry


def propagate_dummy_surface(ray_table, en, ex, surface):
    """
    ray_table = [surf, ray, type, point (type=0) or vector (type=1)]
    en = ray number of previous surface (entrance)
    ex = ray number of this surface (exit)
    """
    ray_table[ex, :, 0, :] = intersect_line_plane_array(surface["origin"], 
                                                        surface["z_axis"],
                                                        ray_table[en, :, 0, :], 
                                                        ray_table[en, :, 1, :])
    ray_table[ex, :, 1, :] = ray_table[en, :, 1, :]


def propagate_conic_surface(ray_table, en, ex, surface):
    """
    ray_table = [surf, ray, type, point (type=0) or vector (type=1)]
    """
    # first step is intersecting the input ray with the curved surface
    sag = np.zeros([ray_table.shape[1]])
    for k in range(10):
        ray_table[ex, :, 0, :] = intersect_line_plane_array(
            surface["origin"] + np.outer(sag, surface["z_axis"]), surface["z_axis"],
            ray_table[en, :, 0, :], ray_table[en, :, 1, :])
        r = distance_point_line_array(
            surface["origin"],
            surface["z_axis"],
            ray_table[ex, :, 0, :])
        sag = sag_conic(r, surface["c"], surface["k"])

    # second step is finding surface normal on the exit surface
    sd = sag_conic_diff(r, surface['c'], surface['k'])
    # d_sag_d_r2(r, surface["c"], surface["k"])
    x = np.dot(ray_table[ex, :, 0, :] - surface["origin"], surface["x_axis"])
    y = np.dot(ray_table[ex, :, 0, :] - surface["origin"], surface["y_axis"])

    # n = normalize(  # array of vectors
    #    np.outer(2*x*sd, surface["x_axis"])
    #    + np.outer(2*y*sd, surface["y_axis"])
    #    - surface["z_axis"])

    n = normalize_vec_array(np.outer(sd*x, surface['x_axis'])
                            + np.outer(sd*y, surface['y_axis'])
                            - surface['z_axis'])

    # third step is finding the direction of the refracted ray
    # https://image.slidesharecdn.com/raytracing-111026083006-phpapp01/95/ray-tracing-23-728.jpg
    i = -ray_table[en, :, 1, :]  # [vectors]
    ci = np.sum(i*n, axis=1)  # [scalars]
    eta = surface["n1"] / surface["n2"]  # scalar
    ct = np.sqrt(1 - eta**2 * (1 - ci**2))  # [scalars]
    vr = (eta * ci - ct)[:, np.newaxis]*n - eta*i  # [vectors]
    ray_table[ex, :, 1, :] = vr

def propagate_cylindrical_surface(ray_table, en, ex, surface):
    """
    ray_table = [surf, ray, type, point (type=0) or vector (type=1)]
    """
    # Intersect rays with curved surface
    saga = np.zeros([ray_table.shape[1]])
    ya = np.zeros_like(saga)
    for k in range(10):
        ray_table[ex, :, 0, :] = intersect_line_plane_array(
            surface["origin"] 
            + np.outer(saga, surface["z_axis"]), 
            surface["z_axis"],
            ray_table[en, :, 0, :], ray_table[en, :, 1, :])
        ya[:] = np.dot(ray_table[ex, :, 0, :] - surface["origin"], surface["y_axis"])
        saga[:] = sag_conic(ya, surface["c"], surface["k"])

    # second step is finding surface normal on the exit surface
    
    x = np.dot(ray_table[ex, :, 0, :] - surface["origin"], surface["x_axis"])
    y = np.dot(ray_table[ex, :, 0, :] - surface["origin"], surface["y_axis"])
    sdx = np.zeros_like(y)
    sdy = sag_conic_diff(y, surface['c'], surface['k'])

    n = normalize_vec_array(np.outer(sdx*x, surface['x_axis'])
                            + np.outer(sdy*y, surface['y_axis'])
                            - surface['z_axis'])

    # third step is finding the direction of the refracted ray
    i = -ray_table[en, :, 1, :]  # [vectors]
    ci = np.sum(i*n, axis=1)  # [scalars]
    eta = surface["n1"] / surface["n2"]  # scalar
    ct = np.sqrt(1 - eta**2 * (1 - ci**2))  # [scalars]
    vr = (eta * ci - ct)[:, np.newaxis]*n - eta*i  # [vectors]
    ray_table[ex, :, 1, :] = vr

def propagate_plane_grating(ray_table, en, ex, surface, ray_properties):
    """
    ray_table = [surf, ray, type, point (type=0) or vector (type=1)]
    
    # https://www.phys.uconn.edu/~eyler/phys4150/R/OSLO%20Optics%20Reference.pdf
    """

    ray_table[ex, :, 0, :] = intersect_line_plane_array(
        surface["origin"], surface["z_axis"],
        ray_table[en, :, 0, :], ray_table[en, :, 1, :])

    # three arrays of scalars in coordinates local to the surface

    exy = (
        surface["n1"] * np.dot(ray_table[en, :, 1, :], surface["y_axis"])
        + surface["m"] * ray_properties['wavl'] / surface["d"]
    ) / surface["n2"]

    exx = surface["n1"] * np.dot(ray_table[en, :, 1, :],
                                 surface["x_axis"]) / surface["n2"]

    exz = np.sqrt(1 - exx*exx - exy*exy)

    # convert to a vector in global coordinates

    ray_table[ex, :, 1, :] = (np.outer(exx, surface["x_axis"])
                              + np.outer(exy, surface["y_axis"])
                              + np.outer(exz, surface["z_axis"]))

def propagate_plane_grating_2(ray_table, en, ex, surface, ray_properties):
    """
    ray_table = [surf, ray, type, point (type=0) or vector (type=1)]
    
    # https://www.phys.uconn.edu/~eyler/phys4150/R/OSLO%20Optics%20Reference.pdf
    """

    ray_table[ex, :, 0, :] = intersect_line_plane_array(
        surface["origin"], surface["z_axis"],
        ray_table[en, :, 0, :], ray_table[en, :, 1, :])

    # three arrays of scalars in coordinates local to the surface

    exx = (
        surface["n1"] * np.dot(ray_table[en, :, 1, :], surface["x_axis"])
        + surface["m"] * ray_properties['wavl'] / surface["d"]
    ) / surface["n2"]

    exy = surface["n1"] * np.dot(ray_table[en, :, 1, :],
                                 surface["y_axis"]) / surface["n2"]

    exz = np.sqrt(1 - exx*exx - exy*exy)

    # convert to a vector in global coordinates

    ray_table[ex, :, 1, :] = (np.outer(exx, surface["x_axis"])
                              + np.outer(exy, surface["y_axis"])
                              + np.outer(exz, surface["z_axis"]))



def propagate_ray(ray_table, ray_properties, geometry):
    """
    ray_table = [surf, ray, type, point (type=0) or vector (type=1)]
    """
    # No tracing to first surface, so we start with index 1
    for i in range(1, len(geometry)):
        g = geometry[i]
        if g["surf"] in ["rotate", "shift", "dummy"]:
            propagate_dummy_surface(ray_table, i-1, i, g)
        elif g["surf"] == "conic":
            propagate_conic_surface(ray_table, i-1, i, g)
        elif g["surf"] == "cylindrical":
            propagate_cylindrical_surface(ray_table, i-1, i, g)
        elif g["surf"] == "plane grating":
            propagate_plane_grating(ray_table, i-1, i, g, ray_properties)
        elif g["surf"] == "plane grating 2":
            propagate_plane_grating_2(ray_table, i-1, i, g, ray_properties)
        else:
            print("Did not recognize surface type", g["surf"])


def new_plot_xy(axes=["x", "y"], **kwargs):
    # axes is an optional list of axes that you actually want to plot
    fig, axs = plt.subplots(len(axes), 1, sharex="col", sharey="row", **kwargs)
    # resolve fact that plt.subplots returns either a scalar axs or a list
    # depending on number of axes. We want always a list for consistency
    # later on.
    if len(axes) == 1:
        return {"fig": fig, "axs": [axs]} | {axes[i]: i for i in range(len(axes))}
    else:
        return {"fig": fig, "axs": axs} | {axes[i]: i for i in range(len(axes))}


def plot_faces(axd, surface_list):
    npts = 50
    x = np.empty((npts, 3))

    for s in surface_list:
        ra = np.linspace(-s["draw_radius"], s["draw_radius"], npts)
        xa = ra + s["draw_dx"]
        ya = ra + s["draw_dy"]
        if "label" in s:
            label = s['label']
        else:
            label = None

        # generate drawing curves in local coordinates of object
        match s["surf"]:
            case ("dummy" | "plane grating" | "plane grating 2" | "rotate" | "shift"):
                zxa = np.zeros_like(ra)
                zya = zxa
            case "conic":
                zxa = sag_conic(xa, s["c"], s["k"])
                zya = sag_conic(ya, s["c"], s["k"])
            case "cylindrical":
                zxa = np.zeros_like(ra)
                zya = sag_conic(ya, s["c"], s["k"])

        # first graph shows the z-x coordinate axis in global space
        if "x" in axd:
            for i in range(npts):
                x[i, :] = s["origin"] + xa[i] * \
                    s["x_axis"] + zxa[i] * s["z_axis"]
            axd["axs"][axd["x"]].plot(x[:, 2], x[:, 0], linewidth=1, label = label)
            axd["axs"][axd["x"]].set_ylabel("X (mm)")

        # second graph shows the z-y coordinate axis in global space
        if "y" in axd:
            for i in range(npts):
                x[i, :] = s["origin"] + ya[i] * \
                    s["y_axis"] + zya[i] * s["z_axis"]
            axd["axs"][axd["y"]].plot(x[:, 2], x[:, 1], linewidth=1, label = label)
            axd["axs"][axd["y"]].set_ylabel("Y (mm)")
        axd["axs"][-1].set_xlabel("Z (mm)")

    [ax.set_aspect("equal") for ax in axd["axs"]]


def print_ray_table(ray_table):
    print('surface ray x y z vx vy vz')
    for s in range(ray_table.shape[0]):
        for r in range(ray_table.shape[1]):
            print(s, r, ray_table[s, r, 0, :], ray_table[s, r, 1, :])


def ray_table_from_fields_vectors(geometry, field_points, vectors):
    ray_table = np.empty([len(geometry), len(field_points)*len(vectors), 2, 3])
    for i in range(len(field_points)):
        for j in range(len(vectors)):
            k = j + i*len(vectors)
            ray_table[0, k, 0, :] = field_points[i]
            ray_table[0, k, 1, :] = vectors[j]
    return ray_table


def ray_table_from_fields_points(geometry, field_points, pupil_points):
    ray_table = np.empty(
        [len(geometry), len(field_points)*len(pupil_points), 2, 3])
    for i in range(len(field_points)):
        for j in range(len(pupil_points)):
            k = j + i*len(pupil_points)
            ray_table[0, k, 0, :] = field_points[i]
            ray_table[0, k, 1, :] = normalize_vec(
                pupil_points[j] - field_points[i])
    return ray_table


def plot_rays(axd, geometry, ray_table, ray_properties, **kwargs):
    keep = [s['draw_radius'] != 0 for s in geometry]

    for i in range(ray_table.shape[1]):
        if "x" in axd:
            axd["axs"][axd["x"]].plot(
                ray_table[:, i, 0, 2][keep], ray_table[:, i, 0, 0][keep], **kwargs)
        if "y" in axd:
            axd["axs"][axd["y"]].plot(
                ray_table[:, i, 0, 2][keep], ray_table[:, i, 0, 1][keep], **kwargs)


def default_ray_properties(ray_table):
    return {'wavl': np.full((ray_table.shape[1]), 633e-6)}
