'''
Optical ray tracing formulas
Francis Deck 10-4-2023
MIT License
'''

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import yaml
import os

epsilon = 1e-3


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
    npd1 = np.dot((x_surf - x_line), n_surf)
    npd2 = np.dot(v_line, n_surf)

    return (x_line + v_line * npd1/npd2)


def intersect_line_plane_array(x_surf, n_surf, x_linea, v_line):
    """
    intersection of line and plane
    x_surf = point on plane
    n_surf = normal to plane
    x_line = [point on line]
    v_line = direction vector of line

    https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    """

    nps = (np.sum((x_surf-x_linea) * n_surf, axis=1) / np.sum(v_line*n_surf, axis=1))[:, np.newaxis]
    return x_linea + v_line * nps


def rotate_matrix(axis, alpha_deg):
    '''
    generic 3d rotation matrix

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


def build_geometry(surface_list):
    """
    Turn surface list into a "geometry" which includes the absolute coordinates and
    orientations of the surfaces, and other properties including defaults to support
    ray tracing. Properties not dealt with here will be simply passed on from the
    surface list.
    """
    def apply_defaults(d, defaults):
        return d | {k: defaults[k] for k in defaults if k not in d}

    geometry_defaults = {"dist": 0,
                         "draw_radius": 10,
                         "draw_dx": 0,
                         "draw_dy": 0,
                         "mirror": False,
                         "material": "vacuum"}

    # This is the coordinate basis that will be cumulatively updated at each surface
    origin = np.array((0, 0, 0))
    x_axis = np.array((1, 0, 0))
    y_axis = np.array((0, 1, 0))
    z_axis = np.array((0, 0, 1))

    # Run through the surfaces, applying the coordinate basis to each surface
    geometry = []
    for s in surface_list:
        g = apply_defaults(s, geometry_defaults)
        if g['material'] == 'vacuum':
            g['glass_coeffs'] = np.array((0, 0, 0, 0, 0, 0))
        else:
            g['glass_coeffs'] = load_glass(g['material'])['coeffs']
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


def help_ray_table():
    """
    ray_table[surf #, ray #, type #, data]
    type = 0
        data = intersection point of ray in space
    type = 1
        data = direction vector of ray
    type = 2
        data = [refractive index, wavelength, <reserved>]
    type = 3
        data = [field #, pupil #, wavelength #]
    """


def fields(ray_table):
    """
    Helper function to return the field numbers of the rays in the ray table
    """
    return ray_table[0, :, 3, 0].astype(int)


def pupils(ray_table):
    """
    Helper function to return the pupil numbers of the rays in the ray table
    """
    return ray_table[0, :, 3, 1].astype(int)


def wavens(ray_table):
    """
    Helper function to return the wavelength numbers of the rays in the ray table
    """
    return ray_table[0, :, 3, 2].astype(int)


def new_ray_table(geometry, fields, pupils, wavls, infinite=False):
    """
    geometry = list generated by build_geometry
    field_points: list of points or angles (if "infinite") for fields
    pupil_points: list of points in space locating the pupils
    wavls: list of wavelengths in mm (sorry)
    infinite: True if fields are vectors instead of positions, for infinite conjugate system
    """
    nrays = len(fields)*len(pupils)*len(wavls)
    ray_table = np.zeros([len(geometry), nrays, 4, 3])
    mg = np.array(np.meshgrid(range(len(fields)), range(
        len(pupils)), range(len(wavls)))).T.reshape(-1, 3)
    # https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays

    # Build the starting ray table

    for i in range(ray_table.shape[1]):  # i = ray
        j, k, m = mg[i]  # j = field, k = pupil, m = wavelength
        if infinite:
            ray_table[0, i, 1, :] = normalize_vec(fields[j])
            ray_table[0, i, 0, :] = pupils[k]
        else:
            ray_table[0, i, 0, :] = fields[j]
            ray_table[0, i, 1, :] = normalize_vec(pupils[k] - fields[j])
        ray_table[0, i, 2, 1] = wavls[m]
        ray_table[0, i, 3, :] = mg[i]

    # Build the refractive index table

    for i in range(ray_table.shape[0]):
        if i > 0:
            ray_table[i, :, 2, 1] = ray_table[i-1, :, 2, 1]
        ray_table[i, :, 2, 0] = glass_index(geometry[i]["glass_coeffs"],
                                            ray_table[i, :, 2, 1])

    return ray_table


def propagate_dummy_surface(ray_table, geometry, surf):
    """
    ray_table = [surf, ray, type, p (type=0), v (type=1), [n w _] (type=2)]
    en = ray number of previous surface (entrance)
    ex = ray number of this surface (exit)
    """
    surface = geometry[surf]
    ray_table[surf, :, 0, :] = intersect_line_plane_array(surface["origin"],
                                                          surface["z_axis"],
                                                          ray_table[surf -
                                                                    1, :, 0, :],
                                                          ray_table[surf-1, :, 1, :])
    ray_table[surf, :, 1, :] = ray_table[surf-1, :, 1, :]


def propagate_conic_surface(ray_table, geometry, surf):
    """
    ray_table = [surf, ray, type, p (type=0), v (type=1), [n w _] (type=2)]
    """
    surface = geometry[surf]
    # first step is intersecting the input ray with the curved surface:
    sag = np.zeros([ray_table.shape[1]])
    # TODO: Handle zero curvature by bypassing this loop
    for k in range(40):
        ray_table[surf, :, 0, :] = intersect_line_plane_array(
            surface["origin"] +
            np.outer(sag, surface["z_axis"]), surface["z_axis"],
            ray_table[surf-1, :, 0, :], ray_table[surf-1, :, 1, :])
        r = distance_point_line_array(
            surface["origin"],
            surface["z_axis"],
            ray_table[surf, :, 0, :])
        sagnew = sag_conic(r, surface["c"], surface["k"])
        sigma = np.max(np.abs(sagnew - sag))
        sag = sagnew
        if sigma < epsilon:
            break

    # second step is finding surface normal on the exit surface
    sd = sag_conic_diff(r, surface['c'], surface['k'])
    x = np.dot(ray_table[surf, :, 0, :] - surface["origin"], surface["x_axis"])
    y = np.dot(ray_table[surf, :, 0, :] - surface["origin"], surface["y_axis"])
    n = normalize_vec_array(np.outer(sd*x, surface['x_axis'])
                            + np.outer(sd*y, surface['y_axis'])
                            - surface['z_axis'])

    # third step is finding the direction of the refracted ray
    # https://image.slidesharecdn.com/raytracing-111026083006-phpapp01/95/ray-tracing-23-728.jpg
    i = -ray_table[surf-1, :, 1, :]  # incident ray vector [vectors]
    ci = np.sum(i*n, axis=1)  # cosine of incident angle [scalars]
    if surface["mirror"]:
        eta = -np.ones(ray_table.shape[1])
    else:
        eta = ray_table[surf-1, :, 2, 0]/ray_table[surf, :, 2, 0]
    # eta = surface["n1"] / surface["n2"]  # scalar
    ct = np.sqrt(1 - eta**2 * (1 - ci**2))  # [scalars]
    vr = ((eta * ci - ct)[:, np.newaxis]*n
          - eta[:, np.newaxis]*i)  # [vectors]
    ray_table[surf, :, 1, :] = vr


def propagate_cylindrical_surface(ray_table, geometry, surf):
    """
    ray_table = [surf, ray, type, p (type=0), v (type=1), [n w _] (type=2)]
    """
    surface = geometry[surf]
    # first step is intersecting the input ray with the curved surface:
    sag = np.zeros([ray_table.shape[1]])
    # TODO: Handle zero curvature by bypassing this loop
    for k in range(40):
        ray_table[surf, :, 0, :] = intersect_line_plane_array(
            surface["origin"] +
            np.outer(sag, surface["z_axis"]), surface["z_axis"],
            ray_table[surf-1, :, 0, :], ray_table[surf-1, :, 1, :])
        ya = np.dot(ray_table[surf, :, 0, :] -
                    surface["origin"], surface["y_axis"])
        sagnew = sag_conic(ya, surface["c"], surface["k"])
        sigma = np.max(np.abs(sagnew - sag))
        sag = sagnew
        if sigma < epsilon:
            break

    # second step is finding surface normal on the exit surface
    sd = sag_conic_diff(ya, surface['c'], surface['k'])
    x = np.dot(ray_table[surf, :, 0, :] - surface["origin"], surface["x_axis"])
    y = np.dot(ray_table[surf, :, 0, :] - surface["origin"], surface["y_axis"])
    n = normalize_vec_array(0*np.outer(sd*x, surface['x_axis'])
                            + np.outer(sd*y, surface['y_axis'])
                            - surface['z_axis'])

    # third step is finding the direction of the refracted ray
    # https://image.slidesharecdn.com/raytracing-111026083006-phpapp01/95/ray-tracing-23-728.jpg
    i = -ray_table[surf-1, :, 1, :]  # incident ray vector [vectors]
    ci = np.sum(i*n, axis=1)  # cosine of incident angle [scalars]
    if surface["mirror"]:
        eta = -np.ones(ray_table.shape[1])
    else:
        eta = ray_table[surf-1, :, 2, 0]/ray_table[surf, :, 2, 0]
    # eta = surface["n1"] / surface["n2"]  # scalar
    ct = np.sqrt(1 - eta**2 * (1 - ci**2))  # [scalars]
    vr = ((eta * ci - ct)[:, np.newaxis]*n
          - eta[:, np.newaxis]*i)  # [vectors]
    ray_table[surf, :, 1, :] = vr


def propagate_plane_grating(ray_table, geometry, surf):
    """
    ray_table = [surf, ray, type, p (type=0), v (type=1), [n w _] (type=2)]

    https://www.phys.uconn.edu/~eyler/phys4150/R/OSLO%20Optics%20Reference.pdf
    """
    surface = geometry[surf]
    en = surf - 1
    ex = surf
    ray_table[surf, :, 0, :] = intersect_line_plane_array(surface["origin"],
                                                          surface["z_axis"],
                                                          ray_table[surf -
                                                                    1, :, 0, :],
                                                          ray_table[surf-1, :, 1, :])

    # three arrays of scalars in coordinates local to the surface

    wavl = ray_table[surf, :, 2, 1]

    if surface["mirror"]:
        n1 = 1
        n2 = -1
    else:
        n1 = ray_table[surf-1, :, 2, 0]
        n2 = ray_table[surf, :, 2, 0]

    exy = (
        n1 * np.dot(ray_table[en, :, 1, :], surface["y_axis"])
        + surface["m"] * wavl / surface["d"]
    ) / n2

    exx = n1 * np.dot(ray_table[en, :, 1, :],
                      surface["x_axis"]) / n2

    exz = np.sqrt(1 - exx*exx - exy*exy)

    # convert to a vector in global coordinates

    ray_table[ex, :, 1, :] = (np.outer(exx, surface["x_axis"])
                              + np.outer(exy, surface["y_axis"])
                              + np.outer(exz, surface["z_axis"]))


def propagate_plane_grating_2(ray_table, geometry, surf):
    """
    ray_table = [surf, ray, type, p (type=0), v (type=1), [n w _] (type=2)]

    https://www.phys.uconn.edu/~eyler/phys4150/R/OSLO%20Optics%20Reference.pdf
    """

    surface = geometry[surf]
    en = surf - 1
    ex = surf
    ray_table[surf, :, 0, :] = intersect_line_plane_array(surface["origin"],
                                                          surface["z_axis"],
                                                          ray_table[surf -
                                                                    1, :, 0, :],
                                                          ray_table[surf-1, :, 1, :])

    # three arrays of scalars in coordinates local to the surface

    wavl = ray_table[surf, :, 2, 1]

    if surface["mirror"]:
        n1 = 1
        n2 = -1
    else:
        n1 = ray_table[surf-1, :, 2, 0]
        n2 = ray_table[surf, :, 2, 0]

    exx = (
        n1 * np.dot(ray_table[en, :, 1, :], surface["x_axis"])
        + surface["m"] * wavl / surface["d"]
    ) / n2

    exy = n1 * np.dot(ray_table[en, :, 1, :],
                      surface["y_axis"]) / n2

    exz = np.sqrt(1 - exx*exx - exy*exy)

    # convert to a vector in global coordinates

    ray_table[ex, :, 1, :] = (np.outer(exx, surface["x_axis"])
                              + np.outer(exy, surface["y_axis"])
                              + np.outer(exz, surface["z_axis"]))


def propagate_ray(ray_table, geometry):
    """
    ray_table = [surf, ray, type, p (type=0), v (type=1), [n w _] (type=2)]
    """
    # No tracing to first surface, so we start with index 1
    for i in range(1, len(geometry)):
        g = geometry[i]
        if g["surf"] in ["rotate", "shift", "dummy"]:
            propagate_dummy_surface(ray_table, geometry, i)
        elif g["surf"] == "conic":
            propagate_conic_surface(ray_table, geometry, i)
        elif g["surf"] == "cylindrical":
            propagate_cylindrical_surface(ray_table, geometry, i)
        elif g["surf"] == "plane grating":
            propagate_plane_grating(ray_table, geometry, i)
        elif g["surf"] == "plane grating 2":
            propagate_plane_grating_2(ray_table, geometry, i)
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


def axdPlot(x, axd, label):
    # first graph shows the z-x coordinate axis in global space
    if "x" in axd:

        axd["axs"][axd["x"]].plot(
            x[:, 2], x[:, 0], linewidth=1, label=label)
        axd["axs"][axd["x"]].set_ylabel("X (mm)")

    # second graph shows the z-y coordinate axis in global space
    if "y" in axd:
        axd["axs"][axd["y"]].plot(
            x[:, 2], x[:, 1], linewidth=1, label=label)
        axd["axs"][axd["y"]].set_ylabel("Y (mm)")
    axd["axs"][-1].set_xlabel("Z (mm)")


def plot_faces(axd, surface_list):
    npts = 50
    x = np.empty((npts, 3))

    for s in surface_list:
        if "label" in s:
            label = s['label']
        else:
            label = None

        # create (for now) circular aperture in local coordinates
        theta_ring = np.linspace(0, 2*np.pi, npts)
        x_ring = s["draw_radius"]*np.sin(theta_ring)
        y_ring = s["draw_radius"]*np.cos(theta_ring)

        # create x or y bar in local coordinates
        bar = np.linspace(-s["draw_radius"], s["draw_radius"], npts)

        # generate drawing curves in local coordinates of object
        match s["surf"]:
            case ("dummy" | "plane grating" | "plane grating 2" | "rotate" | "shift"):
                sag_ring = np.zeros((npts))
                sag_xbar = np.zeros((npts))
                sag_ybar = np.zeros((npts))
            case "conic":
                r_ring = np.ones_like(x_ring)*s["draw_radius"]
                sag_ring = sag_conic(r_ring, s["c"], s["k"])
                sag_xbar = sag_conic(bar, s["c"], s["k"])
                sag_ybar = sag_xbar

            case "cylindrical":
                sag_ring = sag_conic(y_ring, s["c"], s["k"])
                sag_xbar = np.zeros((npts))
                sag_ybar = sag_conic(bar, s["c"], s["k"])

        for i in range(npts):
            x[i, :] = (s["origin"]
                       + x_ring[i]*s["x_axis"]
                       + y_ring[i]*s["y_axis"]
                       + sag_ring[i]*s["z_axis"])

        axdPlot(x, axd, label)

        for i in range(npts):
            x[i, :] = (s["origin"]
                       + bar[i]*s["x_axis"]
                       #+ y_ring[i]*s["y_axis"]
                       + sag_xbar[i]*s["z_axis"])

        axdPlot(x, axd, None)

        for i in range(npts):
            x[i, :] = (s["origin"]
                       + bar[i]*s["y_axis"]
                       #+ y_ring[i]*s["y_axis"]
                       + sag_ybar[i]*s["z_axis"])

        axdPlot(x, axd, None)

    [ax.set_aspect("equal") for ax in axd["axs"]]


def print_ray_table(ray_table):
    print('surface ray x y z vx vy vz n w')
    for s in range(ray_table.shape[0]):
        for r in range(ray_table.shape[1]):
            print(s, r, ray_table[s, r, 0, :],
                  ray_table[s, r, 1, :],
                  ray_table[s, r, 2, 0:2])


def plot_rays(axd, geometry, ray_table, **kwargs):
    keep = [s['draw_radius'] != 0 for s in geometry]

    for i in range(ray_table.shape[1]):
        if "x" in axd:
            axd["axs"][axd["x"]].plot(
                ray_table[:, i, 0, 2][keep], ray_table[:, i, 0, 0][keep], **kwargs)
        if "y" in axd:
            axd["axs"][axd["y"]].plot(
                ray_table[:, i, 0, 2][keep], ray_table[:, i, 0, 1][keep], **kwargs)


"""
Glass data is from this site:

M. N. Polyanskiy, "Refractive index database," https://refractiveindex.info. Accessed on 2023-11-07.
"""
default_glass_root = os.path.dirname(os.path.abspath(__file__)) + '/database/data-nk/'
default_glass_catalog = 'glass/schott/'


def load_glass(glass, root=default_glass_root, catalog=default_glass_catalog):
    if root is None:
        root = default_glass_root
    if catalog is None:
        catalog = default_glass_catalog
    filename = root + catalog + glass + '.yml'
    data = yaml.load(open(filename, 'r'), yaml.Loader)
    return {'coeffs': list(map(float, data['DATA'][0]['coefficients'].split()[1:]))}


def glass_index(coeffs, mm):
    microns = mm*1000
    b1, c1, b2, c2, b3, c3 = coeffs
    l2 = microns**2
    return np.sqrt(1 + b1*l2/(l2 - c1) + b2*l2/(l2 - c2) + b3*l2/(l2 - c3))
