'''
Optical ray tracing formulas
Francis Deck 10-4-2023
MIT License
'''

import numpy as np
import sympy as sp
from time import time
import os
import yaml
import matplotlib.pyplot as plt
import pickle


epsilon = 1e-4
maxiter = 10


def norm_vec(v):
    """
    Norm of a vector
    """
    return np.sqrt(np.dot(v, v))


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

#r, c, k = sp.symbols('r c k')

#sag_conic_sp = c*r**2/(1 + sp.sqrt(1 - (1 + k)*c**2*r**2))
#sag_conic_diff_sp = sp.simplify(sp.diff(sag_conic_sp, r)/r)

#sag_conic = sp.lambdify([r, c, k], sag_conic_sp, 'numpy')
#sag_conic_diff = sp.lambdify([r, c, k], sag_conic_diff_sp, 'numpy')

def sag_conic(r, c, k):
    return c*r**2/(np.sqrt(-c**2*r**2*(k + 1) + 1) + 1)

def sag_conic_diff(r, c, k):
    u = np.sqrt(-c**2*r**2*(k + 1) + 1)
    return c*(c**2*r**2*(k + 1) + 2*u*(u + 1))/(u*(u + 1)**2)


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
    v_line = [direction vector of line]

    https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    """

    dist = (np.sum((x_surf-x_linea) * n_surf, axis=1)
            / np.sum(v_line*n_surf, axis=1))[:, np.newaxis]
    return x_linea + v_line * dist


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
    For convenience and similarity to commercial programs, the surface list specifies the
    positions of the elements in _relative_ terms. This function gives each surface a
    position and orientation within the global coordinate space.

    It also applies defaults to the surfaces, mainly just to keep my programs from blowing
    up if a property is omitted.
    """

    def apply_defaults(d, defaults):
        return d | {k: defaults[k] for k in defaults if k not in d}

    prop_cycle = plt.rcParams['axes.prop_cycle']
    default_colors = prop_cycle.by_key()['color']
    default_catalog = "schott"

    geometry_defaults = {"dist": 0,
                         "draw": True,
                         "draw_radius": 10,
                         "draw_dx": 0,
                         "draw_dy": 0,
                         "mirror": False,
                         "catalog": default_catalog,
                         "material": "vacuum",
                         "c": 0, "k": 0,
                         "R": None}

    # This is the coordinate basis that will be cumulatively updated at each surface
    origin = np.array((0, 0, 0))
    x_axis = np.array((1, 0, 0))
    y_axis = np.array((0, 1, 0))
    z_axis = np.array((0, 0, 1))

    # Run through the surfaces, applying the coordinate basis to each surface
    geometry = []
    for i in range(len(surface_list)):
        s = surface_list[i]
        geometry_defaults['color'] = default_colors[i % len(default_colors)]
        g = apply_defaults(s, geometry_defaults)
        default_catalog = g["catalog"]
        if g["surf"] == "conic":
            if g["R"] is not None:
                g['c'] = 1/g["R"]
        if g['material'] == 'vacuum':
            g['glass_coeffs'] = np.array((0, 0, 0, 0, 0, 0))
        else:
            # For now, just use the 6 coefficient model, which only works for visible wavelengths.
            # g['glass_coeffs'] = load_glass(g['material'])['coeffs']
            g['glass_coeffs'] = glass_list[g["catalog"]][g['material']]['coeffs']
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


help_ray_table = """
    Help for ray table. The table is structured like this:

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


def find_stop(geometry):
    '''
    Given a geometry, find the surface labeled as "stop"
    '''

    istop = [i for i in range(len(geometry)) if 'stop' in geometry[i]]
    if len(istop) != 1:
        print('Must have exactly one stop surface, returning surface 0')
        return 0
    return istop[0]


def pupil_shift(field, aim_consts):
    """
    I accomplish basic ray aiming by applying a linear transformation to each pupil point.
    The apparent "mixing" of the 0 and 1 (x and y) coordinates allows for ray aiming through
    systems that have coordinate transformations. The 3 (Z) coordinate is left alone.

    A single set of constants is applied across the entire ray table, and these constants
    also become the optimization parameters for minimizing the centroid of the ray bundle
    through the pupil surface.

    I don't know how the commercial programs do ray aiming, but I needed to have *something*
    to support basic ray tracing and optimization. This is the first thing I got working

    field = The field point to be aimed
    aim_consts = 6 coefficients of a linear transformation applied to this field point
    returns: A new field point, transformed by aim_consts
    """
    a, b, c, d, e, f = aim_consts
    return np.array([a + b*field[0] + c*field[1],
                     d + e*field[0] + f*field[1],
                     0])


def new_ray_table(geometry, fields, pupils, wavls, infinite,
                  aim_consts=None):
    """
    geometry = list generated by build_geometry
    fields: list of vectors specifying direction of ray at starting surface
    pupils: list of points locating the pupils relative to starting surface
    wavls: list of wavelengths in mm (sorry)
    infinite: True if fields are vectors instead of positions, for infinite conjugate system
    aim_consts: Ray aiming constants, more about this when I get it working
    """
    nrays = len(fields)*len(pupils)*len(wavls)
    ray_table = np.zeros([len(geometry), nrays, 4, 3])

    if aim_consts is None:
        aim_consts = np.zeros((6))

    '''
    Work through every possible combination of a field, pupil, and wavelength
    '''

    i = 0

    for j in range(len(fields)):
        if infinite:
            fields_norm_j = normalize_vec(fields[j])
            pupil_shift_fields_j = pupil_shift(fields[j], aim_consts)
        else:
            pupil_shift_fields_j = pupil_shift(fields[j], aim_consts) - fields[j]
        for k in range(len(pupils)):
            if infinite:
                '''
                When the light source is considered to be at infinity, fields are vectors,
                and pupils are points on the entrance surface, optionally modified by
                ray aiming.
                '''
                ray_table[0, i, 0, :] = pupils[k] + pupil_shift_fields_j
                ray_table[0, i, 1, :] = fields_norm_j
            else:
                '''
                When the light source is at a finite distance, fields are points on the
                source surface, and pupils are points in space.
                '''
                ray_table[0, i, 0, :] = fields[j]
                ray_table[0, i, 1, :] = normalize_vec(pupils[k]
                                                        + pupil_shift_fields_j)
            # for now, only record the wavelength, index is computed below
            for m in range(len(wavls)):
                ray_table[0, i, 2, 1] = wavls[m]
                ray_table[0, i, 3, :] = np.array([j, k, m])
                i = i + 1

    # Build the refractive index table

    for i in range(ray_table.shape[0]):
        if i > 0:
            ray_table[i, :, 2, 1] = ray_table[i-1, :, 2, 1]
        ray_table[i, :, 2, 0] = glass_index(geometry[i]["glass_coeffs"],
                                            ray_table[i, :, 2, 1])

    return ray_table

default_aim_consts = np.zeros([6], dtype = float)

def ray_table_fields_rings(geometry, 
                           max_field_degrees, number_of_fields, 
                           pupil_diameter_mm, number_of_rings, 
                           wavls, aim_consts = default_aim_consts):
    """
    Create a basic circular ray bundle. This is the most common type of ray bundle used for
    general purpose ray tracing, so I created a function for it. 

    geometry = geometry object representing your design
    max_field_degrees = maximum field angle. Fields start from 0 to this angle.
    number_of_fields = number of fields represented, which can be 1
    pupil_diameter_mm = diameter of circular ray bundle in mm
    number_of_rings = number of circular "rings" of rays
    wavls = list of wavelengths, in mm, e.g., 633 nm is 0.000633 mm

    Returns: Initialized ray table
    """
    fields = np.zeros((number_of_fields, 3))
    fields[:, 0] = np.sqrt(np.linspace(0, max_field_degrees**2, number_of_fields))*np.pi/180
    fields[:, 2] = np.ones(number_of_fields)
    field_vectors = fields/np.sum(fields**2, axis = 1)[:, np.newaxis]

    r_ring = np.linspace(0, pupil_diameter_mm/2, number_of_rings)
    pupil_points = [[0, 0, 0]]
    for i in range(1, number_of_rings):
        r = r_ring[i]
        for j in range(i*6):
            theta = j*2*np.pi/i/6
            pupil_points.append([r*np.sin(theta), r*np.cos(theta), 0])

    infinite = True

    return new_ray_table(geometry, field_vectors, pupil_points, wavls, infinite, aim_consts)

def ray_table_fields_1d(geometry, 
                           max_field_degrees, number_of_fields, 
                           pupil_diameter_mm, number_of_pupils, 
                           wavls, aim_consts = default_aim_consts):
    """
    Create a basic circular ray bundle. This is the most common type of ray bundle used for
    general purpose ray tracing, so I created a function for it. 

    geometry = geometry object representing your design
    max_field_degrees = maximum field angle. Fields start from 0 to this angle.
    number_of_fields = number of fields represented, which can be 1
    pupil_diameter_mm = diameter of circular ray bundle in mm
    number_of_rings = number of circular "rings" of rays
    wavls = list of wavelengths, in mm, e.g., 633 nm is 0.000633 mm

    Returns: Initialized ray table
    """
    fields = np.zeros((number_of_fields, 3))
    fields[:, 0] = np.sqrt(np.linspace(0, max_field_degrees**2, number_of_fields))*np.pi/180
    fields[:, 2] = np.ones(number_of_fields)
    field_vectors = fields/np.sum(fields**2, axis = 1)[:, np.newaxis]

    x_pupil = np.linspace(-pupil_diameter_mm/2, pupil_diameter_mm/2, number_of_pupils)

    pupil_points = np.zeros([number_of_pupils, 3])
    pupil_points[:, 0] = x_pupil

    infinite = True

    return new_ray_table(geometry, field_vectors, pupil_points, wavls, infinite, aim_consts)

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
    t0 = time()
    surface = geometry[surf]
    # first step is intersecting the input ray with the curved surface:
    # sag = np.zeros([ray_table.shape[1]])
    sag = np.zeros([ray_table.shape[1]])
    for k in range(40):
        ray_table[surf, :, 0, :] = intersect_line_plane_array(
            surface["origin"] + np.outer(sag, surface["z_axis"]), 
            surface["z_axis"],
            ray_table[surf-1, :, 0, :], # position of previous ray
            ray_table[surf-1, :, 1, :]) # direction of previous ray
        r = distance_point_line_array(
            surface["origin"],
            surface["z_axis"],
            ray_table[surf, :, 0, :])
        sagnew = sag_conic(r, surface["c"], surface["k"])
        sigma = np.max(sagnew - sag)
        sag = sagnew
        if sigma < epsilon:
            break
    times = [time() - t0]
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
    times.append(time() - t0 - times[-1])
    return times

def propagate_cylindrical_surface(ray_table, geometry, surf):
    """
    ray_table = [surf, ray, type, p (type=0), v (type=1), [n w _] (type=2)]
    """
    t0 = time()
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
    times = [time() - t0]

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
    return times


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


def propagate_paraxial(ray_table, geometry, surf):
    surface = geometry[surf]

    # Bring the incoming ray to this surface

    ray_table[surf, :, 0, :] = intersect_line_plane_array(surface["origin"],
                                                          surface["z_axis"],
                                                          ray_table[surf - 1, :, 0, :],
                                                          ray_table[surf-1, :, 1, :])

    # Get directions of chief rays from previous surface

    chief_ray_vectors = ray_table[surf-1, :, 1, :]

    # All chief rays pass through origin of surface, by definition

    chief_ray_points = np.zeros_like(chief_ray_vectors)

    # Find where each chief ray intersects focal plane

    focal_points = intersect_line_plane_array(surface['origin'] + surface["z_axis"]*surface["f"],
                                              surface["z_axis"],
                                              chief_ray_points,
                                              chief_ray_vectors)

    # point each outgoing ray to its corresponding focal point

    ray_table[surf, :, 1, :] = normalize_vec_array(focal_points - ray_table[surf, :, 0, :])


def trace_rays(ray_table, geometry):
    """
    Propagate rays through the system. This is the main ray tracing function.

    The outer loop iterates through the geometry, rather than through the rays, for speed,
    since each propagate function acts on all of the rays at once for a given surface type.
    """

    # No tracing to first surface, so we start with index 1
    
    for i in range(1, len(geometry)):
        times = []
        match geometry[i]["surf"]:
            case ("rotate" | "shift" | "dummy"):
                propagate_dummy_surface(ray_table, geometry, i)
            case "conic":
                times = propagate_conic_surface(ray_table, geometry, i)
            case "cylindrical":
                propagate_cylindrical_surface(ray_table, geometry, i)
            case "plane grating":
                propagate_plane_grating(ray_table, geometry, i)
            case "plane grating 2":
                propagate_plane_grating_2(ray_table, geometry, i)
            case "paraxial":
                propagate_paraxial(ray_table, geometry, i)
            case _:
                print("Did not recognize surface type", geometry[i]["surf"])

        # Propagate ray origin data through the entire table
        ray_table[i, :, 3, :] = ray_table[i - 1, :, 3, :]
        ray_table[i, :, 2, 1] = ray_table[i - 1, :, 2, 1]
        geometry[i]['times'] = times


def print_ray_table(ray_table):
    print('surface ray x y z vx vy vz n w')
    for s in range(ray_table.shape[0]):
        for r in range(ray_table.shape[1]):
            print(s, r, ray_table[s, r, 0, :],
                  ray_table[s, r, 1, :],
                  ray_table[s, r, 2, 0:2])

def create_glass_file():
    """
    The glasses take a long time to load from their original data files, so I turned
    them into a pickle file. I doubt that the glasses are likely to change, and they
    certainly won't unless you download new ones. So it's probably safe to run this
    function only if you know that the data files have changed for some reason,
    or you want to add more catalogs to my list.

    Glass data is from this site:

    M. N. Polyanskiy, "Refractive index database," https://refractiveindex.info. Accessed on 2023-11-07.
    """

    default_glass_root = os.path.dirname(os.path.abspath(__file__)) + '/database/data-nk/'
    # default_glass_catalog = 'glass/schott/'

    glass_list = {}
    root_path = os.path.join(default_glass_root, 'glass')
    glass_catalogs = ['schott', 'hikari', 'hoya', 'ohara', 'sumita']
    for catalog in glass_catalogs:
        print(catalog)
        glass_list[catalog] = {}
        path = os.path.join(root_path, catalog)
        for filename in os.listdir(path):
            if filename.endswith('.yml'):
                glass = filename.split('.')[-2]
                data = yaml.load(open(os.path.join(path, filename), 'r'), yaml.Loader)
                coeffs = list(map(float, data['DATA'][0]['coefficients'].split()[1:]))
                glass_list[catalog][glass] = {'coeffs': coeffs}
    pickle.dump(glass_list, open('glass_list.pickle', 'wb'))

glass_list = pickle.load(open('glass_list.pickle', 'rb'))

def glass_index(coeffs, mm):
    microns = mm*1000
    b1, c1, b2, c2, b3, c3 = coeffs
    l2 = microns**2
    return np.sqrt(1 + b1*l2/(l2 - c1) + b2*l2/(l2 - c2) + b3*l2/(l2 - c3))

print('rayTraceTools loaded')
