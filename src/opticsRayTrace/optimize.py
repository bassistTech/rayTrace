import numpy as np
from scipy.optimize import minimize
import opticsRayTrace.rayTraceTools as rtt


def ray_aiming(geometry, wavelengths, infinite):
    '''
    TODO: Need a better choice of default fields, that are not random.
    '''
    def merit(aim_consts):
        ray_table = rtt.new_ray_table(geometry, fields, pupils,
                                      wavelengths, aim_consts=aim_consts, infinite=True)
        rtt.trace_rays(ray_table, geometry)

        # Compute the merit function based on the squares of X positions
        merit = np.average(ray_table[stop_surface, :, 4, 0]**2) + \
            np.average(ray_table[stop_surface, :, 4, 1]**2)
        return merit
    stop_surface = rtt.find_stop(geometry)
    fields = np.random.normal(size=(10, 3))*0.01
    fields[:, 2] = 1
    pupils = [[0, 0, 0]]
    return minimize(merit, np.zeros([6]))['x']

def get_variable_list(surface_list):
    '''
    Search the surface list for variables, and create a variable list.

    returns [{'surf': surface_number, 'var': parameter_name, 'val': value}]
    '''
    result = []
    for i in range(len(surface_list)):
        surf = surface_list[i]
        if 'variables' in surf.keys():
            for var in surf['variables']:
                result.append({'surf': i, 'var': var, 'val': surf[var]})
    return result

def get_values(variable_list):
    '''
    Return the values of variables in the design, as a straight Numpy array
    '''
    
    return np.array([d['val'] for d in variable_list])

def update_values(surface_list, params):
    '''
    Take a straight numpy array of values and insert them back into the surface
    list. This assumes the surface list has the same structure as the parameter
    list, i.e., you haven't done anything really funky with the surface list.
    '''
    
    variable_list = get_variable_list(surface_list)
    for i in range(len(variable_list)):
        var = variable_list[i]
        surf = surface_list[var['surf']]
        surf[var['var']] = params[i]

def merit_prepare(params, ray_table, surface_list):
    '''
    These are the preliminaries that a merit function always needs before
    performing an analysis of the ray table.

    1. Stuff the parameter values back into the surface list
    2. New parameter values require computing a new geometry
    3. Finally, trace rays through the new geometry

    The ray table and surface list are modified in place.
    '''

    update_values(surface_list, params)
    geometry = rtt.build_geometry(surface_list)
    rtt.trace_rays(ray_table, geometry)
    return geometry

print('optimize loaded')