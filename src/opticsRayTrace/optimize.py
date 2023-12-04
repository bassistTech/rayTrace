import numpy as np
from scipy.optimize import minimize
import opticsRayTrace.rayTraceTools as rtt


def ray_aiming(geometry, wavelengths, infinite):
    def merit(aim_consts):
        ray_table = rtt.new_ray_table(geometry, fields, pupils,
                                      wavelengths, aim_consts=aim_consts, infinite=True)
        rtt.trace_rays(ray_table, geometry)

        # Compute the merit function based on the squares of X positions
        merit = np.average(ray_table[stop_surface, :, 0, 0]**2) + \
            np.average(ray_table[stop_surface, :, 0, 1]**2)
        return merit
    stop_surface = rtt.find_stop(geometry)
    fields = np.random.normal(size=(10, 3))*0.01
    fields[:, 2] = 1
    pupils = [[0, 0, 0]]
    return minimize(merit, np.zeros([6]))['x']
