import numpy as np

import opticsRayTrace.rayTraceTools as rtt

def rms_by_field_xy(ray_table, index = -1):
    '''
    Compute RMS spot size as a function of field number for X ad Y

    index = the surface number of interest, or -1 for the last surface

    Returns: rms[field_number, axis]
    '''
    
    rayx = ray_table[index, :, 0, 0]
    rayy = ray_table[index, :, 0, 1]
    rayfields = ray_table[-1, :, 3, 0].astype(int)
    return np.array([[np.std(rayx[rayfields == f]), np.std(rayy[rayfields == f])] 
                     for f in np.arange(np.max(rayfields) + 1)])

def rms_by_field_radial(ray_table, index = -1):
    '''
    Computes RMS spot size as a function of field number for R

    Returns: rms[field_number]
    '''

    rms_xy = rms_by_field_xy(ray_table, index)
    return np.sqrt(np.average(rms_xy**2, axis = 1))

def average_by_field_xy(ray_table, index = -1):
    '''
    Compute RMS spot size as a function of field number for X ad Y

    index = the surface number of interest, or -1 for the last surface

    Returns: rms[field_number, axis]
    '''
    
    rayx = ray_table[index, :, 0, 0]
    rayy = ray_table[index, :, 0, 1]
    rayfields = ray_table[-1, :, 3, 0].astype(int)
    return np.array([[np.average(rayx[rayfields == f]), np.average(rayy[rayfields == f])] 
                     for f in np.arange(np.max(rayfields) + 1)])

def focal_length(geometry, beam_diameter_mm, wavelength_mm):
    test_angle = 0.01 # in radians
    test_rays = rtt.ray_table_fields_rings(geometry, test_angle*180/np.pi, 2, 
                                           beam_diameter_mm, 3, 
                                           [wavelength_mm])
    rtt.trace_rays(test_rays, geometry)
    spot_positions = average_by_field_xy(test_rays)
    test_position = spot_positions[1, 0]
    focal_length = test_position/test_angle
    return focal_length

print('analyze loaded')