import numpy as np
import sunpy.map

def double_field_of_view(map_):
    """
    Doubles the field of view in both x and y dimensions of a SunPy map, filling
    the new pixels with NaN values.

    Parameters
    ----------
    map_ : sunpy.map.Map
        The input SunPy map to be doubled in field of view.

    Returns
    -------
    sunpy.map.Map
        The new SunPy map with the doubled field of view.
    """

    # Get the dimensions of the map as numerical values
    x_orig = map_.meta['naxis1']
    y_orig = map_.meta['naxis2']

    # Calculate the new dimensions
    x_new = int(x_orig * 2)
    y_new = int(y_orig * 2)

    # Create a new array with NaN values
    data_new = np.full((y_new, x_new), np.nan)

    # Copy the data from the original map to the new array
    data_new[int(.25*x_new)-3:int(.75*x_new)+3, int(.25*y_new)-3:int(.75*y_new)+3] = map_.data

    # Update the meta information of the map
    map_.meta['naxis1'] = x_new
    map_.meta['naxis2'] = y_new
    # map_.meta['cdelt1'] = map_.meta['cdelt1'] / 2
    # map_.meta['cdelt2'] = map_.meta['cdelt2'] / 2
    map_.meta['crpix1'] = 1
    map_.meta['crpix2'] = 1
    map_.meta['crval1'] = -(x_new/2)*map_.meta['cdelt1']
    map_.meta['crval2'] = -(y_new/2)*map_.meta['cdelt2']

    # Create a new map with the new data and meta information
    map_new = sunpy.map.Map(data_new, map_.meta)

    return map_new