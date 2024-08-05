import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from datetime import datetime
from sunpy.coordinates import frames
import sunpy.map
import eispac
import eismaps.utils.find
from matplotlib.colors import LogNorm
from tqdm import tqdm
from sunpy.coordinates.frames import Helioprojective
import sunpy.sun.constants
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.coordinates import Helioprojective, propagate_with_solar_surface, RotatedSunFrame
from sunpy.coordinates import SphericalScreen

def safe_load_map(map_file):
    try:
        map = sunpy.map.Map(map_file)
    except:
        print(f"Error loading map {map_file}. Skipping.")
        map = None
    return map

def make_helioprojective_map(map_files, save_dir, wavelength, measurement, overlap, preserve_limb=None, save_fit=False, save_plot=False, plot_ext='png', plot_dpi=300):
    """
    Make a helioprojective full disk map from a list of maps.
    """

    first_map = safe_load_map(map_files[0])
    if first_map is None:
        return

    # fd_size = 0  # Determine the full disk grid size
    # for file in map_files:
    #     map = sunpy.map.Map(file)
    #     fd_size = max(fd_size, (abs(map.meta['xcen']) + (map.meta['fovx'] / 2)) * 2)
    #     fd_size = max(fd_size, (abs(map.meta['ycen']) + (map.meta['fovy'] / 2)) * 2)
    # fd_size = int(fd_size * 1.05)  # Add some padding
    fd_size = 3500  # Hardcoded

    map_dx = first_map.meta['cdelt1']  # Pixel sizes for the full disk image
    map_dy = first_map.meta['cdelt2']
    fd_width = round(fd_size / map_dx)
    fd_height = round(fd_size / map_dy)

    fd_header = {
        'cdelt1': map_dx,
        'cdelt2': map_dy,
        'crpix1': fd_width / 2,
        'crpix2': fd_height / 2,
        'crval1': 0,
        'crval2': 0,
        'cunit1': 'arcsec',
        'cunit2': 'arcsec',
        'ctype1': 'HPLN-TAN',
        'ctype2': 'HPLT-TAN',
        'date-obs': first_map.meta['date_obs'],
        'dsun_obs': first_map.meta['dsun_obs'],
        'hgln_obs': first_map.meta['hgln_obs'],
        'hglt_obs': first_map.meta['hglt_obs'],
    }

    fd_data = np.full((fd_height, fd_width), np.nan)
    fd_map = sunpy.map.Map(fd_data, fd_header)

    overlap_mask = np.zeros((fd_height, fd_width))
    combined_data = np.full((fd_height, fd_width), np.nan)

    for map_file in map_files:

        map = sunpy.map.Map(map_file)

        if preserve_limb=='drag':  # Manually drag rasters (and don't distort) with a custom point

            def differental_rotate_map_by_drag(map, point):
                # Calculate the time difference between the map and the first map
                map_time = datetime.strptime(map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f')
                first_map_time = datetime.strptime(first_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f')
                duration = map_time - first_map_time
                # Rotate the point by the differential rotation
                diffrot_point = SkyCoord(RotatedSunFrame(base=point, duration=duration))
                transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)
                # Calculate the difference between the original and transformed points
                shift_x = transformed_diffrot_point.Tx - point.Tx
                shift_y = transformed_diffrot_point.Ty - point.Ty
                # Shift the map by the difference
                map = map.shift_reference_coord(shift_x, shift_y)
                return map

            if sunpy.map.is_all_off_disk(map):  # This map is completely off disk, and so can't be shifted
                pass
            elif sunpy.map.contains_limb(map):  # Some of the map is off limb, so need to carefully choose the pixels to drag
                pixel_coords = sunpy.map.all_coordinates_from_map(map)
                distances = np.sqrt(pixel_coords.Tx ** 2 + pixel_coords.Ty ** 2)
                closest = np.argmin(distances)
                closest_coords = pixel_coords[closest]
                map = differental_rotate_map_by_drag(map, closest_coords)
            else: # This map is completely on disk, so can be shifted as normal
                map = differental_rotate_map_by_drag(map, map.center)

        elif preserve_limb=='spherical_screen':

            # Apply differential rotation to the map, compared to the time of the fd_map
            with propagate_with_solar_surface():
                with SphericalScreen(map.observer_coordinate, only_off_disk=True):
                    map = map.reproject_to(fd_map.wcs, algorithm='exact')

        else:

            # Apply differential rotation to the map, compared to the time of the fd_map
            with propagate_with_solar_surface():
                map = map.reproject_to(fd_map.wcs, algorithm='exact')

            if overlap == 'max':
                combined_data = np.where(np.isnan(combined_data), map.data, np.nanmax([combined_data, map.data], axis=0))
            elif overlap == 'mean':
                combined_data = np.where(np.isnan(combined_data), map.data, np.nansum([combined_data, map.data], axis=0))
            elif overlap == 'nan':
                combined_data = np.where(np.isnan(combined_data), map.data, np.nansum([combined_data, map.data], axis=0))

            overlap_mask = np.where(np.isnan(map.data), overlap_mask, overlap_mask + 1)

    if overlap == 'mean':
        fd_map = sunpy.map.Map(combined_data / overlap_mask, fd_map.meta)
    elif overlap == 'mask':
        fd_map = sunpy.map.Map(overlap_mask, fd_map.meta)
    elif overlap == 'nan':
        combined_data = np.where(overlap_mask > 1, np.nan, combined_data)
        fd_map = sunpy.map.Map(combined_data, fd_map.meta)
    else:
        fd_map = sunpy.map.Map(combined_data, fd_map.meta)

    if save_fit:

        map_file_datetime = os.path.basename(map_files[0]).split('.')[0].replace('eis_', '')
        fd_map.save(os.path.join(save_dir, f"eis_{map_file_datetime}.{wavelength}.{measurement}.fd.fits"), overwrite=True)

    if save_plot:
        
        fig = plt.figure()
        ax = plt.subplot(projection=fd_map)
        
        if measurement == 'int':
            fd_map.plot_settings['norm'] = LogNorm(vmin=1e1, vmax=5e3)
            im = fd_map.plot(cmap='gist_heat')
            fd_map.draw_limb(axes=ax, color="k")
        elif measurement == 'vel':
            im = fd_map.plot(cmap='RdBu_r')
            im.set_norm(plt.Normalize(vmin=-10, vmax=10))
        elif measurement == 'ntv':
            im = fd_map.plot(cmap='inferno')
            im.set_norm(plt.Normalize(vmin=0, vmax=40))
        elif measurement == 'mag':
            im = fd_map.plot(cmap='gray')
            im.set_norm(plt.Normalize(vmin=-300, vmax=300))
        elif measurement == 'fip':
            im = fd_map.plot(cmap='CMRmap')
            im.set_norm(plt.Normalize(vmin=0, vmax=3))
        else:
            raise ValueError('Error: measurement is not valid')
        
        im = ax.get_images()
        im_lims = im[0].get_extent()
        ax.set_aspect(abs((im_lims[1]-im_lims[0])/(im_lims[3]-im_lims[2])))
        plt.colorbar(extend='both')
        plt.savefig(os.path.join(save_dir, f"eis_{map_file_datetime}.{wavelength}.{measurement}.fd_hp.{plot_ext}"), dpi=plot_dpi)
        plt.close()

    return fd_map

def make_carrington_map(map_files, save_dir, wavelength, measurement, overlap, deg_per_pix=0.1, save_fit=False, save_plot=False, plot_ext='png', plot_dpi=300):
    """
    Make a Carrington full disk map from a list of maps.
    """

    first_map = safe_load_map(map_files[0])
    if first_map is None:
        return
    
    lon_pixels = int(360 / deg_per_pix)
    lat_pixels = int(180 / deg_per_pix)
    fd_lon = np.linspace(0, 360, lon_pixels) * u.deg
    fd_lat = np.linspace(-90, 90, lat_pixels) * u.deg
    fd_lon, fd_lat = np.meshgrid(fd_lon, fd_lat)

    fd_header = {
        'cunit1': 'deg',
        'cunit2': 'deg',
        'crpix1': lon_pixels/2,
        'crpix2': lat_pixels/2,
        'crval1': 0,
        'crval2': 0,
        'cdelt1': deg_per_pix,
        'cdelt2': deg_per_pix,
        'ctype1': 'CRLN-CEA',
        'ctype2': 'CRLT-CEA',
        'date_obs': first_map.meta['date_obs'],
        'hgln_obs': first_map.meta['hgln_obs'],
        'dsun_obs': first_map.meta['dsun_obs'],
        'hglt_obs': first_map.meta['hglt_obs'],
    }

    fd_data = np.full((lat_pixels, lon_pixels), np.nan)
    fd_map = sunpy.map.Map(fd_data, fd_header)

    overlap_mask = np.zeros((lat_pixels, lon_pixels))

    for map_file in map_files:

        map = sunpy.map.Map(map_file)

        # Change the time of the fd_map to the same as the time of the raster
        fd_map_temp = fd_map
        fd_map_temp.meta['date_obs'] = map.meta['date_obs']

        # Create a WCS object for the target map
        target_wcs = fd_map_temp.wcs

        # Reproject the raster map to the Carrington map
        map_carrington = map.reproject_to(target_wcs)

        if overlap == 'max':
            combined_data = np.where(np.isnan(fd_map.data), map_carrington.data, np.nanmax([fd_map.data, map_carrington.data], axis=0))
        elif overlap == 'mean':
            combined_data = np.where(np.isnan(fd_map.data), map_carrington.data, (fd_map.data + map_carrington.data))
            overlap_mask = np.where(np.isnan(fd_map.data), overlap_mask + 1, overlap_mask)
        elif overlap == 'nan':
            combined_data = np.where(np.isnan(fd_map.data), map_carrington.data, np.nan)

        fd_map = sunpy.map.Map(combined_data, fd_map.meta)

    if overlap == 'mean':
        fd_map.data = fd_map.data / overlap_mask

    if save_fit:

        map_file_datetime = os.path.basename(map_files[0]).split('.')[0].replace('eis_', '')
        fd_map.save(os.path.join(save_dir, f"eis_{map_file_datetime}.{wavelength}.{measurement}.fd.fits"), overwrite=True)

    if save_plot:

        fig = plt.figure()
        ax = plt.subplot(projection=fd_map)

        if measurement == 'int':
            fd_map.plot_settings['norm'] = LogNorm(vmin=1e1, vmax=5e3)
            im = fd_map.plot(cmap='gist_heat')
            fd_map.draw_limb(axes=ax, color="k")
        elif measurement == 'vel':
            im = fd_map.plot(cmap='RdBu_r')
            im.set_norm(plt.Normalize(vmin=-10, vmax=10))
        elif measurement == 'ntv':
            im = fd_map.plot(cmap='inferno')
            im.set_norm(plt.Normalize(vmin=0, vmax=40))
        elif measurement == 'mag':
            im = fd_map.plot(cmap='gray')
            im.set_norm(plt.Normalize(vmin=-300, vmax=300))
        elif measurement == 'fip':
            im = fd_map.plot(cmap='CMRmap')
            im.set_norm(plt.Normalize(vmin=0, vmax=3))
        else:
            raise ValueError('Error: measurement is not valid')

        im = ax.get_images()
        im_lims = im[0].get_extent()
        ax.set_aspect(abs((im_lims[1]-im_lims[0])/(im_lims[3]-im_lims[2])))
        plt.colorbar(extend='both')
        plt.savefig(os.path.join(save_dir, f"eis_{map_file_datetime}.{wavelength}.{measurement}.fd_ca.{plot_ext}"), dpi=plot_dpi)
        plt.close()