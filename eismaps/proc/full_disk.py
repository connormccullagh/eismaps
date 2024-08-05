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

def make_helioprojective_map(map_files, save_dir, wavelength, measurement, overlap, apply_rotation=True, preserve_limb=True, save_fit=False, save_plot=False, plot_ext='png', plot_dpi=300):
    """
    Make a helioprojective full disk map from a list of maps.
    """

    first_map = safe_load_map(map_files[0])
    if first_map is None:
        return

    fd_size = 3500  # Hardcoded to avoid anomolous rasters generating incorrect huge full disk maps and crashing with memory errors

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

        if apply_rotation:

            with propagate_with_solar_surface(rotation_model='howard'):
                if preserve_limb:
                    with SphericalScreen(map.observer_coordinate, only_off_disk=True):
                        map = map.reproject_to(fd_map.wcs, algorithm='exact')
                else:
                    map = map.reproject_to(fd_map.wcs, algorithm='exact')

        else:

            if preserve_limb:
                with SphericalScreen(map.observer_coordinate, only_off_disk=True):
                    map = map.reproject_to(fd_map.wcs, algorithm='exact')
            else:
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
        fd_map.save(os.path.join(save_dir, f"eis_{map_file_datetime}.{wavelength}.{measurement}.fd_hp.fits"), overwrite=True)

    if save_plot:
        
        fig = plt.figure()
        ax = plt.subplot(projection=fd_map)
        
        if measurement == 'int':
            fd_map.plot_settings['norm'] = LogNorm(vmin=1e1, vmax=5e3)
            im = fd_map.plot(cmap='gist_heat')
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
            print(f"Error: plotting information for this measurement is not defined in eismaps. Full disk fits file was saved, but can't plot.")
            return
        
        fd_map.draw_limb(axes=ax, color="k")
        im = ax.get_images()
        im_lims = im[0].get_extent()
        ax.set_aspect(abs((im_lims[1]-im_lims[0])/(im_lims[3]-im_lims[2])))
        plt.colorbar(extend='both')
        plt.savefig(os.path.join(save_dir, f"eis_{map_file_datetime}.{wavelength}.{measurement}.fd_hp.{plot_ext}"), dpi=plot_dpi)
        plt.close()

    return fd_map

def make_carrington_map(map_files, save_dir, wavelength, measurement, overlap, apply_rotation=True, deg_per_pix=0.1, save_fit=False, save_plot=False, plot_ext='png', plot_dpi=300):
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
        if apply_rotation:
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
        fd_map.save(os.path.join(save_dir, f"eis_{map_file_datetime}.{wavelength}.{measurement}.fd_ca.fits"), overwrite=True)

    if save_plot:

        fig = plt.figure()
        ax = plt.subplot(projection=fd_map)

        if measurement == 'int':
            fd_map.plot_settings['norm'] = LogNorm(vmin=1e1, vmax=5e3)
            im = fd_map.plot(cmap='gist_heat')
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
            print(f"Error: plotting information for this measurement is not defined in eismaps. Full disk fits file was saved, but can't plot.")
            return

        im = ax.get_images()
        im_lims = im[0].get_extent()
        ax.set_aspect(abs((im_lims[1]-im_lims[0])/(im_lims[3]-im_lims[2])))
        plt.colorbar(extend='both')
        plt.savefig(os.path.join(save_dir, f"eis_{map_file_datetime}.{wavelength}.{measurement}.fd_ca.{plot_ext}"), dpi=plot_dpi)
        plt.close()

    return fd_map