import glob
import os
import numpy as np
import sunpy.map
from sunpy.time import parse_time
import time
from sunpy.net import Fido, attrs as a
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import pathlib
from eismaps.utils.fov import double_field_of_view
from eismaps.utils.format import change_line_format
from eismaps.utils.logo import print_logo
import eispac


def create_magnetic_flux_maps(base_files, skip_done=False, resample=True):
    """
    Creates magnetic flux maps corresponding to EIS fit files.

    Parameters
    ----------
    base_files : list
        List of paths to EIS fit files.
    """

    for i, base_file in enumerate(base_files):
        
        file_dir = os.path.dirname(base_file)
        file_name = os.path.basename(base_file)
        save_dir = file_dir.replace('data_eis', 'data_hmi')
        save_name = file_name.replace('.fit.h5', '.fits')
        save_path = os.path.join(save_dir, save_name)

        # Check if both FITS and PNG files exist
        fits_file_exists = os.path.exists(save_path)
        png_file_exists = os.path.exists(save_path.replace('.fits', '.png'))

        # Skip processing if both files exist
        if skip_done and fits_file_exists and png_file_exists:
            print(f'Skipping {base_file} as both FITS and PNG files already exist.')
            continue

        print(f'Processing {base_file} ({i+1}/{len(base_files)})')

        fits = eispac.read_fit(base_file)
        map_int = fits.get_map(component=0, measurement='intensity')

        # time_range = a.Time(parse_time(map_int.meta['date_beg']), parse_time(map_int.meta['date_end']))
        # hmi_parameters = time_range & a.Instrument.hmi & a.Physobs.los_magnetic_field

        # hmi_search = Fido.search(hmi_parameters)
        # hmi_result = Fido.fetch(hmi_search[0][int(round(len(hmi_search[0]) / 2, 0))])

        download_success = False
        extend_count = 0
        fail_count = 0
        while not download_success:
            try:
                time_range = a.Time(parse_time(map_int.meta['date_beg']) - extend_count * u.hour, parse_time(map_int.meta['date_end']) + extend_count * u.hour)
                hmi_parameters = time_range & a.Instrument.hmi & a.Physobs.los_magnetic_field
                hmi_search = Fido.search(hmi_parameters)

                if hmi_search.file_num == 0:
                    if extend_count >= 10:
                        print(f"No HMI data found within extended time range for {base_file}. Skipping this file.")
                        break
                    extend_count += 2
                    print(f"No HMI data found, extending the search window by {extend_count} hours.")
                    continue

                hmi_result = Fido.fetch(hmi_search[0][int(len(hmi_search[0]) / 2)])
                map_hmi = sunpy.map.Map(hmi_result[0])
                download_success = True
            except Exception as e:
                fail_count += 1
                if fail_count >= 10:
                    print(f"Failed to download HMI data after 10 attempts for {base_file}. Skipping this file.")
                    break
                print(f"Download failed, retrying... (Attempt {fail_count}/10)")
                time.sleep(60)

        if not download_success:
            # give the user a warning
            print(f"WARNING!!! Failed to download HMI data for {base_file}. Skipping this file.")
            # sleep for 5 seconds
            time.sleep(5)
            continue  # Skip to the next file if data download was not successful

        map_hmi = sunpy.map.Map(hmi_result[0])

        # print(map_hmi.meta) ### WIP: used different meta for map creation
        # exit()

        map_hmi = map_hmi.rotate(order=3)
        map_hmi = double_field_of_view(map_hmi)

        if resample:
            bottom_left = SkyCoord(map_int.meta['crval1'] * u.arcsec, map_int.meta['crval2'] * u.arcsec, frame=map_hmi.coordinate_frame)
            top_right = SkyCoord((map_int.meta['crval1'] + map_int.meta['cdelt1'] * map_int.meta['naxis1']) * u.arcsec, (map_int.meta['crval2'] + map_int.meta['cdelt2'] * map_int.meta['naxis2']) * u.arcsec, frame=map_hmi.coordinate_frame)
            map_hmi = map_hmi.submap(bottom_left, top_right=top_right)
            map_hmi = map_hmi.resample(u.Quantity([map_int.meta['naxis1'], map_int.meta['naxis2']]) * u.pixel)

        pathlib.Path(save_dir).mkdir(parents=True, exist_ok=True)

        map_hmi.save(save_path, overwrite=True)

        fig = plt.figure()
        ax = plt.subplot(projection=map_hmi)
        map_hmi.plot_settings['norm'].vmin = -500
        map_hmi.plot_settings['norm'].vmax = 500
        im = map_hmi.plot()
        im = ax.get_images()
        im_lims =  im[0].get_extent()
        ax.set_aspect(abs((im_lims[1]-im_lims[0])/(im_lims[3]-im_lims[2])))
        plt.title(f'SDO HMI magnetogram {map_int.date_start.strftime("%Y-%m-%d %H:%M:%S")}')
        plt.colorbar(extend='both')
        plt.savefig(save_path.replace('.fits', '.png'), dpi=300)
        plt.close()

    print("Magnetic flux maps creation complete.")

def create_aia_maps(base_files, skip_done=False, resample=True, aia_channel=193):
    """
    Creates AIA maps corresponding to EIS fit files.

    Parameters
    ----------
    base_files : list
        List of paths to EIS fit files.
    """

    for i, base_file in enumerate(base_files):
        
        file_dir = os.path.dirname(base_file)
        file_name = os.path.basename(base_file)
        save_dir = file_dir.replace('data_eis', 'data_aia')
        save_name = file_name.replace('.fit.h5', '.fits')
        save_path = os.path.join(save_dir, save_name)

        # Check if both FITS and PNG files exist
        fits_file_exists = os.path.exists(save_path)
        png_file_exists = os.path.exists(save_path.replace('.fits', '.png'))

        # Skip processing if both files exist
        if skip_done and fits_file_exists and png_file_exists:
            print(f'Skipping {base_file} as both FITS and PNG files already exist.')
            continue

        print(f'Processing {base_file} ({i+1}/{len(base_files)})')

        fits = eispac.read_fit(base_file)
        main_component = fits.fit['main_component']
        map_int = fits.get_map(component=main_component, measurement='intensity')

        download_success = False
        extend_count = 0
        fail_count = 0
        while not download_success:
            try:
                time_range = a.Time(parse_time(map_int.meta['date_beg']) - extend_count * u.hour, parse_time(map_int.meta['date_end']) + extend_count * u.hour)
                aia_parameters = time_range & a.Instrument.aia & a.Physobs.intensity & a.Wavelength(aia_channel * u.angstrom)
                aia_search = Fido.search(aia_parameters)

                if aia_search.file_num == 0:
                    if extend_count >= 10:
                        print(f"No AIA data found within extended time range for {base_file}. Skipping this file.")
                        break
                    extend_count += 2
                    print(f"No AIA data found, extending the search window by {extend_count} hours.")
                    continue
                
                middle_result = int(len(aia_search[0]) / 2)
                aia_result = Fido.fetch(aia_search[0][middle_result])
                map_aia = sunpy.map.Map(aia_result[0])

                if sunpy.map.is_all_off_disk(map_aia):
                    print(f"Off-limb result found from Fido. Looking around for AIA coverage during time window.")

                    # try the result before, then the result after, then the one before the one before, then the one after the one after, etc.
                    # r should be middle_result, then -1, then +1, then -2, then +2, etc.
                    list_length = int(len(aia_search[0]) -2)  # This can be adjusted to change the number of elements
                    number_list = []
                    for j in range(list_length): # Loop through each number needed to fill the list
                        number = j // 2 + 1 # Calculate the base number for each pair
                        if j % 2 == 0: # Determine if the current position should be positive or negative
                            number_list.append(number) # Even index, add the positive number
                        else:
                            number_list.append(-number) # Odd index, add the negative of the number
                    check_every_n = 10

                    for r in number_list:
                        # if r is not a multiple of check_every_n, skip this iteration
                        if abs(r) % check_every_n != 0:
                            continue
                        aia_result = Fido.fetch(aia_search[0][middle_result+r])
                        map_aia = sunpy.map.Map(aia_result[0])
                        if sunpy.map.is_all_off_disk(map_aia):
                            print(f'Still off-limb result found from Fido. Continuing search for AIA coverage during time window.')
                            map_aia = None
                        else:
                            print(f'Found on-disk result from Fido. Continuing with this result.')
                            break
                        
                    if map_aia is None:
                        print(f"No on-disk AIA data found within set time range for {base_file}. Extending search window.")
                        extend_count += 2
                        continue

                download_success = True
            except Exception as e:
                fail_count += 1
                if fail_count >= 10:
                    print(f"Failed to download AIA data after 10 attempts for {base_file}. Skipping this file.")
                    break
                print(f"Download failed, retrying... (Attempt {fail_count}/10)")
                # time.sleep(60)

        if not download_success:
            # give the user a warning
            print(f"WARNING!!! Failed to download AIA data for {base_file}. Skipping this file.")
            # sleep for 5 seconds
            time.sleep(5)
            continue
        
        map_aia = sunpy.map.Map(aia_result[0])

        if resample:
            bottom_left = SkyCoord(map_int.meta['crval1'] * u.arcsec, map_int.meta['crval2'] * u.arcsec, frame=map_aia.coordinate_frame)
            top_right = SkyCoord((map_int.meta['crval1'] + map_int.meta['cdelt1'] * map_int.meta['naxis1']) * u.arcsec, (map_int.meta['crval2'] + map_int.meta['cdelt2'] * map_int.meta['naxis2']) * u.arcsec, frame=map_aia.coordinate_frame)
            map_aia = map_aia.submap(bottom_left, top_right=top_right)
            map_aia = map_aia.resample(u.Quantity([map_int.meta['naxis1'], map_int.meta['naxis2']]) * u.pixel)

        pathlib.Path(save_dir).mkdir(parents=True, exist_ok=True)

        map_aia.save(save_path, overwrite=True)

        fig = plt.figure()
        ax = plt.subplot(projection=map_aia)
        map_aia.plot_settings['norm'].vmin = 0
        map_aia.plot_settings['norm'].vmax = 1e4
        im = map_aia.plot()
        im = ax.get_images()
        im_lims =  im[0].get_extent()
        ax.set_aspect(abs((im_lims[1]-im_lims[0])/(im_lims[3]-im_lims[2])))
        plt.title(f'SDO AIA {aia_channel} Å {map_int.date_start.strftime("%Y-%m-%d %H:%M:%S")}')
        plt.colorbar(extend='both')
        plt.savefig(save_path.replace('.fits', '.png'), dpi=300)
        plt.close()