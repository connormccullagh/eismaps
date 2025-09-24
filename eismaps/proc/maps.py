import eispac
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolor
import sunpy.map
from sunpy.coordinates import frames
import astropy.units as u

CC = 2.9979E+5 # speed of light

def map_sd_filter(map,stds=3,log=False):

    data = map.data
    if log:
        data[data == 0] = np.nan
        data = np.log10(data)
    data_avg = np.nanmean(data)
    data_std = np.nanstd(data)

    map.data[np.abs(data-data_avg) > stds*data_std] = np.nan

    return map
    
def batch(files, measurement=None, clip=False, save_fit=True, save_plot=False,
          output_dir=None, output_dir_tree=False, vel_los_correct=False,
          skip_done=True, mssl_solarb_file_format=False):
    """
    Batch process EIS files to create maps for different measurements.

    If 'err_vel' is requested, it will automatically copy metadata from the 'int' map
    in the same file.
    """
    from sunpy.map import Map
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from sunpy.coordinates import frames
    import astropy.units as u
    import eispac
    if not save_fit and not save_plot:
        print("No output specified. Exiting.")
        return False

    VALID_MEASUREMENTS = ['int', 'vel', 'wid', 'ntv', 'chi2', 'err_int', 'err_vel']
    assert measurement is not None and any(m in VALID_MEASUREMENTS for m in measurement), \
        f"Invalid measurement specified. Must be one of {VALID_MEASUREMENTS}."

    out_files = []

    for file in files:
        file_name = os.path.basename(file)
        file_date = file_name.split('.')[0].split('_')[1]

        if output_dir is None:
            output_dir = os.path.dirname(file)
        if output_dir_tree:
            output_dir = os.path.join(output_dir, file_date[:4], file_date[4:6], file_date[6:8])
            os.makedirs(output_dir, exist_ok=True)

        fit_res = eispac.core.read_fit(file)
        main_component = fit_res.fit['main_component']

        for m in measurement:
            output_file_fit = os.path.join(output_dir, f"{file_name.replace('.fit.h5', f'.{m}.fits')}")
            output_file_png = None
            output_file_gif = None

            if save_fit and skip_done and os.path.exists(output_file_fit):
                print(f"Skipping {m} map for {file_name} because it already exists.")
                continue

            if save_plot:
                output_file_png = os.path.join(output_dir, f"{file_name.replace('.fit.h5', f'.{m}.png')}")
                if mssl_solarb_file_format:
                    output_file_png_filename = os.path.basename(output_file_png)
                    output_file_datetime = f"{output_file_png_filename.split('.')[0].split('_')[1]}_{output_file_png_filename.split('.')[0].split('_')[2]}"
                    output_file_iwin = fit_res.meta['iwin']
                    output_file_lineid = fit_res.meta['line_id']
                    output_file_gif_filename = f"eis_l0_{output_file_datetime}.fits_line_{output_file_iwin}_{output_file_lineid.replace(' ', '_').upper()}.{m}.gif"
                    output_file_gif = os.path.join(output_dir, output_file_gif_filename)

            # -----------------------------
            # Generate the map
            # -----------------------------
            if m == 'ntv':
                m_map = get_ntv_map(fit_res, component=main_component)
            elif m == 'chi2':
                m_map = get_chi2_map(fit_res, component=main_component)
            elif m == 'err_vel':
                data = fit_res.fit['err_vel'][:, :, main_component]

                # Use 'int' measurement in the same file for metadata
                int_map = fit_res.get_map(component=main_component, measurement='int')
                meta = int_map.meta.copy()

                m_map = Map(data, meta)
            else:
                m_map = fit_res.get_map(component=main_component, measurement=m)

            # -----------------------------
            # Clip outliers
            # -----------------------------
            if clip:
                if m in ['int', 'err_int']:
                    m_map = map_sd_filter(m_map, stds=6, log=True)
                elif m in ['vel', 'err_vel']:
                    m_map = map_sd_filter(m_map, stds=6)
                elif m == 'wid':
                    m_map = map_sd_filter(m_map, stds=3)
                elif m == 'ntv':
                    m_map = map_sd_filter(m_map, stds=3)

            # -----------------------------
            # LOS correction
            # -----------------------------
            if vel_los_correct and m in ['vel', 'err_vel']:
                helioprojective_coords = sunpy.map.all_coordinates_from_map(m_map)
                heliographic_coords = helioprojective_coords.transform_to(frames.HeliographicStonyhurst)
                latitude = heliographic_coords.lat.to(u.deg).value
                longitude = heliographic_coords.lon.to(u.deg).value
                observer_lon = m_map.observer_coordinate.lon.to(u.deg).value
                observer_lat = m_map.observer_coordinate.lat.to(u.deg).value

                latitude_rad = np.deg2rad(latitude)
                longitude_rad = np.deg2rad(longitude)
                observer_lon_rad = np.deg2rad(observer_lon)
                observer_lat_rad = np.deg2rad(observer_lat)

                los_factor = (np.sin(latitude_rad) * np.sin(observer_lat_rad) +
                              np.cos(latitude_rad) * np.cos(observer_lat_rad) *
                              np.cos(longitude_rad - observer_lon_rad))

                los_corrected_data = m_map.data / los_factor
                # Replace any NaNs or infinities with a large masked value
                los_corrected_data = np.nan_to_num(los_corrected_data, nan=0.0, posinf=0.0, neginf=0.0)
                m_map = Map(los_corrected_data, m_map.meta)       

            # -----------------------------
            # Save map
            # -----------------------------
            if save_fit:
                m_map.save(output_file_fit, overwrite=True)
                out_files.append(output_file_fit)

            if save_plot:
                plt.figure()
                m_map.plot()
                if mssl_solarb_file_format:
                    plt.savefig(output_file_gif, dpi=100)
                    out_files.append(output_file_gif)
                else:
                    plt.savefig(output_file_png)
                    out_files.append(output_file_png)
                plt.close()

    return out_files

def get_ntv_map(fit_res, component=None):
    import pandas as pd
    from eismaps.utils.roman_numerals import int_to_roman, roman_to_int
    import numpy as np
    import sunpy.map
    import eismaps.utils.width2velocity

    if component is None: component = fit_res.fit['main_component']

    map_wid = fit_res.get_map(component=component,measurement='wid')

    yy = map_wid.data.shape[0]
    xx = map_wid.data.shape[1]

    wavelength_array = np.full((yy,xx),map_wid.wavelength.to_value()) # fixed reference wavelength

    line = map_wid.meta['line_id'] # format e.g. ('line_id': 'Fe XII 195.119')
    line = line.replace('  ',' ') # protecting against e.g. 'S  XXI 135.8'

    element = line.split(' ')[0]
    ion = int(roman_to_int(line.split(' ')[1]))

    t_max, v_therm = eismaps.utils.width2velocity.width2velocity(element, ion)

    obser_fwhm = (map_wid.data*2*np.sqrt(2*np.log(2))) # set the observed fwhm (use the fwhm formula)

    if 'slit_width' in fit_res.meta and fit_res.meta['slit_width'] is not None:
        instr_fwhm = fit_res.meta['slit_width'] # set the instrument fwhm (this is a column vector)
    else:
        print("WARNING: Instrumental width not in metadata and so not used")
        instr_fwhm = np.full([yy,xx],0)

    therm_fwhm = np.sqrt(4*np.log(2))*wavelength_array*v_therm/CC # calculate the thermal fwhm

    # create the empty arrays
    dl_o = np.empty([yy,xx]) # observed fwhm
    dl_i = np.empty([yy,xx]) # instrument fwhm
    dl_t = np.empty([yy,xx]) # thermal fwhm

    dl_o = obser_fwhm # set the observed fwhm
    for j in range(0,xx): dl_i[0:,j] = instr_fwhm # set the instrument fwhm
    dl_t = therm_fwhm # set the thermal fwhm

    dl_nt_2 = (np.square(dl_o) - np.square(dl_i) - np.square(dl_t)) # calculate non-thermal velocity fwhm**2

    dl_nt_2 = np.where( dl_nt_2>0, dl_nt_2, 0 ) # replace negative values with 0

    v_nt = np.sqrt( dl_nt_2 * np.divide(CC**2, (4*np.log(2)*np.square(wavelength_array)), out=np.zeros_like(wavelength_array), where=wavelength_array!=0) ) # calculate non-thermal velocity using v_nt = np.sqrt(dl_nt_2 * (cc**2/(4*np.log(2)*np.square(wavelength)))) but covering the fact that sometimes the wavelength is 0

    map_ntv = sunpy.map.Map(v_nt, map_wid.meta) # create a non-thermal velocity map
    map_wid.meta['measurement'] = 'non-thermal velocity'
    map_wid.meta['bunit'] = 'km/s'

    return map_ntv

# def get_bwa_map(fit_res):
#     # blue wing asymmetry
#     return None

def get_chi2_map(fit_res, component=None):
    if component is None: component = fit_res.fit['main_component']

    map_int = fit_res.get_map(component=component,measurement='int')

    map_chi2 = sunpy.map.Map(fit_res.fit['chi2'], map_int.meta)

    map_chi2.meta['measurement'] = 'chi2'
    map_chi2.meta['bunit'] = 'chi2'
    return map_chi2
