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
    
def batch(files, measurement=None, clip=False, save_fit=True, save_plot=False, output_dir=None, output_dir_tree=False, vel_los_correct=False, skip_done=True):

    if not save_fit and not save_plot: print("No output specified. Exiting."); return False
    VALID_MEASUREMENTS = ['int', 'vel', 'wid', 'ntv']
    assert measurement is not None and any(m in VALID_MEASUREMENTS for m in measurement), f"Invalid measurement specified. Must be one of {VALID_MEASUREMENTS}."

    for file in files:

        file_name = os.path.basename(file)
        file_date = os.path.basename(file).split('.')[0].split('_')[1]

        if output_dir is None: output_dir = os.path.dirname(file)
        if output_dir_tree: output_dir = os.path.join(output_dir, file_date[:4], file_date[4:6], file_date[6:8])

        fit_res = eispac.core.read_fit(file)
        main_component = fit_res.fit['main_component']

        for m in measurement:

            output_file_fit = os.path.join(output_dir, f"{file_name.replace('.fit.h5', f'.{m}.fits')}")

            if save_fit:
                output_file_fit = os.path.join(output_dir, f"{file_name.replace('.fit.h5', f'.{m}.fits')}")
                if skip_done and os.path.exists(output_file_fit): print(f"Skipping {m} map for {file_name} because it already exists."); continue

            if save_plot:
                output_file_png = os.path.join(output_dir, f"{file_name.replace('.fit.h5', f'.{m}.png')}")
                if skip_done and os.path.exists(output_file_fit): print(f"Skipping {m} map for {file_name} because it already exists."); continue

            if m == 'ntv': # not an eispac default, so need eismaps
                m_map = get_ntv_map(fit_res, component=main_component)
            else:
                m_map = fit_res.get_map(component=main_component, measurement=m)

            if clip:
                if m == 'int': m_map = map_sd_filter(m_map, stds=6, log=True)
                if m == 'vel': m_map = map_sd_filter(m_map, stds=6)
                if m == 'wid': m_map = map_sd_filter(m_map, stds=3)
                if m == 'ntv': m_map = map_sd_filter(m_map, stds=3)

            if vel_los_correct and m == 'vel':
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

                los_factor = np.sin(latitude_rad) * np.sin(observer_lat_rad) + \
                            np.cos(latitude_rad) * np.cos(observer_lat_rad) * np.cos(longitude_rad - observer_lon_rad)

                for i in range(m_map.data.shape[0]):
                    for j in range(m_map.data.shape[1]):
                            m_map.data[i,j] = m_map.data[i,j] / los_factor[i,j]

            if save_fit:
                m_map.save(output_file_fit, overwrite=True)

            if save_plot:
                if m == 'int':
                    plt.figure()
                    m_cmap = plt.get_cmap('gist_heat').copy()
                    m_cmap.set_bad(color='gray')
                    m_map.plot_settings['cmap'] = m_cmap
                    m_map.plot_settings['norm'].vmin = 1e-1
                    m_map.plot()
                    plt.colorbar(label='Intensity (DN/s)', extend='max')
                    plt.savefig(output_file_png)
                    plt.close()
                if m == 'vel':
                    plt.figure()
                    m_cmap = plt.get_cmap('RdBu_r').copy()
                    m_cmap.set_bad(color='gray')
                    m_map.plot_settings['cmap'] = m_cmap
                    vmin_percentile = np.nanpercentile(m_map.data, 10) # using percentiles (robust against extreme outliers)
                    vmax_percentile = np.nanpercentile(m_map.data, 90)
                    max_limit = max(abs(vmin_percentile), abs(vmax_percentile)) # set the limits to be symmetrical
                    vmin, vmax = -max_limit, max_limit
                    m_map.plot_settings['norm'].vmin = vmin
                    m_map.plot_settings['norm'].vmax = vmax
                    m_map.plot()
                    plt.colorbar(label='Doppler velocity (km/s)', extend='both')
                    plt.savefig(output_file_png)
                    plt.close()
                if m == 'wid':
                    plt.figure()
                    m_cmap = plt.get_cmap('viridis').copy()
                    m_cmap.set_bad(color='gray')
                    m_map.plot_settings['cmap'] = m_cmap
                    m_map.plot()
                    plt.colorbar(label='Width (Angstrom)', extend='both')
                    plt.savefig(output_file_png)
                    plt.close()
                if m == 'ntv':
                    plt.figure()
                    m_cmap = plt.get_cmap('inferno').copy()
                    m_cmap.set_bad(color='gray')
                    m_map.plot_settings['cmap'] = m_cmap
                    vmax_percentile = np.nanpercentile(m_map.data, 90)
                    m_map.plot_settings['norm'].vmin = 0
                    m_map.plot_settings['norm'].vmax = vmax_percentile
                    m_map.plot()
                    plt.colorbar(label='Non-thermal velocity (km/s)', extend='both')
                    plt.savefig(output_file_png)
                    plt.close()

    return

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

# def get_bwa_map(file):
#     # blue wing asymmetry
#     return None