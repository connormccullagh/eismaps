import eispac
import os
import matplotlib.pyplot as plt
import matplotlib
import time
import numpy as np
import matplotlib.colors as mcolor
import sunpy.map
import sunkit_image.coalignment
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
    
def batch(files, save_fit=False, save_plot=False, measurement=None, clip=False, output_dir=None, los_correct=False, skip_done=True):
    
    if not save_fit and not save_plot:
        print("No output specified. Exiting.")
        return False
    if measurement is None:
        print("No measurement specified. Exiting.")
        return False

    for file in files:

        file_name = os.path.basename(file)
        file_line = file_name.split('.')[1]
        print(f'Processing {file_line} (file: {files.index(file)+1}/{len(files)})')

        file_date = os.path.basename(file).split('.')[0].split('_')[1]

        output_path = os.path.join(output_dir, file_date[:4], file_date[4:6], file_date[6:8]) if output_dir else os.path.dirname(file)
        os.makedirs(output_path, exist_ok=True)
        print(f"Saving to {output_path}")

        fit_res = eispac.core.read_fit(file)

        main_component = fit_res.fit['main_component']

        for m in measurement:
            output_file_png = os.path.join(output_path, file_name.replace('.fit.h5', f'.{m}.png'))
            output_file_fit = os.path.join(output_path, file_name.replace('.fit.h5', f'.{m}.fits'))

            if skip_done:
                if save_fit and os.path.exists(output_file_fit):
                    print(f"Skipping {m} map for {file_name} because it already exists.")
                    continue
                if save_plot and os.path.exists(output_file_png):
                    print(f"Skipping {m} map for {file_name} because it already exists.")
                    continue
            
            if m == 'ntv':
                m_map = get_ntv(file)
            else:
                m_map = fit_res.get_map(component=main_component, measurement=m)

            if m == 'int':
                if clip:
                    m_map = map_sd_filter(m_map, stds=6, log=True)
            if m == 'vel':
                if clip:
                    m_map = map_sd_filter(m_map, stds=6)
            if m == 'wid':
                if clip:
                    m_map = map_sd_filter(m_map, stds=3)
            if m == 'ntv':
                if clip:
                    m_map = map_sd_filter(m_map, stds=3)

            if los_correct:
                helioprojective_coords = sunpy.map.all_coordinates_from_map(m_map)
                heliographic_coords = helioprojective_coords.transform_to(frames.HeliographicStonyhurst)
                latitude = heliographic_coords.lat.to(u.deg).value
                longitude = heliographic_coords.lon.to(u.deg).value

                cos_mu = np.cos(latitude*((2*np.pi)/360)) * np.cos(longitude*((2*np.pi)/360))

                m_map_data = np.zeros(m_map.data.shape)
                for i in range (0, m_map.data.shape[0]):
                    for j in range (0, m_map.data.shape[1]):
                        m_map.data[i][j] = m_map.data[i][j] / cos_mu[i][j]

                # solar_x_len = fd_map.data.shape[1]
                # solar_y_len = fd_map.data.shape[0]

                # # Using np.linspace to ensure the arrays have the exact expected length
                # solar_x = np.linspace(-fd_size/2, fd_size/2, solar_x_len)
                # solar_y = np.linspace(-fd_size/2, fd_size/2, solar_y_len)

                # # create meshgrid for solar_x and solar_y
                # solar_x_array, solar_y_array = np.meshgrid(solar_x, solar_y)

                # # convert solar_x and solar_y to longitude and latitude
                # obs_time = fd_map.meta['date_obs']
                # helioprojective_coordinates = SkyCoord(solar_x_array * u.arcsec, solar_y_array * u.arcsec, frame=frames.Helioprojective(observer="earth", obstime=obs_time))
                # heliographic_coordinates = helioprojective_coordinates.transform_to(frames.HeliographicStonyhurst)
                # latitude = heliographic_coordinates.lat.to(u.deg).value
                # longitude = heliographic_coordinates.lon.to(u.deg).value

                # # apply the doppler velocity correction
                # cos_mu = np.cos(latitude*((2*np.pi)/360)) * np.cos(longitude*((2*np.pi)/360))

                # for i in range (0, fd_map.data.shape[0]):
                #     for j in range (0, fd_map.data.shape[1]):
                #         fd_map.data[i][j] = fd_map.data[i][j] / cos_mu[i][j]

            if save_fit:
                m_map.save(output_file_fit, overwrite=True)

            if save_plot:
                if m == 'int':
                    plt.figure()
                    m_cmap = plt.get_cmap('Oranges_r').copy()
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

                    # m_map.plot_settings['norm'].vmin = -10 # fixed symmetrical
                    # m_map.plot_settings['norm'].vmax = 10
                    
                    # data_mean = np.nanmean(m_map.data) # symmetrical plotting around 3 sigma (unrobust against extreme outliers)
                    # data_std = np.nanstd(m_map.data)
                    # sigma_limit = 3 * data_std
                    # vmin = data_mean - sigma_limit
                    # vmax = data_mean + sigma_limit
                    # max_limit = max(abs(vmin), abs(vmax))
                    # vmin, vmax = -max_limit, max_limit
                    # m_map.plot_settings['norm'].vmin = vmin
                    # m_map.plot_settings['norm'].vmax = vmax

                    vmin_percentile = np.nanpercentile(m_map.data, 10) # using percentiles (robust against extreme outliers)
                    vmax_percentile = np.nanpercentile(m_map.data, 90)
                    max_limit = max(abs(vmin_percentile), abs(vmax_percentile))
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

        # if 'int' in measurement and not skip_int:
        #     map_int = fit_res.get_map(component=main_component,measurement='intensity')
        #     if clip:
        #         map_int = map_sd_filter(map_int, stds=6, log=True)
        #     if save_fit:
        #         map_int.save(file.replace('.fit.h5', '.int.fits'), overwrite=True)
        #     if save_plot:
        #         plt.figure()
        #         cmap_int = plt.get_cmap('Oranges').copy()
        #         cmap_int.set_bad(color='gray')
        #         map_int.plot_settings['cmap'] = cmap_int
        #         map_int.plot_settings['norm'].vmin = 1e-1
        #         map_int.plot()
        #         plt.colorbar(label='Intensity (DN/s)', extend='max')
        #         plt.savefig(os.path.join(output_path, file_name.replace('.fit.h5', '.int.png')))
        #         plt.close()
        # if 'vel' in measurement and not skip_vel:
        #     map_vel = fit_res.get_map(component=main_component,measurement='velocity')
        #     if clip:
        #         map_vel = map_sd_filter(map_vel, stds=6)
        #     if save_fit:
        #         map_vel.save(file.replace('.fit.h5', '.vel.fits'), overwrite=True)
        #     if save_plot:
        #         plt.figure()
        #         cmap_vel = plt.get_cmap('RdBu_r').copy()
        #         cmap_vel.set_bad(color='gray')
        #         map_vel.plot_settings['cmap'] = cmap_vel
        #         map_vel.plot_settings['norm'].vmin = -10 ### TEST WITHOUT THIS. AT LEAST FORCE TO BE SYMMETRICAL
        #         map_vel.plot_settings['norm'].vmax = 10
        #         map_vel.plot()
        #         plt.colorbar(label='Doppler velocity (km/s)', extend='both')
        #         plt.savefig(os.path.join(output_path, file_name.replace('.fit.h5', '.vel.png')))
        #         plt.close()
        # if 'wid' in measurement and not skip_wid:
        #     map_wid = fit_res.get_map(component=main_component,measurement='width')
        #     if clip:
        #         map_wid = map_sd_filter(map_wid, stds=3)
        #     if save_fit:
        #         map_wid.save(file.replace('.fit.h5', '.wid.fits'), overwrite=True)
        #     if save_plot:
        #         plt.figure()
        #         cmap_wid = plt.get_cmap('viridis').copy()
        #         cmap_wid.set_bad(color='gray')
        #         map_wid.plot_settings['cmap'] = cmap_wid
        #         map_wid.plot()
        #         plt.colorbar(label='Width (Angstrom)', extend='both')
        #         plt.savefig(os.path.join(output_path, file_name.replace('.fit.h5', '.wid.png')))
        #         plt.close()
        # if 'ntv' in measurement and not skip_ntv:
        #     map_ntv = get_ntv(file)
        #     if clip:
        #         map_ntv = map_sd_filter(map_ntv, stds=3, log=False)
        #     if save_fit:
        #         map_ntv.save(file.replace('.fit.h5', '.ntv.fits'), overwrite=True)
        #     if save_plot:
        #         plt.figure()
        #         cmap_ntv = plt.get_cmap('inferno').copy()
        #         cmap_ntv.set_bad(color='gray')
        #         map_ntv.plot_settings['cmap'] = cmap_ntv
        #         map_ntv.plot()
        #         plt.colorbar(label='Non-thermal velocity (km/s)', extend='both')
        #         plt.savefig(os.path.join(output_path, file_name.replace('.fit.h5', '.ntv.png')))
        #         plt.close()

    return True

def get_ntv(file):
    import pandas as pd
    from eismaps.utils.roman_numerals import int_to_roman, roman_to_int
    import numpy as np
    import sunpy.map
    import eismaps.utils.width2velocity

    print(f'file: {file}')
    
    fit_res = eispac.core.read_fit(file)
    main_component = fit_res.fit['main_component']
    map_wid = fit_res.get_map(component=main_component,measurement='width')

    yy = map_wid.data.shape[0]
    xx = map_wid.data.shape[1]

    wavelength_array = np.full((yy,xx),map_wid.wavelength.to_value()) # fixed reference wavelength

    file_name = os.path.basename(file)
    file_line = file_name.split('.')[1]
    file_line_parts = file_line.split('_')

    if len(file_line_parts) == 4:
      element = file_line_parts[0].capitalize()
      ion = int(roman_to_int(file_line_parts[1]))
    elif len(file_line_parts) == 5: # protecting against s__xi_111_111 (double _)
      element = file_line_parts[0].capitalize()
      ion = int(roman_to_int(file_line_parts[2]))
    else:
      print(f"Invalid file line format: {file_line}.")
      exit()

    t_max, v_therm = eismaps.utils.width2velocity.width2velocity(element, ion)

    obser_fwhm = (map_wid.data*2*np.sqrt(2*np.log(2))) # set the observed fwhm (use the fwhm formula)

    instr_fwhm = fit_res.meta['slit_width'] # set the instrument fwhm (this is a column vector)

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

# def get_bwa(file):
#     # blue wing asymmetry
#     return None

def coalign(eis_maps, aia_maps, save=False):
    
    shifts = []

    for eis_map, aia_map in zip(eis_maps, aia_maps):
        mc = sunpy.map.Map(
            [
                aia_map,
                eis_map,
            ],
            sequence=True,
        )

        shift = sunkit_image.coalignment.calculate_match_template_shift(mc, layer_index=0)

        shifts.append(shift)

        if save:
            mc_shifted = sunkit_image.coalignment.mapsequence_coalign_by_match_template(mc, shift=shift)

            # save the aligned maps
            map_shifted = mc_shifted[1]
            map_shifted.save(eis_map, overwrite=True)

    # return a list of the shifts needed for each map
    return shifts