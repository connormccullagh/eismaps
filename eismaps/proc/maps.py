import eispac
import os
import matplotlib.pyplot as plt
import matplotlib
import time
import numpy as np
import matplotlib.colors as mcolor
import sunpy.map

CC = 2.9979E+5 # speed of light

def map_sd_filter(map,stds=3,log=True):

    data = map.data
    if log:
        data[data == 0] = np.nan
        data = np.log10(data)
    data_avg = np.nanmean(data)
    data_std = np.nanstd(data)

    map.data[np.abs(data-data_avg) > stds*data_std] = np.nan

    return map

# def batch(files, save_fit=False, save_plot=False, measurement=None, line=None, clip=False, output_dir=None):
def batch(files, save_fit=False, save_plot=False, measurement=None, clip=False, output_dir=None):
## todo: add line as an option

    for file in files:
        file_name = os.path.basename(file)
        file_line = file_name.split('.')[1]

        base_name = file_name.replace('.fit.h5', '') # construct output file paths
        output_path = output_dir if output_dir else os.path.dirname(file)        

        fit_res = eispac.core.read_fit(file)
        main_component = fit_res.fit['main_component']

        if 'int' in measurement:
            map_int = fit_res.get_map(component=main_component,measurement='intensity')
            if clip:
                map_int = map_sd_filter(map_int, stds=6)
            if save_fit:
                map_int.save(os.path.join(output_path, f'{base_name}.int.fits'), overwrite=True)
            if save_plot:
                plt.figure()
                cmap_int = plt.get_cmap('Oranges').copy()
                cmap_int.set_bad(color='gray')
                map_int.plot_settings['cmap'] = cmap_int
                map_int.plot_settings['norm'].vmin = 1e-1
                map_int.plot()
                plt.colorbar(label='Intensity (DN/s)', extend='max')
                plt.savefig(os.path.join(output_path, f'{base_name}.int.png'))
                plt.close()
        if 'vel' in measurement:
            map_vel = fit_res.get_map(component=main_component,measurement='velocity')
            if clip:
                map_vel = map_sd_filter(map_vel, stds=6, log=False)
            if save_fit:
                map_vel.save(os.path.join(output_path, f'{base_name}.vel.fits'), overwrite=True)
            if save_plot:
                plt.figure()
                cmap_vel = plt.get_cmap('RdBu_r').copy()
                cmap_vel.set_bad(color='gray')
                map_vel.plot_settings['cmap'] = cmap_vel
                map_vel.plot_settings['norm'].vmin = -10
                map_vel.plot_settings['norm'].vmax = 10
                # data_max = np.max(np.abs([np.nanmax(map_vel.data), np.nanmin(map_vel.data)]))
                # map_vel.plot_settings['norm'].vmin = -data_max
                # map_vel.plot_settings['norm'].vmax = data_max
                map_vel.plot()
                plt.colorbar(label='Doppler velocity (km/s)', extend='both')
                plt.savefig(os.path.join(output_path, f'{base_name}.vel.png'))
                plt.close()

        if 'wid' in measurement:
            map_wid = fit_res.get_map(component=main_component,measurement='width')
            if clip:
                map_wid = map_sd_filter(map_wid, stds=3, log=False)
            if save_fit:
                map_wid.save(os.path.join(output_path, f'{base_name}.wid.fits'), overwrite=True)
            if save_plot:
                plt.figure()
                cmap_wid = plt.get_cmap('viridis').copy()
                cmap_wid.set_bad(color='gray')
                map_wid.plot_settings['cmap'] = cmap_wid
                map_wid.plot()
                plt.colorbar(label='Width (Angstrom)', extend='both')
                plt.savefig(os.path.join(output_path, f'{base_name}.wid.png'))
                plt.close()
        if 'ntv' in measurement:
            map_ntv = get_ntv(file)
            if clip:
                map_ntv = map_sd_filter(map_ntv, stds=3, log=False)
            if save_fit:
                map_ntv.save(os.path.join(output_path, f'{base_name}.ntv.fits'), overwrite=True)
            if save_plot:
                plt.figure()
                cmap_ntv = plt.get_cmap('inferno').copy()
                cmap_ntv.set_bad(color='gray')
                map_ntv.plot_settings['cmap'] = cmap_ntv
                map_ntv.plot()
                plt.colorbar(label='Non-thermal velocity (km/s)', extend='both')
                plt.savefig(os.path.join(output_path, f'{base_name}.ntv.png'))
                plt.close()

    return True

def get_ntv(file):
    import pandas as pd
    from eismaps.utils.roman_numerals import int_to_roman
    import numpy as np
    import sunpy.map
    import eismaps.utils.width2velocity
    
    fit_res = eispac.core.read_fit(file)
    main_component = fit_res.fit['main_component']
    map_wid = fit_res.get_map(component=main_component,measurement='width')

    yy = map_wid.data.shape[0]
    xx = map_wid.data.shape[1]

    wavelength_array = np.full((yy,xx),map_wid.wavelength.to_value()) # fixed reference wavelength

    file_name = os.path.basename(file)
    file_line = file_name.split('.')[1]
    # file line could be fe_12_195_119 or s__2_125_119
    file_line_parts = file_line.split('_')
    if len(file_line_parts) == 4:
      element = file_line_parts[0].capitalize()
      ion = int(file_line_parts[1])
    elif len(file_line_parts) == 5:
      element = file_line_parts[0].capitalize()
      ion = int(file_line_parts[2])
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