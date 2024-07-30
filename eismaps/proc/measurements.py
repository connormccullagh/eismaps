import eispac

CC = 2.9979E+5 # speed of light

def batch_l2(files,measurement=None):
    import os

    for file in files:
        file_dir = os.path.dirname(file)

        if measurement == 'ntv':
            l2_map = get_ntv(file)
            path = file.replace('.data.h5', f'.ntv.fit')
        elif measurement == 'bwa':
            l2_map = get_bwa(file)
            path = file.replace('.data.h5', f'.bwa.fit')
        eispac.save_fit(l2_map, save_dir=file_dir)

    return path

def get_ntv(file):
    import pandas as pd
    from eismaps.utils.roman_numerals import int_to_roman
    import numpy as np
    import sunpy.map
    
    fit = eispac.read_fit(file)

    map_wid = fit.get_map(measurement='width')

    yy = map_wid.data.shape[0]
    xx = map_wid.data.shape[1]

    wavelength_array = np.full((yy,xx),map_wid.wavelength.to_value()) # fixed reference wavelength

    element = fit.meta['eismaps_element']
    ion = fit.meta['eismaps_ion']

    t_max, v_therm = pd.DataFrame(element,ion)

    obser_fwhm = (map_wid.data*2*np.sqrt(2*np.log(2))) # set the observed fwhm (use the fwhm formula)

    instr_fwhm = fit.meta['slit_width'] # set the instrument fwhm (this is a column vector)

    therm_fwhm = np.sqrt(4*np.log(2))*wavelength_array*v_therm/CC # calculate the thermal fwhm

    dl_o = obser_fwhm # set the observed fwhm
    dl_i = np.full((yy,xx), instr_fwhm)
    dl_t = therm_fwhm

    ###

    # copy the meta data from the original fits file
    ntv_meta = fit.meta.copy()
    ntv_meta['measurement'] = 'non-thermal velocity'
    ntv_meta['bunit'] = 'km/s'
    map_ntv = sunpy.map.Map(v_nt, ntv_meta) # create a non-thermal velocity map

    return map_ntv

def get_bwa(file):
    # blue wing asymmetry
    return None