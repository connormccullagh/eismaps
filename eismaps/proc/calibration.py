"""
calibration.py

This module provides functions for calibrating Hinode/EIS rasters using both the 2014 and 2023 calibration methods.

The calibration methods include:
- Effective area (EA) interpolation for different wavelength ranges.
- Calibration ratios based on pre-flight and in-flight data.
- Conversion from detector units (DN/s) to physical units (ergs/(sr cm^2 s)).
"""

import datetime
import re
import numpy as np
from scipy.interpolate import interp1d
from scipy.io import readsav
import sunpy.map
import os

# Paths to the .sav files used in calibration
CALIBRATION_DATA_PATH = os.path.join(os.path.dirname(__file__), 'calibration_data')
WARREN_2014_FILE = os.path.join(CALIBRATION_DATA_PATH, 'eis_calib_warren_2014.sav')
FIT_EA_2023_FILE = os.path.join(CALIBRATION_DATA_PATH, 'fit_eis_ea_2023-05-04.sav')
PREFLIGHT_SHORT_FILE = os.path.join(CALIBRATION_DATA_PATH, 'preflight_calib_short.sav')
PREFLIGHT_LONG_FILE = os.path.join(CALIBRATION_DATA_PATH, 'preflight_calib_long.sav')

__all__ = ['calib_2014', 'calib_2023', 'interpol_eis_ea', 'eis_ea']

# Utility function to convert time to TAI (Temps Atomique International)
def anytim2tai(time_str):
    """
    Converts a given time string into TAI (Temps Atomique International) format.
    """
    time_str = re.sub(r'[^\w\s.]', '', time_str)
    time_str, _, fractional_part = time_str.partition('.')
    if 'T' in time_str:
        date_str, time_str = time_str.split('T')
        time_str = f"{date_str} {time_str}"
    else:
        time_str = re.sub(r'[^\w\s]', '', time_str)
    dt = datetime.datetime.strptime(time_str, '%Y%m%d %H%M%S')
    seconds_since_epoch = (dt - datetime.datetime(1970, 1, 1)).total_seconds()
    tai_offset = (datetime.datetime(1970, 1, 1) - datetime.datetime(1958, 1, 1)).total_seconds()
    tai_time = seconds_since_epoch + tai_offset + 37
    return tai_time

# Function for interpolating effective area (EA) values for 2023 calibration
def interpol_eis_ea(date, wavelength, short=False, long=False, radcal=False, ea_file=FIT_EA_2023_FILE, quiet=False):
    if np.size(date) != 1:
        raise ValueError('ERROR: please input a single date')
    
    # Ensure wavelength is always treated as an array
    wavelength = np.atleast_1d(wavelength)
    
    in_tai = anytim2tai(date)

    if in_tai < anytim2tai('2006-10-20T10:20:00.000'):
        print('WARNING: Selected date is before the start of normal EIS science operations. Output values may be inaccurate.')

    if not short and not long:
        n_input_wave = np.size(wavelength)
        loc_short = np.where((wavelength >= 165) & (wavelength <= 213))[0]
        loc_long = np.where((wavelength >= 245) & (wavelength <= 292))[0]
        if (len(loc_short) + len(loc_long) < n_input_wave) or (len(loc_short) > 0 and len(loc_long) > 0):
            raise ValueError('ERROR: Invalid wavelength(s). Please only select values in either the short (165 - 213) or long (245 - 292) wavelength bands.')

    if short:
        wavelength = 1
    elif long:
        wavelength = 1000

    fit_ea = readsav(ea_file)['fit_ea']
    fit_dates = fit_ea.date_obs[0].astype(str)
    fit_easw = fit_ea.sw_ea[0]
    fit_ealw = fit_ea.lw_ea[0]
    sw_wave = fit_ea.sw_wave[0]
    lw_wave = fit_ea.lw_wave[0]

    ref_tai = np.array([anytim2tai(date) for date in fit_dates])

    if in_tai < ref_tai[0]:
        if not quiet:
            print(f"WARNING: Selected date is before the first calibrated date on {fit_ea.date_obs[0]}. Returning first fit calibration")
        in_tai = ref_tai[0]

    if in_tai > ref_tai[-1]:
        if not quiet:
            print(f"WARNING: Selected date is after the last calibrated date on {fit_ea.date_obs[-1][-1]}. Returning last fit calibration")
        in_tai = ref_tai[-1]

    if short or (np.size(wavelength) > 0 and np.max(wavelength) < 220):
        ref_ea = fit_easw
        ref_wave = sw_wave
    else:
        ref_ea = fit_ealw
        ref_wave = lw_wave

    n_ref_waves = len(ref_wave)
    new_ea = np.zeros(n_ref_waves)
    for w in range(n_ref_waves):
        ea_values = ref_ea[w, :]
        new_ea[w] = np.interp(in_tai, ref_tai, ea_values)

    if not short and not long:
        out_ea = interp1d(ref_wave, new_ea, kind='cubic')(wavelength)
    else:
        wavelength = ref_wave
        out_ea = new_ea

    if radcal:
        sr_factor = (725.0 / 1.496e8) ** 2
        ergs_to_photons = 6.626e-27 * 2.998e10 * 1.e8
        gain = 6.3
        phot_to_elec = 12398.5 / 3.65
        tau_sensitivity = 1894.0

        print('Returning radcal values for converting [DN/s] to [ergs/(sr cm^2 s)]')
        print('   Note: You may still need to adjust for exposure time and slitsize.')

        radcal = (wavelength * gain) / (out_ea * phot_to_elec)
        radcal = radcal * ergs_to_photons / wavelength / sr_factor
        out_ea = radcal

    return out_ea

# Calibration function for 2023 data
def calib_2023(map):
    """
    Calibrates a given EIS map using the 2023 calibration method.
    """
    match = re.search(r'\d+\.\d+', map.meta['line_id'])
    wvl_value = float(match.group())
    calib_ratio_2023 = eis_ea(wvl_value) / interpol_eis_ea(map.date.value, wvl_value)
    print(f'Calibration ratio for {wvl_value} A: {calib_ratio_2023} using map date {map.date.value} was applied.')
    new_map = sunpy.map.Map(map.data * calib_ratio_2023, map.meta)
    return new_map

# Calibration function for 2014 data
def calib_2014(map):
    """
    Calibrates a given EIS map using the 2014 calibration method.
    """
    match = re.search(r'\d+\.\d+', map.meta['line_id'])
    wvl_value = float(match.group())
    calib_ratio = eis_ea(wvl_value) / eis_ea_nrl(map.date.value, wvl_value)
    print(f'Calibration ratio for {wvl_value} A: {calib_ratio} using map date {map.date.value} was applied.')
    new_map = sunpy.map.Map(map.data * calib_ratio, map.meta)
    return new_map

# Function to calculate effective area (EA) based on preflight or in-flight data
def eis_ea(input_wave, short=False, long=False):
    """
    Returns the effective area (EA) for a given wavelength, either short or long band.
    """
    if short:
        wave, ea = eis_effective_area_read(short=True)
        input_wave = wave
        return ea

    if long:
        wave, ea = eis_effective_area_read(long=True)
        input_wave = wave
        return ea

    if isinstance(input_wave, (int, float)):
        input_wave = np.array([input_wave])

    nWave = len(input_wave)
    ea = np.zeros(nWave)

    for i in range(nWave):
        short, long = is_eis_wavelength(input_wave[i])

        if not short and not long:
            ea[i] = 0.0
        else:
            wave, area = eis_effective_area_read(long=long, short=short)
            ea[i] = np.exp(np.interp(input_wave[i], wave, np.log(area)))

    if nWave == 1:
        ea = ea[0]

    return ea

# Support functions for the 2014 calibration
def eis_ea_nrl(date, wave, short=False, long=False):
    """
    Interpolates the effective area (EA) based on the 2014 calibration using data from Warren et al.
    """
    eis = read_calib_file(WARREN_2014_FILE)
    t = (get_time_tai(date) - get_time_tai(eis['t0'][0].decode('utf-8'))) / (86400 * 365.25)
    ea_knots_SW = eis['a0_sw'][0] * np.exp(-t / eis['tau_sw'][0])
    ea_knots_LW = eis['a0_lw'][0] * np.exp(-t / eis['tau_lw'][0])

    if short:
        wave = eis['wave_area_sw'][0]
    elif long:
        wave = eis['wave_area_lw'][0]

    if isinstance(wave, (int, float)):
        wave = np.array([wave])

    nWave = len(wave)
    ea_out = np.zeros(nWave)

    for i in range(nWave):
        band = eis_get_band(wave[i])
        if band == 'SW':
            w = eis['wave_knots_sw'][0]
            e = np.log(ea_knots_SW)
        elif band == 'LW':
            w = eis['wave_knots_lw'][0]
            e = np.log(ea_knots_LW)
        else:
            print(f"WAVELENGTH OUT OF BOUNDS {wave[i]}")
            continue

        interp_func = interp1d(w, e, kind='linear')
        ea_out[i] = np.exp(interp_func(wave[i]))

    if nWave == 1:
        ea_out = ea_out[0]

    return ea_out

def get_time_tai(date_string):
    """
    Converts a date string to TAI time for use in calibration.
    """
    idl_ref_epoch = datetime.datetime(1979, 1, 1)
    unix_epoch = datetime.datetime(1970, 1, 1)
    epoch_diff = (idl_ref_epoch - unix_epoch).total_seconds()
    date_object = datetime.datetime.fromisoformat(date_string)
    unix_timestamp = date_object.timestamp()
    idl_timestamp = unix_timestamp - epoch_diff + 3600
    return idl_timestamp

def eis_get_band(wave):
    """
    Returns the wavelength band (short or long) for a given wavelength.
    """
    if 165 <= wave <= 213:
        return 'SW'
    elif 245 <= wave <= 292:
        return 'LW'
    else:
        return ''

# Function to read the calibration file for 2014 data
def read_calib_file(file_path=WARREN_2014_FILE):
    """
    Reads the calibration file for the 2014 Warren calibration.
    """
    return readsav(file_path)['eis']

# Function to read the effective area for short and long bands from preflight data
def eis_effective_area_read(short=False, long=False):
    """
    Reads the effective area preflight calibration for the specified band.
    """
    if short:
        preflight = readsav(PREFLIGHT_SHORT_FILE)
    elif long:
        preflight = readsav(PREFLIGHT_LONG_FILE)
    wave = preflight['wave']
    ea = preflight['ea']
    return wave, ea

def is_eis_wavelength(input_wave):
    """
    Determines if a wavelength is in the short (SW) or long (LW) EIS band.
    """
    wave_sw_min = 165
    wave_sw_max = 213
    wave_lw_min = 245
    wave_lw_max = 292

    short = long = False

    if wave_sw_min <= input_wave <= wave_sw_max:
        short = True
    if wave_lw_min <= input_wave <= wave_lw_max:
        long = True

    return short, long