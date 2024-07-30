import pandas as pd
import pkg_resources
import eismaps.utils.roman_numerals as roman_numerals

def width2velocity(element, ion):
    
    element = element.capitalize()  # Ensure the element is in title case
    ion = roman_numerals.int_to_roman(ion)  # Convert the ion number to roman numerals

    dat_path = pkg_resources.resource_filename('eismaps', 'utils/eis_width2velocity.dat')
    ion_equilibriums = pd.read_fwf(dat_path,widths=[7,10,10,10],skiprows=9, header=None)
    ion_equilibriums.columns = ['ELEM', 'ION', 'T_MAX', 'V_THERM']

    # Now access using the exact column names
    t_max = ion_equilibriums.loc[(ion_equilibriums['ELEM'] == element) & (ion_equilibriums['ION'] == ion), 'T_MAX'].values
    v_therm = ion_equilibriums.loc[(ion_equilibriums['ELEM'] == element) & (ion_equilibriums['ION'] == ion), 'V_THERM'].values

    if not t_max or not v_therm:
        raise ValueError(f"No data found for element {element} and ion {ion}")
    
    # turn strings into floats
    t_max = float(t_max[0])
    v_therm = float(v_therm[0])

    return t_max, v_therm