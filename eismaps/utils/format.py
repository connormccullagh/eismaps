import eismaps.utils.roman_numerals as roman_numerals

def change_line_format(line):
    # input: e.g. "Fe XIV 264.700" or "S X 258.375" or "Si V 228.375"
    # output: "fe_14_264_700" or "s__10_258_375" or "si_05_228_375"

    elements, ion, wavelength = line.split()
    formatted_element = elements.lower()

    if len(formatted_element) == 1:
        formatted_element += '_'

    formatted_ion = roman_numerals.roman_to_int(ion)
    formatted_wavelength = wavelength.replace('.', '_')

    # if the formatted ion integer is only one digit, add a zero before
    if formatted_ion < 10:
        formatted_ion = f'0{formatted_ion}'

    return f"{formatted_element}_{formatted_ion}_{formatted_wavelength}"