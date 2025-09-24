import os
import re
from tqdm import tqdm
from astropy.io import fits
import numpy as np
from scipy.io import readsav
import sunpy.map

def load_clean_map(filename):
    """Load a FITS file, replace NaNs with 0, and return a SunPy map."""
    data, header = fits.getdata(filename, header=True)
    data = np.array(data, dtype=float)
    data[np.isnan(data)] = 0
    return sunpy.map.Map((data, header))

def align_fip_eis(
    fip_path='./fip',
    int_vel_ntv_path='./fitted_data',
    density_path='./fitted_data',
    chi2_threshold=3,
    map_types=['int', 'vel', 'ntv', 'fip', 'err_int', 'err_vel',
               'density', 'ntmap', 'mcmc_logt', 'hmi', 'chi2']
):
    def find_files(directory, extension):
        """Recursively find all files with a given extension."""
        return sorted([
            os.path.join(root, file)
            for root, _, files in os.walk(directory)
            for file in files if file.endswith(extension)
        ])

    def make_aligned_filename(file_path, meas, line_tag=None):
        """
        Construct an aligned filename without duplicating line tags.
        Example:
            eis_20151019_131642.fe_15_284_160.int.fits
            -> eis_20151019_131642.fe_15_284_160.aligned_int.fits
        """
        dir_name = os.path.dirname(file_path)
        base_name = os.path.basename(file_path)
    
        # remove known suffixes
        for suffix in [
            '.fits', '.sav',
            '.int', '.vel', '.ntv', '.err_int', '.err_vel', '.chi2', '.density',
            '_aligned_int', '_aligned_vel', '_aligned_ntv',
            '_aligned_err_int', '_aligned_err_vel', '_aligned_density',
            '_aligned_fip', '_aligned_chi2'
        ]:
            if base_name.endswith(suffix):
                base_name = base_name[: -len(suffix)]
    
        return os.path.join(dir_name, f"{base_name}.aligned_{meas}.fits")


    # --- Find all input files ---
    sav_files = find_files(fip_path, '.sav')
    int_files = find_files(int_vel_ntv_path, '.int.fits')
    vel_files = find_files(int_vel_ntv_path, '.vel.fits')
    ntv_files = find_files(int_vel_ntv_path, '.ntv.fits')
    err_int_files = find_files(int_vel_ntv_path, '.err_int.fits')
    err_vel_files = find_files(int_vel_ntv_path, '.err_vel.fits')
    chi2_files = find_files(int_vel_ntv_path, '.chi2.fits') if 'chi2' in map_types else []
    density_files = find_files(density_path, '_density.fits') if 'density' in map_types else []

    print("FIP SAV files:", sav_files)
    print("INT FITS files:", int_files)

    # --- Pair files robustly ---
    paired_files = []
    for sav in sav_files:
        base = os.path.splitext(os.path.basename(sav))[0].replace('_aligned_maps', '').strip()

        # extract timestamp from SAV filename (e.g. 20151018_102719)
        timestamp = base  

        def extract_line_tag(filename, timestamp):
            # works for .int/.vel/.ntv/.err_int/.err_vel/.chi2/.density
            m = re.search(rf"{timestamp}\.(.+?)\.(int|vel|ntv|err_int|err_vel|chi2|density)\.fits$", filename)
            if m:
                return m.group(1)  # spectral line string like fe_12_195_119
            return None

        # collect all line tags
        all_lines = set()
        for flist in [int_files, vel_files, ntv_files, err_int_files, err_vel_files, chi2_files, density_files]:
            for f in flist:
                tag = extract_line_tag(os.path.basename(f), timestamp)
                if tag:
                    all_lines.add(tag)

        # build pairs for each line
        for line_tag in all_lines:
            def match_line(file_list, ext):
                return next((f for f in file_list if f"{timestamp}.{line_tag}.{ext}.fits" in f), None)

            int_file     = match_line(int_files, "int")
            vel_file     = match_line(vel_files, "vel")
            ntv_file     = match_line(ntv_files, "ntv")
            err_int_file = match_line(err_int_files, "err_int")
            err_vel_file = match_line(err_vel_files, "err_vel")
            chi2_file    = match_line(chi2_files, "chi2") if chi2_files else None
            density_file = match_line(density_files, "density") if density_files else None

            skip = False
            for map_type, file_obj in zip(
                ['int','vel','ntv','err_int','err_vel','density','chi2'],
                [int_file, vel_file, ntv_file, err_int_file, err_vel_file, density_file, chi2_file]
            ):
                if map_type in map_types and file_obj is None:
                    skip = True
                    print(f"Skipping {sav} [{line_tag}]: missing {map_type}")
                    break
            if skip:
                continue

            paired_files.append((sav, line_tag, int_file, vel_file, ntv_file,
                                 err_int_file, err_vel_file, density_file, chi2_file))

    created_files = []

    for fip_file, line_tag, int_file, vel_file, ntv_file, err_int_file, err_vel_file, density_file, chi2_file in tqdm(paired_files):
        maps = {}
        if 'int' in map_types and int_file: maps['int'] = load_clean_map(int_file)
        if 'vel' in map_types and vel_file: maps['vel'] = load_clean_map(vel_file)
        if 'ntv' in map_types and ntv_file: maps['ntv'] = load_clean_map(ntv_file)
        if 'err_int' in map_types and err_int_file: maps['err_int'] = load_clean_map(err_int_file)
        if 'err_vel' in map_types and err_vel_file: maps['err_vel'] = load_clean_map(err_vel_file)
        if 'density' in map_types and density_file: maps['density'] = load_clean_map(density_file)
        if 'chi2' in map_types and chi2_file: maps['chi2'] = load_clean_map(chi2_file)

        # --- Apply chi2 threshold ---
        if 'chi2' in maps:
            bad_mask = maps['chi2'].data > chi2_threshold
            for meas, map_obj in maps.items():
                map_obj.data[bad_mask] = 0

        # --- Ensure no NaNs ---
        for meas, map_obj in maps.items():
            map_obj.data[np.isnan(map_obj.data)] = 0

        # --- Load FIP SAV ---
        fip_sav = readsav(fip_file)
        sav_keys = fip_sav.keys()
        fip_key = next(key for key in ['fip_map', 'fipmap'] if key in sav_keys)

        sav_x = fip_sav[fip_key][0][1]
        sav_y = fip_sav[fip_key][0][2]

        first_map = list(maps.values())[0] if maps else None

        # --- Align EIS maps ---
        if first_map:
            x_shift = sav_x - first_map.meta['xcen']
            y_shift = sav_y - first_map.meta['ycen']
            for meas, map_obj in maps.items():
                map_obj.meta['xcen'] = sav_x
                map_obj.meta['ycen'] = sav_y
                map_obj.meta['crval1'] += x_shift
                map_obj.meta['crval2'] += y_shift
                map_obj.meta['measrmnt'] = meas

        # --- Save FIP-derived maps (ntmap, mcmc_logt, fip, hmi) ---
        for fip_measure in ['ntmap', 'mcmc_logt', 'fip', 'hmi']:
            if fip_measure in map_types and fip_measure in fip_sav:
                data_array = np.array(fip_sav[fip_measure][0][0], dtype=float)
                data_array[np.isnan(data_array)] = 0

                if first_map:
                    meta = first_map.meta.copy()
                else:
                    meta = {
                        'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec',
                        'CTYPE1': 'HPLN-TAN', 'CTYPE2': 'HPLT-TAN',
                        'CRPIX1': data_array.shape[1]/2, 'CRPIX2': data_array.shape[0]/2,
                        'CRVAL1': 0, 'CRVAL2': 0
                    }

                meta['xcen'] = sav_x
                meta['ycen'] = sav_y
                meta['measrmnt'] = fip_measure

                fip_map_obj = sunpy.map.Map(data_array, meta)
                out_file = make_aligned_filename(fip_file, fip_measure, line_tag)
                fip_map_obj.save(out_file, overwrite=True)
                created_files.append(out_file)

        # --- Save EIS maps ---
        for meas, map_obj in maps.items():
            file_var = locals()[f"{meas}_file"] if meas != 'density' else density_file
            if file_var is None: continue
            out_file = make_aligned_filename(file_var, meas, line_tag)
            map_obj.save(out_file, overwrite=True)
            created_files.append(out_file)

    return created_files
