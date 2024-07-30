import eispac
import os
import sys
import time
from eismaps.utils.format import change_line_format

def template_summary():
    eispac_dir = os.path.dirname(eispac.__file__)
    template_dir = os.path.join(eispac_dir, 'data/templates')
    template_files = [os.path.join(template_dir, f) for f in os.listdir(template_dir) if f.endswith('.h5')]

    # Dictionary to hold line IDs and corresponding templates
    line_templates = {}

    for template_file in template_files:
        template = eispac.read_template(template_file)
        for line_id in template.template['line_ids']:
            if line_id in line_templates:
                line_templates[line_id].append(template_file)
            else:
                line_templates[line_id] = [template_file]

    # Displaying results
    for line, templates in line_templates.items():
        if len(templates) > 1:
            print(f"Line {line} is duplicated in the following templates:")
            for template in templates:
                print(f"- {os.path.basename(template)}")
            print()  # Blank line for better readability

    # Total lines and unique lines
    total_lines = sum(len(templates) for templates in line_templates.values())
    unique_lines = len(line_templates)
    print(f"Total number of lines: {total_lines}")
    print(f"Number of unique lines: {unique_lines}")

    return line_templates

# def batch(files, ncpu='max', save=True, line=None):

#     for file in files:
#         file_dir = os.path.dirname(file)
#         file_name = os.path.basename(file)

#         file_date = file_name.split('.')[0].split('_')[1] # YYYYMMDD
#         file_time = file_name.split('.')[0].split('_')[2] # HHMMSS

#         wininfo = eispac.read_wininfo(file) # windows observed in the file
#         twin = len(wininfo) # total number of windows

#         template_search_res = eispac.core.match_templates(file) # templates matching each window for the file
#                                                                 # returns a list of lists, one master list for each window

#         # check the number of windows matches the returned list
#         if len(template_search_res) != twin:
#             print(f"Something went wrong when using eispac's template matcher. The number of results returned doesn't match the number of observation windows. Exiting.")
#             sys.exit()

#         iwin = 0 # window counter
#         for templates_for_iwin in template_search_res: # templates_for_iwin is a list of templates for the current window
#             print(f'Fitting window {iwin} of {twin-1}')

#             if len(templates_for_iwin) == 0:
#                 print(f"No templates found for window {iwin}. Skipping.")
#                 iwin += 1
#                 continue

#             for template_path in templates_for_iwin:
#                 print(f"Using template {template_path}")

#                 template_name = os.path.basename(template_path)
#                 template_name_line = template_name.split('.')[0] # get the line name that the template is named after (perhaps not the only line it fits)

#                 template = eispac.read_template(template_path)

#                 # check if any of the lines in the template match the input lines
#                 template_lines = [change_line_format(l) for l in template.template['line_ids']]
#                 if line is not None and not any(tl in line for tl in template_lines):
#                     print(f"Skipping template {template_path} as it doesn't contain the specified emission lines.")
#                     continue

#                 cube = eispac.read_cube(file, iwin) # load the iwin window

#                 fit = eispac.fit_spectra(cube, template, ncpu=ncpu) # fit the window with the chosen template

#                 saved_fits = eispac.save_fit(fit, save_dir=file_dir) # save the fit and get the output file name

#                 # if saved_fits is not a list, make it a one element list
#                 if not isinstance(saved_fits, list):
#                     saved_fits = [saved_fits]

#                 icomponent = 0
#                 for saved_fit in saved_fits:
#                     emission_line = change_line_format(template.template['line_ids'][icomponent]) # the line that was fitted
#                     new_out_name = f'{file_dir}/eis_{file_date}_{file_time}.{emission_line}.{template_name_line}.fit.h5' # the new name of the output file
#                     # if the emission_line is in line, then rename, else delete
#                     if line is not None and emission_line not in line:
#                         os.remove(saved_fit)
#                     else:
#                         os.rename(saved_fit, new_out_name)
#                     icomponent += 1

#                 # the final fit is saved as:
#                 # eis_YYYYMMSS_HHMMSS.{EMISSION_LINE}.{TEMPLATE_NAME (no component)}.fit.h5
#                 # EMISSION_LINE in format: ee_ii_www_www where ee is the element, ii is the ion, www_www is the wavelength range. If ee is single letter, add an extra underscore
                
#             iwin += 1

def fit_specific_line(file, iwin, template, best_match, best_match_index, ncpu='max', save=True, output_dir=None, output_dir_tree=False):

    if output_dir is None:
        print('No output directory specified. Saving to the same directory as the input file.')
        output_dir = os.path.dirname(file)
    else:
        if output_dir_tree:
            file_date = os.path.basename(file).split('.')[0].split('_')[1]
            output_dir = os.path.join(output_dir, file_date[:4], file_date[4:6], file_date[6:8])
        os.makedirs(output_dir, exist_ok=True)

    if save: # if it exists, skip it
        new_filename = os.path.join(output_dir, f"{os.path.basename(file).split('.')[0]}.{best_match.replace(' ', '_').replace('.','_').lower()}.fit.h5")
        if os.path.exists(new_filename):
            print(f"Fit already exists for {best_match}. Skipping.")
            return

        cube = eispac.read_cube(file, iwin)
        template = eispac.read_template(template)
        fit = eispac.fit_spectra(cube, template, ncpu=ncpu)

        saved_fits = eispac.save_fit(fit, save_dir=output_dir)

        # ensure saved_fits is a list even if only one file is returned
        if not isinstance(saved_fits, list):
            saved_fits = [saved_fits]

        # loop over all saved fits, delete the unwanted ones, and rename the one we want to keep
        for index, saved_fit in enumerate(saved_fits):
            if index == best_match_index: # this is the fit file we want to keep
                os.rename(saved_fit, new_filename)
                print(f"Fit saved to {new_filename} (by renaming {saved_fit})")
            else:
                # delete the unwanted fit files
                os.remove(saved_fit)
                print(f"Deleted {saved_fit} as component not needed")

    else:
        cube = eispac.read_cube(file, iwin)
        template = eispac.read_template(template)
        fit = eispac.fit_spectra(cube, template, ncpu=ncpu)
        print(f"Fit for {best_match} completed, but not saved.")

def batch(files, ncpu='max', save=True, output_dir=None, line=None, output_dir_tree=False):
    ## todo: add a line argument to specify the line to fit, which can include lines not from the list of window names (search through all possible templates and line_ids)

    for file in files: # cycle through all the files
        wininfo = eispac.read_wininfo(file)
        templates = eispac.core.match_templates(file) # match templates for the entire file

        for iwin, template_group in enumerate(templates): # cycle through the windows in the file

          wininfo_ion = wininfo[iwin][1].split()[0] + ' ' + wininfo[iwin][1].split()[1]
          wininfo_wavelength = wininfo[iwin][1].split()[2]

          if len(template_group) == 0:
            print(f"No templates found for window {iwin}. Skipping.")
            continue

          # initialize variables to track the best match
          best_match = None
          best_match_index = None
          best_template = None
          smallest_offset = float('inf')
          offset_threadhold = 0.15

          for template in template_group: # cycle through the templates which match the current window
            loaded_template = eispac.read_template(template)

            for line_id_index, line_id in enumerate(loaded_template.template['line_ids']): # cycle through each line in the template

              # if the template line_id contains the "unknown" string, skip it
              if 'unknown' in line_id.lower():
                # skipped as the splitting may not work as expected
                continue

              template_ion = line_id.split()[0] + ' ' + line_id.split()[1]
              template_wavelength = line_id.split()[2]

              if template_ion.lower() == wininfo_ion.lower():
                offset = abs(float(template_wavelength) - float(wininfo_wavelength))
                if offset < smallest_offset:
                  smallest_offset = offset
                  if smallest_offset < offset_threadhold:
                    best_match = line_id
                    best_match_index = line_id_index
                    best_template = template
                    # todo:
                    # if two templates offer good matches, the last one will be chosen
                    # this could be improved by choosing the one with a better fitting template (e.g. if it is closer to the middle of the template region)

          if best_match:
            print(f"Found a match for {wininfo[iwin][1]}: {best_match}")
            fit_specific_line(file, iwin, best_template, best_match, best_match_index, ncpu=ncpu, save=save, output_dir=output_dir, output_dir_tree=output_dir_tree)
          else:
            print(f"No good match found for {wininfo[iwin][1]}")