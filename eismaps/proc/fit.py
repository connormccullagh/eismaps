import eispac
import os
import eismaps.utils.roman_numerals as roman_numerals
from eismaps.utils.format import change_line_format

def fit_specific_line(file, iwin, template, line_label, ncpu='max', save=True, output_dir=None, output_dir_tree=False):
    # Determine the output directory
    if output_dir is None:
        print('No output directory specified. Saving to the same directory as the input file.')
        output_dir = os.path.dirname(file)
    else:
        if output_dir_tree:
            # Create a directory tree based on the file date
            file_date = os.path.basename(file).split('.')[0].split('_')[1]
            output_dir = os.path.join(output_dir, file_date[:4], file_date[4:6], file_date[6:8])
        os.makedirs(output_dir, exist_ok=True)

    # Check if the fit already exists
    new_filename = os.path.join(output_dir, f"{os.path.basename(file).split('.')[0]}.{line_label.replace(' ', '_').replace('.', '_').lower()}.fit.h5")
    if save and os.path.exists(new_filename):
        print(f"Fit already exists for {line_label}. Skipping.")
        return

    # Read the data cube and the template
    cube = eispac.read_cube(file, iwin)
    template = eispac.read_template(template)
    fit = eispac.fit_spectra(cube, template, ncpu=ncpu)

    if save:
        # Save the fit result
        saved_fits = eispac.save_fit(fit, save_dir=output_dir)

        # Ensure saved_fits is a list even if only one file is returned
        if not isinstance(saved_fits, list):
            saved_fits = [saved_fits]

        # Loop over all saved fits, delete the unwanted ones, and rename the one we want to keep
        for saved_fit in saved_fits:
            if line_label in os.path.basename(saved_fit):
                os.rename(saved_fit, new_filename)
                print(f"Fit saved to {new_filename} (by renaming {saved_fit})")
            else:
                os.remove(saved_fit)
                print(f"Deleted {saved_fit} as component not needed")

    else:
        print(f"Fit for {line_label} complete but not saved.")

def batch(files, ncpu='max', save=True, output_dir=None, output_dir_tree=False, lock_to_window=False):
    for file in files:  # Cycle through all the files
        wininfo = eispac.read_wininfo(file)
        templates = eispac.core.match_templates(file)  # Match templates for the entire file

        for iwin, template_group in enumerate(templates):  # Cycle through the windows in the file
            if len(template_group) == 0:
                print(f"No templates found for window {iwin}. Skipping.")
                continue

            templates_to_fit = []

            if lock_to_window:  # If the lock_to_window flag is set, only fit one component per window (named after the windows)

                # If there is a template with 1 component, use that
                if any('1c' in template for template in template_group):
                    for template in template_group:
                        if template.split('.')[-3].replace('c', '') == '1':
                            templates_to_fit.append(template)
                            break
                            ### TODO: If more than one template has 1 component, optimise selection ###

                # Otherwise just choose the first template
                templates_to_fit.append(template_group[0])
                ### TODO: Optimise selection ###

                assert len(templates_to_fit) == 1, f"More than one template selected for window {iwin} with lock_to_window=True"

            else:  # Fit as many lines as possible

                # If there are two templates for the same line, use the one with the highest number of components
                for template in template_group:
                    if template in templates_to_fit:
                        for i, t in enumerate(templates_to_fit):
                            if t == template:
                                if int(t.split('.')[-3].replace('c', '')) < int(template.split('.')[-3].replace('c', '')):  # If the current template has more components than the one already in the list, replace it
                                    templates_to_fit[i] = template
                    else:
                        templates_to_fit.append(template)

            for template_to_fit in templates_to_fit:  # Cycle through the templates to fit
                
                if lock_to_window:
                    # Take the fit name from the window name
                    line_label = change_line_format(wininfo[iwin]['line_id'])
                else:
                    # Take the fit name from the template name
                    line_label = os.path.basename(template_to_fit).split('.')[-4]

                fit_specific_line(file, iwin, template_to_fit, line_label, ncpu=ncpu, save=save, output_dir=output_dir, output_dir_tree=output_dir_tree)