# import os

# def main(ext, dir):
#     collected_files = []
#     for root, dirs, files in os.walk(dir):
#         for file in files:
#             if file.endswith(f"{ext}"):
#                 collected_files.append(os.path.join(root, file))
#     if len(collected_files) == 0:
#         print(f'No {ext} files found in {dir}. Returning None.')
#         return None
#     collected_files.sort()
#     return collected_files

import os
from datetime import datetime

def main(ext, folder, chosen_wavelength = None, unique_times=False, start_str=None, end_str=None):
    collected_files = []
    seen_times = set()

    # Convert time range strings to datetime, if provided
    start_dt = datetime.strptime(start_str, "%Y%m%d_%H%M%S") if start_str else None
    end_dt   = datetime.strptime(end_str, "%Y%m%d_%H%M%S") if end_str else None

    if chosen_wavelength:
        # Split integer and decimal part
        int_part = int(chosen_wavelength)           # 195
        dec_part = int(round((chosen_wavelength % 1) * 1000))  # 119
        
        line_str = f"fe_12_{int_part}_{dec_part}"   # fe_12_195_119

    for root, _, files in os.walk(folder):
        for fname in files:
            if not fname.endswith(ext):
                continue

            # Extract timestamp depending on filename pattern
            ts_str = None
            if fname.startswith("eis_"):  
                # eis_20230101_120000.ext
                ts_str = fname[4:].split('.')[0]
            else:  
                # generic pattern: something_something_TIMESTAMP.ext
                parts = fname.split('_')
                if len(parts) > 2:
                    ts_str = parts[2].split('.')[0]

            file_dt = None
            if ts_str:
                try:
                    file_dt = datetime.strptime(ts_str, "%Y%m%d_%H%M%S")
                except ValueError:
                    pass

            # Apply unique timestamp filter
            if unique_times and ts_str:
                if ts_str in seen_times:
                    continue
                seen_times.add(ts_str)

            # Apply time range filter
            if file_dt:
                if start_dt and file_dt < start_dt:
                    continue
                if end_dt and file_dt > end_dt:
                    continue

            if chosen_wavelength != None:
                if line_str not in fname:
                    continue

            collected_files.append(os.path.join(root, fname))

    collected_files.sort()
    return collected_files
