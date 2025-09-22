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

def main(ext, dir, unique_times=False):
    collected_files = []
    seen_times = set()

    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(f"{ext}"):
                if unique_times:
                    # Extract the timestamp from the filename
                    time_stamp = file.split('_')[2]
                    if time_stamp not in seen_times:
                        collected_files.append(os.path.join(root, file))
                        seen_times.add(time_stamp)
                else:
                    collected_files.append(os.path.join(root, file))

    collected_files.sort()
    return collected_files
    
def find_files_in_range(ext, folder, start_str, end_str):
    from datetime import datetime
    import datetime as dt
    import os
    
    start_dt = datetime.strptime(start_str, "%Y%m%d_%H%M%S")
    end_dt   = datetime.strptime(end_str, "%Y%m%d_%H%M%S")
    files_in_range = []

    for fname in os.listdir(folder):
        if fname.startswith('eis_') and fname.endswith(ext):
            ts_str = fname[4:].split('.')[0]
            try:
                file_dt = datetime.strptime(ts_str, "%Y%m%d_%H%M%S")
            except ValueError:
                continue
            if start_dt <= file_dt <= end_dt:
                files_in_range.append(os.path.join(folder, fname))

    files_in_range.sort()
    return files_in_range
