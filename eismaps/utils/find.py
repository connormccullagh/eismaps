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

    if len(collected_files) == 0:
        print(f'No {ext} files found in {dir}. Returning None.')
        return None

    collected_files.sort()
    return collected_files
