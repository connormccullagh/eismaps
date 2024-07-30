import os

def main(ext, dir):
    collected_files = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(f"{ext}"):
                collected_files.append(os.path.join(root, file))
    if len(collected_files) == 0:
        print(f'No {ext} files found in {dir}. Returning None.')
        return None
    collected_files.sort()
    return collected_files