SOURCE_URLS = { # move this elsewhere as multiple scripts reference it
    'nrl': "https://eis.nrl.navy.mil/level1/hdf5",
    'mssl': "https://vsolar.mssl.ucl.ac.uk/eispac/hdf5"
}

def eispac(file_urls, source='nrl', local_top='data_eis', datetree=False):
    """
    Download files using eispac.download.download_hdf5_data based on a list of URLs.

    :param file_urls: List of URLs to the files to be downloaded.
    """
    from eispac.download import download_hdf5_data

    for file_url in file_urls:
        filename = os.path.basename(file_url)
        print(f"Downloading {filename}...")
        # download_hdf5_data(filename=filename, source=source, datetree=datetree, local_top=local_top)
        download_hdf5_data(filename=filename, datetree=datetree, local_top=local_top)
        print(f"Finished downloading {filename}.")

import os
import requests

def eismaps(file_urls, source='same', local_top='data_eis', datatree=False):
    """
    Download files from a list of URLs and save them to the specified directory.
    Optionally organize the files in a directory structure based on their date.

    :param file_urls: List of URLs to the files to be downloaded.
    :param source: The source from which the files are being downloaded.
    :param local_top: Directory where files should be saved.
    :param datatree: If True, save files in a YYYY/MM/DD directory structure; otherwise, save in the base folder.
    """
    if source != 'same' and source not in SOURCE_URLS:
        raise ValueError(f"Invalid source: {source}")

    if source != 'same':
        source_url = SOURCE_URLS[source]

    for file_url in file_urls:
        filename = file_url.split('/')[-1]

        if source != 'same':
            local_path_parts = file_url.split('/')[-4:]
            local_path = '/'.join(local_path_parts)
            file_url = soirce_url.rstrip('/') + '/' + local_path

        if datatree:
            # Extracts the date from the URL and creates the directory structure
            date_path = '/'.join(file_url.split('/')[-4:-1]) 
            target_dir = os.path.join(local_top, date_path)
        else:
            target_dir = local_top

        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        target_path = os.path.join(target_dir, filename)

        # Check if file already exists
        if not os.path.exists(target_path):
            print(f"Downloading {file_url} to {target_path}...")
            response = requests.get(file_url, stream=True)
            if response.status_code == 200:
                with open(target_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:  # filter out keep-alive new chunks
                            f.write(chunk)
                print(f"Finished downloading {filename}.")
            else:
                print(f"Failed to download {filename}. Status code: {response.status_code}")
        else:
            print(f"{filename} already exists in {target_dir}. Skipping download.")
    return [target_path for file_url in file_urls]
