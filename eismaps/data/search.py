import requests
from bs4 import BeautifulSoup
from datetime import datetime, timedelta

SOURCE_URLS = {
    'nrl': "https://eis.nrl.navy.mil/level1/hdf5",
    'mssl': "https://vsolar.mssl.ucl.ac.uk/eispac/hdf5"
}

def main(start_datetime, end_datetime, source='nrl', base_url=None):
    """
    Fetch and filter files from the EIS website within a given datetime range.

    :param base_url: Base URL of the EIS data source.
    :param start_datetime: Start datetime as a datetime object.
    :param end_datetime: End datetime as a datetime object.
    :return: List of filtered file URLs.
    """

    if source in SOURCE_URLS:
        if base_url is not None:
            print("Overriding custom url for selected source.")
            print("To specify a custom url, use the 'custom' source.")
        base_url = SOURCE_URLS[source]
    elif source == 'custom':
        if base_url is None:
            raise ValueError("Custom source requires a base URL.")
    
    print(f'Searching for files between {start_datetime} and {end_datetime}...')
    print(f'Fetching files from {base_url}...')
    
    filtered_files = []

    # Generate the date range
    current_date = start_datetime.date()
    while current_date <= end_datetime.date():
        date_str = current_date.strftime("%Y/%m/%d")
        url = f"{base_url}/{date_str}/"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')

        # List files for the current date
        for link in soup.find_all('a'):
            href = link.get('href')
            if href and href.endswith('.h5'):
                file_url = f"{url}{href}"
                filename = href
                # Filter by datetime
                datetime_str = '_'.join(filename.split('_')[1:3]).split('.')[0]
                file_datetime = datetime.strptime(datetime_str, "%Y%m%d_%H%M%S")
                if start_datetime <= file_datetime <= end_datetime:
                    filtered_files.append(file_url)

        current_date += timedelta(days=1)

    return filtered_files

def check_source(source, base_url=None):
    """
    Check if an source is functioning correctly. Mostly useful for custom sources.

    :param source: The source to check.
    :return: True if the source is functioning correctly, False otherwise.
    """
    if source in SOURCE_URLS:
        if base_url is not None:
            print("Overriding custom url for selected source.")
            print("To specify a custom url, use the 'custom' source.")
        base_url = SOURCE_URLS[source]

    try:
        response = requests.get(base_url)
        
        # Check if the request was successful
        if response.status_code != 200:
            print(f"Failed to access {base_url}: Status code {response.status_code}")
            return False
      
        print(f"Source at {base_url} appears to be functioning correctly.")
        return True

    except requests.exceptions.RequestException as e:
        print(f"Error checking source status: {e}")
        return False