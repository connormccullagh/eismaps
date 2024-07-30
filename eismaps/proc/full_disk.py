import sunpy.map
import numpy as np
import time
import math
from astropy.coordinates import SkyCoord
from sunpy.coordinates import RotatedSunFrame
import astropy.units as u
from datetime import datetime
import matplotlib.pyplot as plt
import os
from sunpy.coordinates import frames

def check_fd(files):
    # that that the set of files specified all belong to the same full disk observation
    return None

def main(files, high_res=False):

    # get the maximum extent of the full disk
    fd_size = 0
    for file in files:
        map = sunpy.map.Map(file)
        if ((abs(map.meta['xcen'])+(map.meta['fovx']/2))*2 > fd_size):
            fd_size = (abs(map.meta['xcen'])+(map.meta['fovx']/2))*2
        if ((abs(map.meta['ycen'])+(map.meta['fovy']/2))*2 > fd_size):
            fd_size = (abs(map.meta['ycen'])+(map.meta['fovy']/2))*2

    first_map = sunpy.map.Map(files[0])
    last_map = sunpy.map.Map(files[-1])
    
    map_dx = first_map.meta['cdelt1']
    map_dy = first_map.meta['cdelt2']

    measurement = files[0].split('.')[-2] # can be 'int', 'vel', 'wid', 'ntv' for now
    line_id = files[0].split('.')[-3]

    # double check all the maps have the same dx and dy #TODO add checks of other parameters (e.g. the measurement, line)
    for file in files:
        map = sunpy.map.Map(file)
        if (map.meta['cdelt1'] != map_dx) or (map.meta['cdelt2'] != map_dy):
            raise ValueError('Error: maps do not have the same dx and dy')

    fd_width = round(fd_size/map_dx)
    fd_height = round(fd_size/map_dy)

    fd_data = np.full((fd_height, fd_width), 0.0)
    fd_mask = np.full((fd_height, fd_width), 0)
    fd_mask_flat = np.full((fd_height, fd_width), 0)
    fd_overlap_map = np.full((fd_height, fd_width), np.nan)

    # create a list of the times of the maps
    map_times = []
    for file in files:
      map = sunpy.map.Map(file)
      map_times.append([file, map.meta['date_beg'], map.meta['date_end']])
    map_times = np.array(map_times)

    fd_header = {
        'crpix1': 1, # reference pixel (bottom left)
        'crpix2': 1,
        'cdelt1': map_dx,
        'cdelt2': map_dy,
        'cunit1': 'arcsec',
        'cunit2': 'arcsec',
        'crval1': -fd_size/2, # location of reference coordinate
        'crval2': -fd_size/2,
        'ctype1': 'HPLN-TAN',
        'ctype2': 'HPLT-TAN',
        'date_obs': first_map.meta['date_beg'],
        'date_beg': first_map.meta['date_beg'],
        'date_end': last_map.meta['date_end'],
        'dsun_obs': first_map.meta['dsun_obs'],
        'hgln_obs': first_map.meta['hgln_obs'],
        'hglt_obs': first_map.meta['hglt_obs'],
        'instrume': 'EIS',
        'measrmnt': 'int', # TEMP for changing
        'telescop': 'Hinode',
        'xcen': 0,
        'ycen': 0,
        'history': 'Created by eismaps.proc.full_disk.main'
    }

    fd_map = sunpy.map.Map(fd_data, fd_header)

    # insert the maps into the full disk map
    imap = 0
    for file in files:

        map = sunpy.map.Map(file)
        
        map_x_init = map.meta['xcen']
        map_y_init = map.meta['ycen']
        map_rad = np.sqrt(map_x_init**2 + map_y_init**2)

        # if the centre is outside the limb, the sunpy differential rotation can't be applied. We move these maps slightly inside the limb, apply the rotation, then put them back
        if map_rad > 950:
            map_angle = math.atan(map_y_init/map_x_init)
            map_y_init_shifted = 950 * math.sin(map_angle)
            map_x_init_shifted = 950 * math.cos(map_angle)
            map_y_shift = map_y_init - map_y_init_shifted
            map_x_shift = map_x_init - map_x_init_shifted
            map_y_init = map_y_init_shifted
            map_x_init = map_x_init_shifted
        else:
            map_y_shift = 0
            map_x_shift = 0

        # apply the differential rotation
        point = SkyCoord(map_x_init*u.arcsec, map_y_init*u.arcsec, frame=map.coordinate_frame)
        fd_map_beg_time = datetime.strptime(fd_map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
        map_beg_time = datetime.strptime(map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
        diffrot_point = RotatedSunFrame(base=point, duration=((fd_map_beg_time-map_beg_time).seconds)*u.second)
        transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)

        new_x = str(transformed_diffrot_point.Tx)
        new_y = str(transformed_diffrot_point.Ty)
        new_x = new_x.replace('arcsec', '')
        new_y = new_y.replace('arcsec', '')
        new_x = float(new_x)
        new_y = float(new_y)

        # put the map back
        map.meta['xcen'] = new_x + map_x_shift
        map.meta['ycen'] = new_y + map_y_shift

        print(f'this map shifted {new_y-map_y_init} in y and {new_x-map_x_init} in x for time')
        
        # calculate the boundaries of the map in arcsecs and convert to array indices
        map_height, map_width = map.data.shape
        map_left = round( (map.meta['xcen']/map_dx - map_width/2) + fd_width/2 )
        map_right = round( map_left + map_width )
        map_bottom = round( (map.meta['ycen']/map_dy - map_height/2) + fd_height/2 )
        map_top = round( map_bottom + map_height )

        # update the mask to show the number of contributions being made to the particular cells
        fd_mask_flat[map_bottom:map_top, map_left:map_right] = fd_mask_flat[map_bottom:map_top, map_left:map_right] + 1

        # insert the data to the full disk map where there are no overlaps
        fd_data[map_bottom:map_top, map_left:map_right] = np.where(fd_mask_flat[map_bottom:map_top, map_left:map_right]<=1,fd_data[map_bottom:map_top, map_left:map_right] + map.data, fd_data[map_bottom:map_top, map_left:map_right])

        # detect if merging is required
        if np.any(fd_mask_flat[map_bottom:map_top, map_left:map_right] == 2):
            print(f'merging')

            # overlaps are here
            overlap_height_indices = np.where(fd_mask_flat==2)[0]
            overlap_width_indices = np.where(fd_mask_flat==2)[1]

            # add the map to a full disk normalised array
            fd_overlap_map[map_bottom:map_top, map_left:map_right] = map.data

            # overlap for higher ntv and doppler values, and get most extreme magnetogram values
            y_loc = 0
            x_loc = 0
            for j in np.arange(0,len(overlap_height_indices)):
              y_loc = overlap_height_indices[j]
              x_loc = overlap_width_indices[j]
              if fd_overlap_map[y_loc, x_loc] > fd_data[y_loc, x_loc]:
                fd_data[y_loc, x_loc] = fd_overlap_map[y_loc, x_loc]

    # convert values where no raster is present to NaNs so they are not plotted
    fd_data = np.where(fd_mask_flat==0, np.nan, fd_data)

    # convert any pixels where the value is exactly 0 to NaNs so they are not plotted
    fd_data = np.where(fd_data==0, np.nan, fd_data)

    # create the full disk maps
    fd_map = sunpy.map.Map(fd_data, fd_header)

    # if looking at velocity maps, correct for line of sight effects
    if measurement == 'vel':
        solar_x_len = fd_map.data.shape[1]
        solar_y_len = fd_map.data.shape[0]
        # solar_x = np.arange(-fd_size/2, fd_size/2, map_dx)
        # solar_y = np.arange(-fd_size/2, fd_size/2, map_dy)

        # # create an empty numpy array for the pointing coordinates
        # solar_x_array = np.empty([solar_y_len,solar_x_len])
        # solar_y_array = np.empty([solar_y_len,solar_x_len])

        # # fill the solar_x array
        # for i in range(0,solar_y_len):
        #   solar_x_array[i,:] = solar_x

        # # fill the solar_y array
        # for i in range(0,solar_x_len):
        #   solar_y_array[:,i] = solar_y

        # Using np.linspace to ensure the arrays have the exact expected length
        solar_x = np.linspace(-fd_size/2, fd_size/2, solar_x_len)
        solar_y = np.linspace(-fd_size/2, fd_size/2, solar_y_len)

        # create meshgrid for solar_x and solar_y
        solar_x_array, solar_y_array = np.meshgrid(solar_x, solar_y)

        # convert solar_x and solar_y to longitude and latitude
        obs_time = fd_map.meta['date_obs']
        helioprojective_coordinates = SkyCoord(solar_x_array * u.arcsec, solar_y_array * u.arcsec, frame=frames.Helioprojective(observer="earth", obstime=obs_time))
        heliographic_coordinates = helioprojective_coordinates.transform_to(frames.HeliographicStonyhurst)
        latitude = heliographic_coordinates.lat.to(u.deg).value
        longitude = heliographic_coordinates.lon.to(u.deg).value

        # apply the doppler velocity correction
        cos_mu = np.cos(latitude*((2*np.pi)/360)) * np.cos(longitude*((2*np.pi)/360))

        for i in range (0, fd_map.data.shape[0]):
            for j in range (0, fd_map.data.shape[1]):
                fd_map.data[i][j] = fd_map.data[i][j] / cos_mu[i][j]


    # parse the date from first_map.meta["date_beg"] in format YYYY-MM-DDTHH:MM:SS.sss
    fd_map_datetime = datetime.strptime(first_map.meta["date_beg"], "%Y-%m-%dT%H:%M:%S.%f")
    fd_map_datetime = fd_map_datetime.strftime("%Y%m%d_%H%M%S")

    path = './data_eis/full_disks/'
    if not os.path.exists(path):
        os.makedirs(path)

    # save filename should be in format eis_YYYYMMDD_HHMMSS.{emission_line}.{measurement}..fits
    save_filename = f'eis_{fd_map_datetime}.{line_id}.{measurement}.fits'
    fd_map.save(path+save_filename, overwrite=True)

    # plot the full disk non-thermal velocity map
    fig = plt.figure()
    ax = plt.subplot(projection=fd_map)

    if measurement == 'int':
        im = fd_map.plot(cmap='Oranges')
    elif measurement == 'vel':
        im = fd_map.plot(cmap='RdBu_r')
        im.set_norm(plt.Normalize(vmin=-10, vmax=10))
    elif measurement == 'wid':
        im = fd_map.plot(cmap='viridis')
    elif measurement == 'ntv':
        im = fd_map.plot(cmap='inferno')
        im.set_norm(plt.Normalize(vmin=0, vmax=40))
    # fd_map_ntv.draw_limb()
    # fd_map_ntv.draw_grid()
    im = ax.get_images()
    im_lims =  im[0].get_extent()
    ax.set_aspect(abs((im_lims[1]-im_lims[0])/(im_lims[3]-im_lims[2])))
    # plt.title(f'Hinode EIS non-thermal velocity (km/s)\n{fd_map.meta["date_beg"]} - {fd_map.meta["date_end"]}')
    plt.colorbar(extend='both')
    plt.savefig(path+save_filename.replace('.fits', '.png'))
    if high_res:
        plt.savefig(path+save_filename.replace('.fits', '.pdf'), dpi=1000)
    plt.close()

    return True