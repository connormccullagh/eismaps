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
from sunpy.map.header_helper import make_fitswcs_header
import sunkit_image.coalignment
from matplotlib.colors import LogNorm


def check_fd(files):
    # that that the set of files specified all belong to the same full disk observation
    return None

def main(files, high_res=False, base_path='./data_eis/', overlap=None, mode=None, fd_fits=None, los_correct=False):

    if overlap == None:
        raise ValueError('Error: overlap is not specified')

    # if mode is not None but the fd_fits is not specified or does not match the length of files, raise an error
    if mode != None:
        if fd_fits == None:
            raise ValueError('Error: mode is not None but fd_fits is not specified')
        if len(files) != len(fd_fits):
            raise ValueError('Error: mode is not None but fd_fits does not match the length of files')

    if mode == None:
        measurement = files[0].split('.')[-2] # can be 'int', 'vel', 'wid', 'ntv' for now
        line_id = files[0].split('.')[-3]
    else:
        measurement = None
        line_id = None    

    # get the maximum extent of the full disk   
    fd_size = 0
    for file in files:
        map = sunpy.map.Map(file)
        if((abs(map.meta['xcen'])+(map.meta['fovx']/2))*2 > fd_size):
            fd_size = (abs(map.meta['xcen'])+(map.meta['fovx']/2))*2
        if((abs(map.meta['ycen'])+(map.meta['fovy']/2))*2 > fd_size):
            fd_size = (abs(map.meta['ycen'])+(map.meta['fovy']/2))*2

    first_map = sunpy.map.Map(files[0])

    # read in the pixel sizes for this full disk image
    map_dx = first_map.meta['cdelt1']
    map_dy = first_map.meta['cdelt2']

    # set the number of columns (width) and rows (height) of the full disk data array
    fd_width = round(fd_size/map_dx)
    fd_height = round(fd_size/map_dy)

    # create the blank arrays for the full disk image
    fd_data = np.full((fd_height, fd_width), 0.0)
    fd_mask = np.full((fd_height, fd_width), 0)
    fd_mask_flat = np.full((fd_height, fd_width), 0)
    fd_OVERLAP = np.full((fd_height, fd_width), np.nan)

    # create array of fits file names and the associated start times
    map_times = [[files[0], first_map.meta['date_beg'], first_map.meta['date_end']]]
    for file in files[1:]:
      map = sunpy.map.Map(file)
      map_times.append([file, map.meta['date_beg'], map.meta['date_end']])
    map_times = np.array(map_times)

    # create the header file for the full disk map
    fd_header = {
      'bunit': '', # units of array values
      'cdelt1': map_dx, # coordinate increment along axis
      'cdelt2': map_dy,
      'cname1': 'Solar-X', # axis labes
      'cname2': 'Solar-Y',
      'crpix1': 1, # reference pixel (bottom left)
      'crpix2': 1,
      'crval1': -fd_size/2, # coordinate of reference pixel
      'crval2': -fd_size/2,
      'ctype1': 'HPLN-TAN', # !!! fits standard says this should be label of coordiate axis... experiment ###HELP!! WHAT IS THIS???
      'ctype2': 'HPLT-TAN',
      'cunit1': 'arcsec', # unit of axis
      'cunit2': 'arcsec',
      'date_beg': map_times[0,1], # start time of observation
      'date_end': map_times[-1,2], # end time of observation
      'date_obs': map_times[0,1], # start time of observation
      'dsun_obs': first_map.meta['dsun_obs'], # distance to sun centre from spacecraft
      'hgln_obs': first_map.meta['hgln_obs'], # Stonyhurst heliographic longitude of the observer
      'hglt_obs': first_map.meta['hglt_obs'], # Heliographic latitude (Carrington or Stonyhurst) of the observers
      'instrume': '',
      'line_id': line_id,
      'measrmnt': measurement,
      'telescop': '',
      'timesys': 'UTC',
      'xcen': 0.0,
      'ycen': 0.0,
      'history': 'blank fits created'
    }

    # create a blank map with this header for use of parameters below
    fd_map = sunpy.map.Map(fd_data, fd_header)

    # insert raster data into the full disk data array
    iraster_total = len(files)
    for i, file in enumerate(files):
        
        map = sunpy.map.Map(file)

        # if the mode is anything but None, use the data from the HMI fits files to create the full disk maps
        if mode != None:
            map_extra = sunpy.map.Map(fd_fits[i])
            map_data = map_extra.data
        else:
            map_data = map.data

        # get the initial position of the raster
        map_x_init = map.meta['xcen']
        map_y_init = map.meta['ycen']
        map_rad = np.sqrt(map_x_init**2 + map_y_init**2)

        # if the centre of the raster is outside the limb, take the point at the same angle to the horizontal as the centre of the raster, but within the limb so differential rotation can be applied
        if(map_rad > 940): ## TODO impliment accurate limb distance
          map_angle = math.atan(map_y_init/map_x_init)
          map_y_init_shifted = 950 * math.sin(map_angle) # calculate distance from centre of sun to limb at this angle
          map_x_init_shifted = 950 * math.cos(map_angle)
          map_y_SHIFT = map_y_init - map_y_init_shifted # calculate how far the map needs to be shifted to be centred on the limb
          map_x_SHIFT = map_x_init - map_x_init_shifted
          map_y_init = map_y_init_shifted
          map_x_init = map_x_init_shifted
        else:
          map_y_SHIFT = 0
          map_x_SHIFT = 0

        # apply differential rotation to the map using a selected point
        point = SkyCoord(map_x_init*u.arcsec, map_y_init*u.arcsec, frame=map.coordinate_frame) # use the centre of the map as the point to apply differential rotation to (may have been shifted to ensure it is inside the limb)
        fd_map_beg_time = datetime.strptime(fd_map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
        map_beg_time = datetime.strptime(map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
        diffrot_point = RotatedSunFrame(base=point, duration=((fd_map_beg_time-map_beg_time).seconds)*u.second) # set up the rotation
        transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame) # transform the raster centre

        new_x = str(transformed_diffrot_point.Tx)
        new_y = str(transformed_diffrot_point.Ty)
        new_x = new_x.replace('arcsec', '')
        new_y = new_y.replace('arcsec', '')
        new_x = float(new_x)
        new_y = float(new_y)

        # put the map back in the correct position if its centre was outside the limb
        map.meta['xcen'] = new_x + map_x_SHIFT
        map.meta['ycen'] = new_y + map_y_SHIFT

        print(f'this map shifted {new_y-map_y_init} in y and {new_x-map_x_init} in x for time')

        # calculate the boundaries of the map in arcsecs and convert to array indexes
        map_height, map_width = map.data.shape
        map_left = round( (map.meta['xcen']/map_dx - map_width/2) + fd_width/2 )
        map_right = round( map_left + map_width )
        map_bottom = round( (map.meta['ycen']/map_dy - map_height/2) + fd_height/2 )
        map_top = round( map_bottom + map_height )

        # update the mask to show the number of contributions being made to the particular cells
        fd_mask_flat[map_bottom:map_top, map_left:map_right] = fd_mask_flat[map_bottom:map_top, map_left:map_right] + 1

        # insert the data to the full disk map where there are no overlaps        
        fd_data[map_bottom:map_top, map_left:map_right] = np.where(fd_mask_flat[map_bottom:map_top, map_left:map_right]<=1,fd_data[map_bottom:map_top, map_left:map_right] + map_data, fd_data[map_bottom:map_top, map_left:map_right])

        # detect if merging is required
        if np.any(fd_mask_flat[map_bottom:map_top, map_left:map_right] == 2):
          print(f'merging')

          # overlaps are here:
          overlap_height_indices = np.where(fd_mask_flat==2)[0]
          overlap_widths_indices = np.where(fd_mask_flat==2)[1]

          # add the map to a full disk normalised array
          fd_OVERLAP[map_bottom:map_top, map_left:map_right] = map_data

          # treat overlapping pixels
          y_loc = 0
          x_loc = 0
          for j in np.arange(0,len(overlap_height_indices)):
            y_loc = overlap_height_indices[j]
            x_loc = overlap_widths_indices[j]
            if overlap == 'max':
                if abs(fd_OVERLAP[y_loc, x_loc]) > abs(fd_data[y_loc, x_loc]):
                    fd_data[y_loc, x_loc] = fd_OVERLAP[y_loc, x_loc]
            elif overlap == 'mean':
                fd_data[y_loc, x_loc] = (fd_data[y_loc, x_loc] + fd_OVERLAP[y_loc, x_loc])/2
            elif overlap == 'min':
                if abs(fd_OVERLAP[y_loc, x_loc]) < abs(fd_data[y_loc, x_loc]):
                    fd_data[y_loc, x_loc] = fd_OVERLAP[y_loc, x_loc]
            elif overlap == 'sum':
                fd_data[y_loc, x_loc] = fd_data[y_loc, x_loc] + fd_OVERLAP[y_loc, x_loc]
            else:
                raise ValueError('Error: overlap is not valid')

          # flatten the flat mask
          fd_mask_flat = np.where((fd_mask_flat >= 1), 1, 0)

    # convert values where no raster is present to NaNs so they are not plotted
    fd_data = np.where(fd_mask_flat==0, np.nan, fd_data)

    # create the full disk maps
    fd_map = sunpy.map.Map(fd_data, fd_header)

    # # if correcting for LOS effects
    # if los_correct:
    #     solar_x_len = fd_map.data.shape[1]
    #     solar_y_len = fd_map.data.shape[0]

    #     # Using np.linspace to ensure the arrays have the exact expected length
    #     solar_x = np.linspace(-fd_size/2, fd_size/2, solar_x_len)
    #     solar_y = np.linspace(-fd_size/2, fd_size/2, solar_y_len)

    #     # create meshgrid for solar_x and solar_y
    #     solar_x_array, solar_y_array = np.meshgrid(solar_x, solar_y)

    #     # convert solar_x and solar_y to longitude and latitude
    #     obs_time = fd_map.meta['date_obs']
    #     helioprojective_coordinates = SkyCoord(solar_x_array * u.arcsec, solar_y_array * u.arcsec, frame=frames.Helioprojective(observer="earth", obstime=obs_time))
    #     heliographic_coordinates = helioprojective_coordinates.transform_to(frames.HeliographicStonyhurst)
    #     latitude = heliographic_coordinates.lat.to(u.deg).value
    #     longitude = heliographic_coordinates.lon.to(u.deg).value

    #     # apply the doppler velocity correction
    #     cos_mu = np.cos(latitude*((2*np.pi)/360)) * np.cos(longitude*((2*np.pi)/360))

    #     for i in range (0, fd_map.data.shape[0]):
    #         for j in range (0, fd_map.data.shape[1]):
    #             fd_map.data[i][j] = fd_map.data[i][j] / cos_mu[i][j]


    # parse the date from first_map.meta["date_beg"] in format YYYY-MM-DDTHH:MM:SS.sss
    fd_map_datetime = datetime.strptime(first_map.meta["date_beg"], "%Y-%m-%dT%H:%M:%S.%f")
    fd_map_datetime = fd_map_datetime.strftime("%Y%m%d_%H%M%S")
    
    path = base_path + 'full_disks/'
    if not os.path.exists(path):
        os.makedirs(path)

    # save filename should be in format eis_YYYYMMDD_HHMMSS.{emission_line}.{measurement}..fits
    save_filename = f'eis_{fd_map_datetime}.{line_id}.{measurement}.fits'
    fd_map.save(path+save_filename, overwrite=True)

    # plot the full disk non-thermal velocity map
    fig = plt.figure()
    ax = plt.subplot(projection=fd_map)

    if measurement == 'int':
        im = fd_map.plot(cmap='gist_heat')
        im.set_norm(LogNorm(vmin=1e1, vmax=5e3))
    elif measurement == 'vel':
        im = fd_map.plot(cmap='RdBu_r')
        im.set_norm(plt.Normalize(vmin=-20, vmax=20))
    elif measurement == 'wid':
        im = fd_map.plot(cmap='viridis')
    elif measurement == 'ntv':
        im = fd_map.plot(cmap='inferno')
        im.set_norm(plt.Normalize(vmin=0, vmax=40))
    else:
        if mode == 'hmi':
            im = fd_map.plot(cmap='gray')
            im.set_norm(plt.Normalize(vmin=-300, vmax=300))
        elif mode == 'fip':
            im = fd_map.plot(cmap='CMRmap')
            im.set_norm(plt.Normalize(vmin=0, vmax=3))
        else:
            im = fd_map.plot(cmap='Oranges')
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














    
#     map_dx = first_map.meta['cdelt1']
#     map_dy = first_map.meta['cdelt2']

#     measurement = files[0].split('.')[-2] # can be 'int', 'vel', 'wid', 'ntv' for now
#     line_id = files[0].split('.')[-3]

#     # double check all the maps have the same dx and dy #TODO add checks of other parameters (e.g. the measurement, line)
#     for file in files:
#         map = sunpy.map.Map(file)
#         if (map.meta['cdelt1'] != map_dx) or (map.meta['cdelt2'] != map_dy):
#             raise ValueError('Error: maps do not have the same dx and dy')

#     fd_width = round(fd_size/map_dx)
#     fd_height = round(fd_size/map_dy)

#     fd_data = np.full((fd_height, fd_width), 0.0)
#     fd_mask = np.full((fd_height, fd_width), 0)
#     fd_mask_flat = np.full((fd_height, fd_width), 0)
#     fd_overlap_map = np.full((fd_height, fd_width), np.nan)

#     # create a list of the times of the maps
#     map_times = []
#     for file in files:
#       map = sunpy.map.Map(file)
#       map_times.append([file, map.meta['date_beg'], map.meta['date_end']])
#     map_times = np.array(map_times)

#     fd_header = {
#         'crpix1': 1, # reference pixel (bottom left)
#         'crpix2': 1,
#         'cdelt1': map_dx,
#         'cdelt2': map_dy,
#         'cunit1': 'arcsec',
#         'cunit2': 'arcsec',
#         'crval1': -fd_size/2, # location of reference coordinate
#         'crval2': -fd_size/2,
#         'ctype1': 'HPLN-TAN',
#         'ctype2': 'HPLT-TAN',
#         'date_obs': first_map.meta['date_beg'],
#         'date_beg': first_map.meta['date_beg'],
#         'date_end': last_map.meta['date_end'],
#         'dsun_obs': first_map.meta['dsun_obs'],
#         'hgln_obs': first_map.meta['hgln_obs'],
#         'hglt_obs': first_map.meta['hglt_obs'],
#         'instrume': 'EIS',
#         'measrmnt': 'int', # TEMP for changing
#         'telescop': 'Hinode',
#         'xcen': 0,
#         'ycen': 0,
#         'history': 'Created by eismaps.proc.full_disk.main'
#     }

#     fd_map = sunpy.map.Map(fd_data, fd_header)

#     # insert the maps into the full disk map
#     imap = 0
#     for file in files:

#         map = sunpy.map.Map(file)
#         print(file)

#         # find the pixel nearest the centre of the map but which is still on the limb, and apply the differential rotation
#         # use the sunpy functions to do this instead, to make it a more robust approach.
#         ###### TODO ######

#         pixel_coords = sunpy.map.all_coordinates_from_map(map) # make array of all pixel coordinates
#         is_on_disk = sunpy.map.coordinate_is_on_solar_disk(pixel_coords) # make array of booleans for whether each pixel is on the disk

#         # print the map x and y cen
#         print(map.center.Tx.value, map.center.Ty.value)

#         map_xcen = pixel_coords[pixel_coords.shape[0]//2, pixel_coords.shape[1]//2].Tx.value
#         map_ycen = pixel_coords[pixel_coords.shape[0]//2, pixel_coords.shape[1]//2].Ty.value

#         print(map_xcen, map_ycen)
#         exit()

#         if not np.any(is_on_disk): # if the map is completely off the disk...
        
#             # don't apply shifting
#             print('off disk')

#         else:
#             # if the centre of the map is not on the disk, find the nearest pixel that is on the disk
#             if not is_on_disk[int(map.data.shape[0]/2), int(map.data.shape[1]/2)]: # if the centre of the map is not on the disk
#                 # find the nearest pixel that is on the disk
#                 distances = np.zeros((map.data.shape[0], map.data.shape[1]))
#                 for i in range(0, map.data.shape[0]):
#                     for j in range(0, map.data.shape[1]):
#                         if not is_on_disk[i,j]:
#                             distances[i,j] = np.inf
#                         else:
#                             distances[i,j] = np.sqrt((i - int(map.data.shape[0]/2))**2 + (j - int(map.data.shape[1]/2))**2)
#                 nearest_disk_pixel = np.unravel_index(np.argmin(distances, axis=None), distances.shape)
#                 chief_pixel = nearest_disk_pixel

#             else: 
#                 chief_pixel = (int(map.data.shape[0]/2), int(map.data.shape[1]/2))
            
#         # # # # apply the differential rotation
#         # # # point = pixel_coords[chief_pixel[0], chief_pixel[1]]
#         # # # fd_map_beg_time = datetime.strptime(fd_map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
#         # # # map_beg_time = datetime.strptime(map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
#         # # # diffrot_point = RotatedSunFrame(base=point, duration=((fd_map_beg_time-map_beg_time).seconds)*u.second)
#         # # # transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)
#         # # # map = map.shift_reference_coord(transformed_diffrot_point.Tx, transformed_diffrot_point.Ty)

#             point = pixel_coords[chief_pixel[0], chief_pixel[1]]
#             fd_map_beg_time = datetime.strptime(fd_map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
#             map_beg_time = datetime.strptime(map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
#             print(fd_map_beg_time, map_beg_time)
#             print(((fd_map_beg_time-map_beg_time).seconds)*u.second)
#             diffrot_point = SkyCoord(RotatedSunFrame(base=point, duration=((fd_map_beg_time-map_beg_time).seconds)*u.second))
#             transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)
#             shift_to_apply_x = transformed_diffrot_point.Tx - point.Tx
#             shift_to_apply_y = transformed_diffrot_point.Ty - point.Ty
#             print(shift_to_apply_x, shift_to_apply_y,'****************')
#             map = map.shift_reference_coord(shift_to_apply_y, shift_to_apply_x)



#         # transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)
#         # print(f'this map shifted {transformed_diffrot_point.Tx} in x and {transformed_diffrot_point.Ty} in y for time')
#         # map = map.shift_reference_coord(transformed_diffrot_point.Tx, transformed_diffrot_point.Ty)
#         # exit(0)


        
#         # map_x_init = map.center.Tx.value
#         # map_y_init = map.center.Ty.value
#         # map_rad = np.sqrt(map_x_init**2 + map_y_init**2)

#         # # if the centre is outside the limb, the sunpy differential rotation can't be applied. We move these maps slightly inside the limb, apply the rotation, then put them back
#         # if map_rad > 950:
#         #     map_angle = math.atan(map_y_init/map_x_init)
#         #     map_y_init_shifted = 950 * math.sin(map_angle)
#         #     map_x_init_shifted = 950 * math.cos(map_angle)
#         #     map_y_shift = map_y_init - map_y_init_shifted
#         #     map_x_shift = map_x_init - map_x_init_shifted
#         #     map_y_init = map_y_init_shifted
#         #     map_x_init = map_x_init_shifted
#         # else:
#         #     map_y_shift = 0
#         #     map_x_shift = 0

#         # # apply the differential rotation
#         # point = SkyCoord(map_x_init*u.arcsec, map_y_init*u.arcsec, frame=map.coordinate_frame)
#         # fd_map_beg_time = datetime.strptime(fd_map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
#         # map_beg_time = datetime.strptime(map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
#         # diffrot_point = RotatedSunFrame(base=point, duration=((fd_map_beg_time-map_beg_time).seconds)*u.second)
#         # transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)

#         # new_x = str(transformed_diffrot_point.Tx)
#         # new_y = str(transformed_diffrot_point.Ty)
#         # new_x = new_x.replace('arcsec', '')
#         # new_y = new_y.replace('arcsec', '')
#         # new_x = float(new_x)
#         # new_y = float(new_y)

#         # # put the map back
#         # # map.meta['xcen'] = new_x + map_x_shift
#         # # map.meta['ycen'] = new_y + map_y_shift
#         # # convert new_x and new_y to astropy Quantity objects with units of arcseconds

#         # shift_to_apply_x = (new_x-map_x_init) * u.arcsec
#         # shift_to_apply_y = (new_y-map_y_init) * u.arcsec

#         # map = map.shift_reference_coord(shift_to_apply_y, shift_to_apply_x)

#         # print(f'this map shifted {shift_to_apply_y} in y and {shift_to_apply_x} in x for time')

#         pixel_coords = sunpy.map.all_coordinates_from_map(map) # make array of all pixel coordinates
#         map_xcen = pixel_coords[pixel_coords.shape[0]//2, pixel_coords.shape[1]//2].Tx.value
#         map_ycen = pixel_coords[pixel_coords.shape[0]//2, pixel_coords.shape[1]//2].Ty.value

#         print(map_xcen, map_ycen)
        
#         # calculate the boundaries of the map in arcsecs and convert to array indices
#         map_height, map_width = map.data.shape
#         map_left = round( (map_xcen/map_dx - map_width/2) + fd_width/2 )
#         map_right = round( map_left + map_width )
#         map_bottom = round( (map_ycen/map_dy - map_height/2) + fd_height/2 )
#         map_top = round( map_bottom + map_height )

#         print(map_height, map_width)
#         print(map_left, map_right, map_bottom, map_top)

#         continue

#         # update the mask to show the number of contributions being made to the particular cells
#         fd_mask_flat[map_bottom:map_top, map_left:map_right] = fd_mask_flat[map_bottom:map_top, map_left:map_right] + 1

#         print(np.where(fd_mask_flat[map_bottom:map_top, map_left:map_right]<=1))

#         # insert the data to the full disk map where there are no overlaps
#         fd_data[map_bottom:map_top, map_left:map_right] = np.where(fd_mask_flat[map_bottom:map_top, map_left:map_right]<=1,fd_data[map_bottom:map_top, map_left:map_right] + map.data, fd_data[map_bottom:map_top, map_left:map_right])

#         # detect if merging is required
#         if np.any(fd_mask_flat[map_bottom:map_top, map_left:map_right] == 2):
#             print(f'merging')

#             # overlaps are here
#             overlap_height_indices = np.where(fd_mask_flat==2)[0]
#             overlap_width_indices = np.where(fd_mask_flat==2)[1]

#             # add the map to a full disk normalised array
#             fd_overlap_map[map_bottom:map_top, map_left:map_right] = map.data

#             # overlap for higher ntv and doppler values, and get most extreme magnetogram values
#             y_loc = 0
#             x_loc = 0
#             for j in np.arange(0,len(overlap_height_indices)):
#               y_loc = overlap_height_indices[j]
#               x_loc = overlap_width_indices[j]
#               if fd_overlap_map[y_loc, x_loc] > fd_data[y_loc, x_loc]:
#                 fd_data[y_loc, x_loc] = fd_overlap_map[y_loc, x_loc]

#     # convert values where no raster is present to NaNs so they are not plotted
#     fd_data = np.where(fd_mask_flat==0, np.nan, fd_data)

#     # convert any pixels where the value is exactly 0 to NaNs so they are not plotted
#     fd_data = np.where(fd_data==0, np.nan, fd_data)

#     # create the full disk maps
#     fd_map = sunpy.map.Map(fd_data, fd_header)

#     # if looking at velocity maps, correct for line of sight effects
#     if measurement == 'vel':
#         solar_x_len = fd_map.data.shape[1]
#         solar_y_len = fd_map.data.shape[0]
#         # solar_x = np.arange(-fd_size/2, fd_size/2, map_dx)
#         # solar_y = np.arange(-fd_size/2, fd_size/2, map_dy)

#         # # create an empty numpy array for the pointing coordinates
#         # solar_x_array = np.empty([solar_y_len,solar_x_len])
#         # solar_y_array = np.empty([solar_y_len,solar_x_len])

#         # # fill the solar_x array
#         # for i in range(0,solar_y_len):
#         #   solar_x_array[i,:] = solar_x

#         # # fill the solar_y array
#         # for i in range(0,solar_x_len):
#         #   solar_y_array[:,i] = solar_y

#         # Using np.linspace to ensure the arrays have the exact expected length
#         solar_x = np.linspace(-fd_size/2, fd_size/2, solar_x_len)
#         solar_y = np.linspace(-fd_size/2, fd_size/2, solar_y_len)

#         # create meshgrid for solar_x and solar_y
#         solar_x_array, solar_y_array = np.meshgrid(solar_x, solar_y)

#         # convert solar_x and solar_y to longitude and latitude
#         obs_time = fd_map.meta['date_obs']
#         helioprojective_coordinates = SkyCoord(solar_x_array * u.arcsec, solar_y_array * u.arcsec, frame=frames.Helioprojective(observer="earth", obstime=obs_time))
#         heliographic_coordinates = helioprojective_coordinates.transform_to(frames.HeliographicStonyhurst)
#         latitude = heliographic_coordinates.lat.to(u.deg).value
#         longitude = heliographic_coordinates.lon.to(u.deg).value

#         # apply the doppler velocity correction
#         cos_mu = np.cos(latitude*((2*np.pi)/360)) * np.cos(longitude*((2*np.pi)/360))

#         for i in range (0, fd_map.data.shape[0]):
#             for j in range (0, fd_map.data.shape[1]):
#                 fd_map.data[i][j] = fd_map.data[i][j] / cos_mu[i][j]


#     # parse the date from first_map.meta["date_beg"] in format YYYY-MM-DDTHH:MM:SS.sss
#     fd_map_datetime = datetime.strptime(first_map.meta["date_beg"], "%Y-%m-%dT%H:%M:%S.%f")
#     fd_map_datetime = fd_map_datetime.strftime("%Y%m%d_%H%M%S")
    
#     # path = './data_eis/full_disks/'
#     path = base_path + 'full_disks/'
#     if not os.path.exists(path):
#         os.makedirs(path)

#     # save filename should be in format eis_YYYYMMDD_HHMMSS.{emission_line}.{measurement}..fits
#     save_filename = f'eis_{fd_map_datetime}.{line_id}.{measurement}.fits'
#     fd_map.save(path+save_filename, overwrite=True)

#     # plot the full disk non-thermal velocity map
#     fig = plt.figure()
#     ax = plt.subplot(projection=fd_map)

#     if measurement == 'int':
#         im = fd_map.plot(cmap='Oranges')
#     elif measurement == 'vel':
#         im = fd_map.plot(cmap='RdBu_r')
#         im.set_norm(plt.Normalize(vmin=-10, vmax=10))
#     elif measurement == 'wid':
#         im = fd_map.plot(cmap='viridis')
#     elif measurement == 'ntv':
#         im = fd_map.plot(cmap='inferno')
#         im.set_norm(plt.Normalize(vmin=0, vmax=40))
#     # fd_map_ntv.draw_limb()
#     # fd_map_ntv.draw_grid()
#     im = ax.get_images()
#     im_lims =  im[0].get_extent()
#     ax.set_aspect(abs((im_lims[1]-im_lims[0])/(im_lims[3]-im_lims[2])))
#     # plt.title(f'Hinode EIS non-thermal velocity (km/s)\n{fd_map.meta["date_beg"]} - {fd_map.meta["date_end"]}')
#     plt.colorbar(extend='both')
#     plt.savefig(path+save_filename.replace('.fits', '.png'))
#     if high_res:
#         plt.savefig(path+save_filename.replace('.fits', '.pdf'), dpi=1000)
#     plt.close()

#     return True


# # def main(files, high_res=False, base_path='./data_eis/', type='eis'):
    
# #     fd_size = 0
# #     for file in files:
# #         # get the maximum extent of the full disk
# #         map = sunpy.map.Map(file)
# #         map_xcen = map.center.Tx.value
# #         map_ycen = map.center.Ty.value
# #         map_bl = map.bottom_left_coord
# #         map_tr = map.top_right_coord
# #         map_fovx = abs(map_tr.Tx.value - map_bl.Tx.value)
# #         map_fovy = abs(map_tr.Ty.value - map_bl.Ty.value)

# #         if ((abs(map_xcen)+(map_fovx/2))*2 > fd_size):
# #             fd_size = (abs(map_xcen)+(map_fovx/2))*2
# #         if ((abs(map_ycen)+(map_fovy/2))*2 > fd_size):
# #             fd_size = (abs(map_ycen)+(map_fovy/2))*2

# #     map = sunpy.map.Map(files[0])

# #     # get the spatial scale of the files rounded to 0 decimal places
# #     map_dx = map.scale[0]
# #     map_dy = map.scale[1]
# #     map_dx_val = round(map_dx.value, 0)
# #     map_dy_val = round(map_dy.value, 0)
# #     map_dx = map_dx_val * map_dx.unit
# #     map_dy = map_dy_val * map_dy.unit
# #     # set map_dx to the value in arcse/pix
# #     map_dx = map_dx.to(u.arcsec/u.pixel).value
# #     map_dy = map_dy.to(u.arcsec/u.pixel).value

# #     ### could add a check that all the maps have the same dx and dy here

# #     # setup the full disk
# #     fd_width = round(fd_size/map_dx)
# #     fd_height = round(fd_size/map_dy)

# #     fd_data = np.full((fd_height, fd_width), 0.0)
# #     fd_mask = np.full((fd_height, fd_width), 0)
# #     fd_mask_flat = np.full((fd_height, fd_width), 0)
# #     fd_overlap_map = np.full((fd_height, fd_width), np.nan)

# #     # create a list of the times of the maps
# #     map_times = []
# #     for file in files:
# #         map = sunpy.map.Map(file)
# #         if 'date_beg' not in map.meta: # eis uses date_beg/end, hmi uses date-obs
# #             map_times.append([file, map.meta['date-obs'], map.meta['date-obs']])
# #         else:
# #             map_times.append([file, map.meta['date_beg'], map.meta['date_end']])
# #     map_times = np.array(map_times)
# #     print(map_times)

# #     #     fd_header = {
# # #         'crpix1': 1, # reference pixel (bottom left)
# # #         'crpix2': 1,
# # #         'cdelt1': map_dx,
# # #         'cdelt2': map_dy,
# # #         'cunit1': 'arcsec',
# # #         'cunit2': 'arcsec',
# # #         'crval1': -fd_size/2, # location of reference coordinate
# # #         'crval2': -fd_size/2,
# # #         'ctype1': 'HPLN-TAN',
# # #         'ctype2': 'HPLT-TAN',
# # #         'date_obs': first_map.meta['date_beg'],
# # #         'date_beg': first_map.meta['date_beg'],
# # #         'date_end': last_map.meta['date_end'],
# # #         'dsun_obs': first_map.meta['dsun_obs'],
# # #         'hgln_obs': first_map.meta['hgln_obs'],
# # #         'hglt_obs': first_map.meta['hglt_obs'],
# # #         'instrume': 'EIS',
# # #         'measrmnt': 'int', # TEMP for changing
# # #         'telescop': 'Hinode',
# # #         'xcen': 0,
# # #         'ycen': 0,
# # #         'history': 'Created by eismaps.proc.full_disk.main'
# # #     }

# #     if 'date_beg' not in map.meta:
# #       obstime = map.meta['date-obs']
# #     else:
# #       obstime = map.meta['date_beg']
# #     reference_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime = obstime, observer='earth', frame=frames.Helioprojective)
# #     reference_pixel = u.Quantity([int(fd_width/2), int(fd_height/2)], u.pixel)
# #     scale = u.Quantity([map_dx, map_dy], u.arcsec/u.pixel)
# #     fd_header = make_fitswcs_header(fd_data, reference_coord, reference_pixel=reference_pixel, scale=scale)

# #     fd_map = sunpy.map.Map(fd_data, fd_header)

# #     # insert the maps into the full disk map
# #     imap = 0
# #     for file in files:

# #         map = sunpy.map.Map(file)
# #         print(file)


# #         pixel_coords = sunpy.map.all_coordinates_from_map(map) # make array of all pixel coordinates
# #         is_on_disk = sunpy.map.coordinate_is_on_solar_disk(pixel_coords) # make array of booleans for whether each pixel is on the disk

# #         map_xcen = pixel_coords[pixel_coords.shape[0]//2, pixel_coords.shape[1]//2].Tx.value
# #         map_ycen = pixel_coords[pixel_coords.shape[0]//2, pixel_coords.shape[1]//2].Ty.value

# #         if not np.any(is_on_disk): # if the map is completely off the disk...
        
# #             # don't apply shifting
# #             print('off disk')

# #         else:
# #             # if the centre of the map is not on the disk, find the nearest pixel that is on the disk
# #             if not is_on_disk[int(map.data.shape[0]/2), int(map.data.shape[1]/2)]: # if the centre of the map is not on the disk
# #                 # find the nearest pixel that is on the disk
# #                 distances = np.zeros((map.data.shape[0], map.data.shape[1]))
# #                 for i in range(0, map.data.shape[0]):
# #                     for j in range(0, map.data.shape[1]):
# #                         if not is_on_disk[i,j]:
# #                             distances[i,j] = np.inf
# #                         else:
# #                             distances[i,j] = np.sqrt((i - int(map.data.shape[0]/2))**2 + (j - int(map.data.shape[1]/2))**2)
# #                 nearest_disk_pixel = np.unravel_index(np.argmin(distances, axis=None), distances.shape)
# #                 chief_pixel = nearest_disk_pixel

# #             else: 
# #                 chief_pixel = (int(map.data.shape[0]/2), int(map.data.shape[1]/2))
            
# # #         # # # # apply the differential rotation
# # #         # # # point = pixel_coords[chief_pixel[0], chief_pixel[1]]
# # #         # # # fd_map_beg_time = datetime.strptime(fd_map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
# # #         # # # map_beg_time = datetime.strptime(map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
# # #         # # # diffrot_point = RotatedSunFrame(base=point, duration=((fd_map_beg_time-map_beg_time).seconds)*u.second)
# # #         # # # transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)
# # #         # # # map = map.shift_reference_coord(transformed_diffrot_point.Tx, transformed_diffrot_point.Ty)

# #             point = pixel_coords[chief_pixel[0], chief_pixel[1]]
# #             fd_map_beg_time = datetime.strptime(fd_map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
# #             map_beg_time = datetime.strptime(map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
# #             print(fd_map_beg_time, map_beg_time)
# #             print(((fd_map_beg_time-map_beg_time).seconds)*u.second)
# #             diffrot_point = SkyCoord(RotatedSunFrame(base=point, duration=((fd_map_beg_time-map_beg_time).seconds)*u.second))
# #             transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)
# #             shift_to_apply_x = transformed_diffrot_point.Tx - point.Tx
# #             shift_to_apply_y = transformed_diffrot_point.Ty - point.Ty
# #             print(shift_to_apply_x, shift_to_apply_y,'****************')
# #             map = map.shift_reference_coord(shift_to_apply_y, shift_to_apply_x)



# # #         # transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)
# # #         # print(f'this map shifted {transformed_diffrot_point.Tx} in x and {transformed_diffrot_point.Ty} in y for time')
# # #         # map = map.shift_reference_coord(transformed_diffrot_point.Tx, transformed_diffrot_point.Ty)
# # #         # exit(0)


        
# # #         # map_x_init = map.center.Tx.value
# # #         # map_y_init = map.center.Ty.value
# # #         # map_rad = np.sqrt(map_x_init**2 + map_y_init**2)

# # #         # # if the centre is outside the limb, the sunpy differential rotation can't be applied. We move these maps slightly inside the limb, apply the rotation, then put them back
# # #         # if map_rad > 950:
# # #         #     map_angle = math.atan(map_y_init/map_x_init)
# # #         #     map_y_init_shifted = 950 * math.sin(map_angle)
# # #         #     map_x_init_shifted = 950 * math.cos(map_angle)
# # #         #     map_y_shift = map_y_init - map_y_init_shifted
# # #         #     map_x_shift = map_x_init - map_x_init_shifted
# # #         #     map_y_init = map_y_init_shifted
# # #         #     map_x_init = map_x_init_shifted
# # #         # else:
# # #         #     map_y_shift = 0
# # #         #     map_x_shift = 0

# # #         # # apply the differential rotation
# # #         # point = SkyCoord(map_x_init*u.arcsec, map_y_init*u.arcsec, frame=map.coordinate_frame)
# # #         # fd_map_beg_time = datetime.strptime(fd_map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
# # #         # map_beg_time = datetime.strptime(map.meta['date_beg'], "%Y-%m-%dT%H:%M:%S.%f")
# # #         # diffrot_point = RotatedSunFrame(base=point, duration=((fd_map_beg_time-map_beg_time).seconds)*u.second)
# # #         # transformed_diffrot_point = diffrot_point.transform_to(map.coordinate_frame)

# # #         # new_x = str(transformed_diffrot_point.Tx)
# # #         # new_y = str(transformed_diffrot_point.Ty)
# # #         # new_x = new_x.replace('arcsec', '')
# # #         # new_y = new_y.replace('arcsec', '')
# # #         # new_x = float(new_x)
# # #         # new_y = float(new_y)

# # #         # # put the map back
# # #         # # map.meta['xcen'] = new_x + map_x_shift
# # #         # # map.meta['ycen'] = new_y + map_y_shift
# # #         # # convert new_x and new_y to astropy Quantity objects with units of arcseconds

# # #         # shift_to_apply_x = (new_x-map_x_init) * u.arcsec
# # #         # shift_to_apply_y = (new_y-map_y_init) * u.arcsec

# # #         # map = map.shift_reference_coord(shift_to_apply_y, shift_to_apply_x)

# # #         # print(f'this map shifted {shift_to_apply_y} in y and {shift_to_apply_x} in x for time')

# # #         pixel_coords = sunpy.map.all_coordinates_from_map(map) # make array of all pixel coordinates
# # #         map_xcen = pixel_coords[pixel_coords.shape[0]//2, pixel_coords.shape[1]//2].Tx.value
# # #         map_ycen = pixel_coords[pixel_coords.shape[0]//2, pixel_coords.shape[1]//2].Ty.value

# # #         print(map_xcen, map_ycen)
        
# # #         # calculate the boundaries of the map in arcsecs and convert to array indices
# # #         map_height, map_width = map.data.shape
# # #         map_left = round( (map_xcen/map_dx - map_width/2) + fd_width/2 )
# # #         map_right = round( map_left + map_width )
# # #         map_bottom = round( (map_ycen/map_dy - map_height/2) + fd_height/2 )
# # #         map_top = round( map_bottom + map_height )

# # #         print(map_height, map_width)
# # #         print(map_left, map_right, map_bottom, map_top)

# # #         continue

# # #         # update the mask to show the number of contributions being made to the particular cells
# # #         fd_mask_flat[map_bottom:map_top, map_left:map_right] = fd_mask_flat[map_bottom:map_top, map_left:map_right] + 1

# # #         print(np.where(fd_mask_flat[map_bottom:map_top, map_left:map_right]<=1))

# # #         # insert the data to the full disk map where there are no overlaps
# # #         fd_data[map_bottom:map_top, map_left:map_right] = np.where(fd_mask_flat[map_bottom:map_top, map_left:map_right]<=1,fd_data[map_bottom:map_top, map_left:map_right] + map.data, fd_data[map_bottom:map_top, map_left:map_right])

# # #         # detect if merging is required
# # #         if np.any(fd_mask_flat[map_bottom:map_top, map_left:map_right] == 2):
# # #             print(f'merging')

# # #             # overlaps are here
# # #             overlap_height_indices = np.where(fd_mask_flat==2)[0]
# # #             overlap_width_indices = np.where(fd_mask_flat==2)[1]

# # #             # add the map to a full disk normalised array
# # #             fd_overlap_map[map_bottom:map_top, map_left:map_right] = map.data

# # #             # overlap for higher ntv and doppler values, and get most extreme magnetogram values
# # #             y_loc = 0
# # #             x_loc = 0
# # #             for j in np.arange(0,len(overlap_height_indices)):
# # #               y_loc = overlap_height_indices[j]
# # #               x_loc = overlap_width_indices[j]
# # #               if fd_overlap_map[y_loc, x_loc] > fd_data[y_loc, x_loc]:
# # #                 fd_data[y_loc, x_loc] = fd_overlap_map[y_loc, x_loc]

# # #     # convert values where no raster is present to NaNs so they are not plotted
# # #     fd_data = np.where(fd_mask_flat==0, np.nan, fd_data)

# # #     # convert any pixels where the value is exactly 0 to NaNs so they are not plotted
# # #     fd_data = np.where(fd_data==0, np.nan, fd_data)

# # #     # create the full disk maps
# # #     fd_map = sunpy.map.Map(fd_data, fd_header)

# # #     # if looking at velocity maps, correct for line of sight effects
# # #     if measurement == 'vel':
# # #         solar_x_len = fd_map.data.shape[1]
# # #         solar_y_len = fd_map.data.shape[0]
# # #         # solar_x = np.arange(-fd_size/2, fd_size/2, map_dx)
# # #         # solar_y = np.arange(-fd_size/2, fd_size/2, map_dy)

# # #         # # create an empty numpy array for the pointing coordinates
# # #         # solar_x_array = np.empty([solar_y_len,solar_x_len])
# # #         # solar_y_array = np.empty([solar_y_len,solar_x_len])

# # #         # # fill the solar_x array
# # #         # for i in range(0,solar_y_len):
# # #         #   solar_x_array[i,:] = solar_x

# # #         # # fill the solar_y array
# # #         # for i in range(0,solar_x_len):
# # #         #   solar_y_array[:,i] = solar_y

# # #         # Using np.linspace to ensure the arrays have the exact expected length
# # #         solar_x = np.linspace(-fd_size/2, fd_size/2, solar_x_len)
# # #         solar_y = np.linspace(-fd_size/2, fd_size/2, solar_y_len)

# # #         # create meshgrid for solar_x and solar_y
# # #         solar_x_array, solar_y_array = np.meshgrid(solar_x, solar_y)

# # #         # convert solar_x and solar_y to longitude and latitude
# # #         obs_time = fd_map.meta['date_obs']
# # #         helioprojective_coordinates = SkyCoord(solar_x_array * u.arcsec, solar_y_array * u.arcsec, frame=frames.Helioprojective(observer="earth", obstime=obs_time))
# # #         heliographic_coordinates = helioprojective_coordinates.transform_to(frames.HeliographicStonyhurst)
# # #         latitude = heliographic_coordinates.lat.to(u.deg).value
# # #         longitude = heliographic_coordinates.lon.to(u.deg).value

# # #         # apply the doppler velocity correction
# # #         cos_mu = np.cos(latitude*((2*np.pi)/360)) * np.cos(longitude*((2*np.pi)/360))

# # #         for i in range (0, fd_map.data.shape[0]):
# # #             for j in range (0, fd_map.data.shape[1]):
# # #                 fd_map.data[i][j] = fd_map.data[i][j] / cos_mu[i][j]


# # #     # parse the date from first_map.meta["date_beg"] in format YYYY-MM-DDTHH:MM:SS.sss
# # #     fd_map_datetime = datetime.strptime(first_map.meta["date_beg"], "%Y-%m-%dT%H:%M:%S.%f")
# # #     fd_map_datetime = fd_map_datetime.strftime("%Y%m%d_%H%M%S")
    
# # #     # path = './data_eis/full_disks/'
# # #     path = base_path + 'full_disks/'
# # #     if not os.path.exists(path):
# # #         os.makedirs(path)

# # #     # save filename should be in format eis_YYYYMMDD_HHMMSS.{emission_line}.{measurement}..fits
# # #     save_filename = f'eis_{fd_map_datetime}.{line_id}.{measurement}.fits'
# # #     fd_map.save(path+save_filename, overwrite=True)

# # #     # plot the full disk non-thermal velocity map
# # #     fig = plt.figure()
# # #     ax = plt.subplot(projection=fd_map)

# # #     if measurement == 'int':
# # #         im = fd_map.plot(cmap='Oranges')
# # #     elif measurement == 'vel':
# # #         im = fd_map.plot(cmap='RdBu_r')
# # #         im.set_norm(plt.Normalize(vmin=-10, vmax=10))
# # #     elif measurement == 'wid':
# # #         im = fd_map.plot(cmap='viridis')
# # #     elif measurement == 'ntv':
# # #         im = fd_map.plot(cmap='inferno')
# # #         im.set_norm(plt.Normalize(vmin=0, vmax=40))
# # #     # fd_map_ntv.draw_limb()
# # #     # fd_map_ntv.draw_grid()
# # #     im = ax.get_images()
# # #     im_lims =  im[0].get_extent()
# # #     ax.set_aspect(abs((im_lims[1]-im_lims[0])/(im_lims[3]-im_lims[2])))
# # #     # plt.title(f'Hinode EIS non-thermal velocity (km/s)\n{fd_map.meta["date_beg"]} - {fd_map.meta["date_end"]}')
# # #     plt.colorbar(extend='both')
# # #     plt.savefig(path+save_filename.replace('.fits', '.png'))
# # #     if high_res:
# # #         plt.savefig(path+save_filename.replace('.fits', '.pdf'), dpi=1000)
# # #     plt.close()

# # #     return True