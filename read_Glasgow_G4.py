import sys
import uproot as ur
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

pixel_widths = [6.25, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.25]
pixel_borders = np.cumsum(pixel_widths)
bin_widths = [1.25] + [1] * 62 + [1.25]
bin_borders = np.cumsum(bin_widths)

def xy2pixel(x, y):
    """
    returns pixel number (1-64) according to Hamamatsu's convention from x y coordinates (0-7)
    """
    return (7 - y) * 8 + (x + 1)

def G4_coordinates2pixel(x,y):
    """
    takes mm coordinates with (0,0) (x, y) being the center of the MAPMT
    and returns the pixel that's in. Choses the left/lower pixel if the 
    coordinates are exactly in the middle between two pixels.
    """
    x += 24.25 #Shift so that 0,0 is in the lower left corner
    y += 24.25 #Shift so that 0,0 is in the lower left corner
    pixel_col = np.searchsorted(pixel_borders, x) #0 is the leftmost column
    pixel_row = np.searchsorted(pixel_borders, y) #0 is the lowest row
    return xy2pixel(pixel_col, pixel_row)

def get_scan_bin(x,y):
    """returns the 12*12 bin that the light originates from"""
    x += 24.25 #Shift so that 0,0 is in the lower left corner
    y += 24.25 #Shift so that 0,0 is in the lower left corner
    bin_col = np.searchsorted(bin_borders, x) #0 is the leftmost column
    bin_row = np.searchsorted(bin_borders, y) #0 is the lowest row
    return bin_col, bin_row


#load the root file
tfile = ur.open(sys.argv[1])
ttree = tfile.get('PD')
arrays = ttree.arrays(['PD_vid', 'PD_x', 'PD_y', 'PD_vx', 'PD_vy', 'PD_vz'])
vid_array = arrays[b'PD_vid']
x_array = arrays[b'PD_x']
y_array = arrays[b'PD_y']
x_origin = arrays[b'PD_vx']
y_origin = arrays[b'PD_vy']
z_origin = arrays[b'PD_vz']
total_number_events = len(y_array)

# TODO only do this if nrows is not specified
nrows = total_number_events
total_number_events = 1000

#sum the photon counts in each pixel for every event
print("Summing photons...")
photon_counts = np.zeros((total_number_events, 64)) #address by [event, pixel-1]
x_origin_bins = np.zeros((total_number_events, 64)) #address by [event, pixel-1]
y_origin_bins = np.zeros((total_number_events, 64)) #address by [event, pixel-1]
for event in np.arange(total_number_events):
    for x, y, xo, yo, vid in zip(x_array[event], y_array[event], x_origin[event], y_origin[event], vid_array[event]):
        if vid == 2:
            pixel = G4_coordinates2pixel(x,y)
            xbin, ybin = get_scan_bin(xo, yo)
            x_origin_bins[event, pixel-1] = xbin
            y_origin_bins[event, pixel-1] = ybin
            #print(f"xo {xo}, xbin: {get_scan_bin(xo,yo)[0]}, ybin: {get_scan_bin(xo,yo)[1]}, yo:{yo}")
            if  1 <= pixel <= 64:
                photon_counts[event, pixel-1] += 1

# save as a dataframe
print("Making dataframe...")
flat_photon_counts_list = np.asarray([item for sublist in photon_counts for item in sublist], dtype=np.uint16)
flat_xo_bin_list = np.asarray([item for sublist in x_origin_bins for item in sublist], dtype=np.uint8)
flat_yo_bin_list = np.asarray([item for sublist in y_origin_bins for item in sublist], dtype=np.uint8)
flat_pixel_list = np.tile(np.arange(1,65, dtype=np.uint8),total_number_events)
events = np.asarray([int(row / 64) for row in range(nrows*64)], dtype=np.uint16)
dataframe = pd.DataFrame()

dataframe['pno'] = flat_pixel_list[:nrows*64]
dataframe['ph'] = flat_photon_counts_list[:nrows*64]
dataframe['xo'] = flat_xo_bin_list[:nrows*64]
dataframe['yo'] = flat_yo_bin_list[:nrows*64]
dataframe['evtno'] = events