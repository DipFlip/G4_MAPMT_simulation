import sys
import uproot as ur
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

pixel_widths = [6.25, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.25]
pixel_borders = np.cumsum(pixel_widths)
bin_widths = [1.25] + [1] * 46 + [1.25]
bin_borders = np.cumsum(bin_widths)
bin_list_with_zero = np.insert(bin_borders,0,0)
bin_centers = (bin_list_with_zero[1:] + bin_list_with_zero[:-1]) / 2

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
    """returns the bin that the light originates from (there are 48*48 bins)"""
    x += 24.25 #Shift so that 0,0 is in the lower left corner
    y += 24.25 #Shift so that 0,0 is in the lower left corner
    bin_col = np.searchsorted(bin_borders, x) #0 is the leftmost column
    bin_row = np.searchsorted(bin_borders, y) #0 is the lowest row
    return bin_col + 48*bin_row

def load_G4_root_file(file_to_load):
    tfile = ur.open(file_to_load)
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
    total_number_photons = sum(len(l) for l in x_array)
    #sum the photon counts in each pixel for every event
    print("Summing photons...")
    photon_counts = np.zeros((2304, 64)) #address by [event, pixel-1]
    #x_origin_bins = np.zeros((total_number_events, 64)) #address by [event, pixel-1]
    #y_origin_bins = np.zeros((total_number_events, 64)) #address by [event, pixel-1]
    for event in np.arange(total_number_events):
        if x_array[event].size > 0:
            for x, y, xo, yo, vid in zip(x_array[event], y_array[event], x_origin[event], y_origin[event], vid_array[event]):
                # print(f"vid:{vid}, xo {xo}, xbin: {get_scan_bin(xo,yo)[0]}, ybin: {get_scan_bin(xo,yo)[1]}, yo:{yo}")
                if vid == 2:
                    pixel = G4_coordinates2pixel(x,y)
                    bin_number = get_scan_bin(xo, yo)
                    #x_origin_bins[event, pixel-1] = xbin
                    #y_origin_bins[event, pixel-1] = ybin
                    if  1 <= pixel <= 64:
                        # print(f"placing count in pixel {pixel} in bin {bin_number}")
                        photon_counts[bin_number, pixel-1] += 1

    # save as a dataframe
    print("Making dataframe...")
    flat_photon_counts_list = np.asarray([item for sublist in photon_counts for item in sublist], dtype=np.int16)
    flat_xo_list = np.tile(np.repeat(bin_centers,64),48)
    flat_yo_list = np.repeat(bin_centers,64*48)
    # flat_xo_bin_list = np.asarray([item for sublist in x_origin_bins for item in sublist], dtype=np.uint8)
    # flat_yo_bin_list = np.asarray([item for sublist in y_origin_bins for item in sublist], dtype=np.uint8)
    flat_pixel_list = np.tile(np.arange(1,65, dtype=np.uint8),2304)
    # events = np.asarray([int(row / 64) for row in range(nrows*64)], dtype=np.uint16)
    dataframe = pd.DataFrame()

    dataframe['pno'] = flat_pixel_list
    dataframe['ph'] = flat_photon_counts_list
    dataframe['xo'] = flat_xo_list
    dataframe['yo'] = flat_yo_list
    return dataframe


#load the file
print(f"loading file {sys.argv[1]}")
if sys.argv[1].endswith(".root"):
    df = load_G4_root_file(sys.argv[1])
elif sys.argv[1].endswith(".pkl"):
    import pickle
    df = pickle.load( open(sys.argv[1], 'rb'))
else:
    print("load a .root or .pkl file")
# dataframe['evtno'] = events



def plot_at_position(x, y):
    phs = np.zeros(64)
    for pno in np.arange(1,65):
        phs[pno-1] = df.ph[(df.xo.between(x-1,x+1)) & (df.yo.between(y-1,y+1)) & (df.pno == pno)].sum()
    phs = np.reshape(phs, (8,8))
    sns.heatmap(phs, annot=True, fmt='.0f')
    plt.show()

def plot_max_at_each_position():
    maxes = np.zeros((48,48))
    for xid, xbin in enumerate(bin_centers):
        for yid, ybin in enumerate(bin_centers):
            maxes[yid, xid] = df.ph[(df.xo == xbin) & (df.yo == ybin)].max()
    ax = sns.heatmap(maxes)
    ax.invert_yaxis()
    plt.show()
    return maxes