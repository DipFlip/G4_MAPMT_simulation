import sys
import uproot as ur
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

pixel_widths = [6.25, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.25]
pixel_borders = np.cumsum(pixel_widths)
# bin_widths = [1.25] + [1] * 46 + [1.25]
bin_widths = [0.5] * 98 #0.5 mm bins
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
    """returns the bin that the light originates from (there are 98*98 bins)"""
    x += 24.25 #Shift so that 0,0 is in the lower left corner
    y += 24.25 #Shift so that 0,0 is in the lower left corner
    bin_col = np.searchsorted(bin_borders, x) #0 is the leftmost column
    bin_row = np.searchsorted(bin_borders, y) #0 is the lowest row
    # return np.searchsorted(bin_borders, x) + 98*np.searchsorted(bin_borders, y) + 1
    return (bin_col + 98*bin_row)

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
    photon_counts = np.zeros((98*98, 64)) #address by [event, pixel-1]
    #x_origin_bins = np.zeros((total_number_events, 64)) #address by [event, pixel-1]
    #y_origin_bins = np.zeros((total_number_events, 64)) #address by [event, pixel-1]
    for event in np.arange(total_number_events):
        if x_array[event].size > 0:
            for x, y, xo, yo, vid in zip(x_array[event], y_array[event], x_origin[event], y_origin[event], vid_array[event]):
                # print(f"vid:{vid}, xo {xo}, xbin: {get_scan_bin(xo,yo)[0]}, ybin: {get_scan_bin(xo,yo)[1]}, yo:{yo}")
                if vid == 2:
                    pixel = G4_coordinates2pixel(x,y)
                    bin_number = get_scan_bin(xo, yo)
                    # print(f"placing count in pixel {pixel} in bin {bin_number}, from xo {xo}, yo {yo}")
                    #x_origin_bins[event, pixel-1] = xbin
                    #y_origin_bins[event, pixel-1] = ybin
                    if  1 <= pixel <= 64:
                        photon_counts[bin_number, pixel-1] += 1

    # save as a dataframe
    print("Making dataframe...")
    flat_photon_counts_list = np.asarray([item for sublist in photon_counts for item in sublist], dtype=np.int16)
    flat_xo_list = np.tile(np.repeat(bin_centers,64),98)
    flat_yo_list = np.repeat(bin_centers,64*98)
    # flat_xo_bin_list = np.asarray([item for sublist in x_origin_bins for item in sublist], dtype=np.uint8)
    # flat_yo_bin_list = np.asarray([item for sublist in y_origin_bins for item in sublist], dtype=np.uint8)
    flat_pixel_list = np.tile(np.arange(1,65, dtype=np.uint8),98*98)
    # events = np.asarray([int(row / 64) for row in range(nrows*64)], dtype=np.uint16)
    dataframe = pd.DataFrame()

    dataframe['pno'] = flat_pixel_list
    dataframe['ph'] = flat_photon_counts_list
    dataframe['xo'] = flat_xo_list
    dataframe['yo'] = flat_yo_list
    return dataframe


#load the file
if len(sys.argv) > 1:
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
    maxes = np.zeros((98,98))
    for xid, xbin in enumerate(bin_centers):
        for yid, ybin in enumerate(bin_centers):
            maxes[yid, xid] = df.ph[(df.xo == xbin) & (df.yo == ybin)].max()
    ax = sns.heatmap(maxes, vmin = 200)
    ax.invert_yaxis()
    for pos in pixel_borders*2:
        ax.axhline(pos, linestyle='-', color='w') # horizontal lines
        ax.axvline(pos, linestyle='-', color='w') # vertical lines
    plt.show()
    return maxes

def calculate_multiplicity_at_each_position(threshold):
    mults = np.zeros((98,98))
    for xid, xbin in enumerate(bin_centers):
        for yid, ybin in enumerate(bin_centers):
            mults[yid, xid] = df.ph[(df.xo == xbin) & (df.yo == ybin) & (df.ph > threshold)].count()
    return mults

def plot_mults(mults, title=None):
    cmap = plt.cm.terrain
    cmap.set_under(color='white')
    fig, ax = plt.subplots()
    plt.contourf(mults, levels=[-0.5,0.5,1.5,2.5,3.5,4.5], vmin=0.01, extent=[bin_centers.min(),bin_centers.max(),bin_centers.min(),bin_centers.max()], cmap=cmap)
    plt.colorbar(ticks=[0,1,2,3,4], label=f"Multiplicity")
    for pos in pixel_borders:
        ax.axhline(pos, linestyle='-', color='k') # horizontal lines
        ax.axvline(pos, linestyle='-', color='k') # vertical lines
    ax.set_aspect('equal')
    # plt.xlim([34,50])
    # plt.ylim([34,50])
    plt.xlim([16,32.5])
    plt.ylim([16,32.5])
    if title is not None:
        plt.title(title)
    #plt.show()
    return fig, ax

def usable_area():
    usables_grooved = np.zeros(len(thresh))
    usables_ungrooved = np.zeros(len(thresh))
    usables_grooved_ref = np.zeros(len(thresh3))
    usables_ungrooved_ref = np.zeros(len(thresh3))

    for idt, t in enumerate(thresh):
        usables_grooved[idt] = np.sum(grooved_mults[t] == 1) / (24*24)
        usables_ungrooved[idt] = np.sum(ungrooved_mults[t] == 1) / (24*24)
    for idt, t in enumerate(thresh3):
        usables_grooved_ref[idt] = np.sum(grooved_ref_mults[t] == 1) / (24*24)
        usables_ungrooved_ref[idt] = np.sum(ungrooved_ref_mults[t] == 1) / (24*24)
    plt.plot(thresh, usables_ungrooved, '-x', label='Ungrooved')
    plt.plot(thresh, usables_grooved, '-o', label='Grooved')
    plt.plot(thresh3, usables_grooved_ref, '-^', label='Grooved with reflector')
    plt.plot(thresh3, usables_ungrooved_ref, '-v', label='Ungrooved with reflector')
    plt.ylabel('Fraction area with multiplicity 1')
    plt.xlabel('Threshold / arb. unit')
    plt.legend()
    plt.show()
    return usables_grooved, usables_ungrooved

def make_dict_of_mults():
    new_mults_dict = {el:[] for el in thresh3}
    for t in thresh3:
        print(t)
        new_mults_dict[t] = calculate_multiplicity_at_each_position(t)
    return new_mults_dict

def draw_all_mult_maps():
    for t, m in corner_grooved_mults.items():
        fig, ax = plot_mults(m, title=f"Corner ungrooved, threshold {t}")
        fig.savefig(f"output/corner_T{t}_grooved.png")
    #for t, m in ungrooved_mults.items():
        #fig, ax = plot_mults(m, title=f"Ungrooved, threshold {t}")
        #fig.savefig(f"output/corner_T{t}_ungrooved.png")


thresh = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
thresh2 = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000]
thresh3 = np.arange(100,2001,100)
ungrooved_mults = pickle.load( open('pickles/cent_ungrooved_mults.pkl', 'rb'))
grooved_mults = pickle.load( open('pickles/cent_grooved_mults.pkl', 'rb'))
corner_grooved_mults = pickle.load( open('pickles/corn_ungrooved_mults.pkl', 'rb'))
grooved_ref_mults = pickle.load( open('pickles/grooved_ref_mults.pkl', 'rb'))
ungrooved_ref_mults = pickle.load( open('pickles/ungrooved_ref_mults.pkl', 'rb'))