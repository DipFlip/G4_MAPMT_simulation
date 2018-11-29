import sys
import uproot as ur
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

n_bins = 194 #cover all area with 0.25 mm bins
n_bins = 61 #cover 2*2 pixels with 0.2 mm bins
pixel_widths = [6.25, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.25]
pixel_borders = np.cumsum(pixel_widths)
# bin_widths = [1.25] + [1] * 46 + [1.25]
bin_widths = [0.25] * n_bins #0.25 mm bins
bin_widths = [0.2] * n_bins #0.2 mm bins around 2*2 pixels
# bin_widths = [0.5] * 98 #0.5 mm bins
bin_borders = np.cumsum(bin_widths)
bin_borders += 18.15 #to center around center pixels
bin_list_with_zero = np.insert(bin_borders,0, 18.15)
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

def G4_coordinate2_col_row(x,y):
    """returns the column and row number"""
    x += 24.25 #Shift so that 0,0 is in the lower left corner
    y += 24.25 #Shift so that 0,0 is in the lower left corner
    bin_col = np.searchsorted(bin_borders, x) #0 is the leftmost column
    bin_row = np.searchsorted(bin_borders, y) #0 is the lowest row
    return bin_col, bin_row

def get_scan_bin(x,y):
    """returns the bin that the light originates from (there are n_bins*n_bins bins)"""
    x += 24.25 #Shift so that 0,0 is in the lower left corner
    y += 24.25 #Shift so that 0,0 is in the lower left corner
    bin_col = np.searchsorted(bin_borders, x) #0 is the leftmost column
    bin_row = np.searchsorted(bin_borders, y) #0 is the lowest row
    # return np.searchsorted(bin_borders, x) + n_bins*np.searchsorted(bin_borders, y) + 1
    return (bin_col + n_bins*bin_row)

def is_in_groove(xo, yo, zo):
    xo += 24.25 #Shift so that 0,0 is in the lower left corner
    yo += 24.25 #Shift so that 0,0 is in the lower left corner
    zo += 0.5 #Shift so that 0 is the scintillator surface facing MAPMT
    """Returns True if the photon originates in a groove."""
    groove_width = 0.2 #mm
    groove_depth = 0.5 #mm
    groove_depth = 1 #mm
    grooves = [18.25, 24.25, 30.25]
    for groove in grooves:
        if groove - groove_width < xo < groove + groove_width:
            if zo < groove_depth:
                # print(f"xo {xo}, yo {yo}, zo {zo} was in groove")
                return True 
        if groove - groove_width < yo < groove + groove_width:
            if zo < groove_depth:
                # print(f"xo {xo}, yo {yo}, zo {zo} was in groove")
                return True 
    return False

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
    return vid_array, x_array, y_array, x_origin, y_origin, z_origin, total_number_events

def root_to_dataframe(vid_array, x_array, y_array, x_origin, y_origin, z_origin, total_number_events):
    #sum the photon counts in each pixel for every event
    print("Summing photons...")
    photon_counts = np.zeros((n_bins*n_bins, 64)) #address by [event, pixel-1]
    #x_origin_bins = np.zeros((total_number_events, 64)) #address by [event, pixel-1]
    #y_origin_bins = np.zeros((total_number_events, 64)) #address by [event, pixel-1]
    for event in np.arange(total_number_events):
        if x_array[event].size > 0:
            for x, y, xo, yo, zo, vid in zip(x_array[event], y_array[event], x_origin[event], y_origin[event], z_origin[event], vid_array[event]):
                # print(f"vid:{vid}, xo {xo}, xbin: {get_scan_bin(xo,yo)[0]}, ybin: {get_scan_bin(xo,yo)[1]}, yo:{yo}")
                if vid == 2:
                    # if not is_in_groove(xo, yo, zo):
                    pixel = G4_coordinates2pixel(x,y)
                    bin_number = get_scan_bin(xo, yo)
                    if  1 <= pixel <= 64:
                        photon_counts[bin_number, pixel-1] += 1

    # save as a dataframe
    print("Making dataframe...")
    flat_photon_counts_list = np.asarray([item for sublist in photon_counts for item in sublist], dtype=np.int16)
    flat_xo_list = np.tile(np.repeat(bin_centers,64),n_bins)
    flat_yo_list = np.repeat(bin_centers,64*n_bins)
    # flat_xo_bin_list = np.asarray([item for sublist in x_origin_bins for item in sublist], dtype=np.uint8)
    # flat_yo_bin_list = np.asarray([item for sublist in y_origin_bins for item in sublist], dtype=np.uint8)
    flat_pixel_list = np.tile(np.arange(1,65, dtype=np.uint8),n_bins*n_bins)
    # events = np.asarray([int(row / 64) for row in range(nrows*64)], dtype=np.uint16)
    dataframe = pd.DataFrame()

    dataframe['pno'] = flat_pixel_list
    dataframe['ph'] = flat_photon_counts_list
    dataframe['xo'] = flat_xo_list
    dataframe['yo'] = flat_yo_list
    return dataframe

def plot_light_cone_at_pos(x_pos, y_pos):
    delta = 0.5 #mm
    photons_counted = np.zeros((n_bins,n_bins))
    for event in np.arange(total_number_events):
        if x_array[event].size > 0:
            for x, y, xo, yo, vid in zip(x_array[event], y_array[event], x_origin[event], y_origin[event], vid_array[event]):
                # print(f"vid:{vid}, xo {xo}, xbin: {get_scan_bin(xo,yo)[0]}, ybin: {get_scan_bin(xo,yo)[1]}, yo:{yo}")
                if vid == 2:
                    if x_pos-delta < xo + 24.25 < x_pos+delta:
                        if y_pos-delta < yo + 24.25 < y_pos+delta:
                            photons_counted[G4_coordinate2_col_row(x,y)] += 1
                            # print(f"add photon to colrow {G4_coordinate2_col_row(xo,yo)} from xo {xo}, yo {yo}.")
    return photons_counted
#load the file
if len(sys.argv) > 1:
    print(f"loading file {sys.argv[1]}")
    if sys.argv[1].endswith(".root"):
        vid_array, x_array, y_array, x_origin, y_origin, z_origin, total_number_events = load_G4_root_file(sys.argv[1])
    elif sys.argv[1].endswith(".pkl"):
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
    maxes = np.zeros((n_bins,n_bins))
    for xid, xbin in enumerate(bin_centers):
        for yid, ybin in enumerate(bin_centers):
            maxes[yid, xid] = df.ph[(df.xo == xbin) & (df.yo == ybin)].max()
    ax = sns.heatmap(maxes, vmin = 200)
    ax.invert_yaxis()
    for pos in pixel_borders*4:
        ax.axhline(pos, linestyle='-', color='w') # horizontal lines
        ax.axvline(pos, linestyle='-', color='w') # vertical lines
    plt.show()
    return maxes

def calculate_multiplicity_at_each_position(threshold):
    """returns a 2D array of multiplicities for each position"""
    mults = np.zeros((n_bins,n_bins))
    # mults = np.zeros((98,98))
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
    plt.xlim([18.25,30.25])
    plt.ylim([18.25,30.25])
    if title is not None:
        plt.title(title)
    #plt.show()
    return fig, ax

def plot_M1_many(mults_array):
    fig1, ax1 = plt.subplots()
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_choice = list(colors[i] for i in [0, 4, 9, 8, 7])
    thresh = np.transpose(list(mults_array[1].keys()))
    m_counts = np.zeros((len(thresh), len(mults_array)))
    for idm, mults in enumerate(mults_array):
        mults_reduced = mults.copy()
        for idt, t in enumerate(thresh):
            mults_reduced[t] = mults[t][5:35,5:35]
        n_illuminated_bins = mults_reduced[thresh[0]].shape[0]
        for idt, t in enumerate(thresh):
            m_counts[idt, idm] = np.sum(mults_reduced[t] == 1) / (n_illuminated_bins*n_illuminated_bins)
    for multip_line, col in zip(np.transpose(m_counts[:,:]), color_choice):
        ax1.plot(thresh, multip_line, color=col)
    ax1.legend(['Ungrooved','0.5 mm deep grooves,\n facing MAPMT', '0.5 mm deep grooves,\n facing away from MAPMT', 'Grooved all through'])
    from matplotlib.ticker import ScalarFormatter
    [ax.xaxis.set_major_formatter(ScalarFormatter()) for ax in [ax1]]
    ax1.set_ylim([-0, 1])
    [ax.set_xlabel('Threshold / arb. unit') for ax in [ax1]]
    [ax.set_ylabel('$M = 1$ fraction') for ax in [ax1]]
    plt.show()

def plot_0_to_5(mults):
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    thresh = np.transpose(list(mults.keys()))
    m_counts = np.zeros((len(thresh),6))
    mults_reduced = mults.copy()
    for idt, t in enumerate(thresh):
        mults_reduced[t] = mults[t][5:35,5:35]
    n_illuminated_bins = mults_reduced[thresh[0]].shape[0]
    for idt, t in enumerate(thresh):
        for m in [0,1,2,3,4]:
            m_counts[idt, m] = np.sum(mults_reduced[t] == m) / (n_illuminated_bins*n_illuminated_bins)
            # m_counts[idt, m] = np.sum(ungrooved_mults_reduced_big[t] == m) / (24*24)
        # m_counts[idt, 0] = 1 - sum(m_counts[idt,:])
        m_counts[idt, 5] = np.sum(mults_reduced[t] > 4) / (n_illuminated_bins*n_illuminated_bins)
    for multip_line, col in zip(np.transpose(m_counts[:,(0,-1)]), ['k','c']):
        ax1.semilogx(thresh, multip_line, color=col)
    for multip_line in np.transpose(m_counts[:,(1,2,3,4)]):
        ax2.semilogx(thresh, multip_line)
    ax1.legend(['$M = 0$','$M > 4$'])
    ax2.legend(['$M = 1$','$M = 2$','$M = 3$','$M = 4$'])
    from matplotlib.ticker import ScalarFormatter
    [ax.xaxis.set_major_formatter(ScalarFormatter()) for ax in [ax1, ax2]]
    ax2.set_xticks([10,100,500])
    ax1.set_ylim([-0.1, 1.1])
    ax2.set_ylim([-0.0, 1])
    [ax.set_xlabel('Threshold / arb. unit') for ax in [ax1, ax2]]
    [ax.set_ylabel('Fraction of area') for ax in [ax1, ax2]]
    plt.show()

def usable_area():
    usables_grooved = np.zeros(len(thresh))
    usables_ungrooved = np.zeros(len(thresh))
    usables_grooved_ref = np.zeros(len(thresh3_old))
    usables_ungrooved_ref = np.zeros(len(thresh3_old))
    usables_grooved_ref2 = np.zeros(len(thresh3))
    usables_ungrooved_ref2 = np.zeros(len(thresh3_new))

    for idt, t in enumerate(thresh):
        usables_grooved[idt] = np.sum(grooved_mults[t] == 1) / (24*24)
        usables_ungrooved[idt] = np.sum(ungrooved_mults[t] == 1) / (24*24)
    for idt, t in enumerate(thresh3_old):
        usables_grooved_ref[idt] = np.sum(grooved_ref_mults[t] == 1) / (24*24)
        usables_ungrooved_ref[idt] = np.sum(ungrooved_ref_mults[t] == 1) / (24*24)
    for idt, t in enumerate(thresh3):
        usables_grooved_ref2[idt] = np.sum(grooved_ref_mults2[t] == 1) / (48*48)
    for idt, t in enumerate(thresh3_new):
        usables_ungrooved_ref2[idt] = np.sum(ungrooved_ref_mults2[t] == 1) / (48*48)

    plt.plot(thresh, usables_ungrooved, '-x', label='Ungrooved')
    plt.plot(thresh, usables_grooved, '-o', label='Grooved')
    plt.plot(thresh3_old, usables_grooved_ref, '-^', label='Grooved with reflector')
    plt.plot(thresh3_old, usables_ungrooved_ref, '-v', label='Ungrooved with reflector, 0.5 mm bins')
    plt.plot(thresh3, usables_grooved_ref2, '->', label='Grooved with reflector, 0.25 mm bins')
    plt.plot(thresh3_new, usables_ungrooved_ref2, '-<', label='Ungrooved with reflector, 0.25 mm bins')
    plt.ylabel('Fraction area with multiplicity 1')
    plt.xlabel('Threshold / arb. unit')
    plt.legend()
    plt.show()
    return usables_grooved, usables_ungrooved

def make_dict_of_mults():
    new_mults_dict = {el:[] for el in t_tot2}
    for t in t_tot2:
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
thresh3 = np.arange(50,501,50)
thresh3_new = np.arange(50,501,50)/4
thresh3_old = np.arange(100,2001,100)
ungrooved_mults = pickle.load( open('pickles/cent_ungrooved_mults.pkl', 'rb'))
grooved_mapmt_side_half = pickle.load( open('pickles/new_halfgroove_mults_0.25.pkl', 'rb'))
grooved_mapmt_side_all = pickle.load( open('pickles/new_allgroove_mults_0.25.pkl', 'rb'))
filtered_grooved_mapmt_side_all_mults = pickle.load( open('pickles/filtered_allgroove_mults_0.20.pkl', 'rb'))
filtered_grooved_mapmt_side_half_mults = pickle.load( open('pickles/filtered_halfgroove_mults_0.20.pkl', 'rb'))
filtered_nogrooved_mapmt_rachel_mults = pickle.load( open('pickles/filtered_nogroove_rachel_mults_0.20.pkl', 'rb'))
filtered_nogrooved_mapmt_laura_mults = pickle.load( open('pickles/filtered_nogroove_laura_mults_0.20.pkl', 'rb'))
filtered_nogrooved_mapmt_mults = pickle.load( open('pickles/filtered_nogroove_mults_0.20.pkl', 'rb'))
filtered_topgrooved_mapmt_mults = pickle.load( open('pickles/filtered_topgroove_rachel_mults_0.20.pkl', 'rb'))
grooved_mults = pickle.load( open('pickles/cent_grooved_mults.pkl', 'rb'))
corner_grooved_mults = pickle.load( open('pickles/corn_ungrooved_mults.pkl', 'rb'))
grooved_ref_mults = pickle.load( open('pickles/grooved_ref_mults.pkl', 'rb'))
ungrooved_ref_mults = pickle.load( open('pickles/ungrooved_ref_mults.pkl', 'rb'))
grooved_ref_mults2 = pickle.load( open('pickles/grooved_ti_0.25_100_1000.pkl', 'rb'))
ungrooved_ref_mults2 = pickle.load( open('pickles/ungrooved_ti_0.25_100_1000.pkl', 'rb'))


# df = pickle.load( open('pickles/grooved_ref_df0.25.pkl', 'rb'))
# df = ungrooved_mults
#creating big run like in the alpha paper
# t1 = np.arange(10, 601, 10)
# t_low = np.arange(0, 126, 5)
# t_lowlow = np.append(np.arange(0,10,1),np.arange(10,31,2))
# t_tot = np.append(t_low, t1[12:])
# t_tot2 = np.append(t_lowlow, t_tot[7:])
# t_tot2 = t_tot2[::5]

t1 = np.arange(10, 601, 20)
t_low = np.arange(0, 126, 5)
t_lowlow = np.append(np.arange(0,10,1),np.arange(10,31,2))
t_tot = np.append(t_low, t1[6:])
t_tot2 = np.append(t_lowlow, t_tot[7:])
t_n = np.append(np.arange(0,10,1),np.arange(10,31,2))
t_n2 = np.append(t_n, np.arange(33, 55, 3))
t_n3 = np.append(t_n2, np.arange(56, 147, 6))
t_n4 = np.append(t_n3, np.arange(150, 300, 20))
t_tot2 = np.append(t_n4, np.arange(300, 601, 30))
# t_tot2 = t_tot2[::3]
np.append(np.append(np.append(np.arange(0,10,1),np.arange(11,31,2)), np.append(np.arange(33, 55, 3), np.arange(56, 147, 6)))
    , np.append(np.arange(150, 300, 20), np.arange(300, 600, 30)))