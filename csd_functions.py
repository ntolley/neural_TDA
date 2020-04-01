import numpy as np
import scipy
from scipy.signal import decimate
from scipy import interpolate
import h5py
import load, load_hnn
import pandas as pd
import icsd
import neo
import quantities as pq
import os
from os.path import isfile, join
from os import listdir
    
#Adjust scalars to center of each contour range, inserts nodes at edges that span range
def resample_tree(node_data, edge_data, bounds):
    return

#Returns the minimum distance between two nodes
def min_node_dist(node_data):
    return


#Takes in a .param dictionary, parameter sweeps defined by an array of values, generates all combinations 
def dict_expand(array_dict):
    return

def csd_interp(data_path, ds_step):
    # calculate pixel (um) positions of model objects
    top_l5 = 1466
    soma_l5 = 177.5
    bottom_l5 = -72

    top_l2 = 1466
    soma_l2 = 1000
    bottom_l2 = 769

    num_contacts_hnn = 20
    hnn_values = [top_l2, soma_l2, bottom_l2, soma_l5, bottom_l5 ]
    perlayer = (top_l5 - bottom_l5)/(num_contacts_hnn - 5.0)
    spacing_um_hnn = perlayer

    # print("contact spacing (microns)",spacing_um_hnn)
    first_contact = top_l5+2*spacing_um_hnn
    last_contact = bottom_l5-2*spacing_um_hnn

    spacing_hnn = []
    for i,x in enumerate(np.linspace(first_contact,last_contact,num_contacts_hnn)):
        # print("contact %d: %.2f" % (i+1, x))
        spacing_hnn.append(x)


    #load data
    dir = data_path
    sampr_hnn,LFP_hnn,dt_hnn, tt_hnn, CSD_hnn, maxlfp_hnn, ntrial_hnn = load_hnn.loadHNNdir(dir,spacing_um_hnn)

    #Down sample fixed interval 
    # ds_step = 10
    tt_hnn = tt_hnn[::ds_step]

    # Average iCSD from HNN
    z_data_hnn = np.linspace(spacing_um_hnn*1E-6, 2300E-6, num_contacts_hnn) * pq.m  # [m]
    diam_hnn = 500E-6 * pq.m                              # [m]
    h_hnn = spacing_um_hnn * 1E-6 * pq.m                                 # [m]

    sigma_hnn = 0.3 * pq.S / pq.m                         # [S/m] or [1/(ohm*m)]
    sigma_top_hnn = 0.3 * pq.S / pq.m                     # [S/m] or [1/(ohm*m)]

    # Create dictionary with iCSD calculations
    icsd_hnn = {}
    for key in LFP_hnn:

        lfp_data_hnn = LFP_hnn[key] * 1E-6 * pq.V        # [uV] -> [V]
        lfp_data_hnn = lfp_data_hnn[:,::ds_step]

        # Input dictionaries for monkey data
        delta_input_hnn = {
            'lfp' : lfp_data_hnn,
            'coord_electrode' : z_data_hnn,
            'diam' : diam_hnn,          # source diameter
            'sigma' : sigma_hnn,        # extracellular conductivity
            'sigma_top' : sigma_hnn,    # conductivity on top of cortex
            'f_type' : 'gaussian',  # gaussian filter
            'f_order' : (3, 1),     # 3-point filter, sigma = 1.
        }

        delta_icsd_hnn = icsd.DeltaiCSD(**delta_input_hnn)
        icsd_hnn[key] = delta_icsd_hnn.get_csd()

    avgiCSD_hnn = load_hnn.getAvgERP(icsd_hnn, sampr_hnn, tt_hnn, maxlfp_hnn, ntrial_hnn)
    # avgiCSD_hnn = load_hnn.downsample(avgiCSD_hnn, sampr_hnn, 2000) #Down sample to LFP range

    #Prepare HNN data
    # set timerange from 30 to 80 ms
    # ewindowms_hnn = 35
    # tmin_hnn = 37
    # tmax_hnn = tmin_hnn + ewindowms_hnn
    # (idx_min_hnn, idx_max_hnn) = (np.where(tt_hnn==tmin_hnn)[0][0], np.where(tt_hnn==tmax_hnn)[0][0])
    # X_hnn = tt_hnn[idx_min_hnn:idx_max_hnn]
    X_hnn = tt_hnn

    # mask = np.zeros(len(tt_hnn), dtype=bool)
    # mask[idx_min_hnn:idx_max_hnn] = True
    # avgCSD_trim_X_hnn = avgiCSD_hnn[:,mask]
    avgCSD_trim_X_hnn = avgiCSD_hnn
    



    Y_hnn = range(avgiCSD_hnn.shape[0])
    # X_hnn = range(avgiCSD_hnn.shape[1])
    # CSD_spline_hnn=scipy.interpolate.RectBivariateSpline(Y_hnn, X_hnn, avgCSD_trim_X_hnn)
    CSD_spline_hnn=scipy.interpolate.RectBivariateSpline(Y_hnn, X_hnn, avgCSD_trim_X_hnn)


    # trim channels
    (idx_min_hnn, idx_max_hnn) = (0, 18)
    HNN_Y_plot = np.linspace(idx_min_hnn,idx_max_hnn,num=1000)
    # print("HNN channels: %d-%d"%(HNN_Y_plot.min(), HNN_Y_plot.max()))
    # mask = np.zeros(len(Y_hnn), dtype=bool)
    # mask[idx_min_hnn:idx_max_hnn] = True
    # avgCSD_trim_XY_hnn = avgCSD_trim_X_hnn[mask,:]

    # HNN_ymax = (spacing_hnn[int(HNN_Y_plot.min())]-spacing_hnn[int(HNN_Y_plot.max())])
    # HNN_ymin = 0

    Z_hnn = CSD_spline_hnn(HNN_Y_plot, X_hnn)

    # normalize to abs(Z.min)
    Z_hnn = Z_hnn/abs(Z_hnn.min())

    return Z_hnn


#Converts 2D array of to long array indexed by x,y,z coordinates
def grid2points(csd_grid):
    csd_points = [[r, c, csd_grid[r,c]] for r in range(csd_grid.shape[0]) for c in range(csd_grid.shape[1])]
    return csd_points

#Load tree data from predefined directory
def load_tree(prefix):

    #Preparing directory
    s_dir = 'D:/Jones_Lab/hnn_params/' + prefix + '/' +  prefix + '_skeleton/'
    # d_dir = 'D:/Jones_Lab/hnn_params/' + prefix + '/' +  prefix + '_data/'
    names = [f for f in listdir(s_dir) if isfile(join(s_dir, f))]
    names = [f_new.replace('_arcs.csv','') for f_new in names]
    names = [f_new.replace('_nodes.csv','') for f_new in names]
    names = np.unique(names)

    #Choose file to plot
    plot_name = names[6]


    # csd_surface_df = pd.read_csv(d_dir + plot_name + '.csv', sep=',')
    csd_nodes_df = pd.read_csv(s_dir + plot_name + '_nodes.csv', sep=',')
    csd_connectivity_df = pd.read_csv(s_dir + plot_name + '_arcs.csv', sep=',')


    # surface_points = np.array(csd_surface_df[['Points:0','Points:1','Points:2']]) #Use if extracting points from paraview
    # surface_points = np.array(csd_surface_df) #Use if using original surface file

    node_points = np.array(csd_nodes_df[['Points:0','Points:1','Points:2']])
    node_connectivity = np.array(csd_connectivity_df[['upNodeId','downNodeId']])
    node_color = np.array(csd_nodes_df[['CriticalType']])

    return node_points, node_connectivity, node_color