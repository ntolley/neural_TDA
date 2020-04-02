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
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
    

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

def plot_graph(G):
    new_points = np.array([G.nodes[node_idx]['Position'] for node_idx in list(G.nodes)])
    new_connectivity = np.array([list(e) for e in list(G.edges())])

    fig = plt.figure(figsize = (8,6))
    ax = plt.axes(projection='3d')

    num_pairs = new_connectivity.shape[0]
    for pair in range(num_pairs):
        pairID = new_connectivity[pair]
        pos_start, pos_end = G.nodes[pairID[0]]['Position'], G.nodes[pairID[1]]['Position']
        xdata, ydata, zdata = [pos_start[0], pos_end[0]], [pos_start[1], pos_end[1]], [pos_start[2], pos_end[2]]

        ax.plot(xdata,ydata,zdata, 'k', linewidth=0.2)

    ax.scatter(new_points[:, 0], new_points[:, 1], new_points[:, 2], 'b')

    plt.show()
    return

# adds or removes nodes based on interval
def edge_edit(B, start_node, end_node, interval_points):
    # print(str([start_node,end_node]))

    #Height of current node
    start_val = B.nodes[start_node]['Position'][2]
    end_val = B.nodes[end_node]['Position'][2]

    #Upper bound index for the start and ending node heights
    start_bound_idx, end_bound_idx = np.sum(interval_points < start_val), np.sum(interval_points < end_val)

    #Update node positions to center of interval
    B.nodes[start_node]['Position'][2] == np.mean(interval_points[start_bound_idx-1:start_bound_idx]) #*** Speed up with computing interval width ***
    B.nodes[end_node]['Position'][2] == np.mean(interval_points[end_bound_idx-1:end_bound_idx])

    #Represents how many intervals the edge spans
    bound_diff = abs(start_bound_idx - end_bound_idx)

    #Edge spans exactly 1 bound as desired, update node position and return
    if bound_diff == 1:
        print(str([start_node,end_node]), 'skip')

        return None

    #Connected nodes in same interval, merge together and update attribute dictionary
    elif bound_diff == 0:
        print(str([start_node,end_node]), 'merge')

        #Update attribute to indicate merging
        B.nodes[start_node]['Merged'].append(end_node)

        #Join end node neighbors to start node 
        merged_neighbors = list(B.neighbors(end_node))
        merged_neighbors.remove(start_node)
        B.remove_node(end_node)

        new_edges = [[start_node, new_neighbor] for new_neighbor in merged_neighbors]
        B.add_edges_from(new_edges)

        #update neighbors graph_search by reference

        return merged_neighbors

    #Edge spans more than one interval, subdivide
    elif bound_diff > 1:
        print(str([start_node,end_node]), 'insert')

        B.remove_edge(start_node, end_node)
        
        #Create new nodes with unique id's
        num_insertions = bound_diff-1
        max_node = max(list(B.nodes())) + 1 
        new_nodes = [start_node] + list(range(max_node + 1, max_node + num_insertions)) + [end_node]
        new_edges = [[new_nodes[i], new_nodes[i+1]] for i in range(len(new_nodes)-1)]

        #Update positions for new nodes
        new_height = [np.mean(interval_points[start_bound_idx-1 + i : start_bound_idx + i ]) for i in range(bound_diff+1)] 
        new_x = np.linspace(B.nodes[start_node]['Position'][0], B.nodes[end_node]['Position'][0], bound_diff+1)
        new_y = np.linspace(B.nodes[start_node]['Position'][1], B.nodes[end_node]['Position'][1], bound_diff+1)

        new_attributes = {new_nodes[idx] : {'Position' : [new_x[idx], new_y[idx],new_height[idx]], 'Visited' : 1, 'Merged':[],'Inserted':[]} for idx in range(1,len(new_nodes)-1)}

        #Insert connected nodes at center of interval
        B.add_edges_from(new_edges)
        nx.set_node_attributes(B,new_attributes)

        for node_idx in range(len(new_nodes)):
            B.node[new_nodes[node_idx]]['Inserted'].append(new_nodes)
     

        return None

#Recursive function that goes through all edges 
def graph_search(A, start, interval_points):
    if A.nodes[start]['Visited'] == 1:
        return
    else:
        A.nodes[start]['Visited'] = 1
        neighbors = list(A.neighbors(start))

        #Insert or remove nodes based on interval position
        for neighbor_id in neighbors:
            #If an edge was merged, the neighbor list updates to reflect new neighbors
            neighbor_update = edge_edit(A, start, neighbor_id, interval_points)
            if neighbor_update != None:
                neighbors.extend(neighbor_update)
                continue #Go on to next neighbor

            else:
                graph_search(A,neighbor_id,interval_points)
            

        
        return


   