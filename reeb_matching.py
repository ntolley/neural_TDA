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
import cv2
import os
sns.set()
    

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

def save_graph(G, current_node,save_name,count,save_dir='D:/Jones_Lab/MRG_rotating/'):
    plt.ioff()
    new_points = np.array([G.nodes[node_idx]['Position'] for node_idx in list(G.nodes)])
    new_connectivity = np.array([list(e) for e in list(G.edges())])



    fig = plt.figure(figsize = (8,6))
    ax = plt.axes(projection='3d')
    plt.axis('off')

    ax.set_xlim(0,1100)
    ax.set_ylim(0,600)
    ax.set_zlim(-1,4)
    ax.view_init(25, count*0.5)

    num_pairs = new_connectivity.shape[0]
    for pair in range(num_pairs):
        pairID = new_connectivity[pair]
        pos_start, pos_end = G.nodes[pairID[0]]['Position'], G.nodes[pairID[1]]['Position']
        xdata, ydata, zdata = [pos_start[0], pos_end[0]], [pos_start[1], pos_end[1]], [pos_start[2], pos_end[2]]

        ax.plot(xdata,ydata,zdata, 'k', linewidth=0.2)

    ax.scatter(new_points[:, 0], new_points[:, 1], new_points[:, 2], 'b')

    current_point = G.nodes[current_node]['Position']
    ax.scatter(current_point[0],current_point[1],current_point[2],'r', s=100)

    fig.savefig(save_dir + save_name + '.png')
    

    return


# adds or removes nodes based on interval
def edge_edit(B, start_node, end_node, interval_dict):
    # print(str([start_node,end_node]))
    interval_points = interval_dict['interval_points']

    half_width = interval_dict['half_width']

    #Height of current node
    start_val = B.nodes[start_node]['Position'][2]
    end_val = B.nodes[end_node]['Position'][2]

    #Upper bound index for the start and ending node heights
    start_bound_idx, end_bound_idx = np.sum(interval_points < start_val), np.sum(interval_points < end_val)

    #Update node positions to center of interval
    B.nodes[start_node]['Position'][2] = interval_points[start_bound_idx] - half_width
    B.nodes[end_node]['Position'][2] = interval_points[end_bound_idx] - half_width

    #Represents how many intervals the edge spans
    bound_diff = abs(start_bound_idx - end_bound_idx)


    #Edge spans exactly 1 bound as desired, update node position and return
    if bound_diff == 1:
        # print(str([start_node,end_node]), 'skip')

        return None

    #Connected nodes in same interval, merge together and update attribute dictionary
    elif bound_diff == 0:
        # print(str([start_node,end_node]), 'merge')

        #Update attribute to indicate merging
        B.nodes[start_node]['Merged'].append(end_node)
        B.nodes[start_node]['New_Merge'].append(end_node)

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
        # print(str([start_node,end_node]), 'insert')

        B.remove_edge(start_node, end_node)
        
        #Create new nodes with unique id's
        num_insertions = bound_diff-1
        max_node = max(list(B.nodes())) + 1 
        new_nodes = [start_node] + list(range(max_node, max_node + num_insertions)) + [end_node]
        new_edges = [[new_nodes[i], new_nodes[i+1]] for i in range(len(new_nodes)-1)]

        #Update positions for new nodes
        # new_height = [np.mean(interval_points[start_bound_idx-1 + i : start_bound_idx + i ]) for i in range(bound_diff+1)] 
        # new_height = [interval_points[start_bound_idx + i] - half_width for i in range(bound_diff+1)] 
        new_height = np.linspace(B.nodes[start_node]['Position'][2], B.nodes[end_node]['Position'][2], bound_diff+1)
        new_x = np.linspace(B.nodes[start_node]['Position'][0], B.nodes[end_node]['Position'][0], bound_diff+1)
        new_y = np.linspace(B.nodes[start_node]['Position'][1], B.nodes[end_node]['Position'][1], bound_diff+1)

        new_attributes = {new_nodes[idx] : {'Position' : [new_x[idx], new_y[idx],new_height[idx]], 'Visited' : 1, 'Merged':[],'Inserted':[],'New_Merge':[]} for idx in range(1,len(new_nodes)-1)}

        #Insert connected nodes at center of interval
        B.add_edges_from(new_edges)
        nx.set_node_attributes(B,new_attributes)

        for node_idx in range(len(new_nodes)):
            B.node[new_nodes[node_idx]]['Inserted'].append(new_nodes)
     

        return None

#Recursive function that goes through all edges 
def graph_search(A, start, interval_dict, save, count=[0]):


    if A.nodes[start]['Visited'] == 1:
        return
    else:
        A.nodes[start]['Visited'] = 1
        neighbors = list(A.neighbors(start))

        #Insert or remove nodes based on interval position
        for neighbor_id in neighbors:

            #save graphs for movie
            if save == 1:
                count[0] = count[0] + 1 #Track number of calls, update by reference
                save_graph(A, start, 'g_' + str(count[0]), count[0])



            #If an edge was merged, the neighbor list updates to reflect new neighbors
            neighbor_update = edge_edit(A, start, neighbor_id, interval_dict)
            if neighbor_update != None:
                neighbors.extend(neighbor_update)
                continue #Go on to next neighbor

            else:
                graph_search(A,neighbor_id,interval_dict,save,count)
            

        
        return

def compute_intervals(node_points, num_intervals):

    min_elevation, max_elevation = min(node_points[:,2]), max(node_points[:,2])

    interval_width = (max_elevation-min_elevation)/(num_intervals-1)
    low, high = min_elevation - (interval_width/2), max_elevation + (interval_width/2)

    #Centers the minimum and maximum nodes at the min and max intervals
    # interval_bound = [[low + interval_width*pos, low + interval_width*(pos+1)] for pos in range(num_intervals)]
    interval_points = np.array([low + interval_width*pos for pos in range(num_intervals+1)])
    interval_dict = {'interval_points': interval_points, 'half_width':(interval_width/2)}

    return interval_dict


def MRG_clear_visits(A, first_pass):

    for node in list(A.nodes):
        A.nodes[node]['Visited'] = 0
        A.nodes[node]['New_Merge'] = []

        if first_pass == 1:
            A.nodes[node]['Merged'] = []
            A.nodes[node]['Inserted'] = []

    return

def make_movie(image_folder, save_dir, images):

    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape

    video = cv2.VideoWriter(save_dir, 0, 30, (width,height))

    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

    cv2.destroyAllWindows()
    video.release()