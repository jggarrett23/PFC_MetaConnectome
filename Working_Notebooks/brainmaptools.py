import pickle
import numpy as np
import networkx as nx
import pandas as pd
import operator
import scipy
import math
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
from collections import OrderedDict
#from brainx import util, detect_modules, modularity


#Utils for workspace files

def build_key_codes_from_workspaces(workspaces, datadir):
    """Return a list of unique Brainmap ids from a list of workspace files. 
    Function reads the individual files and builds a unique key from column one and five
    requires pandas imported as pd 
    
    
    Parameters:
    ------------
    workspaces: a python list
                list of filenames for workspaces
    datadir: string
                location of where all the workspace csvs are
    
    
    Returns:
    ----------
    key_codes: a python nested list of numeric key codes representing brainmap study+contrast unique ids
               every index is the list of key codes in the corresponding workspace. 
               length of list = number of workspaces
               length of each sublist= number of keys in workspace
    
    """
    key_codes=[]
    for x in workspaces:
        file =datadir+x[:-1]
        df = pd.read_csv(file, names = ['zero','one','two', 'three','four','five','six','seven','eight', 'nine'])
        true_only= df[df['zero']>0]
        one_col = true_only['one']
        five_col = true_only['five']
        key_codes.append(one_col.map(str).values + five_col.map(str).values)
    return key_codes

def build_keycodes_from_excel_csv(excel_csv_file):
    df=pd.read_csv(excel_csv_file)
    regionlist=df.keys()
    relabel_dict={idx:x[:] for idx, x in enumerate(regionlist)}
    
    keycodes=[]
    for key in df.keys():
        studies_in_region=df[key][pd.notnull(df[key])]
        nstudies=len(studies_in_region)
        regions_in_study_list=[]
        [regions_in_study_list.append(studies) for studies in studies_in_region]
        keycodes.append(regions_in_study_list)
    return keycodes


def build_jaccard(key_codes):
    array_dims=len(key_codes), len(key_codes)
    jaccard=np.zeros(shape=array_dims)
    for col_idx, x in enumerate(key_codes):
        for row_idx, y in enumerate(key_codes):
            intersect=float(len(set(x) & set(y)))
            union=float(len(list(set(x) | set(y))))
            if union:
                jaccard[col_idx, row_idx]=intersect/union
    return jaccard
    

def build_n_coactives_array(key_codes):
    array_dims=len(key_codes), len(key_codes)
    n_coactives_array =np.zeros(shape=array_dims)
    
    for col_idx, x in enumerate(key_codes):
        for row_idx, y in enumerate(key_codes):
            n_coactives_array[col_idx, row_idx]=len(set(x) & set(y))
    return n_coactives_array
    
def normalize_n_coactives_array(n_coactives_array):
    #Use sklearn's normalize function w L2 remove diagonal first
    #This is not a preferred method
    norm_coactives_array=n_coactives_array.copy()
    np.fill_diagonal(norm_coactives_array, 0)
    
    return norm_coactives_array

def build_region_labels_dict(regionlist):
    relabel_dict ={idx:x[:-5] for idx, x in enumerate(regionlist)}
    return relabel_dict

#Different methods for making control graphs from BrainMap 
def select_n_random_keycodes(keycodes, nstudies):
    unique_keycodes=[]
    for x in keycodes:
        unique_keycodes.extend(x)
    
    unique_keycodes=set(unique_keycodes)
    rands=np.random.rand(len(unique_keycodes))
    random_keys_tuple=zip(rands, unique_keycodes)
    selected_tuple=sorted(random_keys_tuple, reverse=1)[0:nstudies]
    selected_keys=zip(*selected_tuple)[1]
    selected_keycodes=[set(selected_keys) & set(x) for x in keycodes]
    return selected_keycodes

def build_average_graph_from_random_keycodes(keycodes, nstudies, niters):
    jaccards=np.ndarray([len(keycodes),len(keycodes), niters])
    for x in range(niters):
        select_n_keys=select_n_random_keycodes(keycodes, nstudies)
        jaccard=build_jaccard(build_n_coactives_array(select_n_keys))
        jaccards[:,:,x]=jaccard
    G=nx.from_numpy_matrix(np.mean(jaccards, 2))
    return G
 
def number_of_contrasts(key_codes):
    contrast_list=[]
    for x in key_codes:
        contrast_list.extend(x)
    ncontrasts=len(set(contrast_list))
    return ncontrasts
     
def significant_connection_threshold(in_array,total_contrasts,threshold):
    """Function thresholds the array to keep edges that are statistically significant
    
    Parameters:
    -----------
    in_array: coactivation matrix
    array_lenght: dimension of the array
    total_contrasts: total number of contrasts reported by all papers taken from database
    threshold: threshold for significance
    
    Return:
    -------
    Returns thresh_array with significant connections kept and nonsignificant connections made to 0"""
    
    thresh_array=in_array.copy()
    for x in range(np.shape(thresh_array)[0]):
        for y in range(np.shape(thresh_array)[0]):
            if x != y: 
                # p=m/N m=independent activations of region X; N=total number of contrasts
                p=thresh_array[x][x]/total_contrasts 
                # null hypothesis is Binomial distribution of (k;n,p)* Binomial distribution of (m-k;N-n,p)
                # k=Coactivation of region X & Y, n=independent activation of region Y
                null=(scipy.stats.binom.pmf(thresh_array[x][y], thresh_array[y][y], p))*(scipy.stats.binom.pmf((thresh_array[x][x]-thresh_array[x][y]), (total_contrasts-thresh_array[y][y]), p))
                # dependence between activations between both regions defined by p_one and p_zero
                p_one=(thresh_array[x][y])/(thresh_array[y][y]) 
                p_zero=(thresh_array[x][x]-thresh_array[x][y])/(total_contrasts-thresh_array[y][y])
                # likelihood regions are functionally connected
                alternate=(scipy.stats.binom.pmf(thresh_array[x][y], thresh_array[y][y], p_one))*(scipy.stats.binom.pmf((thresh_array[x][x]-thresh_array[x][y]), (total_contrasts-thresh_array[y][y]), p_zero))
                # calculation of p value
                p_val=(-2*(math.log10(null/alternate)))
                # setting connection between region X and Y to zero if insignificant
                if p_val > threshold:
                    thresh_array[x][y]=0
    return thresh_array


# Working with networkx graphs
def applycost_to_g(G,cost):
    """Threshold graph to achieve cost.

    Return the graph after the given cost has been applied. Leaves weights intact.

    Parameters
    ----------
    G: input networkx Graph

    cost: float
        Fraction or proportion of all possible undirected edges desired in the
        output graph.


    Returns
    -------
    G_cost: networkx Graph
        Networkx Graph with cost applied

    threshold: float
        Correlations below this value have been set to 0 in
        thresholded_corr_mat."""
        
    Complete_G=nx.complete_graph(G.number_of_nodes())
    n_edges_keep=int(Complete_G.number_of_edges()*cost)
    weights=[(G[x][y]['weight']) for x,y in G.edges_iter()]
    sorted_weights=sorted(weights, reverse=1)
    thresh=sorted_weights[n_edges_keep+1]
    remove_edgelist=[(x,y) for x,y in G.edges_iter() if G[x][y]['weight']<thresh]
    G_cost=G.copy()
    G_cost.remove_edges_from(remove_edgelist)
    return G_cost

def remove_edges_by_weight(G, max_weight):
    G_removed=G.copy()
    for x,y in G.edges_iter():
        if G_removed[x][y]['weight']<max_weight:
            G_removed.remove_edge(x,y)
            
    for x,y in G.degree_iter():
        if y<2:
            G_removed.remove_node(x)
    return G_removed 

def remove_weight_edge_attribute(G):
    G_binary=G.copy()
    for x,y in G_binary.edges_iter():
        del G_binary[x][y]['weight']
    return G_binary


def remove_edgeless_nodes(G):
    to_remove=[]
    degree = G.degree()
    for x in degree:
        if degree[x]==0:
            to_remove.append(x)

    G.remove_nodes_from(to_remove)
    return G
    
    
def plot_weight_histogram(G):
    histo=[G[x][y]['weight'] for x,y in G.edges()]
    return plt.hist(histo)

def build_binarized_graph(G):
    """Takes graph converts it to a np array, binarizes it, builds networkx graph with binary edges"""
    
    binary_array=nx.to_numpy_matrix(G)
    binary_array=np.where(binary_array>0, 1,0)
    binary_G=nx.from_numpy_matrix(binary_array)
    #if nx.is_connected(binary_G)==0:
     #   print "Graph is not connected, removing nodes"
      #  binary_G=_remove_edgeless_nodes(binary_G)
    #else:
#        print "Graph is connected"
    return binary_G
    

# Domain filtering
def domain_filter_keycodes(key_codes, studies_filtered_by_domain, domain):
    """ filters a keycodes list by behavioral domain
        
        Parameters
        ------------
    
        key_codes: list of lists
            list of key_codes lists by region
    
        studies_filtered_by_domain: dict
            keys = behavioral domain, values are key_codes that correspond to that behavioral domain 
        
        domain: selected behavioral domain string
        Includes: 'Memory', 'Working Memory', 'Emotion', 'Attention', 'Language', 'Vision', 'Audition'

    Returns
    ------------
    domain_filtered_codes: list of lists
                    list of key_codes lists by region filtered by a particular domain
    
    """
    domainlist=studies_filtered_by_domain[domain]
    domain_filtered_codes=[set(x) & set(domainlist) for x in key_codes]
    return domain_filtered_codes
    

#Analyses
def run_basic_metrics(G, top_n=5):
    """runs a bunch of basic metrics and returns a dict"""
    basic_metrics=dict()
    basic_metrics['degrees']=nx.degree(G)
    basic_metrics['cpl']=[nx.average_shortest_path_length(g) for g in nx.connected_component_subgraphs(G) if g.number_of_nodes()>1]
    basic_metrics['ccoeff']=nx.clustering(G)
    basic_metrics['degree_cent']=nx.degree_centrality(G)
    basic_metrics['between_cent']=nx.betweenness_centrality(G) #nodes w/ high betweeness is intersection of many short paths & controls information flow
    basic_metrics['eigenvector_cent']=nx.eigenvector_centrality_numpy(G) #takes into account interactions of differnt lengths & dispersions
    for x in ['degrees','degree_cent','between_cent','ccoeff', 'eigenvector_cent']:
        sorted_x = sorted(basic_metrics[x].items(), key=operator.itemgetter(1))
        tops=[]
        for y in range(top_n):
            tops.append(sorted_x[-top_n:][y][0])
        basic_metrics['top'+x]=tops
    
    return basic_metrics

def run_weighted_metrics(G, top_n=5, comm=True):
    """runs a bunch of basic metrics and returns a dict"""
    weighted_metrics=dict()
    weighted_metrics['degrees']=nx.degree(G, weight='weight')
    weighted_metrics['cpl']=[nx.average_shortest_path_length(g, weight='weight') for g in nx.connected_component_subgraphs(G) if g.number_of_nodes()>1]
    weighted_metrics['ccoeff']=nx.clustering(G, weight='weight')
    weighted_metrics['degree_cent']=nx.degree_centrality(G)
    weighted_metrics['between_cent']=nx.betweenness_centrality(G, weight='weight')
    weighted_metrics['eigenvector_cent']=nx.eigenvector_centrality_numpy(G) #takes into account interactions of differnt lengths & dispersions
    if comm==True:
        weighted_metrics['communicability']=nx.communicability(G) #(Estrada & Hatano, 2008)(Crofts & Higham, 2009)
    for x in ['degrees','degree_cent','between_cent','ccoeff','eigenvector_cent']:
        sorted_x = sorted(weighted_metrics[x].items(), key=operator.itemgetter(1))
        tops=[]
        for y in range(top_n):
            tops.append(sorted_x[-top_n:][y][0])
        weighted_metrics['top'+x]=tops
    
    return weighted_metrics
    
    
def build_influence_matrix(n_coactives_array):
    n_coactive_mat=n_coactives_array.copy()
    diagonal=n_coactive_mat.diagonal()
    a=n_coactive_mat/diagonal[:, np.newaxis] #dividing by row
    b=n_coactive_mat/diagonal #dividing by column
    influence_mat=a-b # positive: rows influence column (A infl. B) , negative: col influence row (B inf. A) --> only in the upper triangle 
    #influence_mat=np.triu(influence_mat)
    return influence_mat
    
def build_influence_digraph(n_coactives_array):
    influence_mat=build_influence_matrix(n_coactives_array)
    influence_di_mat=influence_mat*(influence_mat>0)
    influence_diG=nx.DiGraph(influence_di_mat)
    return influence_diG
    

def make_brainx_style_partition(community_part_dict):
    bx_part=[]
    for x in list(set(community_part_dict.values())):
        sub_part=[y for y in community_part_dict.keys() if community_part_dict[y]==x ]
        #I could make sub_part a set
        bx_part.append(sub_part)
    return bx_part

from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection

#Basic Metrics Radar Plots
def radar_factory(num_vars, frame='circle'):
    """Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle' | 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)
    # rotate theta such that the first axis is at the top
    theta += np.pi/2

    def draw_poly_patch(self):
        verts = unit_poly_verts(theta)
        return plt.Polygon(verts, closed=True, edgecolor='k')

    def draw_circle_patch(self):
        # unit circle centered on (0.5, 0.5)
        return plt.Circle((0.5, 0.5), 0.5)

    patch_dict = {'polygon': draw_poly_patch, 'circle': draw_circle_patch}
    if frame not in patch_dict:
        raise ValueError('unknown value for `frame`: %s' % frame)

    class RadarAxes(PolarAxes):

        name = 'radar'
        # use 1 line segment to connect specified points
        RESOLUTION = 1
        # define draw_frame method
        draw_patch = patch_dict[frame]

        def fill(self, *args, **kwargs):
            """Override fill so that line is closed by default"""
            closed = kwargs.pop('closed', True)
            return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super(RadarAxes, self).plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            return self.draw_patch()

        def _gen_axes_spines(self):
            if frame == 'circle':
                return PolarAxes._gen_axes_spines(self)
            # The following is a hack to get the spines (i.e. the axes frame)
            # to draw correctly for a polygon frame.

            # spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
            spine_type = 'circle'
            verts = unit_poly_verts(theta)
            # close off polygon by repeating first vertex
            verts.append(verts[0])
            path = Path(verts)

            spine = Spine(self, spine_type, path)
            spine.set_transform(self.transAxes)
            return {'polar': spine}

    register_projection(RadarAxes)
    return theta

def unit_poly_verts(theta):
    """Return vertices of polygon for subplot axes.

    This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
    """
    x0, y0, r = [0.5] * 3
    verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
    return verts

def make_fingerprint_data(region_list, domain_list, metric_dict, metric_tuple):
    data=[[] for l in range(len(metric_tuple)+1)]
    #data[0]=['Emotion','Language','Working Memory','Attention','Rest','Audition','Base','Memory','Vision']
    for x in domain_list:
        for metric in metric_tuple:
            for region in region_list:
                if region not in metric_dict[x][metric].keys():
                    metric_dict[x][metric].update({region: 0})
    measure_dict={metric: {region: {x: metric_dict[x][metric][region] for x in domain_list} for region in region_list} for metric in metric_tuple}
    data_raw=[[OrderedDict(measure_dict[metric][region]).values() for region in region_list] for metric in metric_tuple]
    data_tuple=tuple((metric, data_raw[n]) for n in range(len(data_raw)) for metric in metric_tuple)
    data[0]= OrderedDict(measure_dict[metric][region]).keys()
    data[1]=data_tuple[0]
    data[2]=data_tuple[6]
    data[3]=data_tuple[12]
    data[4]=data_tuple[18]
    data[5]=data_tuple[24]
    return data

#Plotting
def plot_pretty_adj_matrix(G, nodelist, tick_space=5):
    mat=nx.to_numpy_matrix(G, nodelist)
    tick_range=range(0,G.number_of_nodes(), tick_space)
    selected_labels=[nodelist[y] for x,y in enumerate(tick_range)]
    
    fig = plt.figure()
    a = fig.add_subplot(1, 1, 1)
    a.imshow(mat, interpolation='nearest')
    a.set_xticks(tick_range)
    a.set_xticklabels(selected_labels, rotation='vertical')
    a.set_yticks(tick_range)
    a.set_yticklabels(selected_labels)
    fig.colorbar(a.imshow(mat, interpolation='nearest'))
    return a
    
def plot_weight_histogram(G):
    histo=[G[x][y]['weight'] for x,y in G.edges()]
    return plt.hist(histo)
