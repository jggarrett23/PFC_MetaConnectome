import pickle
import numpy as np
import networkx as nx
import pandas as pd
import brainmaptools as brainmap
import matplotlib.pyplot as plt
import community
import brainx.modularity as mod
import operator
import csv
import scipy as scipy
import scipy.io as io
import math as math
import bct


#BrainNet Viewer functions
def node_file(partition, y, filename, path='/Users/owner/Documents/MATLAB/Add_Ons/BrainNetViewer_20171031/Data/NodeFiles/CSV/Label_coordinates.csv'):
    
    """Function creates csv node file for BrainNet Viewer
   
    Parameters:
    -----------
    partition=Louvain style partition
    y= dict of metrics with regions as keys and metric as value. Sizes nodes according to value.
    filename= File will by default be sent to a prespecified directory as a .csv file with inputted name.
    Change variable outdir to make it a directory on computer.
    
    Return:
    -------
    df= Dataframe representing outputted .csv node file."""
    
    df=pd.read_csv(path, header=None, names=['x','y','z','regions'], sep='\t')
    df.insert(4,'module',0,True)
    df.insert(5,'size',0,True)
    df.insert(6,'nodes',0,True)
    df=df[['x','y','z','module','size','regions','nodes']]
    
    for region in partition.keys():
        for n in range(len(df)):
            if region == df.loc[n,'regions']:
                df.loc[n,'nodes']=region
                df.loc[n,'module']=partition[region]
                df.loc[n, 'size']=y[region]
    df=df[['x','y','z','module','size','nodes']]
    df=df[df.nodes != 0]
    df.index=range(len(df))
    outdir='/Users/owner/Documents/MATLAB/Add_Ons/BrainNetViewer_20171031/Data/NodeFiles/CSV/'
    df.to_csv(outdir+filename, header=False, index=False, index_label=False)
    
    return df

def sim_node_file(sim_partition, y, filename, path='/Users/owner/Documents/MATLAB/Add-Ons/BrainNetViewer_20171031/Data/NodeFiles/CSV/Label_coordinates.csv'):
    
    """"Function creates csv node file for BrainNet Viewer
    
    Parameters:
    -----------
    sim_partition= Simulated annealing OR Newman partition already indexed as node names.
    y= dict of metrics with regions as keys and metric as value. Sizes nodes according to value.
    filename= File will by default be sent to a prespecified directory as a .csv file with inputted name.
    Change variable outdir to make it a directory on computer.
    
    Return:
    -------
    df= Dataframe representing outputted .csv node file."""
    
    df=pd.read_csv(path, header=None, names=['x','y','z','regions'])
    df.insert(4,'module',1,True)
    df.insert(5,'size',1,True)
    df.insert(6,'nodes',0,True)
    df=df[['x','y','z','module','size','regions','nodes']]
    
    for mods in range(len(sim_partition)):
        for region in range(len(sim_partition[mods])):
            for n in range(len(fd)):
                if sim_partition[mods][region] == fd.loc[n, 'regions']:
                    df.loc[n, 'nodes']=sim_partition[mods][region]
                    df.loc[n, 'module']=mods
                    df.loc[n, 'size']=y[(sim_partition[mods][region])]
    df=df[['x','y','z','module','size','nodes']]
    df=df[df.nodes != 0]
    df.index=range(len(fd))
    outdir='/Users/owner/Documents/MATLAB/Add-Ons/BrainNetViewer_20171031/Data/NodeFiles/CSV/'
    df.to_csv(outdir+filename, header=False, index=False, index_label=False)
    
    return df
    
#Modularity Functions
    
def random_modular_matrix(G, iterations=1, sim=False):
    
    """Function creates a random modularity matrix to be used as a threshold for actual modularity
    matrix. Based off of Newman/Girvan Null Model Graph.
    
    Parameters:
    -----------
    G: Network X Graph
    iterations: number of iterations to run
    sim: Whether or not you want to run simmulated annealing 
    
    Returns:
    --------
    random_matrix= Array where T(ij) is the number of times two nodes show up in the same partition"""
    
    #Note: this function is very slow due to simulated annealing, especially if being run
    #on multiple graphs.
    
    rand_G={}
    for n in range(iterations):
        #create random graph with equal amount of nodes and edges as input
        rand_G[n]=nx.gnm_random_graph(G.number_of_nodes(), G.number_of_edges())
    #run newman modularity on graph (uses newman modularity function from bct)   
    rand_new={n: bct.modularity_und(np.squeeze(np.asarray(nx.to_numpy_matrix(G)))) for n in range(iterations)}
    #create newman modularity partition
    rand_new_part={n: dict(zip(G.nodes(), rand_new[n][0])) for n in range(iterations)}
    #converts newman modularity partition into one similar to louvain
    rand_new_partition={n: brainmap.make_brainx_style_partition(rand_new_part[n]) for n in range(iterations)}
    #runs louvain modularity on random graph
    rand_lou={n: community.best_partition(rand_G[n]) for n in range(iterations)}
    
    if sim == True:
        rand_sim={n: mod.simulated_annealing(rand_G[n]) for n in range(iterations)}
        
    #preallocate random modularity matrices to speed up processing
    random_matrixs={n: np.zeros([rand_G[n].number_of_nodes(),rand_G[n].number_of_nodes()]) for n in range(iterations)}
    #preallocate sum matrix counting the number of times random two nodes are in the same module
    sum_matrix=np.zeros([rand_G[n].number_of_nodes(),rand_G[n].number_of_nodes()])
    for n in range(iterations):
        for r,k in enumerate (rand_G[n].nodes()):
            for d,b in enumerate (rand_G[n].nodes()):
                if k != b:
                    if rand_lou[n][k] == rand_lou[n][b]:
                        random_matrixs[n][r][d] += 1
                    for v in range (len(rand_new_partition[n])):
                        if k in rand_new_partition[n][v]:
                            if b in rand_new_partition[n][v]:
                                random_matrixs[n][r][d] += 1
                    if sim == True:
                        for rpart in range(len(rand_sim[n][0])):
                            if k in rand_sim[n][0].index_as_node_names()[rpart]:
                                if b in rand_sim[n][0].index_as_node_names()[rpart]:
                                    random_matrixs[n][r][d] += 1
        sum_matrix=random_matrixs[n]+sum_matrix
    avg_rand_matrix=sum_matrix/iterations
    return sum_matrix, avg_rand_matrix, np.max

def consensus_modular_matrix(G, louvain_partition, newman_partition, random_matrix, simpart=False, sima_partition=None):
    
    """Function creates consensus modularity matrix for Graph that has been ran through multiple
    modularity detection algorithims.
    
    Parameters:
    -----------
    G: Network X Graph
    louvain_partition: Louvain parition of G
    newman_partition: Newman partition of G indexed as node names
    sima_partition: Simulated Annealing partition
    random_matrix: Random modular matrix of Null Graph equivalent to number of nodes and 
    edges of G.
    
    Returns:
    --------
    consensus_matrix: Consensus modularity matrix of G across multiple modularity methods"""
    
    consensus_matrix=np.zeros([G.number_of_nodes(), G.number_of_nodes()])
    for i,t in enumerate (G.nodes()):
        for j,z in enumerate (G.nodes()):
            if t != z:
                if louvain_partition[t] == louvain_partition[z]:
                    consensus_matrix[i][j] += 1
                for module in range (len(newman_partition)):
                    if t in newman_partition[module]:
                        if z in newman_partition[module]:
                            consensus_matrix[i][j] += 1
                if simpart == True:
                    for part in range(len(sima_partition[0])):
                        if t in sima_partition[0].index_as_node_names()[part]:
                            if z in sima_partition[0].index_as_node_names()[part]:
                                consensus_matrix[i][j] += 1
    
    #Threshold graph according to random null model matrix largest value                       
    if consensus_matrix[i][j] == np.amax(random_matrix):
        consensus_matrix[i][j]=0
    return consensus_matrix
    
def random_g_mod(G):
    """"Function generates random graphs of input and runs newman modularity on it
    
    Parameters:
    -----------
    G: UNWEIGHTED network X graph, necessary for newman modularity
    
    Return:
    -------
    newman modularity of random graph"""
    
    n=G.number_of_nodes()
    m=G.number_of_edges()
    rand_g=nx.gnm_random_graph(n,m)
    rand_g_mod=mod.newman_partition(rand_g)
    return rand_g_mod

def mutual_info(B,N):
    """Run mutual information on both random graphs newman partitions
    
    Parameters:
    ----------
    B: first random network x graph newman partition
    N: second random network x graph newman partition
    
    Return:
    -------
    mutual information"""
    x=mod.mutual_information(B.index, N.index)
    return x

def significant_regions(metric_dict, domain_list, metric_key, node_label_dict, z_percentage):
    """Calculates z-scores for each region based off of given graph metric.
    Returns dict of regions that have high z-scores
    for one domain (specialization).
    
    Parameters:
    -----------
    metric_dict=dictionary of region and its graph theory metric value.
    domain_list=list of domains representing each type of network.
    metric_key=metric you want to get significant regions for.
    node_label_dict=dict of all possible nodes in a graph
    z_percentage= confidence interval value based off of z percentage (i.e 95% = 1.96)
    
    Return:
    -------
    specialized_regions=dictionary of regions statistically significant for one or more domains{domain: [node, z-score, metric-score, pop mean, std}
    z_score_domain=dict of z-scores for each region for one domain {domain: node: z-score}"""
    
    z_score_domain={}
    specialized_regions={}
    for x in domain_list:
        mean_cent={x: np.mean(metric_dict[x][metric_key].values()) for x in domain_list}
        standdev={x: np.std(metric_dict[x][metric_key].values()) for x in domain_list}
        sem={x: scipy.stats.sem(metric_dict[x][metric_key].values()) for x in domain_list}
        z_score_domain[x]={node: float(metric_dict[x][metric_key][node] - mean_cent[x])/standdev[x] for node in metric_dict[x][metric_key].keys()}
        specialized_regions[x]=[(node, "Z = {}".format(z_score_domain[x][node]), metric_key+" = {}".format(metric_dict[x][metric_key][node]), "M = {}".format(mean_cent[x]), "STD = {}".format(standdev[x]), "SEM = {}".format(sem[x])) for node in z_score_domain[x].keys() if z_score_domain[x][node] > z_percentage]
    return specialized_regions, z_score_domain