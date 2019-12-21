#Modularity Functions
    
def random_modular_matrix(G):
    
    """Function creates a random modularity matrix to be used as a threshold for actual modularity
    matrix. Based off of Newman/Girvan Null Model Graph.
    
    Parameters:
    -----------
    G: Network X Graph
    
    Returns:
    --------
    random_matrix= Array where T(ij) is the number of times two nodes show up in the same partition"""
    
    #Note: this function is very slow due to simulated annealing, especially if being run
    #on multiple graphs.
    
    rand_G=nx.gnm_random_graph(G.number_of_nodes(), G.number_of_edges())
    rand_new=mod.newman_partition(rand_G)
    rand_newman=rand_new.index_as_node_names()
    rand_lou=community.best_partition(rand_G)
    random_matrix=np.zeros([rand_G.number_of_nodes(),rand_G.number_of_nodes()])
    rand_sim=mod.simulated_annealing(rand_G)
    for r,k in enumerate (rand_G.nodes()):
        for d,b in enumerate (rand_G.nodes()):
            if r != d:
                if rand_lou[k] == rand_lou[b]:
                    random_matrix[r][d] += 1
                for v in range (len(rand_newman)):
                    if t in rand_newman[v]:
                        if z in rand_newman[v]:
                            random_matrix[r][d] += 1
                for rpart in range(len(rand_sim[0])):
                    if t in rand_sim[0].index_as_node_names()[rpart]:
                        if z in rand_sim[0].index_as_node_names()[rpart]:
                            random_matrix[r][d] += 1
    return random_matrix

def consensus_modular_matrix(G, louvain_partition, newman_partition, sima_partition, random_matrix):
    
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
                for part in range(len(sima_partition[0])):
                    if t in sima_partition[0].index_as_node_names()[part]:
                        if z in sima_partition[0].index_as_node_names()[part]:
                            consensus_matrix[i][j] += 1
    
    #Threshold graph according to random null model matrix largest value                       
    if consensus_matrix[i][j] == np.amax(random_matrix):
        consensus_matrix[i][j]=0
    return consensus_matrix
