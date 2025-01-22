#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 19:06:20 2024

@author: gabrielgibberd
"""

import pandas as pd
import networkx as nx
import numpy as np
import random


# Load data from CSV into a DataFrame
csv_file_path = '/Users/gabrielgibberd/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCL 2024 Y3 T2/Principles of scientific computing/connectome_lets_go/traced-roi-connections.csv' 

def create_gLR_graph(entire_connectome_csv):
    """
    Create a directed graph representing connections in the connectome.

    Parameters:
    - entire_connectome_csv (str): Path to the CSV file containing the connectome data.

    Returns:
    - G (networkx.DiGraph): Directed graph representing connections in the connectome.
    """
    
    df = pd.read_csv(csv_file_path)
    edges_in_roi = df[(df['roi'] == 'gL(R)')]
    G = nx.from_pandas_edgelist(edges_in_roi, source='bodyId_pre', target='bodyId_post', edge_attr='weight', create_using=nx.DiGraph())
    return G

def subset_graph_by_weight (G, threshold_weight):
    """
    Split the input graph into two graphs based on edge weights.

    Parameters:
    - G (networkx.DiGraph): Input graph.
    - threshold_weight (float): Threshold value for edge weight.

    Returns:
    - G_more (networkx.DiGraph): Graph containing edges with weight above threshold.
    - G_less (networkx.DiGraph): Graph containing edges with weight equal to or below threshold.
    """
    
    G_more = nx.DiGraph()
    G_less = nx.DiGraph()
    for edge in G.edges(data = True):
        if edge[2]['weight'] > threshold_weight:  # Check if weight is more than 1
            G_more.add_edge(edge[0], edge[1], weight=edge[2]['weight']) 
        if edge[2]['weight'] <= 1:  # Check if weight is more than 1
            G_less.add_edge(edge[0], edge[1], weight=edge[2]['weight'])   
    return G_more, G_less
 

def initialise_motif_graphs():
    
    """
    Initialize a list of motif graphs.

    Returns:
    - motif_graph_list (list): List of networkx.DiGraph objects representing motif graphs.
    """
    
    motif_matricies = [np.array([[0,0,0],[1,0,0],[1,0,0]]), np.array([[0,0,0],[0,0,1],[1,0,0]]),
             np.array([[0,0,0],[0,0,0],[1,1,0]]), np.array([[0,0,1],[1,0,0],[1,0,0]]),
             np.array([[0,0,0],[1,0,1],[1,0,0]]), np.array([[0,1,1],[0,0,0],[1,0,0]]),
             np.array([[0,1,0],[0,0,1],[1,0,0]]), np.array([[0,0,1],[1,0,1],[1,0,0]]), 
             np.array([[0,1,1],[1,0,0],[1,0,0]]), np.array([[0,0,1],[1,0,0],[1,1,0]]), 
             np.array([[0,0,0],[1,0,1],[1,1,0]]), np.array([[0,1,1],[1,0,0],[1,1,0]]), 
             np.array([[0,1,1],[1,0,1],[1,1,0]])]

    motif_graph_list = []

    for i in motif_matricies:
        motif_graph_list.append(nx.from_numpy_array(i, create_using=nx.DiGraph))
        
    return motif_graph_list

#initialise the motif list once to save time
motif_graph_list = initialise_motif_graphs()

    
def motif_classifier(motif_graph_list, subgraph, motif_count):
    
    """
    Classify a subgraph based on motifs and update motif count.

    Parameters:
    - motif_graph_list (list): List of motif graphs.
    - subgraph (networkx.Graph): Subgraph to classify.
    - motif_count (numpy.ndarray): Array to store motif counts.

    Returns:
    - motif_count (numpy.ndarray): Updated motif count array.
    """
    
    for i in range(len(motif_graph_list)):
        if nx.is_isomorphic(subgraph, motif_graph_list[i]):
            motif_count[i] += 1
            
    return motif_count

def random_connected_nodes(G):
    
    """
    Generate a random connected subgraph from the input graph.

    Parameters:
    - G (networkx.DiGraph): Input graph.

    Returns:
    - subgraph (networkx.Graph): Random connected subgraph.
    """
    
    num_nodes = 3
    random_nodes = random.sample(list(G.nodes()), num_nodes)
    
    random_nodes = set(random_nodes)
        
    subgraph = nx.subgraph(G, random_nodes)
    
    if nx.is_weakly_connected(subgraph):
        return subgraph
   
    
def random_sampling(G, motif_graph_list, iterations):
    
    """
    Perform random sampling of motifs in the graph.

    Parameters:
    - G (networkx.DiGraph): Input graph.
    - motif_graph_list (list): List of motif graphs.
    - iterations (int): Number of iterations for sampling.

    Returns:
    - motif_count (numpy.ndarray): Array containing motif counts.
    """
    
    count = 0
    total_count = 0
    motif_count = np.zeros(13)
    explored_nodes = []
    
    while count < iterations:
        subgraph = random_connected_nodes(G)
        total_count += 1
        
        if subgraph is not None:
            sub_nodes_list = sorted(subgraph.nodes())
            if sub_nodes_list not in explored_nodes:
                motif_classifier(motif_graph_list, subgraph, motif_count)
                count += 1
                explored_nodes.append(sub_nodes_list)
    
    return motif_count



def random_graph_from_base(G):
    """
    Create a random graph with the same number of nodes and edges as the input graph.

    Parameters:
    - G (networkx.DiGraph): Input graph.

    Returns:
    - random_G (networkx.DiGraph): Random graph with same size as input graph.
    """
    
    random_G = nx.gnm_random_graph(len(G.nodes), len(G.edges), directed=True)
    
    return random_G


def repeated_sampling(G, iterations, sample_size):
    """
    Perform repeated random sampling of motifs in the graph.

    Parameters:
    - G (networkx.DiGraph): Input graph.
    - iterations (int): Number of iterations for sampling.
    - sample_size (int): Size of each sample.

    Returns:
    - motif_counts_array (numpy.ndarray): Array containing motif counts for each iteration.
    """
    
    motif_counts_array = np.zeros((iterations, 13))
    for i in range(iterations):
        motif_count = random_sampling(G, motif_graph_list, sample_size)
        for j in range(len(motif_count)):
            motif_counts_array[i][j] = motif_count[j]
    return motif_counts_array

def repeated_sampling_multiple_randoms(G, iterations, sample_size):
    """
    Perform repeated random sampling of motifs in multiple random graphs.

    Parameters:
    - G (networkx.DiGraph): Input graph.
    - iterations (int): Number of iterations for sampling.
    - sample_size (int): Size of each sample.

    Returns:
    - motif_counts_array (numpy.ndarray): Array containing motif counts for each iteration.
    """
    motif_counts_array = np.zeros((iterations, 13))
    for i in range(iterations):
        new_random_G = random_graph_from_base(G)
        motif_count = random_sampling(new_random_G, motif_graph_list, sample_size)
        for j in range(len(motif_count)):
            motif_counts_array[i][j] = motif_count[j]
    return motif_counts_array
    

def p_values(model_array, random_array):
    
    """
    Compute p-values for motifs based on model and random arrays. The p-value calculation is described in the investigation.

    Parameters:
    - model_array (numpy.ndarray): Array containing motif counts for the main graph.
    - random_array (numpy.ndarray): Array containing motif counts for random graphs.

    Returns:
    - None
    """
    
    model_means = np.mean(model_array, axis=0)
    counts = np.zeros(13)
    p_values_motifs= np.zeros(13)
    for i in range(model_array.shape[1]):
        for j in range(model_array.shape[0]):
            if random_array[j][i] > model_means[i]:
                counts[i] += 1
    for i in range(len(counts)):
        p_values_motifs[i] = counts[i]/model_array.shape[0]
        
    return p_values_motifs
    
        




