#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 11:22:58 2024

@author: gabrielgibberd
"""
import pandas as pd
import networkx as nx



def remove_ROI_insufficient_connections(csv_file_path):
    
    """
    Remove ROIs with under 5000 connections from the DataFrame.

    Parameters:
    - csv_file_path (str): Path to the CSV file containing the connectome data.

    Returns:
    - roi_list_filtered (list): List of ROIs with sufficient connections.
    - filtered_df (pandas.DataFrame): DataFrame with ROIs filtered for sufficient connections.
    """
    
    df = pd.read_csv(csv_file_path)
    
    roi_lengths = df.groupby('roi').size()
    roi_lengths_filtered = roi_lengths[roi_lengths >= 5000]
    roi_list_filtered = roi_lengths_filtered.index.tolist()
    
    filtered_df = df[df['roi'].isin(roi_list_filtered)]


    return roi_list_filtered, filtered_df

def calculate_edge_to_node_ratio(roi_list_filtered, filtered_df):
    
    """
    Calculate the edge-to-node ratio for each ROI and find the ROI with the maximum ratio.

    Parameters:
    - roi_list_filtered (list): List of ROIs with sufficient connections.
    - filtered_df (pandas.DataFrame): DataFrame with ROIs filtered for sufficient connections.

    Returns:
    - max_ratio_ROI (str): ROI with the maximum edge-to-node ratio.
    """
    
    graph_list= []
    
    for roi in roi_list_filtered:
        # Filter the DataFrame for the current ROI
        roi_df = filtered_df[filtered_df['roi'] == roi]
    
        # Create a graph
        G = nx.from_pandas_edgelist(roi_df, source='bodyId_pre', target='bodyId_post', create_using=nx.DiGraph())
        graph_list.append(G)
        
    ratios = []
    ratio_dict = {}
    for i in range(len(graph_list)):
        
        edge_to_node_ratio = graph_list[i].number_of_edges() / graph_list[i].number_of_nodes()
        ratios.append(edge_to_node_ratio)
        ratio_dict[edge_to_node_ratio] =  roi_list_filtered[i]
    
    max_ratio = max(ratios)
    return ratio_dict[max_ratio]

       
    
    

#%%

