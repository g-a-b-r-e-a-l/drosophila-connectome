#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 16:50:57 2024

@author: gabrielgibberd
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_from_array(motif_counts_array, colour, title): 
    """
    Plot motif distributions from an array of motif counts.

    Parameters:
    - motif_counts_array (numpy.ndarray): Array containing motif counts.
    - colour (str): Colour for the plot.
    - title (str): Title of the plot.

    Returns:
    - None
    """
    indices = range(1, motif_counts_array.shape[1] + 1)
    
    normalized_array = motif_counts_array / motif_counts_array.sum(axis=1, keepdims=True)

    # Calculate mean and standard deviation of proportions
    means = np.mean(normalized_array, axis=0)
    std_devs = np.std(normalized_array, axis=0)
        
    plt.figure(figsize=(12, 6))
    plt.bar(indices, means, yerr=std_devs, capsize=5, width=0.8, color=colour)
    plt.title(title)
    plt.xlabel('Motif Type')
    plt.ylabel('Proportion of Motifs Present')
    plt.xticks(indices)

    plt.show()
    
def plot_from_arrays(motif_counts_array, random_counts_array, colours, title): 
    
    """
    Plot motif distributions from arrays of motif counts for comparison.

    Parameters:
    - motif_counts_array (numpy.ndarray): Array containing motif counts for the main graph.
    - random_counts_array (numpy.ndarray): Array containing motif counts for random graphs.
    - colours (list): List of colors for the plots.
    - title (str): Title of the plot.

    Returns:
    - None
    """
    
    indices = list(range(1, motif_counts_array.shape[1] + 1))  # Convert range to list
    bar_width = 0.35
    
    normalized_model_array = motif_counts_array / motif_counts_array.sum(axis=1, keepdims=True)
    normalized_random_array = random_counts_array / random_counts_array.sum(axis=1, keepdims=True)
    
    # Calculate mean and standard deviation of proportions
    model_means = np.mean(normalized_model_array, axis=0)
    model_std_devs = np.std(normalized_model_array, axis=0)
    
    random_means = np.mean(normalized_random_array, axis=0)
    random_std_devs = np.std(normalized_random_array, axis=0)
        
    plt.figure(figsize=(12, 6))
    plt.bar(indices, model_means, yerr=model_std_devs, capsize=5, width=bar_width, label='Entire ROI Graph', color=colours[0])
    plt.bar([ind + bar_width for ind in indices], random_means, yerr=random_std_devs, capsize=5, width=bar_width, label='Random Graph', color=colours[1])  # Add bar_width to each index

    plt.title(title)
    plt.xlabel('Motif Type')
    plt.ylabel('Proportion of Motifs Present')
    plt.xticks(indices)

    plt.show()