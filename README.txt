#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 12:21:46 2024

@author: gabrielgibberd
"""

# Drosophila Connectome Analysis

This project investigates the connectome of the Drosophila brain using network analysis techniques. The project is divided into several steps, each implemented as a Python script.

## Files

- **user_file.py**: This script is the main user-facing file that orchestrates the analysis process. It imports functions from other modules and executes them in sequence to perform the investigation.

- **ROI_selection_criteria.py**: This module contains functions related to selecting regions of interest (ROIs) from the connectome based on specific criteria, such as edge-to-node ratios.

- **functions_file.py**: This module contains various utility functions used in the analysis, such as creating graphs from CSV data, conducting motif analysis, and calculating p-values.

- **graph_plotting.py**: This module provides functions to plot motif distributions and compare them between different graph datasets.

- **drosophila_connectome.csv**: This CSV file contains the connectome data for the Drosophila brain.

## Usage

1. Ensure all Python files and the CSV file are in the same directory.

2. Run the `user_file.py` script. This will execute the analysis process in the following steps:

   - Select ROIs from the connectome based on edge-to-node ratios.
   - Create a graph from the entire ROI.
   - Conduct motif analysis on the entire ROI graph and a benchmark set of random graphs.
   - Subset the ROI graph into two graphs with different edge weights.
   - Conduct motif analysis on the less and more weighted graphs and compare the results.
   - Plot separate motif counts for the analysis of results.

## Motif Analysis Plots

The analysis generates several plots to visualize motif distributions and compare them between different graph datasets. The plots include:

- **Comparison of motif spectra**: Compares motif distributions between the entire ROI graph and a benchmark set of random graphs. Also calculates p-values for each motif count.
  
- **Comparison of motif spectra from less and more weighted graphs**: Compares motif distributions between graphs with edge weights below and above a threshold. Also calculates p-values for each motif count.

- **Motif distribution of less graph**: Plots the motif distribution for the less weighted graph.

- **Motif distribution of more graph**: Plots the motif distribution for the more weighted graph.

