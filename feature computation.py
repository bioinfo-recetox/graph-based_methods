#loading libraries
import numpy as np
import networkx as nx
import os
import json

# Define the path to the folder containing graphs and where the features will be saved
input_folder = "<path to folder with graphs>"
output_folder = "<path where will be saved features>"

# function for computing features
def compute_measures(G):
    # initialize dictionary
    measures = {}
    
    # Connectedness measures
    measures['average_degree_'] = round(sum(dict(G.degree()).values()) / len(G), 3)
    measures['average_cluster_coefficient'] = round(nx.average_clustering(G), 3)
    measures['number_of_triangles'] = sum(nx.triangles(G).values())/3
    measures['number_of_nodes'] = G.number_of_nodes() 
    measures['number_of_edges'] = G.number_of_edges()     
    measures['maximum_degree'] = max(dict(G.degree()).values()) 

    # Eccentricity
    eccentricities = nx.eccentricity(G)
    measures['radius'] = int(np.min(list(eccentricities.values())))
    measures['diameter'] = int(np.max(list(eccentricities.values())))
    measures['average_eccentricity'] = round(np.mean(list(eccentricities.values())), 3)

    # Weight-related measures
    weights = [val['weight'] for u, v, val in G.edges(data=True)]
    measures['total_weight'] = round(sum(weights), 3)
    measures['mean_weight'] = round(np.mean(weights), 3)
    measures['max_weight'] = round(max(weights), 3)
    measures['min_weight'] = round(min(weights), 3)
    measures['max_min_weight_ratio'] = round(max(weights) / min(weights), 3)
    measures['standard_deviation_weight'] = round(np.std(weights), 3)

    # Eigenvalues and energy
    adj_matrix = nx.to_numpy_array(G)
    eigenvalues = np.linalg.eigvalsh(adj_matrix)
    measures['largest_eigenvalue'] = round(max(eigenvalues), 3)
    measures['second_largest_eigenvalue'] = round(sorted(eigenvalues)[-2], 3)
    return measures

# Iterate through each file in the input folder
for file in os.listdir(input_folder):
    # Read the weighted edgelist into a graph
    graph = nx.read_weighted_edgelist(os.path.join(input_folder, file))

    # Compute features of the graph
    features = compute_measures(graph)

    # Extract the file name without extension
    file_name_without_extension = os.path.splitext(file)[0]

    # Write the computed features to a JSON file
    with open(output_folder + file_name_without_extension + ".json", "w") as json_file:
        json.dump(features, json_file)
        json_file.close()

