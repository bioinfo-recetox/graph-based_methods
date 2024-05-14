import numpy as np
import networkx as nx
import os
from scipy.spatial import Delaunay

# Define the path to the folder containing point coordinates and where the Delaunay and Gabriel graphs will be saved
input_path = '<path to folder with files with point coordinates>'
output_path_delaunay = '<path to folder where will be Delaunay graphs saved as weighted edgelists>'
output_path_gabriel = '<path to folder where will be Gabriel graphs saved as weighted edgelists>'

# Function to calculate the Euclidean distance between two points
def calculate_distance(point1, point2):
    return np.linalg.norm(point1 - point2)

# Function for generating a Delaunay graph from a set of points
def delaunayGraph(np_points):
    # Compute Delaunay triangulation
    triangulation = Delaunay(np_points)
    edges = []
    edge_weights = {}

    # Generate edges from Delaunay triangulation simplices
    for simplex in triangulation.simplices:
        edges.extend([(simplex[0], simplex[1]), (simplex[1], simplex[2]), (simplex[2], simplex[0])])

    # Get unique nodes from the edges
    nodes = np.unique(edges)

    # Create an empty undirected graph
    delaunayNetworkxGraph = nx.Graph()
    delaunayNetworkxGraph.add_nodes_from(nodes)

    # Calculate distances and assign weights to edges
    for edge in edges:
        point1, point2 = edge
        distance = calculate_distance(np_points[point1], np_points[point2])
        edge_weights[edge] = distance

    # Add edges with weights to the graph
    delaunayNetworkxGraph.add_edges_from(edges)
    nx.set_edge_attributes(delaunayNetworkxGraph, edge_weights, 'weight')
    return triangulation, delaunayNetworkxGraph

# Function for generating a Gabriel graph from a Delaunay graph
def gabriel_and_delaunay(np_points):
    # Compute Delaunay triangulation and Delaunay graph
    triangulation, delaunayNetworkxGraph = delaunayGraph(np_points)

    # Get the number of points
    n=len(np_points)

    # Create a helper graph for Gabriel graph construction
    help_graph = np.eye(n, dtype=bool)

    # Initialize Gabriel graph
    gabriel_graph = np.zeros((n, n), dtype=bool)

    # Initialize lists for edges and edge weights
    edges=[]
    edge_weights = {}

    # Iterate through Delaunay triangulation simplices
    for simplex in triangulation.simplices:
        for i in range(3):
            for j in range(i + 1, 3):
                point_i, point_j = simplex[i], simplex[j]

                # Check if the edge between point_i and point_j is already considered
                if help_graph[point_i, point_j] == True:
                    break

                # Calculate radius of the circumscribed circle
                radius = calculate_distance(np_points[point_i], np_points[point_j]) / 2

                # Calculate middle point of the edge
                middle_point = ((np_points[point_i][0] + np_points[point_j][0]) / 2, (np_points[point_i][1] + np_points[point_j][1]) / 2)

                # Check if the edge is a Gabriel edge
                is_gabriel_edge = True
                for k in range(n):
                    if k != point_i and k != point_j:
                        if np_points[k][0] > middle_point[0]-radius and np_points[k][0] < middle_point[0]+radius:
                            if np_points[k][1] > middle_point[1]-radius and np_points[k][1] < middle_point[1]+radius:
                                distance = calculate_distance(middle_point, np_points[k])
                                if  distance < radius:
                                    is_gabriel_edge = False
                                    break
                
                # If the edge is a Gabriel edge, add it to the Gabriel graph
                if is_gabriel_edge:
                    edges.extend([(point_i, point_j)])
                    gabriel_graph[point_i, point_j] = gabriel_graph[point_j, point_i] = True

                # Mark the edge as considered in the helper graph
                help_graph[point_i, point_j] = help_graph[point_j, point_i] = True


    # Create an empty undirected graph for the Gabriel graph
    gabrielNetworkxGraph = nx.Graph()

    # Get unique nodes from the edges
    nodes = np.unique(edges)
    gabrielNetworkxGraph.add_nodes_from(nodes)

    # Add edges with weights to the Gabriel graph
    for edge in edges:
        point1, point2 = edge
        distance = calculate_distance(np_points[point1], np_points[point2])
        edge_weights[edge] = distance
    gabrielNetworkxGraph.add_edges_from(edges)
    nx.set_edge_attributes(gabrielNetworkxGraph, edge_weights, 'weight')
    return gabrielNetworkxGraph, delaunayNetworkxGraph

# Iterate through each file in the input folder
for file in os.listdir(input_path):
    file_path = os.path.join(input_path, file)

    # Load data points from the file
    data_points = np.loadtxt(file_path)
    
    # Extract the file name without extension
    file_name_without_extension = os.path.splitext(file)[0]

    # Generate Gabriel graph from Delaunay graph
    gabriel_final, delaunay_final = gabriel_and_delaunay(data_points)

    # Write the graph as a weighted edgelist to a file
    nx.write_weighted_edgelist(delaunay_final, output_path_delaunay + 'delaunay_' + file_name_without_extension + '.txt')
    nx.write_weighted_edgelist(gabriel_final, output_path_gabriel + 'gabriel_' + file_name_without_extension + '.txt')
