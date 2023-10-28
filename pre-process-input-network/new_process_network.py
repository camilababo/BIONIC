import argparse
import math
import os
import pandas as pd
import networkx as nx
from sklearn.neighbors import radius_neighbors_graph, NearestNeighbors


def pre_process_data(expression_data_file, sep=',', index_col=0, header=0):
    """
    This function reads a csv file and pre-processes the data for the conversion to an edgelist file.
    The pre-processing includes deleting non-relevant probes (ATM, ATC, control probes, etc.), deleting probes with
    empty values and performing a log2 transformation on the data.

    :param expression_data_file: The expression data file
    :param sep: The separator used in the csv file. (Default: ,)
    :param index_col: The index column for the csv file. (Default: 0)
    :param header: The header row for the csv file. (Default: 0)
    :return: The updated dataframe
    """
    print("Function pre_process_data() is running...")

    # Read the file
    df = pd.read_csv(expression_data_file, sep=sep, index_col=index_col, header=header)

    # Filter non-relevant probes
    # Delete rows that have at least on substring (probe aliases) that start with ATM, ATC or AFFX-
    indices_to_drop = []
    for index in df.index:
        for value in index.split(';'):
            if value.startswith(("ATM", "ATC", "AFFX")):
                indices_to_drop.append(index)
                break

    df.drop(indices_to_drop, axis=0, inplace=True)

    # Delete probes with empty values
    filtered_df = df.dropna(axis=1, how='all')

    # Log2 transformation
    transformed_df = filtered_df.applymap(lambda x: math.log2(x + 1))

    # Save the updated dataframe
    directory_path = os.path.dirname(expression_data_file)
    file_name = os.path.basename(expression_data_file)
    updated_file_path = os.path.join(directory_path, file_name.replace('.csv', '_preprocessed.csv'))

    transformed_df.to_csv(updated_file_path, sep=sep)

    print("Function pre_process_data() has finished running.")


def data_to_edgelist(expression_data_file,
                     edgelist_file,
                     sep=',',
                     index_col=0,
                     header=0,
                     radius=0.5,
                     node_list_file=None):
    """
    This function builds an edgelist file from expression data. The network is built using
    sklearn.neighbors.radius_neighbors_graph, converting the network to a NetworkX object by nx.from_scipy_sparse_array
     and saving it as a weighted edgelist through nx.write_weighted_edgelist.

    :param index_col: number of the column to be used as index
    :param header: number of the row to be used as header
    :param expression_data_file: The expression data file
    :param sep: separator for the expression data file
    :param edgelist_file: The edgelist file
    :param radius: The radius parameter for radius_neighbors_graph
    :param node_list_file: txt file containing an entries that should be added to the graph independently of previous
    calculations.
    :return: None
    """
    if expression_data_file is None:
        raise ValueError("expression_data_file cannot be None.")

    if edgelist_file is None:
        raise ValueError("edgelist_file cannot be None.")

    if radius <= 0:
        raise ValueError("radius must be greater than 0.")

    print("Function data_to_edgelist() is running...")

    # Read the expression data file and transform it to a numpy array
    dataframe = pd.read_csv(expression_data_file, sep=sep, index_col=index_col, header=header)
    data_array = dataframe.values.astype(float)

    # Compute the graph through a nearest neighbors approach
    graph = radius_neighbors_graph(data_array, radius, mode='distance', metric='euclidean')

    # Convert the graph to a NetworkX object
    G_object = nx.from_scipy_sparse_matrix(graph)

    # Remove all isolate nodes from the graph
    G_object.remove_nodes_from(list(nx.isolates(G_object)))

    if node_list_file is not None:
        with open(node_list_file, 'r') as file:
            node_list = [line.strip() for line in file]

        print("Number of nodes in relabeled_graph before adding nodes:", len(G_object.nodes))

        # Create an instance of NearestNeighbors with the same parameters as used in the radius_neighbors_graph method
        nearest_neighbors = NearestNeighbors(radius=radius, metric='euclidean')
        nearest_neighbors.fit(data_array)

        for node in node_list:
            if node not in G_object.nodes:
                if node in dataframe.index:
                    print("Adding node:", node)
                    index = dataframe.index.get_loc(node)
                    node_value = data_array[index]
                    G_object.add_node(node)

                    # Find the neighbors of the new node using the radius_neighbors method
                    distances, indices = nearest_neighbors.radius_neighbors(node_value.reshape(1, -1))

                    # Connect the new node to its neighbors in the relabeled_graph
                    for neighbor_index, distance in zip(indices[0], distances[0]):
                        if dataframe.index[neighbor_index] in G_object.nodes:
                            neighbor_node = dataframe.index[neighbor_index]
                            G_object.add_edge(node, neighbor_node, weight=distance)

                else:
                    print(f"Node {node} not found in the expression data file.")

        print("Number of nodes in relabeled_graph after adding nodes:", G_object.number_of_nodes())

    num_nodes = G_object.number_of_nodes()
    average_degree = sum(dict(G_object.degree()).values()) / num_nodes

    print(f"The network has {num_nodes} nodes, {G_object.number_of_edges()} edges and an average degree of "
          f"{average_degree:.1f}.")

    mapping = {i: str(index) for i, index in enumerate(dataframe.index)}

    relabeled_graph = nx.relabel_nodes(G_object, mapping)

    # Save the graph as a GML file
    nx.write_weighted_edgelist(relabeled_graph, edgelist_file)

    print("Function data_to_edgelist() finished.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This function builds an edgelist file from expression data.')
    parser.add_argument('-d', '--expression_data_file', type=str, help='The file containing the gene expression data.')
    parser.add_argument('-e', '--edgelist_file', type=str, help='The file to save the edgelist.')
    parser.add_argument('-s', '--sep', type=str, default=',', help='The separator for the expression data file.')
    parser.add_argument('-i', '--index_col', type=int, default=0, help='The index column for the expression data file.')
    parser.add_argument('-a', '--header', type=int, default=0, help='The header row for the expression data file.')
    parser.add_argument('-r', '--radius', type=float, default=0.5,
                        help='The radius parameter for radius_neighbors_graph.')
    parser.add_argument('-n', '--node', type=str, help='The file containing nodes that should be added independently '
                                                       'of computation.')
    parser.add_argument('-f', '--function', type=str, default='data_to_edgelist')

    args = parser.parse_args()

    if args.function == 'pre_process_data':
        pre_process_data(args.expression_data_file, args.sep, args.index_col, args.header)

    elif args.function == 'data_to_edgelist':
        data_to_edgelist(args.expression_data_file, args.edgelist_file, args.sep, args.index_col, args.header,
                         args.radius, args.node)
