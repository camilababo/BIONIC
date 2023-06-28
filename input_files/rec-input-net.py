import argparse

import numpy as np
import pandas as pd
import networkx as nx
from sklearn.neighbors import radius_neighbors_graph


def data_to_edgelist(expression_data_file, edgelist_file, sep=',', index_col=0, header=0, radius=0.5):
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
    :return: None
    """

    # Read the expression data file and transform it to a numpy array
    dataframe = pd.read_csv(expression_data_file, sep=sep, index_col=index_col, header=header)
    data_array = dataframe.values.astype(float)

    # Compute the graph
    graph = radius_neighbors_graph(data_array, radius, mode='distance', metric='euclidean')

    # Convert the graph to a NetworkX object
    G_object = nx.from_scipy_sparse_matrix(graph)

    mapping = {i: str(index) for i, index in enumerate(dataframe.index)}

    relabeled_graph = nx.relabel_nodes(G_object, mapping)

    # Save the graph as a GML file
    nx.write_weighted_edgelist(relabeled_graph, edgelist_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This function builds an edgelist file from expression data.')
    parser.add_argument('-d', '--expression_data_file', type=str, help='The file containing the gene expression data.')
    parser.add_argument('-o', '--edgelist_file', type=str, help='The file to save the edgelist.')
    parser.add_argument('-s', '--sep', type=str, default=',', help='The separator for the expression data file.')
    parser.add_argument('-i', '--index_col', type=int, default=0, help='The index column for the expression data file.')
    parser.add_argument('-e', '--header', type=int, default=0, help='The header row for the expression data file.')
    parser.add_argument('-r', '--radius', type=float, default=0.5,
                        help='The radius parameter for radius_neighbors_graph.')
    parser.add_argument('-f', '--function', type=str, default='data_to_edgelist')

    args = parser.parse_args()

    if args.function == 'data_to_edgelist':
        data_to_edgelist(args.expression_data_file, args.edgelist_file, args.sep, args.index_col, args.header,
                         args.radius)