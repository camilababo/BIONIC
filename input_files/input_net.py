# import itertools
# from itertools import combinations
# from scipy.stats import pearsonr
import pandas as pd
from tqdm import tqdm
import numpy as np
import random


def data_to_edgelist(expression_data_file, edgelist_file):
    """
    This function takes the expression data and creates an edgelist file for the network.

    :param expression_data_file: The expression data file
    :param edgelist_file: The edgelist file
    :return: an edgelist file
    """
    print("The function data_to_edgelist() is starting...")

    # Read the expression data
    df = pd.read_csv(expression_data_file, sep=',', index_col=0, header=0)

    # Get a list of all the genes
    all_genes = list(df.T.columns)

    # Remove genes of chloroplast and mitochondria organs (ATC and ATM)
    not_organ_genes = [x for x in all_genes if not x.startswith('ATC') and not x.startswith('ATM')]
    expression_data_df = df.loc[not_organ_genes]
    filtered_genes = list(expression_data_df.T.columns)

    # Calculate the correlation coefficient matrix
    corr_matrix = np.corrcoef(expression_data_df)

    # Get the upper triangle of the correlation coefficient matrix
    upper_triangle_indices = np.triu_indices(corr_matrix.shape[0], k=1)
    corr_values = corr_matrix[upper_triangle_indices]
    absolute_corr_values = np.abs(corr_values)

    # Calculate the threshold value for the 99.5 percentile
    threshold = np.percentile(absolute_corr_values, 99.5)
    print(f"The threshold is {threshold}.")

    # Create a list of edges and their weights
    edges = []
    num_edges = len(upper_triangle_indices[0])
    with tqdm(total=num_edges) as pbar:
        for i in range(num_edges):
            src_index = upper_triangle_indices[0][i]
            trg_index = upper_triangle_indices[1][i]
            src = filtered_genes[src_index]
            trg = filtered_genes[trg_index]
            weight = 1 - absolute_corr_values[i]

            # Add the edge and weight to the list if the weight is above the threshold
            if weight > threshold:
                edges.append((src, trg, weight))

            # Update the progress bar every 1000 edges
            if i % 1000 == 0:
                pbar.update(1000)

        # Write the list of edges to the output file
        with open(edgelist_file, 'w') as f:
            for src, trg, weight in edges:
                f.write(f"{src} {trg} {weight}\n")

    # Update the progress bar to 100% after the loop has finished
    pbar.update(num_edges - pbar.n)
    print("The function data_to_edgelist() is complete...")


def data_to_edgelist_positive(expression_data_file, edgelist_file):
    """
    This function takes the expression data and creates an edgelist file for the network with weights based on
    positive correlation values and not absolute values.
    :param expression_data_file:
    :param edgelist_file:
    :return: an edgelist file
    """

    print("The function data_to_edgelist() is starting...")

    # Read the expression data
    df = pd.read_csv(expression_data_file, sep=',', index_col=0, header=0)

    # Get a list of all the genes
    all_genes = list(df.T.columns)

    # Remove genes of chloroplast and mitochondria organs (ATC and ATM)
    not_organ_genes = [x for x in all_genes if not x.startswith('ATC') and not x.startswith('ATM')]
    expression_data_df = df.loc[not_organ_genes]
    filtered_genes = list(expression_data_df.T.columns)

    # Calculate the correlation coefficient matrix
    corr_matrix = np.corrcoef(expression_data_df)

    # Get the upper triangle of the correlation coefficient matrix
    upper_triangle_indices = np.triu_indices(corr_matrix.shape[0], k=1)
    corr_values = corr_matrix[upper_triangle_indices]

    # Calculate the threshold value for the 99.5 percentile
    threshold = np.percentile(corr_values, 99.5)
    print(f"The threshold is {threshold}.")

    # Create a list of edges and their weights
    edges = []
    num_edges = len(upper_triangle_indices[0])
    with tqdm(total=num_edges) as pbar:
        for i in range(num_edges):
            src_index = upper_triangle_indices[0][i]
            trg_index = upper_triangle_indices[1][i]
            src = filtered_genes[src_index]
            trg = filtered_genes[trg_index]
            weight = corr_values[i]

            # Add the edge and weight to the list if the weight is above the threshold
            if weight > threshold:
                edges.append((src, trg, weight))

            # Update the progress bar every 1000 edges
            if i % 1000 == 0:
                pbar.update(1000)

            # Write the list of edges to the output file
        with open(edgelist_file, 'w') as f:
            for src, trg, weight in edges:
                f.write(f"{src} {trg} {weight}\n")

            # Update the progress bar to 100% after the loop has finished
        pbar.update(num_edges - pbar.n)
        print("The function data_to_edgelist_positive() is complete...")


def generate_binary_edgelist(expression_data_file, edgelist_file):
    """
    This function takes the expression data and creates a binary edgelist file for the network.

    :param expression_data_file: The expression data file
    :param edgelist_file: The edgelist file
    :return: a binary edgelist file
    """
    print("The function generate_binary_edgelist() is starting...")

    # Read the expression data
    df = pd.read_csv(expression_data_file, sep=',', index_col=0, header=0)

    # Get a list of all the genes
    all_genes = list(df.T.columns)

    # Remove genes of chloroplast and mitochondria organs (ATC and ATM)
    not_organ_genes = [x for x in all_genes if not x.startswith('ATC') and not x.startswith('ATM')]
    expression_data_df = df.loc[not_organ_genes]
    filtered_genes = list(expression_data_df.T.columns)

    # Calculate the correlation coefficient matrix
    corr_matrix = np.corrcoef(expression_data_df)

    # Get the upper triangle of the correlation coefficient matrix
    upper_triangle_indices = np.triu_indices(corr_matrix.shape[0], k=1)

    # Create a list of edges
    edges = []
    num_edges = len(upper_triangle_indices[0])
    with tqdm(total=num_edges) as pbar:
        for i in range(num_edges):
            src_index = upper_triangle_indices[0][i]
            trg_index = upper_triangle_indices[1][i]
            src = filtered_genes[src_index]
            trg = filtered_genes[trg_index]

            # Add the edge to the list
            edges.append((src, trg))

            # Update the progress bar every 1000 edges
            if i % 1000 == 0:
                pbar.update(1000)

        # Write the list of edges to the output file
        with open(edgelist_file, 'w') as f:
            for src, trg in edges:
                f.write(f"{src} {trg}\n")

    # Update the progress bar to 100% after the loop has finished
    pbar.update(num_edges - pbar.n)
    print("The function generate_binary_edgelist() is complete...")


def generate_random_edgelist(expression_data_file, edgelist_file):
    """
    This function takes the expression data and creates a random edgelist file for the network.

    :param expression_data_file: The expression data file
    :param edgelist_file: The edgelist file
    :return: None
    """
    print("Creating a random edgelist from expression data...")

    # Read the expression data
    df = pd.read_csv(expression_data_file, sep=',', index_col=0, header=0)

    # Get a list of all the genes
    all_genes = list(df.T.columns)

    # Remove genes of chloroplast and mitochondria organs (ATC and ATM)
    not_organ_genes = [x for x in all_genes if not x.startswith('ATC') and not x.startswith('ATM')]
    expression_data_df = df.loc[not_organ_genes]
    filtered_genes = list(expression_data_df.T.columns)

    # Get the number of nodes
    num_nodes = len(filtered_genes)

    # Generate a list of all possible edges
    edges = []
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            edges.append((i, j))

    # Randomly shuffle the list of edges
    random.shuffle(edges)

    # Generate random weights for the edges
    weights = np.random.rand(len(edges))

    # Write the list of edges and weights to the output file
    with open(edgelist_file, 'w') as f:
        for (src_index, trg_index), weight in zip(edges, weights):
            src = filtered_genes[src_index]
            trg = filtered_genes[trg_index]
            f.write(f"{src} {trg} {weight}\n")

    print("Random edgelist creation complete.")


if __name__ == '__main__':
    # data_to_edgelist("../outputs/separate_conditions/16mad_flowering.csv", "../outputs/separate_conditions"
    #                                                                            "/16mad_weighted.txt")
    # data_to_edgelist("../outputs/separate_conditions/25mad_flowering.csv", "../outputs/separate_conditions"
    #                                                                            "/25mad_weighted.txt")
    generate_binary_edgelist("../outputs/separate_conditions/mad/nx-edgelist/16mad_flowering.csv",
                             "../outputs/separate_conditions/mad/16mad_binary.txt")
    # generate_binary_edgelist("../outputs/separate_conditions/25mad_flowering.csv", "../outputs"
    #                                                                                    "/separate_conditions"
    #                                                                                    "/25mad_binary.txt")
    # data_to_edgelist_positive('../outputs/separate_conditions/bens_rec/16mad_ben_flowering.csv', "../outputs"
    #                                                                                     "/separate_conditions"
    #                                                                                     "/16mad_positive.txt")
    # data_to_edgelist_positive('../outputs/separate_conditions/bens_rec/25mad_ben_flowering.csv', "../outputs"
    #                                                                                     "/separate_conditions"
    #                                                                                     "/25mad_positive.txt")
    # generate_binary_edgelist('../outputs/separate_conditions/scaled_twenty_five_mad_probe.csv', '25mad_binary.txt')
    # generate_random_edgelist('../outputs/separate_conditions/scaled_sixteen_mad_probe.csv',
    #                          '../bionic/random_run/16mad_random.txt')
    # generate_random_edgelist('../outputs/separate_conditions/scaled_twenty_five_mad_probe.csv',
    #                          '../bionic/random_run/25mad_random.txt')
