# import itertools
# from itertools import combinations
# from scipy.stats import pearsonr
import pandas as pd
from tqdm import tqdm
import numpy as np


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
            weight = 1 - abs(corr_matrix[src_index, trg_index])

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


if __name__ == '__main__':
    data_to_edgelist('../outputs/previous_attemps/locus_tag_processed_scaled.csv', 'edgelist_corrected.txt')
